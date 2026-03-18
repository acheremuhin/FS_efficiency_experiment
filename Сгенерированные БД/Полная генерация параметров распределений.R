library(philentropy)
library(writexl)
library(lhs)
library(scales)
set.seed(0)

# --- Глобальные параметры ---
pop_size <- 10
n_generations <- 30
mutation_prob <- 0.2
x_list <- list(
  beta      = seq(1e-6, 1 - 1e-6, length.out = 1000),
  gamma     = seq(1e-4, 20, length.out = 1000),
  lnorm     = seq(1e-4, 20, length.out = 1000),
  invgamma  = seq(1e-4, 20, length.out = 1000),
  betaprime = seq(1e-4, 20, length.out = 1000),
  t         = seq(-1000, 1000, length.out = 1000)
)

# --- JSD-дивергенция ---
js_divergence <- function(pdf_matrix) {
  if (nrow(pdf_matrix) < 2) return(0)
  pairs <- combn(nrow(pdf_matrix), 2)
  js_values <- apply(pairs, 2, function(idx) {
    i <- idx[1]; j <- idx[2]
    JSD(rbind(pdf_matrix[i, ], pdf_matrix[j, ]))
  })
  sum(js_values[is.finite(js_values)])
}

# --- Генетический алгоритм для распределений ---
run_ga <- function(name, pdf_fun, init_fun, mutate_fun, x, n_params, arg_names) {
  initialize_population <- function() replicate(pop_size, init_fun(), simplify = FALSE)
  tournament_selection <- function(pop, fitnesses, k = 3) {
    idx <- sample(1:length(pop), k)
    pop[[idx[which.max(fitnesses[idx])]]]
  }
  
  fitness_function <- function(param_set) {
    pdf_matrix <- t(sapply(param_set, function(pars) {
      args <- as.list(pars)
      names(args) <- arg_names[seq_along(args)]
      dens <- tryCatch(do.call(pdf_fun, c(list(x = x), args)), error = function(e) rep(0, length(x)))
      dens <- dens / sum(dens)
      if (any(is.na(dens)) || any(!is.finite(dens))) return(rep(0, length(x)))
      return(dens)
    }))
    js_divergence(pdf_matrix)
  }
  
  mutate_population <- function(parent) mutate_fun(parent)
  
  population <- initialize_population()
  best_fitness <- -Inf
  best_population <- NULL
  
  for (gen in 1:n_generations) {
    fitnesses <- sapply(population, fitness_function)
    current_fitness <- js_divergence(t(sapply(population, function(pars) {
      args <- as.list(pars)
      names(args) <- arg_names[seq_along(args)]
      dens <- tryCatch(do.call(pdf_fun, c(list(x = x), args)), error = function(e) rep(0, length(x)))
      dens <- dens / sum(dens)
      if (any(is.na(dens)) || any(!is.finite(dens))) return(rep(0, length(x)))
      return(dens)
    })))
    
    if (current_fitness > best_fitness) {
      best_fitness <- current_fitness
      best_population <- population
    }
    
    population <- lapply(1:pop_size, function(i) {
      parent1 <- tournament_selection(population, fitnesses)
      parent2 <- tournament_selection(population, fitnesses)
      if (runif(1) < 0.7) {
        weights <- runif(length(parent1))
        mapply(function(a, b, w) w * a + (1 - w) * b, parent1, parent2, weights, SIMPLIFY = TRUE)
      } else {
        mutate_population(parent1)
      }
    })
  }
  
  mat <- matrix(unlist(best_population), ncol = n_params, byrow = TRUE)
  colnames(mat) <- arg_names
  as.data.frame(mat)
}

# --- Вторичный ГА: 6 параметров (n, m_sig, m_nonsig и т.д.) ---
generate_lhs_population <- function(n) {
  param_ranges <- list(
    n_rows = c(10^1.0, 10^2.0),
    num_sig = c(1, 50),
    fac_sig = 0,
    num_nonsig = c(1, 50),
    fac_nonsig = 0,
    r_target = c(-0.99, 0.99)
  )
  sample <- randomLHS(n, length(param_ranges))
  result <- data.frame(matrix(nrow = n, ncol = 0))
  i <- 1
  for (name in names(param_ranges)) {
    lo <- param_ranges[[name]][1]; hi <- param_ranges[[name]][2]
    if (name == "n_rows") {
      result[[name]] <- round(lo * ((hi / lo)^sample[, i]))
    } else if (name == "r_target") {
      result[[name]] <- round(lo + (hi - lo) * sample[, i], 3)
    } else {
      result[[name]] <- round(lo + (hi - lo) * sample[, i])
    }
    i <- i + 1
  }
  result
}

meta_fitness <- function(df) {
  norm <- as.data.frame(lapply(df, rescale))
  dist_sum <- sum(dist(norm))
  dist_sum
}

run_secondary_ga <- function(n_gen = 30, pop_size = 10) {
  population <- generate_lhs_population(pop_size)
  best <- population
  best_score <- meta_fitness(best)
  
  for (g in 1:n_gen) {
    new_pop <- population
    for (i in 1:nrow(population)) {
      for (col in colnames(population)) {
        if (runif(1) < 0.2) {
          lo <- switch(col,
                       n_rows = 10^1, num_sig = 1, fac_sig = 0, num_nonsig = 1, fac_nonsig = 0, r_target = -0.99
          )
          hi <- switch(col,
                       n_rows = 10^2, num_sig = 50, fac_sig = 0, num_nonsig = 50, fac_nonsig = 0, r_target = 0.99
          )
          if (col == "n_rows") {
            new_pop[i, col] <- round(lo * ((hi / lo)^runif(1)))
          } else if (col == "r_target") {
            new_pop[i, col] <- round(runif(1, lo, hi), 3)
          } else {
            new_pop[i, col] <- round(runif(1, lo, hi))
          }
        }
      }
    }
    score <- meta_fitness(new_pop)
    if (score > best_score) {
      best <- new_pop
      best_score <- score
    }
    population <- new_pop
  }
  best
}

# --- Пользовательские плотности ---
dinvgamma <- function(x, alpha, beta) {
  coef <- (beta^alpha) / gamma(alpha)
  coef * x^(-alpha - 1) * exp(-beta / x)
}
dbetaprime <- function(x, alpha, beta) {
  coef <- 1 / beta(alpha, beta)
  coef * x^(alpha - 1) / (1 + x)^(alpha + beta)
}

# --- Генерация всех листов ---
results <- list()

add_block <- function(name, pdf_fun, init_fun, mutate_fun, x, n_params, arg_names) {
  dist_params <- run_ga(name, pdf_fun, init_fun, mutate_fun, x, n_params, arg_names)
  meta_params <- run_secondary_ga(n_gen = n_generations, pop_size = pop_size)
  results[[name]] <<- cbind(dist_params, meta_params)
}

add_block("I_Beta", dbeta, function() runif(2, 0.5, 10), function(p) lapply(p, function(x) if (runif(1) < mutation_prob) runif(1, 0.5, 10) else x), x_list$beta, 2, c("shape1", "shape2"))
add_block("II_SymBeta", function(x, alpha) dbeta(x, alpha, alpha), function() runif(1, 0.1, 10), function(p) sapply(p, function(x) if (runif(1) < mutation_prob) runif(1, 0.1, 10) else x), x_list$beta, 1, c("alpha"))
add_block("III_Gamma", dgamma, function() runif(2, 0.5, 10), function(p) lapply(p, function(x) if (runif(1) < mutation_prob) runif(1, 0.5, 10) else x), x_list$gamma, 2, c("shape", "rate"))
add_block("IV_LogNormal", dlnorm, function() c(runif(1, -2, 2), runif(1, 0.1, 2)), function(p) list(if (runif(1) < mutation_prob) runif(1, -2, 2) else p[[1]], if (runif(1) < mutation_prob) runif(1, 0.1, 2) else p[[2]]), x_list$lnorm, 2, c("meanlog", "sdlog"))
add_block("V_InvGamma", dinvgamma, function() runif(2, 1.5, 10), function(p) lapply(p, function(x) if (runif(1) < mutation_prob) runif(1, 1.5, 10) else x), x_list$invgamma, 2, c("alpha", "beta"))
add_block("VI_BetaPrime", dbetaprime, function() runif(2, 0.5, 10), function(p) lapply(p, function(x) if (runif(1) < mutation_prob) runif(1, 0.5, 10) else x), x_list$betaprime, 2, c("alpha", "beta"))
add_block("VII_t", function(x, df) dt(x, df), function() runif(1, 1, 100), function(p) sapply(p, function(x) if (runif(1) < mutation_prob) runif(1, 1, 100) else x), x_list$t, 1, c("df"))

write_xlsx(results, "ga_distributions_with_meta.xlsx")

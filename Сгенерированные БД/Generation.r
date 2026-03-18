library(readxl)
library(MASS)

set.seed(123)
path <- "ga_distributions_with_meta.xlsx"
sheets <- excel_sheets(path)

generate_dataset <- function(params, dist_name, index) {
  n <- as.integer(params[["n_rows"]])
  num_sig <- as.integer(params[["num_sig"]])
  fac_sig <- 0
  num_nonsig <- as.integer(params[["num_nonsig"]])
  fac_nonsig <- 0
  r_target <- as.numeric(params[["r_target"]])
  
  # Извлечение параметров распределения
  dist_params <- as.numeric(params[1:(which(names(params) == "n_rows") - 1)])
  
  # Значимые количественные переменные
  X_num <- round(matrix(rnorm(n * num_sig), nrow = n),4)
  beta_num <- runif(num_sig, 0.5, 2)
  
  
  # Линейная комбинация
  lin_comb <- X_num %*% beta_num
  lin_comb <- scale(lin_comb)
  
  # Зависимая переменная с нужным распределением
  y_raw <- switch(dist_name,
                  "I_Beta" = rbeta(n, dist_params[1], dist_params[2]),
                  "II_SymBeta" = 2 * rbeta(n, dist_params[1], dist_params[1]) - 1,
                  "III_Gamma" = rgamma(n, shape = dist_params[1], rate = dist_params[2]),
                  "IV_LogNormal" = rlnorm(n, meanlog = dist_params[1], sdlog = dist_params[2]),
                  "V_InvGamma" = 1 / rgamma(n, shape = dist_params[1], rate = dist_params[2]),
                  "VI_BetaPrime" = rbeta(n, dist_params[1], dist_params[2]) / 
                    rbeta(n, dist_params[2], dist_params[1]),
                  "VII_t" = rt(n, df = dist_params[1]),
                  stop("Unknown distribution")
  )
  y_raw <- scale(y_raw)
  
  # Формирование зависимой переменной с заданной корреляцией
  y_final <- r_target * lin_comb + sqrt(1 - r_target^2) * y_raw
  y_final <- round(as.numeric(scale(y_final, center = TRUE, scale = TRUE)),4)
  
  # Шумовые количественные переменные
  X_num_noise <- round(matrix(rnorm(n * num_nonsig), nrow = n),4)
  
  # Собираем датафрейм
  df <- data.frame(target = y_final)
  for (i in 1:num_sig) df[[paste0("NS", i)]] <- X_num[, i]
  for (i in 1:num_nonsig) df[[paste0("NNS", i)]] <- X_num_noise[, i]
  
  # Имя файла
  filename <- paste0(dist_name, "_", index, ".csv")
  write.csv(df, file = filename, row.names = FALSE)
  return(filename)
}

# Обработка всех листов
all_files <- c()
for (sheet in sheets) {
  dist_name <- sheet
  df <- read_excel(path, sheet = sheet)
  for (i in 1:nrow(df)) {  # строки со 2-й
    row_params <- as.list(df[i, ])
    tryCatch({
      file <- generate_dataset(row_params, dist_name, i)
      all_files <- c(all_files, file)
    }, error = function(e) {
      message(paste("Error in", dist_name, "row", i, ":", e$message))
    })
  }
}

cat("Генерация завершена. Всего файлов:", length(all_files), "\n")

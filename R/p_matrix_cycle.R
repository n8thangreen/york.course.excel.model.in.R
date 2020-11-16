
#'
p_matrix_cycle <- function(p_matrix, age, cycle) {

  tpProg <- 0.01
  tpDcm <- 0.15
  tpDn_lookup <-
    c("(34,44]" = 0.0017,
      "(44,54]" = 0.0044,
      "(54,64]" = 0.0138,
      "(64,74]" = 0.0379,
      "(74,84]" = 0.0912,
      "(84,100]" = 0.1958)
  Effect <- 0.5

  # matrix containing transition probabilities for without_drug

  age_grp <- cut(age, breaks = c(34,44,54,64,74,84,100))

  tpDn <- tpDn_lookup[age_grp]

  p_matrix["without_drug",
           "Asymptomatic_disease", "Progressive_disease"] <- tpProg*cycle

  p_matrix["without_drug",
           "Asymptomatic_disease", "Dead"] <- tpDn

  p_matrix["without_drug",
           "Asymptomatic_disease", "Asymptomatic_disease"] <- 1 - tpProg*cycle - tpDn

  p_matrix["without_drug",
           "Progressive_disease", "Dead"] <- tpDcm + tpDn

  p_matrix["without_drug",
           "Progressive_disease", "Progressive_disease"] <- 1 - tpDcm - tpDn

  p_matrix["without_drug", "Dead", "Dead"] <- 1

  # matrix containing transition probabilities for with_drug

  p_matrix["with_drug",
           "Asymptomatic_disease", "Progressive_disease"] <- tpProg*(1 - Effect)*cycle

  p_matrix["with_drug",
           "Asymptomatic_disease", "Dead"] <- tpDn

  p_matrix["with_drug",
           "Asymptomatic_disease", "Asymptomatic_disease"] <-
    1 - tpProg*(1 - Effect)*cycle - tpDn

  p_matrix["with_drug",
           "Progressive_disease", "Dead"] <- tpDcm + tpDn

  p_matrix["with_drug",
           "Progressive_disease", "Progressive_disease"] <- 1 - tpDcm - tpDn

  p_matrix["with_drug", "Dead", "Dead"] <- 1

  return(p_matrix)
}


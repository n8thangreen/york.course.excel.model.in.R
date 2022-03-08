
#
p_matrix_cycle <- function(p_matrix, age, cycle,
                           tpProg = 0.01,
                           tpDcm = 0.15,
                           effect = 0.5) {

  tpDn_lookup <-
    c("(34,44]" = 0.0017,
      "(44,54]" = 0.0044,
      "(54,64]" = 0.0138,
      "(64,74]" = 0.0379,
      "(74,84]" = 0.0912,
      "(84,100]" = 0.1958)

  age_grp <- cut(age, breaks = c(34,44,54,64,74,84,100))

  tpDn <- tpDn_lookup[age_grp]

  # Matrix containing transition probabilities for without_drug

  p_matrix["Asymptomatic_disease", "Progressive_disease", "without_drug"] <- tpProg*cycle

  p_matrix["Asymptomatic_disease", "Dead", "without_drug"] <- tpDn

  p_matrix["Asymptomatic_disease", "Asymptomatic_disease", "without_drug"] <- 1 - tpProg*cycle - tpDn

  p_matrix["Progressive_disease", "Dead", "without_drug"] <- tpDcm + tpDn

  p_matrix["Progressive_disease", "Progressive_disease", "without_drug"] <- 1 - tpDcm - tpDn

  p_matrix["Dead", "Dead", "without_drug"] <- 1

  # Matrix containing transition probabilities for with_drug

  p_matrix["Asymptomatic_disease", "Progressive_disease", "with_drug"] <- tpProg*(1 - effect)*cycle

  p_matrix["Asymptomatic_disease", "Dead", "with_drug"] <- tpDn

  p_matrix["Asymptomatic_disease", "Asymptomatic_disease", "with_drug"] <-
    1 - tpProg*(1 - effect)*cycle - tpDn

  p_matrix["Progressive_disease", "Dead", "with_drug"] <- tpDcm + tpDn

  p_matrix["Progressive_disease", "Progressive_disease", "with_drug"] <- 1 - tpDcm - tpDn

  p_matrix["Dead", "Dead", "with_drug"] <- 1

  return(p_matrix)
}

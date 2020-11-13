
######################
# basic Markov model #
######################

## model set-up ----

# define the number of treatments and their names
t_names <- c("without_drug", "with_drug")
n_treatments <- 2 #length(t_names)

# define the number of states and their names
s_names  <- c("Asymptomatic_disease", "Progressive_disease", "Dead")
n_states <- 3 #length(s_names)

# specify number of cycles
n_cycles <- 5 #44

Initial_age <- 55

cAsymp <- 500
cDeath <- 1000
cDrug <- 1000
cProg <- 3000
uAsymp <- 0.95
uProg <- 0.75
oDr <- 0.06
cDr <- 0.06

trans_c_matrix <-
  matrix(c(0, 0, cDeath,
           0, 0, cDeath),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

state_c_matrix <-
  matrix(c(cAsymp, cProg, 0,
           cAsymp + cDrug, cProg, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

state_q_matrix <-
  matrix(c(uAsymp, uProg, 0,
           uAsymp, uProg, 0),
         byrow = TRUE,
         nrow = n_treatments,
         dimnames = list(t_names,
                         s_names))

state_q_matrix
state_c_matrix
trans_c_matrix

# transition probabilities
p_matrix <- array(data = 0,
                  dim = c(n_treatments, n_states, n_states),
                  dimnames = list(t_names,
                                  from = s_names,
                                  to = s_names))

# store population output for each cycle

pop <- array(data = NA,
             dim = c(n_treatments, n_cycles, n_states),
             dimnames = list(treatment = t_names,
                             cycle = NULL,
                             state = s_names))

pop[, cycle = 1, "Asymptomatic_disease"] <- 1000
pop[, cycle = 1, "Progressive_disease"] <- 0
pop[, cycle = 1, "Dead"] <- 0

trans <- array(data = NA,
               dim = c(n_treatments, n_cycles, n_states),
               dimnames = list(treatment = t_names,
                               cycle = NULL,
                               state = s_names))

trans[, cycle = 1, "Asymptomatic_disease"] <- 0
trans[, cycle = 1, "Progressive_disease"] <- 0
trans[, cycle = 1, "Dead"] <- 0

# add up costs and QALYs for each cycle at a time for each drug

cycle_costs <- array(NA,
                     dim = c(n_treatments, n_cycles),
                     dimnames = list(treatment = t_names,
                                     cycle = NULL))

cycle_QALYs <- array(NA,
                     dim = c(n_treatments, n_cycles),
                     dimnames = list(treatment = t_names,
                                     cycle = NULL))

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)


## run model ----

for (i in 1:n_treatments) {

  age <- Initial_age

  for (j in 2:n_cycles) {

    p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)

    p_matrix_trans <- p_matrix[treatment = i, , ]
    diag(p_matrix_trans) <- 0

    # matrix multiplication
    pop[treatment = i, cycle = j, ] <-
      pop[treatment = i, cycle = j - 1, ] %*% p_matrix[treatment = i, , ]

    trans[treatment = i, cycle = j, ] <-
      pop[treatment = i, cycle = j - 1, ] %*% p_matrix_trans[, ]

    age <- age + 1
  }

  cycle_QALYs[i, ] <-
    pop[treatment = i, , ] %*% state_q_matrix[treatment = i, ]

  cycle_costs[i, ] <-
    pop[treatment = i, , ] %*% state_c_matrix[treatment = i, ] +
    trans[treatment = i, , ] %*% trans_c_matrix[treatment = i, ]

  # discounted
  cycle_QALYs[i, ] <- cycle_QALYs[i, ] * 1/(1 + oDr)^(1:n_cycles)
  cycle_costs[i, ] <- cycle_costs[i, ] * 1/(1 + cDr)^(1:n_cycles)

  total_costs[i] <- sum(cycle_costs[treatment = i, ])
  total_QALYs[i] <- sum(cycle_QALYs[treatment = i, ])
}


## results ----

# incremental costs and QALYs of with_drug vs to without_drug
c_incr <- total_costs["with_drug"] - total_costs["without_drug"]
q_incr <- total_QALYs["with_drug"] - total_QALYs["without_drug"]

# incremental cost effectiveness ratio
ICER <- c_incr/q_incr


plot(x = q_incr, y = c_incr,
     xlim = c(0,20),
     ylim = c(0, 10e6),
     pch = 16, cex = 1.5,
     xlab = "QALY differential",
     ylab = "Cost differential (£)",
     frame.plot = FALSE)
abline(a = 0, b = 30000) #wtp



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
# starts at 1 not 0 as in spreadsheet model
n_cycles <- 46

Initial_age <- 55
n_pop <- 1000

cAsymp <- 500
cDeath <- 1000
cDrug <- 1000
cProg <- 3000
uAsymp <- 0.95
uProg <- 0.75
oDr <- 0.06
cDr <- 0.06
tpDcm <- 0.15

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

# state_q_matrix
# state_c_matrix
# trans_c_matrix

# transition probabilities
p_matrix <- array(data = 0,
                  dim = c(n_treatments, n_states, n_states),
                  dimnames = list(t_names,
                                  from = s_names,
                                  to = s_names))

# p_matrix[1, , ]

# only include disease -> death transition cost due to disease
p_matrix_trans <-
  matrix(c(0, 0, 0,
           0, 0, tpDcm,
           0, 0, 0),
         byrow = TRUE,
         nrow = n_states,
         dimnames = list(s_names,
                         s_names))

# store population output for each cycle

pop <- array(data = NA,
             dim = c(n_treatments, n_cycles, n_states),
             dimnames = list(treatment = t_names,
                             cycle = NULL,
                             state = s_names))

pop[, cycle = 1, "Asymptomatic_disease"] <- n_pop
pop[, cycle = 1, "Progressive_disease"] <- 0
pop[, cycle = 1, "Dead"] <- 0

# head(pop[1, , ])

trans <- array(data = NA,
               dim = c(n_treatments, n_cycles, n_states),
               dimnames = list(treatment = t_names,
                               cycle = NULL,
                               state = s_names))

trans[, cycle = 1, "Asymptomatic_disease"] <- 0
trans[, cycle = 1, "Progressive_disease"] <- 0
trans[, cycle = 1, "Dead"] <- 0

# add up costs and QALYs for each cycle at a time for each drug

cycle_empty_array <-
  array(NA,
        dim = c(n_treatments, n_cycles),
        dimnames = list(treatment = t_names,
                        cycle = NULL))

cycle_state_costs <- cycle_trans_costs <- cycle_empty_array
cycle_costs <- cycle_QALYs <- cycle_empty_array
LE <- LYs <- cycle_empty_array
cycle_QALE <- cycle_empty_array

total_costs <- setNames(c(NA, NA), t_names)
total_QALYs <- setNames(c(NA, NA), t_names)


## run model ----

for (i in 1:n_treatments) {

  age <- Initial_age

  for (j in 2:n_cycles) {

    p_matrix <- p_matrix_cycle(p_matrix, age, j - 1)

    # matrix multiplication
    pop[treatment = i, cycle = j, ] <-
      pop[treatment = i, cycle = j - 1, ] %*% p_matrix[treatment = i, , ]

    trans[treatment = i, cycle = j, ] <-
      pop[treatment = i, cycle = j - 1, ] %*% p_matrix_trans

    age <- age + 1
  }

  cycle_state_costs[i, ] <-
    (pop[treatment = i, , ] %*% state_c_matrix[treatment = i, ]) * 1/(1 + cDr)^(1:n_cycles - 1)

  cycle_trans_costs[i, ] <-
    (trans[treatment = i, , ] %*% trans_c_matrix[treatment = i, ]) * 1/(1 + cDr)^(1:n_cycles - 2)

  cycle_costs[i, ] <- cycle_state_costs[i, ] + cycle_trans_costs[i, ]

  LE[i, ] <- pop[treatment = i, , ] %*% c(1,1,0)

  LYs[i, ] <- LE[i, ] * 1/(1 + oDr)^(1:n_cycles - 1)

  cycle_QALE[i, ] <-
    pop[treatment = i, , ] %*% state_q_matrix[treatment = i, ]

  cycle_QALYs[i, ] <- cycle_QALE[i, ] * 1/(1 + oDr)^(1:n_cycles - 1)

  # include base year (jump at start or end of cycle)?

  total_costs[i] <- sum(cycle_costs[treatment = i, -1])
  total_QALYs[i] <- sum(cycle_QALYs[treatment = i, -1])
}


## results ----

# incremental costs and QALYs of with_drug vs to without_drug
c_incr <- total_costs["with_drug"] - total_costs["without_drug"]
q_incr <- total_QALYs["with_drug"] - total_QALYs["without_drug"]

# incremental cost effectiveness ratio
ICER <- c_incr/q_incr

plot(x = q_incr/n_pop, y = c_incr/n_pop,
     xlim = c(0, 1100/n_pop),
     ylim = c(0, 10e6/n_pop),
     pch = 16, cex = 1.5,
     xlab = "QALY differential",
     ylab = "Cost differential (£)",
     frame.plot = FALSE)
abline(a = 0, b = 30000)  # willingness-to-pay


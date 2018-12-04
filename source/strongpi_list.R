# write functions for all situations in scott et al
slist = list(function(xx) -1 + 9 * (xx[, 1] - 0.5)^2 - 5 * abs(xx[, 2]),
             function(xx) 20 * (xx[ ,1] - 0.75))

tdlist = NULL
tdlist[[1]] = list(atoms = c(0.5),
                   probs = c(1),
                   variances = c(0))
tdlist[[2]] = list(atoms = c(0.5, 1, 1.5),
                   probs = c(1/2, 1/3, 1/6),
                   variances = c(0, 0.1, 1))







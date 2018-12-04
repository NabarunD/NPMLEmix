# write functions for all situations in scott et al
slist = list(function(xx) -3 + 1.5 * xx[ ,1] + 1.5 * xx[ ,2],
             function(xx) -3.25 + 3.5 * xx[ ,1]^2 - 3.5 * xx[ ,2]^2,
             function(xx) -1.5 * (xx[ ,1] - 0.5)^2 - 5 * abs(xx[ ,2]),
             function(xx) -4.25 + 2 * xx[ ,1]^2 + 2 * xx[ ,2]^2 - 2 * xx[ ,1] * xx[ ,2],
             function(xx) rep(-3, nrow(xx)),
             function(xx) -1.5 * xx[ ,1],
             function(xx) 20 * (xx[ ,1] - 0.5),
             function(xx) 20 * (xx[ ,1] - 0.75))

tdlist = NULL
tdlist[[1]] = list(atoms = c(-2, 0, 2),
                   probs = c(0.48, 0.04, 0.48),
                   variances = c(1, 16, 1))
tdlist[[2]] = list(atoms = c(-1.25, 0, 1.25),
                   probs = c(0.4, 0.2, 0.4),
                   variances = c(2, 4, 2))
tdlist[[3]] = list(atoms = c(0, 0, 0),
                   probs = c(0.3, 0.4, 0.3),
                   variances = c(0.1, 1, 9))
tdlist[[4]] = list(atoms = c(-3, -1.5, 1.5, 3),
                   probs = c(0.2, 0.3, 0.3, 0.2),
                   variances = c(0.01, 0.01, 0.01, 0.01))
tdlist[[5]] = list(atoms = c(3),
                   probs = c(1),
                   variances = c(1))
tdlist[[6]] = list(atoms = c(0.5),
                   probs = c(1),
                   variances = c(0))
tdlist[[7]] = list(atoms = c(1),
                   probs = c(1),
                   variances = c(0))

# these guassian convolution define the (inner) mixing distribution of theta
# f1 is a further gaussian convolution of these, can be thought of as adding 1 to each variance

start_seed_matrix = matrix(c(2135, 3153, 6666, 8534,
                             2919, 331, 8869, 3303,
                             8519, 1417, 3879, 4221,
                             5854, 6302, 8429, 9839,
                             5558, 1408, 161, 5579), byrow = TRUE, ncol = 4)

covariate_seed_matrix = matrix(c(2851, 367, 2982, 2867,
                                 4089, 7295, 9140, 6274,
                                 7179, 3869, 9955, 5056,
                                 1357, 1191, 2585, 5700,
                                 2781, 6532, 6939, 4732), byrow = TRUE, ncol = 4)







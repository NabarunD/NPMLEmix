slist = list(function(xx) -2 + 3.5 * xx[ ,1]^2 - 3.5 * xx[ ,2]^2,
             function(xx) -3 + 1.5 * xx[ ,1] + 1.5 * xx[ ,2],
             function(xx) -1 + 9 * (xx[, 1] - 0.5)^2 - 5 * abs(xx[, 2]),
             function(xx) 20 * (xx[ ,1] - 0.75))

#slist = list(function(xx) rep(-3,nrow(xx)),
#            function(xx) rep(3,nrow(xx)),
#             function(xx) rep(-5,nrow(xx)),
#            function(xx) rep(0,nrow(xx)))


tdlist = NULL
tdlist[[1]] = list(atoms = c(-1.25, 0, 1.25),
                   probs = c(0.4, 0.2, 0.4),
                   variances = c(2, 4, 2))
tdlist[[2]] = list(atoms = c(0, 0, 0),
                   probs = c(0.3, 0.4, 0.3),
                   variances = c(0.1, 1, 9))
tdlist[[3]] = list(atoms = c(0.5, 1, 1.5),
                   probs = c(1/2, 1/3, 1/6),
                   variances = c(0, 0.1, 1))
tdlist[[4]] = list(atoms = c(-2, 0, 2),
                   probs = c(0.48, 0.04, 0.48),
                   variances = c(1, 16, 1))




# slist = list(function(xx) -3 + 1.5 * xx[ ,1] + 1.5 * xx[ ,2],
#              function(xx) -2 + 3.5 * xx[ ,1]^2 - 3.5 * xx[ ,2]^2,
#              function(xx) -2 + 2 * xx[ ,1]^2 + 2 * xx[ ,2]^2 - 2 * xx[ ,1] * xx[ ,2],
#              function(xx) rep(-3, nrow(xx)),
#              function(xx) -1.5 * (xx[ ,1] - 0.5)^2 - 5 * abs(xx[ ,2]),
#              function(xx) -1 + 9 * (xx[, 1] - 0.5)^2 - 5 * abs(xx[, 2]),
#              function(xx) 20 * (xx[ ,1] - 0.75))
# 
# tdlist = NULL
# tdlist[[1]] = list(atoms = c(-2, 0, 2),
#                    probs = c(0.48, 0.04, 0.48),
#                    variances = c(1, 16, 1))
# tdlist[[2]] = list(atoms = c(-1.25, 0, 1.25),
#                    probs = c(0.4, 0.2, 0.4),
#                    variances = c(2, 4, 2))
# tdlist[[3]] = list(atoms = c(0, 0, 0),
#                    probs = c(0.3, 0.4, 0.3),
#                    variances = c(0.1, 1, 9))
# tdlist[[4]] = list(atoms = c(0.5),
#                    probs = c(1),
#                    variances = c(0))
# tdlist[[5]] = list(atoms = c(0.5, 1, 1.5),
#                    probs = c(1/2, 1/3, 1/6),
#                    variances = c(0, 0.1, 1))
# 
# 
# 
# 
# 
# 





# odds_ratio
od <- function(x,n,n.size){
  # both Wald and Rao Score
  # x <- rdirichlet(1,c(1,1,1,1))
  #alt_joint_distrb <- x

  # alt_joint_distrb <- c(0.2166556, 0.4410802, 0.2019513, 0.1403129)
  #  construct prob( |z| >= 1.96 | when the alternative hypothesis is true)
  # Capital letters are reserved for the true alternative distribtuion
  cell_counts <- rmultinom(n, n.size, prob= c(x[1,1], x[1,2], x[1,3], x[1,4]))
  A <- cell_counts[1,1]/n.size
  B <- cell_counts[2,1]/n.size
  C <- cell_counts[3,1]/n.size
  D <- cell_counts[4,1]/n.size

  OddsRatio <- A*(1-A-B-C)/(B*C)
  lnOR <- log(OddsRatio)
  AsyVarOR <- (1/A + 1/B + 1/C + 1/(1 - A - B - C))
  # Variance under Null for Odds ratio
  Var_OR_Null <- n.size*(1/((A+B)*(A+C)) + 1/((A+B)*(B+D)) + 1/((A+C)*(C+D)) + 1/((B+D)*(C+D)))


   # z Wald of odds
  z_wald_lnOR <- lnOR/sqrt(AsyVarOR)
  count_wald_OR <- ifelse(abs(z_wald_lnOR) > 1.96,1,0)
  z_rao_OR <- lnOR/sqrt(Var_OR_Null)
  count_rao_OR <- ifelse(abs(z_rao_OR) > 1.96,1,0)

  output <- list(A=A,B=B,C=C,D=D, lnOR = lnOR,
                 z_wald_lnOR =  z_wald_lnOR,
                 count_wald_OR = count_wald_OR,
                 z_rao_OR = z_rao_OR,
                 count_rao_OR = count_rao_OR)
  return(output)
}

y <-od(rdirichlet(1,c(1,1,1,1)),n=100,n.size = 100)







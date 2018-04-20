# 04202018
library(gtools)

phi <- function(x,n,n.size){
  # both Wald and Rao Score
  # x <- rdirichlet(1,c(1,1,1,1))
  alt_joint_distrb <- x

  # alt_joint_distrb <- c(0.2166556, 0.4410802, 0.2019513, 0.1403129)
  #  construct prob( |z| >= 1.96 | when the alternative hypothesis is true)
  # Capital letters are reserved for the true alternative distribtuion
  cell_counts <- rmultinom(n, n.size, prob= c(x[1,1], x[1,2], x[1,3], x[1,4]))
  A <- cell_counts[1,i]/n.size
  B <- cell_counts[2,i]/n.size
  C <- cell_counts[3,i]/n.size
  D <- cell_counts[4,i]/n.size

  var_A <- A*(1-A)
  var_B<- B*(1-B)
  var_C <- C*(1-C)
  var_D <- D*(1-D)

  cov_AB <- -A*B
  cov_AC <- -A*C
  cov_AD <- -A*D
  cov_BC <- -B*C
  cov_BD <- -B*D
  cov_CD <- -C*D


  row_1 <- A + B;row_2 <- C + D
  col_1 <- A + C;col_2 <- B + D

  pop_phi <- ((A*D)-(B*C))/sqrt(row_1*row_2*col_1*col_2)

  G <- (A+B)*(C+D)*(A+C)*(B+D)
  DphiA <- D*G^(-0.5) - 0.5*pop_phi*G^(-1)*(C+D)*(B+D)*(2*A+B+C)
  DphiD <- A*G^(-0.5) - 0.5*pop_phi*G^(-1)*(A+B)*(A+C)*(2*D+B+C)
  DphiB <- -C*G^(-0.5) - 0.5*pop_phi*G^(-1)*(A+C)*(C+D)*(2*B+A+D)
  DphiC <- -B*G^(-0.5) - 0.5*pop_phi*G^(-1)*(A+B)*(B+D)*(2*C+A+D)

  pop_asy_var_phi <- (var_A*DphiA^2 + var_B*DphiB^2 + var_C*DphiC^2 + var_D*DphiD^2 +
                        2*cov_AB*DphiA*DphiB + 2*cov_AC*DphiA*DphiC + 2*cov_AD*DphiA*DphiD +
                        2*cov_BC*DphiB*DphiC + 2*cov_BD*DphiB*DphiD + 2*cov_CD*DphiC*DphiD)/n

  # z Wald of phi
  z_wald <- pop_phi/sqrt(pop_asy_var_phi)
  count_wald_phi <- ifelse(abs(z_wald) > 1.96,1,0)
  z_rao <- sqrt(n)*pop_phi
  count_rao_phi <- ifelse(abs(z_rao) > 1.96,1,0)

  output <- list(A=A,B=B,C=C,D=D,phi= pop_phi,pop_var=pop_asy_var_phi,
                 z_wald_phi=z_wald,cnt_phi_wald = count_wald_phi,z_rao_phi = z_rao,
                 cnt_phi_rao = count_rao_phi)
  return(output)
}
y <-phi(rdirichlet(1,c(1,1,1,1)),n=100,n.size = 100)
y <- as.data.frame(y)
head(y)

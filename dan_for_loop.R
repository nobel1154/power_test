# dan_for_loop
total_output <- data.frame()
for(j in 1:10){
  diroutput <- rdirichlet(1,c(1,1,1,1))

  for (i in 1:10){
    y <- phi(diroutput,n=100,n.size=100)
    total_output <- rbind(total_output,y)
  }
}
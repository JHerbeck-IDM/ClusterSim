
y <- c(0, 1)
z <- NULL

for (i in c(1:10)){
  
  z <- c(z, i)
  
  x <- runif(1, min = 0, max = 1)
  
  if (x > 0.9){
    
    y <- c(y, x+x)
    
  }
}

WriteMat = function(mat,filename)
{
    cat(dim(mat)[1], " ", dim(mat)[2],"\n",file=filename)
    write.table(mat, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
}
WriteVec = function(vec,filename)
{
    cat(length(vec),"\n",file=filename)
    write.table(vec, file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
}
ComputeTtest = function(i, datG1, datG2)
{
  t = t.test(datG1[,i], datG2[,i], var.equal=FALSE)
  x = datG1[,i]
  y = datG2[,i]
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  ny <- length(y)
  my <- mean(y)
  vy <- var(y)
  stderrx <- sqrt(vx/nx)
  stderry <- sqrt(vy/ny)
  stderr <- sqrt(stderrx^2 + stderry^2)
  df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
  return(list(effectsize = my-mx, effectsd = stderr, statistic = t$statistic, 
              p.value = t$p.value))
}


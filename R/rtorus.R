rtorus <-
function (n, r, R, ct = c(0, 0, 0)) 
{
    s <- runif(n, 0, 2 * pi)
    t <- runif(n, 0, 2 * pi)
    Rg <- runif(n, r, R)
    rg <- runif(n, 0, r)
    x <- cos(s) * (R + rg * cos(t))
    y <- sin(s) * (R + rg * cos(t))
    z <- r * sin(t)
    return(cbind(ct[1] + x, ct[2] + y, ct[3] + z))
}

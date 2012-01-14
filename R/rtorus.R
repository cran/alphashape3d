rtorus <-
function (n, r, R, ct = c(0, 0, 0), rotx = NULL) 
{
    s <- runif(n, 0, 2 * pi)
    t <- runif(n, 0, 2 * pi)
    Rg <- runif(n, r, R)
    rg <- runif(n, 0, r)
    x <- ct[1] + cos(s) * (R + rg * cos(t))
    y <- ct[2] + sin(s) * (R + rg * cos(t))
    z <- ct[3] + r * sin(t)
    if (!is.null(rotx)) {
        aux1 <- cos(rotx) * y - sin(rotx) * z
        aux2 <- sin(rotx) * y + cos(rotx) * z
        y <- aux1
        z <- aux2
    }
    return(cbind(x, y, z))
}

inashape3d <-
function (as3d, indexAlpha = 1, points) 
{
    points = matrix(as.numeric(points), ncol = 3)
    triangles <- as3d$triang
    x <- as3d$x
    if (class(indexAlpha) == "character" & (indexAlpha == "ALL" | 
        indexAlpha == "all")) 
        indexAlpha = 1:length(as3d$alpha)
    if (any(indexAlpha > length(as3d$alpha)) | any(indexAlpha <= 
        0)) {
        if (max(indexAlpha) > length(as3d$alpha)) 
            error = max(indexAlpha)
        else error = min(indexAlpha)
        stop(paste("indexAlpha out of bound : valid range = 1:", 
            length(as3d$alpha), ", problematic value = ", error, 
            sep = ""), call. = TRUE)
    }
    insideVect = NULL
    for (iAlpha in indexAlpha) {
        tr <- triangles[triangles[, 8 + iAlpha] == 2, c("tr1", 
            "tr2", "tr3")]
        nbIntersect = intersectiontriangle(tr, x, points)
        insideVect <- c(insideVect, list(nbIntersect%%2 != 0))
    }
    if (length(indexAlpha) == 1) 
        insideVect <- insideVect[[1]]
    return(insideVect)
}
intersectiontriangle <-
function (triangles, coord, points) 
{
    points = as.matrix(points, nc = 3)
    inside <- integer(dim(points)[1])
    tempdir = as.numeric(2 * runif(3) - 1)
    retour <- .C("pointinashape", as.integer(triangles), dim(triangles)[1], 
        as.numeric(coord), dim(coord)[1], as.numeric(points), 
        dim(points)[1], as.integer(inside), tempdir)
    return(retour[[7]])
}

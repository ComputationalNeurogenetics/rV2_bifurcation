findCircleCenter <- function(p1, p2, p3) {
  # Calculate the midpoints of the chords
  mid1 <- c((p1[1] + p2[1])/2, (p1[2] + p2[2])/2)
  mid2 <- c((p2[1] + p3[1])/2, (p2[2] + p3[2])/2)
  
  # Calculate the slopes of the lines (p1p2 and p2p3)
  slope1 <- (p2[2] - p1[2]) / (p2[1] - p1[1])
  slope2 <- (p3[2] - p2[2]) / (p3[1] - p2[1])
  
  # Calculate the slopes of the perpendicular bisectors
  perpSlope1 <- -1 / slope1
  perpSlope2 <- -1 / slope2
  
  # Calculate the center (intersection of the perpendicular bisectors)
  # Using the formula y = mx + c, where c is y - mx
  c1 <- mid1[2] - perpSlope1 * mid1[1]
  c2 <- mid2[2] - perpSlope2 * mid2[1]
  
  # Solving the two linear equations for x and y
  center_x <- (c2 - c1) / (perpSlope1 - perpSlope2)
  center_y <- perpSlope1 * center_x + c1
  
  return(c(center_x, center_y))
}

distance2D <- function(point1, point2) {
  # Extracting coordinates
  x1 <- point1[1]
  y1 <- point1[2]
  x2 <- point2[1]
  y2 <- point2[2]
  
  # Calculating distance
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  
  return(distance)
}

# Function to find centroids of polygons in a spatstat 'owin' object
findCentroids <- function(window) {
  if (!inherits(window, "owin")) stop("Input must be an owin object")
  
  # If the window is a single polygon
  if (window$type == "polygonal") {
    x <- mean(window$bdry[[1]]$x)
    y <- mean(window$bdry[[1]]$y)
    centroids <- data.frame(X = x, Y = y)
  } 
  
  return(centroids)
}

# convert spatstat objects to sp classes

owin2Polygons <- function(x, id="1") {
  stopifnot(is.owin(x))
  x <- as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

tess2SP <- function(x) {
  stopifnot(is.tess(x))
  y <- tiles(x)
  nam <- names(y)
  z <- list()
  for(i in seq(y))
    z[[i]] <- owin2Polygons(y[[i]], nam[i])
  return(SpatialPolygons(z))
}
""
owin2SP <- function(x) {
  stopifnot(is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}
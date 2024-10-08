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
  distance <- (sqrt((x2 - x1)^2 + (y2 - y1)^2))
  #distance <- (sqrt((x2 - x1)^2 + (y2 - y1)^2)*(211.30/2048))
  return(distance)
}

distanceToCircle <- function(p1, p2, p3, x, y) {
  # Step 1: Find the center of the circle
  center <- findCircleCenter(p1, p2, p3)
  
  # Step 2: Calculate the radius of the circle
  radius <- distance2D(center, p1)
  
  # Step 3: Calculate the distance from the given point to the center of the circle
  point <- c(x, y)
  dist_to_center <- distance2D(center, point)
  
  # Step 4: Calculate the distance from the point to the circle's circumference
  distance_to_circumference <- abs(dist_to_center - radius)
  
  return(distance_to_circumference)
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

owin2SP <- function(x) {
  stopifnot(is.owin(x))
  y <- owin2Polygons(x)
  z <- SpatialPolygons(list(y))
  return(z)
}

# Function to calculate the average diameter of cells based on pairwise distances
average_cell_diameter <- function(polygon) {
  coords <- polygon$coords
  
  if (nrow(coords) < 2) {
    stop("Not enough points to calculate a diameter")
  }
  
  total_dist <- 0
  count <- 0
  
  for (i in 1:(nrow(coords)-1)) {
    for (j in (i+1):nrow(coords)) {
      dist <- sqrt((coords[i,1] - coords[j,1])^2 + (coords[i,2] - coords[j,2])^2)
      total_dist <- total_dist + dist
      count <- count + 1
    }
  }
  
  # Calculate the average diameter
  average_diameter <- total_dist / count
  
  return(average_diameter)
}

# Define the function
aggregate_image_values <- function(rois.ss, im.1, im.2, im.3, cores = 1) {
  combined.agg <- mclapply(1:length(rois.ss), function(cell.i) {
    cell <- rois.ss[[cell.i]]
    
    # Process im.1
    im.1.values <- im.1[i = cell]
    im.1.values <- tibble(f.int = im.1.values[im.1.values > IQR(im.1.values, na.rm = TRUE) * 5], 
                          ch = 1, 
                          cell.i = cell.i)
    
    # Process im.2
    im.2.values <- im.2[i = cell]
    im.2.values <- tibble(f.int = im.2.values[im.2.values > IQR(im.2.values, na.rm = TRUE) * 5], 
                          ch = 2, 
                          cell.i = cell.i)
    
    # Process im.3
    im.3.values <- im.3[i = cell]
    im.3.values <- tibble(f.int = im.3.values[im.3.values > IQR(im.3.values, na.rm = TRUE) * 5], 
                          ch = 3, 
                          cell.i = cell.i)
    
    # Combine the processed data
    data.tb <- rbind(im.1.values, im.2.values, im.3.values)
    return(data.tb)
  }, mc.cores = cores)
  
  # Combine all the aggregated data into a single tibble
  comb.agg.tb <- do.call(rbind, combined.agg)
  
  return(comb.agg.tb)
}

# Function to process the combined aggregated data
process_combined_agg <- function(combined.agg, replicate_name) {
  
  # Summarize by cell and channel
  count.per.cell.channel <- combined.agg %>%
    group_by(cell.i, ch) %>%
    summarize(count.dots = mean(f.int, na.rm = TRUE), .groups = 'drop')
  
  # Normalize counts per cell
  count.per.cell.channel <- count.per.cell.channel %>%
    group_by(cell.i) %>%
    mutate(norm_count_dots = count.dots / sum(count.dots)) %>%
    ungroup()  # Ungroup after mutation
  
  # Add replicate information
  count.per.cell.channel$replicate <- replicate_name
  
  return(count.per.cell.channel)
}

# Define the function
process_cell_distances <- function(rois.ss, rois, three.points, count.per.cell.channel) {
  
  # Step 1: Calculate centroids of cells
  cell.centroids <- lapply(rois.ss, findCentroids)
  
  # Step 2: Calculate average cell size
  cell.size.avg <- mean(sapply(rois, average_cell_diameter))
  
  # Step 3: Compute distances of cell centroids to a reference circle
  cell.distances <- sapply(cell.centroids, function(cent) {
    distanceToCircle(
      p1 = c(three.points$x[1], three.points$y[1]),
      p2 = c(three.points$x[2], three.points$y[2]),
      p3 = c(three.points$x[3], three.points$y[3]),
      x = cent$X, y = cent$Y
    )
  })
  
  # Step 4: Normalize distances by average cell size
  cell.distances <- cell.distances / cell.size.avg
  
  # Step 5: Create a tibble for cell distances
  dist.tb <- tibble(cell.i = 1:length(cell.distances), cell.distances = cell.distances)
  
  # Step 6: Merge the distances with count data and arrange by cell distances
  count.per.cell.channel <- left_join(count.per.cell.channel, dist.tb, by = "cell.i") %>%
    arrange(cell.distances)
  
  return(count.per.cell.channel)
}

# Define the function
process_and_filter_cell_data <- function(rois.ss, rois, three.points, count.per.cell.channel, quantile_threshold = 0.25) {
  
  # Step 1: Calculate centroids of cells
  cell.centroids <- lapply(rois.ss, findCentroids)
  
  # Step 2: Calculate average cell size
  cell.size.avg <- mean(sapply(rois, average_cell_diameter))
  
  # Step 3: Compute distances of cell centroids to a reference circle
  cell.distances <- sapply(cell.centroids, function(cent) {
    distanceToCircle(
      p1 = c(three.points$x[1], three.points$y[1]),
      p2 = c(three.points$x[2], three.points$y[2]),
      p3 = c(three.points$x[3], three.points$y[3]),
      x = cent$X, y = cent$Y
    )
  })
  
  # Step 4: Normalize distances by average cell size
  cell.distances <- cell.distances / cell.size.avg
  
  # Step 5: Create a tibble for cell distances and merge with count data
  dist.tb <- tibble(cell.i = 1:length(cell.distances), cell.distances = cell.distances)
  count.per.cell.channel <- left_join(count.per.cell.channel, dist.tb, by = "cell.i") %>%
    arrange(cell.distances)
  
  # Step 6: Compute the 25th percentile threshold of avg dots per cell
  min.thr <- quantile(
    count.per.cell.channel %>% 
      group_by(cell.i) %>% 
      summarise(avg.dots = mean(count.dots, na.rm = TRUE)) %>% 
      pull(avg.dots), 
    quantile_threshold, na.rm = TRUE
  )
  
  # Step 7: Filter cells that have an average count of dots greater than the threshold
  cells.incl <- count.per.cell.channel %>%
    group_by(cell.i) %>%
    summarise(avg.dots = mean(count.dots, na.rm = TRUE)) %>%
    dplyr::filter(avg.dots > min.thr) %>%
    pull(cell.i)
  
  # Step 8: Filter the original data to include only the selected cells
  count.per.cell.channel.filt <- dplyr::filter(count.per.cell.channel, cell.i %in% cells.incl)
  
  return(count.per.cell.channel.filt)
}
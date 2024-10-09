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
# Define the function to aggregate image values
aggregate_image_values <- function(rois.ss, im.1, im.2, im.3, cores = 1) {
  combined.agg <- mclapply(1:length(rois.ss), function(cell.i) {
    cell <- rois.ss[[cell.i]]
    
    # Subset data once for all channels
    im1_subset <- im.1[i = cell]
    im2_subset <- im.2[i = cell]
    im3_subset <- im.3[i = cell]
    
    # Precompute IQR thresholds
    im1_threshold <- IQR(im1_subset, na.rm = TRUE) * 5
    im2_threshold <- IQR(im2_subset, na.rm = TRUE) * 5
    im3_threshold <- IQR(im3_subset, na.rm = TRUE) * 5
    
    # Filter values based on thresholds
    im1_filtered <- im1_subset[im1_subset > im1_threshold]
    im2_filtered <- im2_subset[im2_subset > im2_threshold]
    im3_filtered <- im3_subset[im3_subset > im3_threshold]
    
    # If all filtered values are empty, return an empty tibble
    if (length(im1_filtered) == 0 && length(im2_filtered) == 0 && length(im3_filtered) == 0) {
      return(tibble(f.int = numeric(0), ch = integer(0), cell.i = integer(0)))
    }
    
    # Combine the filtered values into one tibble
    f.int <- c(im1_filtered, im2_filtered, im3_filtered)
    ch <- rep(1:3, times = c(length(im1_filtered), length(im2_filtered), length(im3_filtered)))
    
    data.tb <- tibble(f.int = f.int, ch = ch, cell.i = cell.i)
    
    return(data.tb)
  }, mc.cores = cores)
  
  # Combine all the aggregated data into a single tibble
  comb.agg.tb <- do.call(rbind, combined.agg)
  
  # Check for mismatch in cell count
  if (length(unique(comb.agg.tb$cell.i)) != length(rois.ss)) {
    warning("Some cells from rois.ss are missing in comb.agg.tb")
  }
  return(comb.agg.tb)
}

process_combined_agg <- function(combined.agg, replicate_name, rois.ss) {
  
  # Step 1: Summarize by cell and channel (ensure unique combinations of cell.i and ch)
  count.per.cell.channel <- combined.agg %>%
    group_by(cell.i, ch) %>%
    summarize(count.dots = mean(f.int, na.rm = TRUE), .groups = 'drop') %>%
    distinct(cell.i, ch, .keep_all = TRUE)  # Ensure each (cell.i, ch) pair is unique
  
  # Step 2: Ensure that every cell from `rois.ss` is represented, even if missing in `combined.agg`
  all_cells <- tibble(cell.i = 1:length(rois.ss))
  count.per.cell.channel <- right_join(all_cells, count.per.cell.channel, by = "cell.i") %>%
    mutate(count.dots = replace_na(count.dots, 0))  # Set missing dots to 0
  
  # Step 3: Normalize counts per cell (ensure that sums across channels equal 1)
  count.per.cell.channel <- count.per.cell.channel %>%
    group_by(cell.i) %>%
    mutate(total_count_dots = sum(count.dots, na.rm = TRUE)) %>%
    mutate(norm_count_dots = ifelse(total_count_dots > 0, count.dots / total_count_dots, 0)) %>%
    ungroup()  # Ungroup after mutation
  
  # Step 4: Add replicate information
  count.per.cell.channel$replicate <- replicate_name
  
  return(count.per.cell.channel)
}

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
  
  # Ensure that every cell is represented, even if some are missing
  count.per.cell.channel <- left_join(count.per.cell.channel, dist.tb, by = "cell.i") %>%
    arrange(cell.distances)
  
  return(count.per.cell.channel)
}

process_and_filter_cell_data <- function(rois.ss, rois, three.points, count.per.cell.channel, quantile_threshold = 0.25) {
  
  # Step 1: Calculate centroids of cells
  cell.centroids <- lapply(rois.ss, findCentroids)
  
  # Step 2: Calculate average cell size
  cell.size.avg <- mean(sapply(rois, average_cell_diameter))
  
  # Step 3: Compute distances of cell centroids to a reference circle
  cell.distances <- sapply(cell.centroids, function(cent) {
    # Check if centroid is valid before calling distanceToCircle
    if (is.null(cent$X) || is.null(cent$Y) || any(is.na(c(cent$X, cent$Y)))) {
      return(NA)  # Return NA for invalid centroids
    }
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
  
  # Step 6: Merge the distances with count data
  count.per.cell.channel <- left_join(count.per.cell.channel, dist.tb, by = "cell.i") %>%
    arrange(cell.distances)
  
  # Step 7: Handle NA distances (optional: remove or fill with a placeholder)
  count.per.cell.channel <- count.per.cell.channel %>%
    filter(!is.na(cell.distances))  # Remove rows where distances could not be calculated
  
  # Handle cases where count.dots is missing (use 0 for missing cells)
  count.per.cell.channel$count.dots <- replace_na(count.per.cell.channel$count.dots, 0)
  
  # Step 8: Compute the 25th percentile threshold of avg dots per cell
  min.thr <- quantile(
    count.per.cell.channel %>% 
      group_by(cell.i) %>% 
      summarise(avg.dots = mean(count.dots, na.rm = TRUE)) %>% 
      pull(avg.dots), 
    quantile_threshold, na.rm = TRUE
  )
  
  # Step 9: Filter cells that have an average count of dots greater than the threshold
  cells.incl <- count.per.cell.channel %>%
    group_by(cell.i) %>%
    summarise(avg.dots = mean(count.dots, na.rm = TRUE)) %>%
    dplyr::filter(avg.dots > min.thr) %>%
    pull(cell.i)
  
  # Step 10: Filter the original data to include only the selected cells
  count.per.cell.channel.filt <- dplyr::filter(count.per.cell.channel, cell.i %in% cells.incl)
  
  return(count.per.cell.channel.filt)
}

plot_sliding_mean <- function(count.per.cell.channel.filt, window_size = 8) {
  # Function to compute sliding mean for a given channel
  compute_sliding_mean <- function(channel, data, window_size) {
    filtered_data <- data %>%
      filter(ch == channel) %>%
      arrange(cell.distances)
    
    tibble(
      ch = channel,
      count.dots = slide_dbl(filtered_data %>% pull(count.dots), ~ mean(.x), .before = window_size),
      cell.distances = filtered_data %>% pull(cell.distances)
    )
  }
  
  # Compute sliding means for channels 1, 2, and 3
  ch1.slid.mean <- compute_sliding_mean(1, count.per.cell.channel.filt, window_size)
  ch2.slid.mean <- compute_sliding_mean(2, count.per.cell.channel.filt, window_size)
  ch3.slid.mean <- compute_sliding_mean(3, count.per.cell.channel.filt, window_size)
  
  # Combine the results
  slided.tb <- bind_rows(ch1.slid.mean, ch2.slid.mean, ch3.slid.mean)
  
  # Create the ggplot
  p <- ggplot(slided.tb, aes(x = cell.distances, y = count.dots)) + 
    geom_point(aes(colour = factor(ch))) + 
    scale_color_manual(values = c("magenta2", "dodgerblue", "forestgreen"), labels = c("Tal1", "Sox4", "Insm1")) + 
    labs(color = "Channel") + 
    theme_minimal() + 
    ggtitle(paste("Sliding mean (", window_size, ") values", sep = "")) + 
    geom_smooth(aes(colour = factor(ch)), method = "loess") + 
    xlab("Average cell diameters from VZ") + 
    ylab("Mean of Intensity") + 
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Rotate x-axis labels for readability
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    ) + 
    scale_x_continuous(breaks = seq(0, max(slided.tb$cell.distances), by = 5)) +  # Increase space between ticks
    ylim(0, NA)
  
  return(p)
}

plot_sliding_mean_replicates <- function(data, window_size = 8) {
  # Function to compute sliding mean for a given channel and replicate
  compute_sliding_mean <- function(channel, replicate_name, data, window_size) {
    filtered_data <- data %>%
      filter(ch == channel, replicate == replicate_name) %>%
      arrange(cell.distances)
    
    tibble(
      ch = channel,
      count.dots = slide_dbl(filtered_data %>% pull(count.dots), ~ mean(.x), .before = window_size),
      cell.distances = filtered_data %>% pull(cell.distances),
      replicate = replicate_name
    )
  }
  
  # Get the unique replicates from the dataset
  replicates <- unique(data$replicate)
  
  # Compute sliding means for each channel and replicate
  ch1.slid.mean <- bind_rows(lapply(replicates, function(rep) compute_sliding_mean(1, rep, data, window_size)))
  ch2.slid.mean <- bind_rows(lapply(replicates, function(rep) compute_sliding_mean(2, rep, data, window_size)))
  ch3.slid.mean <- bind_rows(lapply(replicates, function(rep) compute_sliding_mean(3, rep, data, window_size)))
  
  # Combine the results for all channels and replicates
  slided.tb <- bind_rows(ch1.slid.mean, ch2.slid.mean, ch3.slid.mean)
  
  # Create the ggplot with different symbols for each replicate and a fitted line per channel
  p <- ggplot(slided.tb, aes(x = cell.distances, y = count.dots, shape = replicate)) + 
    geom_point(aes(colour = factor(ch)), size = 2) +  # Increase dot size
    scale_shape_manual(values = c(16, 17, 18)) +  # Set different symbols for each replicate
    scale_color_manual(values = c("magenta2", "dodgerblue", "forestgreen"), labels = c("Tal1", "Sox4", "Insm1")) + 
    labs(color = "Channel", shape = "Replicate") + 
    theme_minimal() + 
    ggtitle(paste("Sliding mean (", window_size, ") values across replicates", sep = "")) + 
    geom_smooth(
      aes(x = cell.distances, y = count.dots, colour = factor(ch)),  
      method = "loess", se = TRUE, level = 0.95,  # Add confidence interval with 95% level
      inherit.aes = FALSE, 
      data = slided.tb %>% reframe(cell.distances = cell.distances, count.dots = count.dots, ch = ch),
      fill = "grey", alpha = 0.6  # Grey area for confidence interval
    ) + 
    xlab("Average cell diameters from VZ") + 
    ylab("Mean of Intensity") + 
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    ) + 
    scale_x_continuous(breaks = seq(0, max(slided.tb$cell.distances), by = 1)) + 
    ylim(0, NA)
  
  return(p)
}

# 1. Apply Symmetrical Sliding Window Approach on count.dots_norm
apply_sliding_window <- function(data, window_size = 6) {
  # Helper function to apply sliding window per replicate and channel
  compute_sliding_mean <- function(channel, replicate_name, data, window_size) {
    filtered_data <- data %>%
      filter(ch == channel, replicate == replicate_name) %>%
      arrange(cell.distances)
    
    tibble(
      ch = channel,
      # Apply symmetrical sliding window to count.dots_norm
      count.dots_norm_slid = slide_dbl(filtered_data %>% pull(count.dots_norm), ~ mean(.x), .before = window_size, .after = window_size),
      cell.distances = filtered_data %>% pull(cell.distances),
      replicate = replicate_name
    )
  }
  
  # Get unique replicates and channels
  replicates <- unique(data$replicate)
  
  # Apply sliding window for each channel and replicate
  ch1.slid.mean <- bind_rows(lapply(replicates, function(rep) compute_sliding_mean(1, rep, data, window_size)))
  ch2.slid.mean <- bind_rows(lapply(replicates, function(rep) compute_sliding_mean(2, rep, data, window_size)))
  ch3.slid.mean <- bind_rows(lapply(replicates, function(rep) compute_sliding_mean(3, rep, data, window_size)))
  
  # Combine the results for all channels and replicates
  slided.tb <- bind_rows(ch1.slid.mean, ch2.slid.mean, ch3.slid.mean)
  
  return(slided.tb)
}

# Wrapper function to process the entire pipeline
process_image_data <- function(rois.ss, im.1, im.2, im.3, rois, three.points, cores = 1, replicate_name, quantile_threshold = 0.25) {
  
  # Step 1: Aggregate image values (Find maximum intensity dots per cell per channel)
  combined.agg <- aggregate_image_values(rois.ss, im.1, im.2, im.3, cores = cores)
  
  # Step 2: Count dots per cell per channel
  count.per.cell.channel <- process_combined_agg(combined.agg, replicate_name = replicate_name, rois.ss)
  
  # Step 3: Calculate distances and filter cells
  count.per.cell.channel.filt <- process_and_filter_cell_data(rois.ss, rois, three.points, count.per.cell.channel, quantile_threshold = quantile_threshold)
  
  # Return the filtered data
  return(count.per.cell.channel.filt)
}
library(Seurat)
library(dplyr)
library(data.table)
library(sf)
library(terra)
library(progress)

options(future.globals.maxSize = 250010886400)
setwd("C:/Users/s215357/Desktop/STAR Protocols R Data")
convert_parquet_to_csv <- function(parquet_path, chunk_size = 1e6) {

  # Define the output CSV file path
  csv_output <- gsub('\\.parquet$', '.csv', parquet_path)

  # Read in the Parquet file as a Table (not DataFrame)
  parquet_file <- arrow::read_parquet(parquet_path, as_data_frame = FALSE)

  # Initialize start position
  start <- 0

  # Loop through the Parquet file and write each chunk to CSV
  while (start < parquet_file$num_rows) {

    # Define the end position of the current chunk
    end <- min(start + chunk_size, parquet_file$num_rows)

    # Slice the Parquet file to get the current chunk and convert it to a data frame
    chunk <- as.data.frame(parquet_file$Slice(start, end - start))

    # Write the chunk to CSV, appending if it's not the first chunk
    fwrite(chunk, csv_output, append = start != 0)

    # Move to the next chunk
    start <- end
  }

  # Return the output file path for reference
  return(csv_output)
}

convert_parquet_to_csv("C:/Users/s215357/Desktop/STAR Protocols R Data/004365__983/transcripts.parquet")

ReadXenium <- function (data.dir, outs = c("matrix", "microns"), type = "centroids",
                        mols.qv.threshold = 20)
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"),
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"),
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir,
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir,
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid,
                                     y = cell_info$y_centroid, cell = cell_info$cell_id,
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir,
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir,
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location,
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}

create_and_save_seurat <- function(data_dir, seurat_object_name, rds_file_name) {
  data <- ReadXenium(data.dir = data_dir, outs = c("matrix", "microns"), type = c("centroids", "segmentations"))
  
  segmentations.data <- list(
    centroids = CreateCentroids(data$centroids),
    segmentation = CreateSegmentation(data$segmentations)
  )
  
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = "Xenium"
  )
  
  seurat_obj <- CreateSeuratObject(
    counts = data$matrix[["Gene Expression"]],
    assay = "Xenium"
  )
  
  seurat_obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  seurat_obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  seurat_obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  seurat_obj[["fov"]] <- coords
  
  assign(seurat_object_name, seurat_obj, envir = .GlobalEnv)
  
  saveRDS(seurat_obj, file = rds_file_name)
  
}

create_and_save_seurat("C:/Users/s215357/Desktop/STAR Protocols R Data/0043635__983", 
                       "xenium_obj", "xenium_obj.rds")

# ---- IHC MCT1 (Opal 780) ----
# Set working directory
setwd("C:/Users/s215357/Desktop/External_Annotation_Integration_QuPath_Project/Exported_Channels")

# Load the TIFF image (assuming a single-channel grayscale image)
tiff_file <- "Channel_6.tiff"  # Update with the actual TIFF file name
image_raster <- rast(tiff_file)  # Use rast() instead of raster()

# Get raster dimensions
width <- ncol(image_raster)  # Number of columns = X dimension
height <- nrow(image_raster) # Number of rows = Y dimension

# Convert the raster to a matrix
image_matrix <- as.matrix(image_raster, wide=TRUE)  # Explicitly convert to matrix

# Create a data table of all pixel coordinates and intensity values
pixel_data <- data.table(
  x = rep(1:width, each = height),   # X-coordinates
  y = rep(height:1, times = width),  # Y-coordinates (flip to match image origin)
  intensity = as.vector(image_raster)  # Flatten the raster into a vector
)

setwd("C:/Users/s215357/Desktop/STAR Protocols R Data/0043635__983")
# Load the cell boundary data
boundary_data <- fread("cell_boundaries.csv")

# Define the conversion factor (adjust as needed)
micrometers_to_pixels <- 1 / 0.2125  

# Convert boundary data into spatial polygons
boundary_polygons <- boundary_data %>%
  mutate(
    vertex_x = vertex_x * micrometers_to_pixels,
    vertex_y = vertex_y * micrometers_to_pixels
  ) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_sfc(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
  st_as_sf(crs = NA)

# Set keys for fast subsetting
setkey(pixel_data, x, y)

# Initialize progress bar
pb <- progress_bar$new(
  format = "  Processing cell :current/:total [:bar] :percent in :elapsedfull",
  total = nrow(boundary_polygons), clear = FALSE, width = 60
)

# Initialize data frame to store results
mean_intensities <- data.frame(cell_id = character(), mean_intensity = numeric(), stringsAsFactors = FALSE)

# Process each cell polygon
for (i in seq_len(nrow(boundary_polygons))) {
  pb$tick()  # Update progress bar
  
  # Extract cell ID and polygon
  cell_id <- boundary_polygons$cell_id[i]
  cell_polygon <- boundary_polygons[i, ]
  
  # Get the bounding box
  bbox <- st_bbox(cell_polygon)
  xmin <- round(as.numeric(bbox["xmin"]))
  ymin <- round(as.numeric(bbox["ymin"]))
  xmax <- round(as.numeric(bbox["xmax"]))
  ymax <- round(as.numeric(bbox["ymax"]))
  
  # Subset pixels in the bounding box
  pixel_subset <- pixel_data[x %between% c(xmin, xmax) & y %between% c(ymin, ymax)]
  
  # Convert to sf object for spatial intersection
  pixel_points <- st_as_sf(pixel_subset, coords = c("x", "y"), crs = NA)
  
  # Check which pixels are inside the polygon
  pixels_in_polygon <- st_intersects(pixel_points, cell_polygon, sparse = FALSE)
  
  # Keep only the pixels inside the cell
  pixel_subset <- pixel_subset[as.logical(pixels_in_polygon)]
  
  # Compute mean intensity
  mean_intensity <- mean(pixel_subset$intensity, na.rm = TRUE)
  
  # Store the result
  mean_intensities <- rbind(mean_intensities, data.frame(cell_id = cell_id, mean_intensity = mean_intensity))
}

setwd("C:/Users/s215357/Desktop/STAR Protocols R Data")
# Save the final MFI per cell
write.csv(mean_intensities, "MFI_MCT1_xenium_obj.csv", row.names = FALSE)

# ---- IHC CD31 (Opal 520) ----
# Set working directory
setwd("C:/Users/s215357/Desktop/External_Annotation_Integration_QuPath_Project/Exported_Channels")

# Load the TIFF image (assuming a single-channel grayscale image)
tiff_file <- "Channel_7.tiff"  # Update with the actual TIFF file name
image_raster <- rast(tiff_file)  # Use rast() instead of raster()

# Get raster dimensions
width <- ncol(image_raster)  # Number of columns = X dimension
height <- nrow(image_raster) # Number of rows = Y dimension

# Convert the raster to a matrix
image_matrix <- as.matrix(image_raster, wide=TRUE)  # Explicitly convert to matrix

# Create a data table of all pixel coordinates and intensity values
pixel_data <- data.table(
  x = rep(1:width, each = height),   # X-coordinates
  y = rep(height:1, times = width),  # Y-coordinates (flip to match image origin)
  intensity = as.vector(image_raster)  # Flatten the raster into a vector
)

setwd("C:/Users/s215357/Desktop/STAR Protocols R Data/0043635__983")
# Load the cell boundary data
boundary_data <- fread("cell_boundaries.csv")

# Define the conversion factor (adjust as needed)
micrometers_to_pixels <- 1 / 0.2125  

# Convert boundary data into spatial polygons
boundary_polygons <- boundary_data %>%
  mutate(
    vertex_x = vertex_x * micrometers_to_pixels,
    vertex_y = vertex_y * micrometers_to_pixels
  ) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_sfc(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
  st_as_sf(crs = NA)

# Set keys for fast subsetting
setkey(pixel_data, x, y)

# Initialize progress bar
pb <- progress_bar$new(
  format = "  Processing cell :current/:total [:bar] :percent in :elapsedfull",
  total = nrow(boundary_polygons), clear = FALSE, width = 60
)

# Initialize data frame to store results
mean_intensities <- data.frame(cell_id = character(), mean_intensity = numeric(), stringsAsFactors = FALSE)

# Process each cell polygon
for (i in seq_len(nrow(boundary_polygons))) {
  pb$tick()  # Update progress bar
  
  # Extract cell ID and polygon
  cell_id <- boundary_polygons$cell_id[i]
  cell_polygon <- boundary_polygons[i, ]
  
  # Get the bounding box
  bbox <- st_bbox(cell_polygon)
  xmin <- round(as.numeric(bbox["xmin"]))
  ymin <- round(as.numeric(bbox["ymin"]))
  xmax <- round(as.numeric(bbox["xmax"]))
  ymax <- round(as.numeric(bbox["ymax"]))
  
  # Subset pixels in the bounding box
  pixel_subset <- pixel_data[x %between% c(xmin, xmax) & y %between% c(ymin, ymax)]
  
  # Convert to sf object for spatial intersection
  pixel_points <- st_as_sf(pixel_subset, coords = c("x", "y"), crs = NA)
  
  # Check which pixels are inside the polygon
  pixels_in_polygon <- st_intersects(pixel_points, cell_polygon, sparse = FALSE)
  
  # Keep only the pixels inside the cell
  pixel_subset <- pixel_subset[as.logical(pixels_in_polygon)]
  
  # Compute mean intensity
  mean_intensity <- mean(pixel_subset$intensity, na.rm = TRUE)
  
  # Store the result
  mean_intensities <- rbind(mean_intensities, data.frame(cell_id = cell_id, mean_intensity = mean_intensity))
}

setwd("C:/Users/s215357/Desktop/STAR Protocols R Data")
# Save the final MFI per cell
write.csv(mean_intensities, "MFI_CD31_xenium_obj.csv", row.names = FALSE)

# ---- H&E Binary Mask ----
# Set working directory
setwd("C:/Users/s215357/Desktop/project testing/635_1050_02T_High/annotations/export")

# Load the TIFF image (assuming a single-channel grayscale image)
tiff_file <- "morphology_focus_0000_binary_mask.tif"  # Update with the actual TIFF file name
image_raster <- rast(tiff_file)  # Use rast() instead of raster()

# Get raster dimensions
width <- ncol(image_raster)  # Number of columns = X dimension
height <- nrow(image_raster) # Number of rows = Y dimension

# Convert the raster to a matrix
image_matrix <- as.matrix(image_raster, wide=TRUE)  # Explicitly convert to matrix

# Create a data table of all pixel coordinates and intensity values
pixel_data <- data.table(
  x = rep(1:width, each = height),   # X-coordinates
  y = rep(height:1, times = width),  # Y-coordinates (flip to match image origin)
  intensity = as.vector(image_raster)  # Flatten the raster into a vector
)

# # Save the extracted pixel intensities to a CSV (optional)
# write.csv(pixel_data, "extracted_pixel_intensities.csv", row.names = FALSE)

setwd("C:/Users/s215357/Desktop/project testing/635_1050_02T_High")
# Load the cell boundary data
boundary_data <- fread("cell_boundaries.csv")

# Define the conversion factor (adjust as needed)
micrometers_to_pixels <- 1 / 0.2125  

# Convert boundary data into spatial polygons
boundary_polygons <- boundary_data %>%
  mutate(
    vertex_x = vertex_x * micrometers_to_pixels,
    vertex_y = vertex_y * micrometers_to_pixels
  ) %>%
  group_by(cell_id) %>%
  summarise(geometry = st_sfc(st_polygon(list(cbind(vertex_x, vertex_y))))) %>%
  st_as_sf(crs = NA)

# Set keys for fast subsetting
setkey(pixel_data, x, y)

# Initialize progress bar
pb <- progress_bar$new(
  format = "  Processing cell :current/:total [:bar] :percent in :elapsed",
  total = nrow(boundary_polygons), clear = FALSE, width = 60
)

# Initialize data frame to store results
mean_intensities <- data.frame(cell_id = character(), mean_intensity = numeric(), stringsAsFactors = FALSE)

# Process each cell polygon
for (i in seq_len(nrow(boundary_polygons))) {
  pb$tick()  # Update progress bar
  
  # Extract cell ID and polygon
  cell_id <- boundary_polygons$cell_id[i]
  cell_polygon <- boundary_polygons[i, ]
  
  # Get the bounding box
  bbox <- st_bbox(cell_polygon)
  xmin <- round(as.numeric(bbox["xmin"]))
  ymin <- round(as.numeric(bbox["ymin"]))
  xmax <- round(as.numeric(bbox["xmax"]))
  ymax <- round(as.numeric(bbox["ymax"]))
  
  # Subset pixels in the bounding box
  pixel_subset <- pixel_data[x %between% c(xmin, xmax) & y %between% c(ymin, ymax)]
  
  # Convert to sf object for spatial intersection
  pixel_points <- st_as_sf(pixel_subset, coords = c("x", "y"), crs = NA)
  
  # Check which pixels are inside the polygon
  pixels_in_polygon <- st_intersects(pixel_points, cell_polygon, sparse = FALSE)
  
  # Keep only the pixels inside the cell
  pixel_subset <- pixel_subset[as.logical(pixels_in_polygon)]
  
  # Compute mean intensity
  mean_intensity <- mean(pixel_subset$intensity, na.rm = TRUE)
  
  # Store the result
  mean_intensities <- rbind(mean_intensities, data.frame(cell_id = cell_id, mean_intensity = mean_intensity))
}

# Save the final MFI per cell
write.csv(mean_intensities, "MFI_binarymask_635_1050_02T_High.csv", row.names = FALSE)
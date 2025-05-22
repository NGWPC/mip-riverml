# Define path to FEMA BLE data and output directory
fema <- '/Volumes/MyBook/mip_full_collection/'
outdir <- glue::glue("{fema}/riverML2")
dir.create(outdir, showWarnings = FALSE)

# Helper function to compute channel area (CA) from a transect
.findCA <- function(df, depth) {
  Y <- NULL
  t <- filter(df, Y <= depth)
  
  suppressWarnings({
    x <- pmax(
      0, 
      DescTools::AUC(x = t$x, y = rep(depth, nrow(t)), absolutearea = FALSE) - 
        DescTools::AUC(x = t$x, y = t$Y, absolutearea = FALSE)
    )
  })
  
  ifelse(is.na(x), 0, x)
}

# Get all BLE subdirectories excluding previous riverML runs
ble <- list.dirs(fema, recursive = FALSE)
ble <- ble[!grepl('riverML', ble)]

# Reference fabric path
ref_path <- "/Users/mikejohnson/hydrofabric/v3.0/reference_fabric.gpkg"

# Loop through each BLE HUC directory
for (b in 464:length(ble)) {
  
  dir <- grep('submodels', list.dirs(ble[b], recursive = FALSE), value = TRUE)
  
  # Find all GPKG submodel files and extract metadata
  subs <- list.files(dir, recursive = TRUE, pattern = ".gpkg$", full.names = TRUE) |> 
    as.data.frame() |> 
    setNames("file") |> 
    mutate(
      reach = gsub('.*/', '', file),
      reach = gsub('.gpkg', '', reach),
      name = gsub('/submodels', "", gsub(fema, "", dir))
    ) 
  
  outdir_here <- glue::glue("{outdir}/{subs$name[1]}.gpkg")
  
  if (file.exists(outdir_here)) {
    message("\tAlready processed ", basename(subs$name[1]), " - skipping")
  } else {
    message("Processing ", basename(ble[b]), " (", b ," in ", length(ble), ")")
    subs_data <- list()
    
    for (v in 1:nrow(subs)) {
      message("\tProcessing ", basename(subs$file[v]), " (", v ," in ", nrow(subs), ")")
      
      transects <- read_sf(subs$file[v], 'XS') |> 
        st_transform(5070) 
      
      ll <- list()
      
      for (j in 1:nrow(transects)) {
        # Clean and parse station-elevation point strings
        cleaned <- gsub("\\[|\\]|\\(|\\)", "", transects$station_elevation_points[j])
        cleaned <- strsplit(cleaned, ", ")[[1]]
        
        df <- as.data.frame(matrix(as.numeric(cleaned), ncol = 2, byrow = TRUE))
        names(df) <- c("x", "Y")
        
        # Parse left and right bank station locations
        pins <- transects$bank_stations[j] %>% 
          gsub("\\[|\\]|\\'", "", .) |> 
          strsplit(",\\s*") |> 
          unlist() |> 
          as.numeric()
        
        # Subset and clean elevation data
        result <- dplyr::filter(df, dplyr::between(x, pins[1], pins[2])) 
        result$Y <- clean_elev(result$Y)
        
        if (nrow(result) <= 2 | diff(range(result$Y)) < .25) {
          warning("No channel in transect ", j, " for ", basename(subs$file[v]))
        } else {
          result$Ym <- max(result$Y) - min(result$Y)
          result$TW <- max(result$x) - min(result$x)
          result$flowpath_id <- subs$reach[v]
          result$river_station <- transects$river_station[j] 
          result$A <- .findCA(result, max(result$Y))
          result$r <- result$A / ((result$Ym * result$TW) - result$A)
          result$domain <- subs$name[v]
          
          ll[[j]] <- dplyr::distinct(dplyr::select(result, -x, -Y)) |> 
            slice(1) |> 
            left_join(
              select(transects[j,], 
                     river_station, river_reach_rs, 
                     source_river, source_reach, source_river_station),
              by = c('river_station')
            ) |> 
            st_as_sf()
        }
      }
      subs_data[[v]] <- bind_rows(ll)
    }
    
    huc_xs <- data.table::rbindlist(subs_data)
    
    if (nrow(huc_xs) == 0) {
      warning("No channels in submodel ", v, " for ", subs$reach[v])
    } else {
      huc_xs <- st_as_sf(huc_xs)
      
      # Compute representative XS features per flowpath
      representive_features <- huc_xs |> 
        tidyr::drop_na(flowpath_id) |> 
        dplyr::group_by(flowpath_id) |>
        arrange(river_station) |> 
        dplyr::summarise(
          r = mean(r[is.finite(r)]),
          TW = mean(TW),
          Y = mean(Ym),
          geom = geom[ceiling(n()/2)]
        )
      
      # Write output layers
      write_sf(huc_xs, outdir_here, layer = "XS")
      write_sf(representive_features, outdir_here, layer = "representative_xs")
      
      read_sf(
        ref_path, "flowpaths",
        wkt_filter = st_as_text(st_as_sfc(st_bbox(st_union(huc_xs))))
      ) |> 
        write_sf(outdir_here, layer = "reference_fabric")
    }
  }
}

# Load and export final dataset

xs <- map(list.files(outdir, full.names = TRUE), 
          ~read_sf(.x, 'representative_xs'))

y <- nhdplusTools::get_vaa(c('ftype', 'streamorde')) |> 
  mutate(comid = as.character(comid))

out_xs <- bind_rows(xs) |> 
  st_drop_geometry() |> 
  left_join(y, by = c('flowpath_id' = 'comid'))

arrow::write_parquet(out_xs, "riverML_ripple_beta.parquet")


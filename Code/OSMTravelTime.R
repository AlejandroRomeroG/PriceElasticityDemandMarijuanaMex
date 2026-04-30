###############################################################################
#  OSMTravelTime.R
#  Helpers for OpenStreetMap/OSRM driving-time instruments
###############################################################################

required_osm_pkgs <- c("jsonlite", "dplyr", "sf", "stringr", "tibble")
missing_osm_pkgs <- required_osm_pkgs[
  !vapply(required_osm_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_osm_pkgs) > 0) {
  stop("Missing required package(s) for OSM driving-time instruments: ",
       paste(missing_osm_pkgs, collapse = ", "), call. = FALSE)
}

osrm_server <- function() {
  Sys.getenv("OSRM_SERVER", "https://router.project-osrm.org")
}

osrm_profile <- function() {
  Sys.getenv("OSRM_PROFILE", "driving")
}

osrm_sleep_seconds <- function() {
  as.numeric(Sys.getenv("OSRM_SLEEP_SECONDS", "0.2"))
}

osrm_allow_network <- function() {
  Sys.getenv("OSRM_ALLOW_NETWORK", Sys.getenv("ALLOW_PUBLIC_OSRM", "0")) == "1"
}

osrm_src_chunk <- function() {
  value <- suppressWarnings(as.integer(Sys.getenv("OSRM_SRC_CHUNK", "20")))
  if (is.na(value) || value < 1) 20L else value
}

osrm_dst_chunk <- function() {
  value <- suppressWarnings(as.integer(Sys.getenv("OSRM_DST_CHUNK", "20")))
  if (is.na(value) || value < 1) 20L else value
}

municipality_centroids_osm <- function(shapefile_path) {
  shp <- sf::st_read(shapefile_path, quiet = TRUE) |>
    sf::st_transform(6372)

  cent_proj <- sf::st_centroid(sf::st_geometry(shp))
  xy <- sf::st_coordinates(cent_proj)

  cent_ll <- sf::st_as_sf(
    data.frame(cvegeo = shp$CVEGEO),
    geometry = cent_proj,
    crs = 6372
  ) |>
    sf::st_transform(4326)
  ll <- sf::st_coordinates(cent_ll)

  tibble::tibble(
    cvegeo = shp$CVEGEO,
    x_m = xy[, 1],
    y_m = xy[, 2],
    lon = ll[, 1],
    lat = ll[, 2]
  )
}

hub_points_from_states <- function(state_shapefile_path, hub_state_codes) {
  shp <- sf::st_read(state_shapefile_path, quiet = TRUE) |>
    sf::st_transform(6372) |>
    dplyr::filter(CVE_ENT %in% hub_state_codes)

  cent_proj <- sf::st_centroid(sf::st_geometry(shp))
  xy <- sf::st_coordinates(cent_proj)

  cent_ll <- sf::st_as_sf(
    data.frame(state_code = shp$CVE_ENT),
    geometry = cent_proj,
    crs = 6372
  ) |>
    sf::st_transform(4326)
  ll <- sf::st_coordinates(cent_ll)

  tibble::tibble(
    state_code = shp$CVE_ENT,
    x_m = xy[, 1],
    y_m = xy[, 2],
    lon = ll[, 1],
    lat = ll[, 2]
  ) |>
    dplyr::mutate(
      hub_id = state_code
    ) |>
    dplyr::select(hub_id, state_code, x_m, y_m, lon, lat)
}

chunk_indices <- function(n, chunk_size) {
  split(seq_len(n), ceiling(seq_len(n) / chunk_size))
}

osrm_table_url <- function(src, dst, server = osrm_server(),
                           profile = osrm_profile()) {
  coords <- paste(
    c(paste(src$lon, src$lat, sep = ","),
      paste(dst$lon, dst$lat, sep = ",")),
    collapse = ";"
  )
  sources <- paste(seq_len(nrow(src)) - 1L, collapse = ";")
  destinations <- paste(nrow(src) + seq_len(nrow(dst)) - 1L, collapse = ";")

  paste0(
    sub("/$", "", server), "/table/v1/", profile, "/", coords,
    "?sources=", sources,
    "&destinations=", destinations,
    "&annotations=duration"
  )
}

duration_block_minutes <- function(src, dst, server = osrm_server(),
                                   profile = osrm_profile(),
                                   max_tries = 3,
                                   min_split_size = 10) {
  url <- osrm_table_url(src, dst, server, profile)
  out <- NULL

  for (attempt in seq_len(max_tries)) {
    out <- tryCatch(
      jsonlite::fromJSON(url, simplifyVector = FALSE),
      error = function(e) e
    )

    if (!inherits(out, "error") && identical(out$code, "Ok")) {
      mat <- do.call(
        rbind,
        lapply(out$durations, function(row) {
          vapply(row, function(x) {
            if (is.null(x)) NA_real_ else as.numeric(x)
          }, numeric(1))
        })
      )
      return(mat / 60)
    }

    if (attempt < max_tries) {
      Sys.sleep(1 + attempt)
    }
  }

  msg <- if (inherits(out, "error")) {
    conditionMessage(out)
  } else if (!is.null(out$message)) {
    out$message
  } else {
    out$code
  }

  if ((nrow(src) > 1 || nrow(dst) > 1) &&
      (nrow(src) + nrow(dst)) > min_split_size) {
    message(
      "OSRM block failed for ", nrow(src), " x ", nrow(dst),
      "; retrying as smaller blocks."
    )

    if (nrow(src) >= nrow(dst) && nrow(src) > 1) {
      mid <- ceiling(nrow(src) / 2)
      cuts <- list(seq_len(mid), seq.int(mid + 1, nrow(src)))
      return(do.call(rbind, lapply(cuts, function(ii) {
        duration_block_minutes(src[ii, , drop = FALSE], dst,
                               server = server, profile = profile,
                               max_tries = max_tries,
                               min_split_size = min_split_size)
      })))
    }

    if (nrow(dst) > 1) {
      mid <- ceiling(nrow(dst) / 2)
      cuts <- list(seq_len(mid), seq.int(mid + 1, nrow(dst)))
      return(do.call(cbind, lapply(cuts, function(jj) {
        duration_block_minutes(src, dst[jj, , drop = FALSE],
                               server = server, profile = profile,
                               max_tries = max_tries,
                               min_split_size = min_split_size)
      })))
    }
  }

  stop("OSRM table request failed: ", msg, call. = FALSE)
}

osrm_duration_matrix <- function(src, dst, cache_path,
                                 src_id = "cvegeo", dst_id = "cvegeo",
                                 src_chunk = osrm_src_chunk(),
                                 dst_chunk = osrm_dst_chunk(),
                                 server = osrm_server(),
                                 profile = osrm_profile(),
                                 allow_network = osrm_allow_network()) {
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)

  src <- src |>
    dplyr::filter(!is.na(.data[[src_id]]), !is.na(lon), !is.na(lat)) |>
    dplyr::distinct(.data[[src_id]], .keep_all = TRUE)
  dst <- dst |>
    dplyr::filter(!is.na(.data[[dst_id]]), !is.na(lon), !is.na(lat)) |>
    dplyr::distinct(.data[[dst_id]], .keep_all = TRUE)

  src_names <- as.character(src[[src_id]])
  dst_names <- as.character(dst[[dst_id]])

  if (file.exists(cache_path)) {
    cached_mat <- readRDS(cache_path)
    if (identical(rownames(cached_mat), src_names) &&
        identical(colnames(cached_mat), dst_names) &&
        all(!is.na(cached_mat))) {
      message("Using cached OSRM matrix: ", cache_path)
      return(cached_mat)
    }

    mat <- matrix(NA_real_, nrow = nrow(src), ncol = nrow(dst),
                  dimnames = list(src_names, dst_names))
    overlap_src <- intersect(rownames(cached_mat), src_names)
    overlap_dst <- intersect(colnames(cached_mat), dst_names)
    if (length(overlap_src) > 0 && length(overlap_dst) > 0) {
      mat[overlap_src, overlap_dst] <- cached_mat[overlap_src, overlap_dst]
    }
  } else {
    mat <- matrix(NA_real_, nrow = nrow(src), ncol = nrow(dst),
                  dimnames = list(src_names, dst_names))
  }

  if (!isTRUE(allow_network)) {
    stop(
      "Cached OSRM matrix is missing or incomplete: ", cache_path,
      ". Provide a local cache, or explicitly set OSRM_ALLOW_NETWORK=1 ",
      "to query the configured OSRM server.",
      call. = FALSE
    )
  }

  message("Building OSRM matrix via ", server,
          " (", nrow(src), " x ", nrow(dst), " pairs).")

  src_chunks <- chunk_indices(nrow(src), src_chunk)
  dst_chunks <- chunk_indices(nrow(dst), dst_chunk)

  for (ii in src_chunks) {
    for (jj in dst_chunks) {
      if (all(!is.na(mat[ii, jj, drop = FALSE]))) next

      block <- duration_block_minutes(
        src[ii, , drop = FALSE],
        dst[jj, , drop = FALSE],
        server = server,
        profile = profile
      )

      mat[ii, jj] <- block
      saveRDS(mat, cache_path)
      Sys.sleep(osrm_sleep_seconds())
    }
  }

  if (anyNA(mat)) {
    stop("OSRM matrix contains missing driving times: ", cache_path, call. = FALSE)
  }

  mat
}

nearest_hub_travel_time <- function(purchase_centroids, hub_points,
                                    cache_path) {
  mat <- osrm_duration_matrix(
    src = hub_points,
    dst = purchase_centroids,
    cache_path = cache_path,
    src_id = "hub_id",
    dst_id = "cvegeo"
  )

  nearest_idx <- max.col(-t(mat), ties.method = "first")

  tibble::tibble(
    cvegeo = colnames(mat),
    travel_time_nearest_hub_min = apply(mat, 2, min, na.rm = TRUE),
    nearest_hub_state = rownames(mat)[nearest_idx]
  ) |>
    dplyr::mutate(
      ln_travel_time_hub = log(pmax(travel_time_nearest_hub_min, 1))
    )
}

gravity_from_travel_time <- function(time_mat, seizure_muni_month,
                                     purchase_months,
                                     kg_col,
                                     out_col,
                                     lag_months = 0,
                                     missing_lag = c("na", "zero"),
                                     month_col = NULL) {
  missing_lag <- match.arg(missing_lag)
  if (is.null(month_col)) {
    month_col <- if ("month_id" %in% names(seizure_muni_month)) "month_id" else "mes"
  }
  if (!month_col %in% names(seizure_muni_month)) {
    stop("Month column not found in seizure_muni_month: ", month_col, call. = FALSE)
  }

  weight_mat <- 1 / (1 + time_mat)

  n_s <- ncol(weight_mat)
  out <- vector("list", length(purchase_months))
  names(out) <- as.character(purchase_months)

  for (t in purchase_months) {
    seizure_month <- t - lag_months

    if (seizure_month < min(seizure_muni_month[[month_col]], na.rm = TRUE)) {
      lag_value <- if (missing_lag == "zero") 0 else NA_real_
      out[[as.character(t)]] <- tibble::tibble(
        cvegeo = rownames(weight_mat),
        !!month_col := t,
        value = lag_value
      )
      next
    }

    sz_t <- seizure_muni_month |>
      dplyr::filter(.data[[month_col]] == seizure_month, .data[[kg_col]] > 0) |>
      dplyr::filter(cvegeo_seizure %in% colnames(weight_mat))

    if (nrow(sz_t) == 0) {
      out[[as.character(t)]] <- tibble::tibble(
        cvegeo = rownames(weight_mat),
        !!month_col := t,
        value = 0
      )
      next
    }

    kg_vec <- rep(0, n_s)
    names(kg_vec) <- colnames(weight_mat)
    kg_vec[sz_t$cvegeo_seizure] <- sz_t[[kg_col]]

    out[[as.character(t)]] <- tibble::tibble(
      cvegeo = rownames(weight_mat),
      !!month_col := t,
      value = as.numeric(weight_mat %*% kg_vec)
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::rename(!!out_col := value)
}

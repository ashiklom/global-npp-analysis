library(drake)
library(ggplot2)
import::from(magrittr, "%>%")

get_glcf_data <- function() {
  all_files <- list.files("glcf_data", pattern = "\\.tif$", full.names = TRUE)
  raster_list <- Map(raster::raster, all_files)
  time_series <- raster::brick(raster_list)
  ts_all_values <- raster::getValues(time_series)
  ts_rownums <- seq_len(nrow(ts_all_values))
  ts_all_coords <- raster::xyFromCell(time_series, ts_rownums)
  colnames(ts_all_coords) <- c("lon", "lat")
  list(values = ts_all_values, coords = ts_all_coords)
}

geom_hl <- function(start, end, ...) {
  annotate(geom = "rect", xmin = start, xmax = end, ymin = -Inf, ymax = Inf, ...)
}

pca_var_plot <- function(pca_result) {
  pca_sd <- pca_result$sdev
  df <- tibble::tibble(
    x = c(0, seq_along(pca_sd)),
    y = c(0, cumsum(pca_sd/sum(pca_sd)))
  )
  ggplot(df) +
    aes(x = x, y = y) +
    geom_step() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw() +
    labs(
      x = "PCA Component",
      y = "Frac. variance explained"
    )
}

pkgconfig::set_config("drake::strings_in_dots" = "literals")

plan <- drake_plan(
  glcf_all = get_glcf_data(),
  cumulative_npp = rowSums(glcf_all[["values"]]),
  threshold = 0,
  above_threshold = which(cumulative_npp >= threshold),
  glcf_values = glcf_all[["values"]][above_threshold, ],
  glcf_coords = glcf_all[["coords"]][above_threshold, ],
  glcf_2000 = cbind(
    data.table::as.data.table(glcf_coords),
    glcf_values[, 20]
  ),
  glcf_2000_plot = ggplot(glcf_2000) +
    aes(x = lon, y = lat, fill = V2) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "green4") +
    ## coord_map("mollweide") +
    guides(fill = guide_legend(title = "NPP")) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      legend.position = "bottom"
    ),
  glcf_2000_plot_usa = glcf_2000_plot +
    ## coord_map("albers", c(-115, -85)) +
    xlim(-130, -70) + ylim(25, 50),
  nino3 = readr::read_table("nino3_hadlsst1.dat", col_names = c("year", 1:12)) %>%
    dplyr::filter(year >= 1980, year <= 2000) %>%
    tidyr::gather("month", "value", -year) %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(nino3 = mean(value)),
  nino3_plot = ggplot(nino3) +
    aes(x = year, y = nino3) +
    geom_line() +
    xlab("Year") + ylab("NINO3 index"),
  npp_diff = (glcf_values[, -1] - glcf_values[, 1]),
  pca_result = princomp(npp_diff),
  pca_sd_plot = pca_var_plot(pca_result),
  loadings = cbind(
    data.table::data.table(year = 1982:2000),
    matrix(c(pca_result$loadings), 19, 19)
  ),
  loadings_tidy = data.table::melt(loadings, id.vars = "year") %>%
    dplyr::left_join(nino3) %>%
    dplyr::mutate(variable = as.integer(gsub("V", "", variable))) %>%
    dplyr::rename(component = variable),
  loadings_ts = loadings_tidy %>%
    dplyr::filter(component <= 10) %>%
    ggplot() +
    aes(x = year, y = value) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ component, scales = "free_y"),
  scores = cbind(
    data.table::as.data.table(glcf_coords),
    pca_result$scores
  ),
  enso_corr = loadings_tidy %>%
    dplyr::group_by(component) %>%
    dplyr::summarize(R = cor(value, nino3)),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_format = "pdf_document",
    output_file = file_out("report.pdf")
  )
)
plan_config <- drake_config(plan)
make(plan)

if (FALSE) {

  readd(loadings_tidy) %>%
    dplyr::filter(variable %in% paste0("V", 1:6)) %>%
    tidyr::gather("var2", "value", -year, -variable) %>%
    ggplot() +
    aes(x = year, y = value) +
    geom_line() +
    ## geom_hline(yintercept = 0, linetype = "dashed") +
    ## geom_hl(1997, 1998, fill = "red", alpha = 0.2) +
    facet_grid(var2 ~ variable, scales = "free_y")

  library(raster)
  library(data.table)

  # Subtract first year NPP to get difference
  npp_diff <- ts_dt - ts_dt[, 1]

  # Use linear algebra to manually calculate regression slopes and
  # intercepts for each cell
  npp_diff_t_all <- t(npp_diff)
  npp_diff_sums <- colSums(npp_diff_t_all)
  hist(npp_diff_sums)
  npp_diff_t <- npp_diff_t_all[, colSums(npp_diff_t_all) > 5]

  y <- seq_len(ncol(npp_diff))
  N <- nrow(npp_diff)
  xy <- npp_diff_t * y
  sxy <- sum(xy)
  sx <- colSums(npp_diff_t)
  sy <- sum(y)
  sx2 <- colSums(npp_diff_t ^ 2)
  s2x <- sx ^ 2
  m <- (N * sxy - sx * sy) / (N * sx2 - s2x)
  hist(m)


  pca_out <- princomp(npp_diff)
  ## pca_out <- princomp(ts_dt, cor = TRUE)

  scores <- pca_out$scores
  scores_dt <- as.data.table(ts_coords)
  scores_dt <- cbind(scores_dt, scores)

  ld <- matrix(c(pca_out$loadings), 20, 20)
  colnames(ld) <- paste0("comp_", 1:20)
  loadings <- cbind(
    data.table(year = 1981:2000),
    ld
  )

  pca_sd <- pca_out$sdev
  plot(
    0:20,
    c(0, cumsum(pca_sd/sum(pca_sd))),
    type = "l",
    xlab = "PCA Component",
    ylab = "Frac. variance explained"
  )
  abline(h = 1, lty = "dashed")

  loadings_tidy <- melt(loadings[, 1:8], id.vars = "year")

  ggplot(loadings_tidy) +
    aes(x = year, y = value) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~variable) +
    theme(axis.title = element_blank())

  ggplot(scores_dt) +
    aes(x = lon, y = lat, fill = Comp.1) +
    geom_tile() +
    scale_fill_gradient2()

  ggplot(scores_dt) +
    aes(x = lon, y = lat, fill = Comp.2) +
    geom_tile() +
    scale_fill_gradient2()
}

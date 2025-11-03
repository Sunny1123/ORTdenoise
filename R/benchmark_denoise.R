#' Benchmark denoising on first channel with Gaussian noise
#'
#' @param img input image (matrix/array, or integer specifying size)
#' @param sigma noise standard deviation (in 0–1 scale)
#' @param iterations number of noisy realizations (default 500)
#' @param minlossimprove parameter for smoothimage
#' @param smband parameter for smoothimage
#'
#' @return list with mean_rmse, sd_rmse, all_rmse
#' @export
benchmark_denoise <- function(img, sigma = 0.05, iterations = 500,
                              minlossimprove = 0.0001, smband = 10) {
  # Case 1: integer -> generate synthetic triangular image
  if (length(img) == 1 && is.numeric(img)) {
    size <- as.integer(img)
    mat <- matrix(0, nrow = size, ncol = size)

    # Fill below diagonal
    for (i in 1:size) {
      for (j in 1:size) {
        if (i - j > 0) {
          mat[i, j] <- 1
        }
      }
    }
    # Fill lower-right triangle
    for (i in 1:size) {
      for (j in 1:size) {
        if (i + j > size) {
          mat[i, j] <- 1
        }
      }
    }
    img <- mat
  }

  # Case 2: if img is 3D array, use first channel
  if (length(dim(img)) == 3) {
    img <- img[,,1]
  }
  # print(dim(img))
  # image(img)
  # Ensure numeric
  # img <- as.numeric(img)
  # img <- matrix(img, nrow = nrow(img), ncol = ncol(img))

  # Rescale if necessary (0–255 -> 0–1)
  if (max(img, na.rm = TRUE) > 1) {
    img <- img / 255
  }

  # Ground truth
  gt <- as.vector(img)

  # Store RMSEs
  rmses <- numeric(iterations)

  for (i in seq_len(iterations)) {
    # Add Gaussian noise
    noisy <- img + rnorm(length(gt), mean = 0, sd = sigma)
    noisy <- matrix(noisy, nrow = nrow(img), ncol = ncol(img))

    # Denoise
    res <- smoothimage(noisy,
                       minlossimprove = minlossimprove,
                       smband = smband)$image

    # RMSE
    rmses[i] <- sqrt(mean((as.vector(res) - gt)^2))
    cat("Processing validation image", i, "/", iterations, "\n")
    flush.console()
  }

  invisible(list(
    mean_rmse = mean(rmses),
    sd_rmse   = sd(rmses),
    all_rmse  = rmses
  ))
}

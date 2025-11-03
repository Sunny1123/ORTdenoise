# run_sidd_srgb<- function() {
#
#     # ---------------------- Denoiser ----------------------
#     denoise_rgb <- function(img, minlossimprove=1e-4, smband=5) {
#       out <- array(0, dim=dim(img))
#       for (c in 1:dim(img)[3]) {
#         out[,,c] <- smoothimage(img[,,c], minlossimprove, smband)
#       }
#       return(out)
#     }
#
#     # ---------------------- Base64 helpers ----------------------
#     array_to_base64string <- function(x) {
#       raw <- writeBin(as.vector(x), raw(), size=4) # float32 assumption
#       base64encode(raw)
#     }
#
#     # ---------------------- Paths and download ----------------------
#     url <- "https://competitions.codalab.org/my/datasets/download/0d8a1e68-155d-4301-a8cd-9b829030d719"
#     data_dir <- "data"
#     input_file <- file.path(data_dir, "BenchmarkNoisyBlocksSrgb.mat")
#
#     if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
#
#     if (!file.exists(input_file)) {
#       cat("Downloading input file...\n")
#       tryCatch({
#         download.file(url, input_file, mode="wb")
#         cat("Downloaded successfully to", input_file, "\n")
#       }, error=function(e) {
#         cat("Automatic download failed. Please download manually from:\n", url, "\n")
#         stop("Missing required input file.")
#       })
#     } else cat("Found existing file:", input_file, "\n")
#
#     # ---------------------- Load .mat file ----------------------
#     key <- "BenchmarkNoisyBlocksSrgb"
#     inputs <- readMat(input_file)[[key]] # N x slice x H x W x C
#     dims <- dim(inputs)
#     cat("inputs dims:", paste(dims, collapse=" x "), "\n")
#     N <- dims[1]
#     slices <- dims[2]
#
#     # ---------------------- Run denoising ----------------------
#     output_blocks <- list()
#     n_blocks <- 0
#
#     for (i in 1:N) {          # loop over images
#       for (j in 1:slices) {   # loop over slices
#         in_block <- inputs[i,j,,,,drop=TRUE]  # H x W x C
#         if (length(dim(in_block))==2) in_block <- array(in_block, dim=c(dim(in_block),1))
#         out_block <- denoise_rgb(in_block)
#         stopifnot(all(dim(in_block)==dim(out_block)))
#
#         n_blocks <- n_blocks + 1
#         output_blocks[[n_blocks]] <- array_to_base64string(out_block)
#         cat(sprintf("Processed image %d / %d, slice %d / %d, dims: %s\n",
#                     i, N, j, slices, paste(dim(in_block), collapse=" x ")))
#       }
#
#     }
#
#     # ---------------------- Save submission CSV ----------------------
#     df <- data.frame(
#       ID = 0:(n_blocks-1),
#       BLOCK = unlist(output_blocks),
#       stringsAsFactors = FALSE
#     )
#
#     out_file <- file.path(data_dir, "SubmitSrgb.csv")
#     write.csv(df, out_file, row.names=FALSE)
#     cat("Saved submission to", out_file, "\n")
#     cat("Number of blocks =", n_blocks, "\n")
#     cat("Now submit SubmitSrgb.csv at: kaggle.com/competitions/sidd-benchmark-srgb-psnr\n")
#   }
run_sidd_srgb <- function() {

  library(R.matlab)
  library(base64enc)

  # ---------------------- Denoiser ----------------------
  denoise_rgb <- function(img, minlossimprove=1e-4, smband=5) {
    out <- array(0, dim=dim(img))
    for (c in 1:dim(img)[3]) {
      out[,,c] <- smoothimage(img[,,c], minlossimprove, smband)
    }
    return(out)
  }

  # ---------------------- Base64 helper ----------------------
  array_to_base64string <- function(x) {
    raw <- writeBin(as.vector(x), raw(), size=4) # float32 assumption
    base64encode(raw)
  }

  # ---------------------- Paths and download ----------------------
  url <- "https://competitions.codalab.org/my/datasets/download/0d8a1e68-155d-4301-a8cd-9b829030d719"
  data_dir <- "data"
  input_file <- file.path(data_dir, "BenchmarkNoisyBlocksSrgb.mat")

  if (!dir.exists(data_dir)) dir.create(data_dir, recursive=TRUE)

  if (!file.exists(input_file)) {
    cat("Downloading input file...\n")
    download.file(url, input_file, mode="wb")
    cat("Downloaded successfully to", input_file, "\n")
  } else cat("Found existing file:", input_file, "\n")

  # ---------------------- Load .mat file ----------------------
  key <- "BenchmarkNoisyBlocksSrgb"
  inputs <- readMat(input_file)[[key]] # N x slice x H x W x C
  dims <- dim(inputs)
  cat("inputs dims:", paste(dims, collapse=" x "), "\n")
  N <- dims[1]; slices <- dims[2]; H <- dims[3]; W <- dims[4]; C <- dims[5]

  # ---------------------- Prepare CSV ----------------------
  out_file <- file.path(data_dir, "SubmitSrgb.csv")
  cat("ID,BLOCK\n", file=out_file)  # write header

  # ---------------------- Process and write incrementally ----------------------
  block_id <- 0
  for (i in 1:N) {
    for (j in 1:slices) {
      in_block <- inputs[i,j,,,,drop=TRUE]  # H x W x C
      if (length(dim(in_block))==2) in_block <- array(in_block, dim=c(dim(in_block),1))

      out_block <- denoise_rgb(in_block)
      stopifnot(all(dim(in_block)==dim(out_block)))

      # Convert to base64 and write immediately
      b64 <- array_to_base64string(out_block)
      cat(sprintf("%d,%s\n", block_id, b64), file=out_file, append=TRUE)

      # Free memory
      rm(in_block, out_block, b64)
      gc()

      block_id <- block_id + 1
      cat(sprintf("Processed image %d/%d, slice %d/%d\n", i, N, j, slices))
    }
  }

  cat("Saved submission to", out_file, "\n")
  cat("Number of blocks =", block_id, "\n")
  cat("Submit the CSV at: kaggle.com/competitions/sidd-benchmark-srgb-psnr\n")
}

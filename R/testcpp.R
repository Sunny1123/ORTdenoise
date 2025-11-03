main =function(size=256,
               sigma = 0.1,
               smband=10,tolmult=0.001){

  tol=tolmult*sigma^2
  minlossimprove= tol#+sigma^2/20
  #img= peppers/max(peppers)
  #img = brain/max(brain)
  #Image type -1
  img=matrix(0,nrow=size,ncol=size)

  for(i in 1:size){
    for(j in 1:size){
      if(i-j>0 ){
        img[i,j]=1
      }
    }
  }
  for(i in 1:size){
    for(j in 1:size){
      if(i+j>size ){
        img[i,j]=1
      }
    }
  }
  dim(img)
  #img=readPNG("cat.png")[,,1]
  nimg=img+matrix(rnorm(nrow(img)*ncol(img),0,sigma),nrow=nrow(img),ncol=ncol(img)) # Noisy image

  t=Sys.time()
  res=smoothimage(nimg,minlossimprove,smband)
  print(Sys.time()-t)
  b= array(0,dim=dim(img))
  l=length(res$tree)
  print(l)
  for(i in 1:l)
  {
    loc = res$tree[[i]]$coord[,2:3]
    b[loc]=i
  }
  par(mfrow=c(1,3))
  image(nimg,useRaster=TRUE,axes=FALSE,col = gray(0:256 / 256))
  image(res$image,useRaster=TRUE,axes=FALSE,col = gray(0:256 / 256))
  image(b,useRaster=TRUE,axes=FALSE,col = gray(0:(2*l) / (2*l)))
  print(mean((nimg-img)^2)^0.5)
  print(mean((res$image-img)^2)^0.5)
  print(ssim(nimg,img))
  print(ssim(res$image,img))

}



main3d<- function(size=c(128,128,128),
                  sigma=0.1,
                  smband=4,
                  tolmult=0.01,
                  margin=0.05,
                  slice_idx=NULL,
                  invert=FALSE) {

  tol <- tolmult * sigma^2
  minlossimprove <- tol
  if(is.null(slice_idx)) slice_idx <- round(seq(size[3]*0.1, size[3]*0.9, length.out=5))

  nx <- size[1]; ny <- size[2]; nz <- size[3]
  xs <- seq(-1, 1, length.out=nx)
  ys <- seq(-1, 1, length.out=ny)
  zs <- seq(-1, 1, length.out=nz)

  img <- array(0, dim=c(nx,ny,nz))

  # Tetrahedron vertices
  # s <- 1 - margin
  # h <- sqrt(3)/2 * s
  # V1 <- c(-s/2, -h/3, -1 + margin)
  # V2 <- c(s/2, -h/3, -1 + margin)
  # V3 <- c(0, 2*h/3, -1 + margin)
  # V4 <- c(0, 0, 1 - margin)

  V1 <- c(-1 + margin, -1 + margin, -1 + margin)
  V2 <- c( 1 - margin, -1 + margin, -1 + margin)
  V3 <- c( 0,          1 - margin, -1 + margin)
  V4 <- c( 0,          0,           1 - margin)

  # -------------------------------------------------
  # Plane coefficients: ax + by + cz + d = 0
  # -------------------------------------------------
  # Base plane: V1,V2,V3
  a1 <- 0; b1 <- 0; c1 <- 1; d1 <- -(-1 + margin)

  # Side 1: V1,V2,V4
  a2 <- (V2[2]-V1[2])*(V4[3]-V1[3]) - (V2[3]-V1[3])*(V4[2]-V1[2])
  b2 <- (V2[3]-V1[3])*(V4[1]-V1[1]) - (V2[1]-V1[1])*(V4[3]-V1[3])
  c2 <- (V2[1]-V1[1])*(V4[2]-V1[2]) - (V2[2]-V1[2])*(V4[1]-V1[1])
  d2 <- - (a2*V1[1] + b2*V1[2] + c2*V1[3])

  # Side 2: V2,V3,V4
  a3 <- (V3[2]-V2[2])*(V4[3]-V2[3]) - (V3[3]-V2[3])*(V4[2]-V2[2])
  b3 <- (V3[3]-V2[3])*(V4[1]-V2[1]) - (V3[1]-V2[1])*(V4[3]-V2[3])
  c3 <- (V3[1]-V2[1])*(V4[2]-V2[2]) - (V3[2]-V2[2])*(V4[1]-V2[1])
  d3 <- - (a3*V2[1] + b3*V2[2] + c3*V2[3])

  # Side 3: V3,V1,V4
  a4 <- (V1[2]-V3[2])*(V4[3]-V3[3]) - (V1[3]-V3[3])*(V4[2]-V3[2])
  b4 <- (V1[3]-V3[3])*(V4[1]-V3[1]) - (V1[1]-V3[1])*(V4[3]-V3[3])
  c4 <- (V1[1]-V3[1])*(V4[2]-V3[2]) - (V1[2]-V3[2])*(V4[1]-V3[1])
  d4 <- - (a4*V3[1] + b4*V3[2] + c4*V3[3])

  planes <- list(
    list(a=a1,b=b1,c=c1,d=d1),
    list(a=a2,b=b2,c=c2,d=d2),
    list(a=a3,b=b3,c=c3,d=d3),
    list(a=a4,b=b4,c=c4,d=d4)
  )

  # Compute origin signs
  origin <- c(0,0,0)
  origin_signs <- sapply(planes, function(p) sign(p$a*origin[1] + p$b*origin[2] + p$c*origin[3] + p$d))

  # -------------------------------------------------
  # Fill tetrahedron voxels
  # -------------------------------------------------
  for(i in 1:nx){
    for(j in 1:ny){
      for(k in 1:nz){
        x <- xs[i]; y <- ys[j]; z <- zs[k]
        vals <- sapply(planes, function(p) sign(p$a*x + p$b*y + p$c*z + p$d))
        if(all(vals == origin_signs)) img[i,j,k] <- 1
      }
    }
  }

  # cat("Voxels inside tetrahedron:", sum(img > 0), "\n")

  # -------------------------------------------------
  # Noise + smoothing
  # -------------------------------------------------
  nimg <- img + array(rnorm(prod(dim(img)), 0, sigma), dim=dim(img))
  t=Sys.time()
  res   <- smoothimage(nimg, minlossimprove, smband)$image
  print(Sys.time()-t)
  # print(dim(res))
  # print(dim(nimg))
  # # -------------------------------------------------
  # # Slice visualization
  # # -------------------------------------------------
  # # -------------------------------------------------
  # # Slice visualization with consistent intensity scale
  # # -------------------------------------------------
  # if(is.null(slice_idx)) slice_idx <- round(c(nz*0.25, nz*0.5, nz*0.75))
  #
  # # Compute global min/max across all slices (noisy + smooth)
  # val_min <- min(c(nimg, res))
  # val_max <- max(c(nimg, res))
  #
  # plots_noisy <- list()
  # plots_smooth <- list()
  #
  # fill_scale <- if (invert) {
  #   scale_fill_gradient(low="black", high="white", limits=c(val_min, val_max), guide="none")
  # } else {
  #   scale_fill_gradient(low="white", high="black", limits=c(val_min, val_max), guide="none")
  # }
  #
  # for (i in seq_along(slice_idx)) {
  #   k <- slice_idx[i]
  #
  #   # Prepare data frames
  #   df_noisy <- expand.grid(x=1:nx, y=1:ny)
  #   df_noisy$val <- as.vector(nimg[,,k])
  #
  #   df_smooth <- expand.grid(x=1:nx, y=1:ny)
  #   df_smooth$val <- as.vector(res[,,k])
  #
  #   # Noisy plot
  #   plots_noisy[[i]] <- ggplot(df_noisy, aes(x, y, fill=val)) +
  #     geom_raster() +
  #     fill_scale +
  #     coord_fixed() +
  #     theme_void() +
  #     ggtitle(paste("Noisy Z =", k))
  #
  #   # Smooth plot
  #   plots_smooth[[i]] <- ggplot(df_smooth, aes(x, y, fill=val)) +
  #     geom_raster() +
  #     fill_scale +
  #     coord_fixed() +
  #     theme_void() +
  #     ggtitle(paste("Smooth Z =", k))
  # }
  #
  # # Arrange plots: noisy on top, smooth below
  # grid.arrange(grobs=c(plots_noisy, plots_smooth), nrow=2)
  # inside main3d, replace the plotting section with this:
  # --- paste/replace this plotting section inside main3d after `res` is computed ---

  # Ensure we have exactly 5 columns (if user supplied a different length, pick/adjust to 5)
  if (is.null(slice_idx) || length(slice_idx) < 5) {
    slice_idx <- round(seq(nz*0.1, nz*0.9, length.out = 5))
  } else if (length(slice_idx) > 5) {
    slice_idx <- slice_idx[round(seq(1, length(slice_idx), length.out = 5))]
  }


  # -------------------------------------------------
  # Compact slice visualization
  # -------------------------------------------------
  val_min <- min(c(nimg, res))
  val_max <- max(c(nimg, res))

  fill_scale <- if (invert) {
    ggplot2::scale_fill_gradient(low="black", high="white",
                                 limits=c(val_min, val_max), guide="none")
  } else {
    ggplot2::scale_fill_gradient(low="white", high="black",
                                 limits=c(val_min, val_max), guide="none")
  }

  plots_noisy <- vector("list", length(slice_idx))
  plots_smooth <- vector("list", length(slice_idx))

  shrink_theme <- ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "pt"))

  for (idx in seq_along(slice_idx)) {
    k <- slice_idx[idx]

    df_noisy <- expand.grid(x=1:nx, y=1:ny)
    df_noisy$val <- as.vector(nimg[,,k])

    df_smooth <- expand.grid(x=1:nx, y=1:ny)
    df_smooth$val <- as.vector(res[,,k])

    plots_noisy[[idx]] <- ggplot2::ggplot(df_noisy, ggplot2::aes(x, y, fill=val)) +
      ggplot2::geom_raster() + fill_scale +
      ggplot2::coord_fixed() + ggplot2::theme_void() + shrink_theme

    plots_smooth[[idx]] <- ggplot2::ggplot(df_smooth, ggplot2::aes(x, y, fill=val)) +
      ggplot2::geom_raster() + fill_scale +
      ggplot2::coord_fixed() + ggplot2::theme_void() + shrink_theme
  }

  # Column titles
  col_grobs <- lapply(slice_idx, function(k) {
    grid::textGrob(label=paste0("Z = ", k),
                   gp=grid::gpar(fontface="bold", fontsize=11))
  })

  # Row labels
  left_noisy  <- grid::textGrob("Noisy",  rot=90,
                                gp=grid::gpar(fontface="bold", fontsize=12))
  left_smooth <- grid::textGrob("Smooth", rot=90,
                                gp=grid::gpar(fontface="bold", fontsize=12))

  # Assemble grobs
  grobs <- c(
    list(grid::nullGrob()),      # 1 top-left empty
    col_grobs,                   # 2..6
    list(left_noisy, left_smooth), # 7,8
    plots_noisy,                 # 9..13
    plots_smooth                 # 14..18
  )

  layout_matrix <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(7, 9,10,11,12,13),
    c(8,14,15,16,17,18)
  )

  # Arrange with minimal gaps
  gridExtra::grid.arrange(
    grobs = grobs,
    layout_matrix = layout_matrix,
    heights = c(0.05, 0.475, 0.475),
    widths  = c(0.05, rep(0.95/5, 5)),
    padding = grid::unit(0, "pt")
  )

  # return as before
  invisible(list(noisy = nimg, smooth = res))
  # --- end of plotting section ---

}
main3diter <- function(size=c(128,128,128),
                   sigma=0.1,
                   smband=4,
                   tolmult=0.01,
                   margin=0.05,
                   slice_idx=NULL,
                   invert=FALSE,
                   n_iter=10,
                   plot_instance=1) {

  tol <- tolmult * sigma^2
  minlossimprove <- tol
  nx <- size[1]; ny <- size[2]; nz <- size[3]
  xs <- seq(-1, 1, length.out=nx)
  ys <- seq(-1, 1, length.out=ny)
  zs <- seq(-1, 1, length.out=nz)

  # -------------------------------------------------
  # Generate ground truth tetrahedron
  # -------------------------------------------------
  img <- array(0, dim=c(nx,ny,nz))

  V1 <- c(-1 + margin, -1 + margin, -1 + margin)
  V2 <- c( 1 - margin, -1 + margin, -1 + margin)
  V3 <- c( 0,          1 - margin, -1 + margin)
  V4 <- c( 0,          0,           1 - margin)

  # Plane equations
  a1 <- 0; b1 <- 0; c1 <- 1; d1 <- -(-1 + margin)
  a2 <- (V2[2]-V1[2])*(V4[3]-V1[3]) - (V2[3]-V1[3])*(V4[2]-V1[2])
  b2 <- (V2[3]-V1[3])*(V4[1]-V1[1]) - (V2[1]-V1[1])*(V4[3]-V1[3])
  c2 <- (V2[1]-V1[1])*(V4[2]-V1[2]) - (V2[2]-V1[2])*(V4[1]-V1[1])
  d2 <- - (a2*V1[1] + b2*V1[2] + c2*V1[3])

  a3 <- (V3[2]-V2[2])*(V4[3]-V2[3]) - (V3[3]-V2[3])*(V4[2]-V2[2])
  b3 <- (V3[3]-V2[3])*(V4[1]-V2[1]) - (V3[1]-V2[1])*(V4[3]-V2[3])
  c3 <- (V3[1]-V2[1])*(V4[2]-V2[2]) - (V3[2]-V2[2])*(V4[1]-V2[1])
  d3 <- - (a3*V2[1] + b3*V2[2] + c3*V2[3])

  a4 <- (V1[2]-V3[2])*(V4[3]-V3[3]) - (V1[3]-V3[3])*(V4[2]-V3[2])
  b4 <- (V1[3]-V3[3])*(V4[1]-V3[1]) - (V1[1]-V3[1])*(V4[3]-V3[3])
  c4 <- (V1[1]-V3[1])*(V4[2]-V3[2]) - (V1[2]-V3[2])*(V4[1]-V3[1])
  d4 <- - (a4*V3[1] + b4*V3[2] + c4*V3[3])

  planes <- list(
    list(a=a1,b=b1,c=c1,d=d1),
    list(a=a2,b=b2,c=c2,d=d2),
    list(a=a3,b=b3,c=c3,d=d3),
    list(a=a4,b=b4,c=c4,d=d4)
  )

  origin <- c(0,0,0)
  origin_signs <- sapply(planes, function(p) sign(p$a*origin[1] + p$b*origin[2] + p$c*origin[3] + p$d))

  for(i in 1:nx){
    for(j in 1:ny){
      for(k in 1:nz){
        x <- xs[i]; y <- ys[j]; z <- zs[k]
        vals <- sapply(planes, function(p) sign(p$a*x + p$b*y + p$c*z + p$d))
        if(all(vals == origin_signs)) img[i,j,k] <- 1
      }
    }
  }

  # -------------------------------------------------
  # Monte Carlo iterations for RMSE
  # -------------------------------------------------
  rmse_noisy  <- numeric(n_iter)
  rmse_smooth <- numeric(n_iter)
  saved_noisy <- NULL
  saved_smooth <- NULL

  for (it in 1:n_iter) {
    nimg <- img + array(rnorm(prod(dim(img)), 0, sigma), dim=dim(img))
    res  <- smoothimage(nimg, minlossimprove, smband)$image

    rmse_noisy[it]  <- sqrt(mean((nimg - img)^2))
    rmse_smooth[it] <- sqrt(mean((res  - img)^2))

    if (it == plot_instance) {
      saved_noisy  <- nimg
      saved_smooth <- res
    }
  }

  cat(sprintf("REMSE over %d runs:\n", n_iter))
  cat(sprintf("  Noisy REMSE  = %.5f\n", sqrt(mean(rmse_noisy^2))))
  cat(sprintf("  Noisy REMSE sd  = %.5f\n", sd(rmse_noisy)))
  cat(sprintf("  Smooth REMSE = %.5f\n", sqrt(mean(rmse_smooth^2))))
  cat(sprintf("  Smooth REMSE sd  = %.5f\n", sd(rmse_smooth)))

  # -------------------------------------------------
  # Plot single instance (same as before)
  # -------------------------------------------------
  if (is.null(slice_idx) || length(slice_idx) < 5) {
    slice_idx <- round(seq(nz*0.1, nz*0.9, length.out = 5))
  } else if (length(slice_idx) > 5) {
    slice_idx <- slice_idx[round(seq(1, length(slice_idx), length.out = 5))]
  }

  val_min <- min(c(saved_noisy, saved_smooth))
  val_max <- max(c(saved_noisy, saved_smooth))

  fill_scale <- if (invert) {
    ggplot2::scale_fill_gradient(low="black", high="white",
                                 limits=c(val_min, val_max), guide="none")
  } else {
    ggplot2::scale_fill_gradient(low="white", high="black",
                                 limits=c(val_min, val_max), guide="none")
  }

  plots_noisy <- vector("list", length(slice_idx))
  plots_smooth <- vector("list", length(slice_idx))
  shrink_theme <- ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "pt"))

  for (idx in seq_along(slice_idx)) {
    k <- slice_idx[idx]

    df_noisy <- expand.grid(x=1:nx, y=1:ny)
    df_noisy$val <- as.vector(saved_noisy[,,k])

    df_smooth <- expand.grid(x=1:nx, y=1:ny)
    df_smooth$val <- as.vector(saved_smooth[,,k])

    plots_noisy[[idx]] <- ggplot2::ggplot(df_noisy, ggplot2::aes(x, y, fill=val)) +
      ggplot2::geom_raster() + fill_scale +
      ggplot2::coord_fixed() + ggplot2::theme_void() + shrink_theme

    plots_smooth[[idx]] <- ggplot2::ggplot(df_smooth, ggplot2::aes(x, y, fill=val)) +
      ggplot2::geom_raster() + fill_scale +
      ggplot2::coord_fixed() + ggplot2::theme_void() + shrink_theme
  }

  col_grobs <- lapply(slice_idx, function(k) {
    grid::textGrob(label=paste0("Z = ", k),
                   gp=grid::gpar(fontface="bold", fontsize=11))
  })

  left_noisy  <- grid::textGrob("Noisy",  rot=90,
                                gp=grid::gpar(fontface="bold", fontsize=12))
  left_smooth <- grid::textGrob("Smooth", rot=90,
                                gp=grid::gpar(fontface="bold", fontsize=12))

  grobs <- c(
    list(grid::nullGrob()), col_grobs,
    list(left_noisy, left_smooth),
    plots_noisy, plots_smooth
  )

  layout_matrix <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(7, 9,10,11,12,13),
    c(8,14,15,16,17,18)
  )

  gridExtra::grid.arrange(
    grobs = grobs,
    layout_matrix = layout_matrix,
    heights = c(0.05, 0.475, 0.475),
    widths  = c(0.05, rep(0.95/5, 5)),
    padding = grid::unit(0, "pt")
  )

  invisible(list(RMSE=list(noisy=rmse_noisy, smooth=rmse_smooth),
                 example=list(noisy=saved_noisy, smooth=saved_smooth)))
}
main3diter2 <- function(image=128,
                       sigma=0.1,
                       smband=4,
                       tolmult=0.01,
                       margin=0.05,
                       slice_idx=NULL,
                       invert=FALSE,
                       n_iter=10,
                       plot_instance=1) {

  # -------------------------------------------------
  # Step 1: Load image or create tetrahedron
  # -------------------------------------------------
  if (is.numeric(image) && length(image) == 1) {
    # Case A: integer input -> synthetic tetrahedron
    size <- rep(image, 3)
    nx <- size[1]; ny <- size[2]; nz <- size[3]
    xs <- seq(-1, 1, length.out=nx)
    ys <- seq(-1, 1, length.out=ny)
    zs <- seq(-1, 1, length.out=nz)
    img <- array(0, dim=c(nx,ny,nz))

    # Define tetrahedron
    V1 <- c(-1 + margin, -1 + margin, -1 + margin)
    V2 <- c( 1 - margin, -1 + margin, -1 + margin)
    V3 <- c( 0,          1 - margin, -1 + margin)
    V4 <- c( 0,          0,           1 - margin)

    # Plane equations
    a1 <- 0; b1 <- 0; c1 <- 1; d1 <- -(-1 + margin)
    a2 <- (V2[2]-V1[2])*(V4[3]-V1[3]) - (V2[3]-V1[3])*(V4[2]-V1[2])
    b2 <- (V2[3]-V1[3])*(V4[1]-V1[1]) - (V2[1]-V1[1])*(V4[3]-V1[3])
    c2 <- (V2[1]-V1[1])*(V4[2]-V1[2]) - (V2[2]-V1[2])*(V4[1]-V1[1])
    d2 <- - (a2*V1[1] + b2*V1[2] + c2*V1[3])

    a3 <- (V3[2]-V2[2])*(V4[3]-V2[3]) - (V3[3]-V2[3])*(V4[2]-V2[2])
    b3 <- (V3[3]-V2[3])*(V4[1]-V2[1]) - (V3[1]-V2[1])*(V4[3]-V2[3])
    c3 <- (V3[1]-V2[1])*(V4[2]-V2[2]) - (V3[2]-V2[2])*(V4[1]-V2[1])
    d3 <- - (a3*V2[1] + b3*V2[2] + c3*V2[3])

    a4 <- (V1[2]-V3[2])*(V4[3]-V3[3]) - (V1[3]-V3[3])*(V4[2]-V3[2])
    b4 <- (V1[3]-V3[3])*(V4[1]-V3[1]) - (V1[1]-V3[1])*(V4[3]-V3[3])
    c4 <- (V1[1]-V3[1])*(V4[2]-V3[2]) - (V1[2]-V3[2])*(V4[1]-V3[1])
    d4 <- - (a4*V3[1] + b4*V3[2] + c4*V3[3])

    planes <- list(
      list(a=a1,b=b1,c=c1,d=d1),
      list(a=a2,b=b2,c=c2,d=d2),
      list(a=a3,b=b3,c=c3,d=d3),
      list(a=a4,b=b4,c=c4,d=d4)
    )

    origin <- c(0,0,0)
    origin_signs <- sapply(planes, function(p) sign(p$a*origin[1] + p$b*origin[2] + p$c*origin[3] + p$d))

    for(i in 1:nx){
      for(j in 1:ny){
        for(k in 1:nz){
          x <- xs[i]; y <- ys[j]; z <- zs[k]
          vals <- sapply(planes, function(p) sign(p$a*x + p$b*y + p$c*z + p$d))
          if(all(vals == origin_signs)) img[i,j,k] <- 1
        }
      }
    }

  } else if (is.character(image) && grepl("\\.nii(\\.gz)?$", image)) {
    # Case B: NIfTI input
    if (!requireNamespace("RNifti", quietly=TRUE)) {
      stop("Please install the 'RNifti' package to load .nii files.")
    }
    if (!requireNamespace("RNifti", quietly=TRUE)) {
      stop("Please install the 'RNifti' package to load .nii files.")
    }
    img <- RNifti::readNifti(image)
    if(max(img)>1) img = (img-min(img))/(max(img)-min(img))
    # Original dimensions
    orig_dim <- dim(img)

    # Target size
    target_dim <- c(128, 128, 128)

    # Build target grid
    xout <- seq(1, orig_dim[1], length.out = target_dim[1])
    yout <- seq(1, orig_dim[2], length.out = target_dim[2])
    zout <- seq(1, orig_dim[3], length.out = target_dim[3])

    # Allocate resampled array
    img_resampled <- array(0, dim = target_dim)

    # Step 1: interpolate along x for each (y,z)
    temp_x <- array(0, dim = c(target_dim[1], orig_dim[2], orig_dim[3]))
    for (y in seq_len(orig_dim[2])) {
      for (z in seq_len(orig_dim[3])) {
        temp_x[,y,z] <- approx(seq_len(orig_dim[1]), img[,y,z], xout, method="linear", rule=2)$y
      }
    }

    # Step 2: interpolate along y for each (x,z)
    temp_y <- array(0, dim = c(target_dim[1], target_dim[2], orig_dim[3]))
    for (x in seq_len(target_dim[1])) {
      for (z in seq_len(orig_dim[3])) {
        temp_y[x,,z] <- approx(seq_len(orig_dim[2]), temp_x[x,,z], yout, method="linear", rule=2)$y
      }
    }

    # Step 3: interpolate along z for each (x,y)
    for (x in seq_len(target_dim[1])) {
      for (y in seq_len(target_dim[2])) {
        img_resampled[x,y,] <- approx(seq_len(orig_dim[3]), temp_y[x,y,], zout, method="linear", rule=2)$y
      }
    }

    # Replace original image with resampled version
    img <- img_resampled
    nx <- dim(img)[1]; ny <- dim(img)[2]; nz <- dim(img)[3]


  } else {
    stop("`image` must be either an integer (cube size) or a .nii/.nii.gz filename")
  }

  # -------------------------------------------------
  # Monte Carlo iterations
  # -------------------------------------------------
  tol <- tolmult * sigma^2
  minlossimprove <- tol
  rmse_noisy  <- numeric(n_iter)
  rmse_smooth <- numeric(n_iter)
  saved_noisy <- NULL
  saved_smooth <- NULL

  for (it in 1:n_iter) {
    nimg <- img + array(rnorm(prod(dim(img)), 0, sigma), dim=dim(img))
    res  <- smoothimage(nimg, minlossimprove, smband)$image

    rmse_noisy[it]  <- sqrt(mean((nimg - img)^2))
    rmse_smooth[it] <- sqrt(mean((res  - img)^2))

    if (it == plot_instance) {
      saved_noisy  <- nimg
      saved_smooth <- res
    }
  }

  cat(sprintf("REMSE over %d runs:\n", n_iter))
  cat(sprintf("  Noisy REMSE   = %.5f ± %.5f\n", mean(rmse_noisy), sd(rmse_noisy)))
  cat(sprintf("  Smooth REMSE  = %.5f ± %.5f\n", mean(rmse_smooth), sd(rmse_smooth)))

  # -------------------------------------------------
  # Slice plotting (unchanged)
  # -------------------------------------------------
  if (is.null(slice_idx) || length(slice_idx) < 5) {
    slice_idx <- round(seq(nz*0.1, nz*0.9, length.out = 5))
  } else if (length(slice_idx) > 5) {
    slice_idx <- slice_idx[round(seq(1, length(slice_idx), length.out = 5))]
  }

  val_min <- min(c(saved_noisy, saved_smooth))
  val_max <- max(c(saved_noisy, saved_smooth))

  fill_scale <- if (invert) {
    ggplot2::scale_fill_gradient(low="black", high="white",
                                 limits=c(0, 1), guide="none")
  } else {
    ggplot2::scale_fill_gradient(low="white", high="black",
                                 limits=c(0, 1), guide="none")
  }

  plots_noisy <- vector("list", length(slice_idx))
  plots_smooth <- vector("list", length(slice_idx))
  shrink_theme <- ggplot2::theme(plot.margin = grid::unit(c(0,0,0,0), "pt"))

  for (idx in seq_along(slice_idx)) {
    k <- slice_idx[idx]

    df_noisy <- expand.grid(x=1:nx, y=1:ny)
    df_noisy$val <- as.vector(saved_noisy[,,k])
    df_noisy$val <- pmin(pmax(df_noisy$val, 0), 1)
    df_smooth <- expand.grid(x=1:nx, y=1:ny)
    df_smooth$val <- as.vector(saved_smooth[,,k])
    df_smooth$val <- pmin(pmax(df_smooth$val, 0), 1)

    plots_noisy[[idx]] <- ggplot2::ggplot(df_noisy, ggplot2::aes(x, y, fill=val)) +
      ggplot2::geom_raster() + fill_scale +
      ggplot2::coord_fixed() + ggplot2::theme_void() + shrink_theme

    plots_smooth[[idx]] <- ggplot2::ggplot(df_smooth, ggplot2::aes(x, y, fill=val)) +
      ggplot2::geom_raster() + fill_scale +
      ggplot2::coord_fixed() + ggplot2::theme_void() + shrink_theme
  }

  col_grobs <- lapply(slice_idx, function(k) {
    grid::textGrob(label=paste0("Z = ", k),
                   gp=grid::gpar(fontface="bold", fontsize=11))
  })

  left_noisy  <- grid::textGrob("Noisy",  rot=90,
                                gp=grid::gpar(fontface="bold", fontsize=12))
  left_smooth <- grid::textGrob("Smooth", rot=90,
                                gp=grid::gpar(fontface="bold", fontsize=12))

  grobs <- c(
    list(grid::nullGrob()), col_grobs,
    list(left_noisy, left_smooth),
    plots_noisy, plots_smooth
  )

  layout_matrix <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(7, 9,10,11,12,13),
    c(8,14,15,16,17,18)
  )

  gridExtra::grid.arrange(
    grobs = grobs,
    layout_matrix = layout_matrix,
    heights = c(0.05, 0.475, 0.475),
    widths  = c(0.05, rep(0.95/5, 5)),
    padding = grid::unit(0, "pt")
  )

  invisible(list(RMSE=list(noisy=rmse_noisy, smooth=rmse_smooth),
                 example=list(noisy=saved_noisy, smooth=saved_smooth)))
}

#' @export
smoothimage <- function(img, minlossimprove = 0.0001, smband = 4) {
  if(max(img>10)) scale_factor=255 else scale_factor=1
  img=img/scale_factor
  # Flatten image
  y <- as.vector(img)
  dims <- dim(img)
  p <- length(dims)

  # Build coordinate matrix with dummy column 1 using expand.grid
  dim_seq <- lapply(dims, seq_len)
  coords_df <- expand.grid(dim_seq)      # each row is a coordinate
  X <- cbind(1, as.matrix(coords_df))    # add dummy column

  # Create root node
  root_intensity <- y
  root_coord <- X

  # Build the tree using C++ function
  tree <- gettree_cpp(root_intensity, root_coord,
                      minlossimprove = minlossimprove,
                      minsize = 1)

  # Apply Treetomat smoothing
  res <- TreetomatPixelSIMD(tree, img,dim(img), smband)

  # Reshape back to original image dimensions
  res_array <- array(res, dim = dims)*scale_factor
  return(list(image=res_array,tree=tree))
}

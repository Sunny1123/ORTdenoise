# benchmark_denoise

## Benchmark denoising on first channel with Gaussian noise

**Usage:**
```

benchmark_denoise(
  img,
  sigma = 0.05,
  iterations = 500,
  minlossimprove = 1e-04,
  smband = 10
)

```

### Arguments

- **img**: input image (matrix/array, or integer specifying size)

- **sigma**: noise standard deviation (in 0–1 scale)

- **iterations**: number of noisy realizations (default 500)

- **minlossimprove**: parameter for smoothimage

- **smband**: parameter for smoothimage

### Value

list with mean_rmse, sd_rmse, all_rmse

**Description:** 
Benchmark denoising on first channel with Gaussian noise

---

# demo2d

Demonstrates denoising performance on simple 2d test image

**Usage:**
```

demo2d(size = 256, sigma = 0.1, smband = 10, tolmult = 0.001)

```

### Arguments

- **size**: resolution of the test image to be taken

- **sigma**: standard deviation of the noise to add

- **smband**: smoothing band

- **tolmult**: Multiplier to calculate cutoff for branching

**Description:** 
Title demo2d
Demonstrates denoising performance on simple 2d test image

### Examples
```r

demo2d(size=256,sigma=0.1,smband=10,tolmult=0.001)

```

---

# demo3d

Demo denoising on 3D images in 0-1 scale

**Usage:**
```

demo3d(
  size = c(128, 128, 128),
  sigma = 0.1,
  smband = 4,
  tolmult = 0.01,
  margin = 0.05,
  slice_idx = NULL,
  invert = FALSE
)

```

### Arguments

- **size**: dimensions of the image

- **sigma**: standard deviation of the noise to add

- **smband**: smoothing band

- **tolmult**: multiplier of tolerance, default is set to 0.01, determines the cutoff for branching

- **margin**: space to leave between the demo image and boundary

- **slice_idx**: slices to plot

- **invert**: Invert image color

### Value

list of noisy image and smoothed image

**Description:** 
Title demo3d
Demo denoising on 3D images in 0-1 scale

### Examples
```r

demo3d(size=c(128,128,128),sigma=0.1,smband=4,tolmult=0.01,margin=0.05,slice_idx=NULL,invert=FALSE)

```

---

# demo3diter

Multiple iterations of 3D demo

**Usage:**
```

demo3diter(
  image = 128,
  sigma = 0.1,
  smband = 4,
  tolmult = 0.01,
  margin = 0.05,
  slice_idx = NULL,
  invert = FALSE,
  n_iter = 10,
  plot_instance = 1
)

```

### Arguments

- **image**: Integer value represents dimensions of the images, location of a .nii file can be used as a list of predefined images

- **sigma**: standard deviation of the noise to add

- **smband**: smoothing band

- **tolmult**: multiplier of tolerance, default is set to 0.01, determines the cutoff for branching

- **margin**: space to leave between the demo image and boundary

- **slice_idx**: slices to plot

- **invert**: Invert image color

- **n_iter**: number of iterations

- **plot_instance**: how many instances to plot

### Value

rmse of noisy and smoothed image, sample instances of noisy and smoothed

**Description:** 
Title demo3diter
Multiple iterations of 3D demo

### Examples
```r

demo3diter(image=128,sigma=0.1,smband=4,tolmult=0.01,margin=0.05,slice_idx=NULL,invert=FALSE,n_iter=10,plot_instance=1)

```

---

# smoothimage

**Usage:**
```

smoothimage(img, minlossimprove = 1e-04, smband = 4, mode =1)

```

### Arguments

- **img**: Image to be denoised

- **minlossimprove**: cutoff for branching

- **smband**: smoothing band

- **mode**: mode 1 uses the c implementation which is faster but use approximations, 2 uses the R implemetation but does not work on 3d mages
### Value

smoothed image, and decision tree

**Description:** 
Title smoothimage

### Examples
```r

smoothimage(img,0,0001,4,1)

```

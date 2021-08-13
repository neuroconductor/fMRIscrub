test_that("pscrub works", {
  psx <- testthat::expect_warning(fMRIscrub:::pscrub_multi(
    Dat1,
    projection = "all"
  ))
  testthat::expect_warning(fMRIscrub:::plot.pscrub_multi(psx))

  psx <- testthat::expect_warning(fMRIscrub:::pscrub_multi(
    Dat2,
    projection = c("PCATF", "ICA_kurt"),
    kurt_quantile = .90,
    cutoff = 5,
    verbose = TRUE
  ))
  fMRIscrub:::plot.pscrub_multi(psx, title = "My Plot")

  psx <- testthat::expect_warning(scrub(Dat1))
  print(psx)

  psx <- pscrub(
    Dat2[,Dat2[1,]!=0], "PCA", 1, center=FALSE, PESEL=FALSE, kurt_quantile=.8,
    full_PCA = TRUE, cutoff=5, verbose=TRUE
  )
  plot(psx)

  psx <- testthat::expect_warning(pscrub(
    matrix(rnorm(10000), ncol=50)
  ))

  psx <- testthat::expect_warning(pscrub(
    matrix(rnorm(10000), nrow=100) + 100, nuisance=dct_bases(100, 2)
  ))

  psx <- testthat::expect_warning(pscrub(
    Dat2, projection="PCATF", nuisance=cbind(1, dct_bases(nrow(Dat2), 12)),
    comps_mean_dt=2, comps_var_dt=2, get_dirs=TRUE, get_outliers=FALSE
  ))
  plot(psx)
})

test_that("DVARS works", {
  dv <- DVARS(Dat1)
  plot(dv)

  dv <- DVARS(scale(Dat2[,Dat2[1,]!=0]), normalize=FALSE, cutoff_DPD=3, verbose=TRUE)
  plot(dv)
})

test_that("ciftiTools-related functions work", {
  if (is.null(ciftiTools.getOption("wb_path"))) {
    skip("Connectome Workbench is not available.")
  }

  cii_fname <- "../Data/rfMRI_REST1_LR_Atlas_MSMAll_6k.dtseries.nii"
  if (!file.exists(cii_fname)) {
    skip("Could not find the CIFTI file.")
  }

  psx <- pscrub(ciftiTools::read_xifti(cii_fname, brainstructures="left"))
  plot(psx)
  psx2 <- pscrub(t(as.matrix(
    ciftiTools::read_xifti(cii_fname, brainstructures="left")
  )))
  testthat::expect_equal(psx$measure, psx2$measure)

  dv <- DVARS(ciftiTools::read_xifti(cii_fname, brainstructures="right"))
  plot(dv)
  dv <- scrub_xifti(cii_fname, "DVARS", c("left", "right"))
  print(dv)
})

test_that("Miscellaneous functions work", {
  testthat::expect_warning(summary(pscrub(fsl_bptf(Dat2))))
})

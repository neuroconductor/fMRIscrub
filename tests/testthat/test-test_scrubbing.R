test_that("clever works", {
  clev <- testthat::expect_warning(fMRIscrub:::clever_multi(
    Dat1,
    projection = "all"
  ))
  testthat::expect_warning(fMRIscrub:::plot.clever_multi(clev))

  clev <- testthat::expect_warning(fMRIscrub:::clever_multi(
    Dat2,
    projection = c("PCATF", "ICA_kurt"),
    kurt_quantile = .90,
    cutoff = 5,
    verbose = TRUE
  ))
  fMRIscrub:::plot.clever_multi(clev, title = "My Plot")

  clev <- testthat::expect_warning(scrub(Dat1))
  print(clev)

  clev <- clever(
    Dat2[,Dat2[1,]!=0], "PCA", 1, center=FALSE, PESEL=FALSE, kurt_quantile=.8,
    full_PCA = TRUE, cutoff=5, verbose=TRUE
  )
  plot(clev)

  clev <- testthat::expect_warning(clever(
    matrix(rnorm(10000), ncol=50)
  ))

  clev <- testthat::expect_warning(clever(
    matrix(rnorm(10000), nrow=100) + 100, nuisance=dct_bases(100, 2)
  ))

  clev <- testthat::expect_warning(clever(
    Dat2, projection="PCATF", nuisance=cbind(1, dct_bases(nrow(Dat2), 12)),
    comps_mean_dt=2, comps_var_dt=2, get_dirs=TRUE, get_outliers=FALSE
  ))
  plot(clev)
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

  clev <- clever(ciftiTools::read_xifti(cii_fname, brainstructures="left"))
  plot(clev)
  clev2 <- clever(t(as.matrix(
    ciftiTools::read_xifti(cii_fname, brainstructures="left")
  )))
  testthat::expect_equal(clev$measure, clev2$measure)

  dv <- DVARS(ciftiTools::read_xifti(cii_fname, brainstructures="right"))
  plot(dv)
  dv <- scrub_xifti(cii_fname, "DVARS", c("left", "right"))
  print(dv)
})

test_that("Miscellaneous functions work", {
  testthat::expect_warning(summary(clever(fsl_bptf(Dat2))))
})

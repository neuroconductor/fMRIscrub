test_that("clever works with all methods and default settings", {
  data(Dat1)
  clev <- clever(
    Dat1, 
    projection = "all", 
    out_meas = "all"
  )
  myplot <- plot(clev, "all")
})

test_that("clever works with some custom parameter settings", {
  data(Dat1)
  clev <- clever(
    Dat1, 
    projection = "all", 
    out_meas = "all",
    DVARS = FALSE,
    PCATF_kwargs = list(lambda = 0.1, niter_max = 1100, TOL = 2e-8),
    kurt_quantile = .95,
    lev_cutoff = 5,
    rbd_cutoff = .9,
    lev_images = FALSE,
    verbose = TRUE
  )
  myplot <- plot(clev, "all", plot_title = "My Plot")
})

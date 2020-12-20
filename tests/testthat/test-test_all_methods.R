test_that("clever works with all methods and default settings", {
  data(Dat1)
  clev <- clever_multi(
    Dat1, 
    measures = "all", 
    projections = "all"
  )
  myplot <- plot(clev, "all")
})

test_that("clever works with some custom parameter settings", {
  data(Dat2)
  clev <- clever_multi(
    Dat2,
    projections = c("PCATF", "ICA_kurt"), 
    measures = "all",
    kurt_quantile = .95,
    outlier_cutoffs = list(leverage=5, robdist=.99),
    verbose = TRUE
  )
  myplot <- plot(clev, "all", title = "My Plot")
})

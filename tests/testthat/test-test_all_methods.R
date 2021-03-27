test_that("clever works with all methods and default settings", {
  data(Dat1)
  clev <- clever:::clever_multi(
    Dat1,
    projection = "all"
  )
  #myplot <- plot(clev, "all")
})

test_that("clever works with some custom parameter settings", {
  data(Dat2)
  clev <- clever:::clever_multi(
    Dat2,
    projection = c("PCATF", "ICA_kurt"),
    kurt_quantile = .90,
    cutoff = 5,
    verbose = TRUE
  )
  #myplot <- plot(clev, "all", title = "My Plot")
})
test_that("clever works with all methods", {
  data(Dat1)
  clev <- clever(Dat1, 
    projection_methods="all", 
    outlyingness_methods="all")
})

# Build --> Install and Restart
# [Edit this] path to the Workbench for your computer.
my_wb <- "../workbench"

library(testthat)
library(ciftiTools)
if (interactive()) { ciftiTools.setOption("wb_path", my_wb) }

library(fMRIscrub)

library(ggplot2)
library(cowplot)
library(ica)

test_check("fMRIscrub")

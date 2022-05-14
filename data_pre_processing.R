setwd("/Users/atiteyk2/Documents/matlab_manuscript")
getwd()

############## Pre-processing lung cancer mutation binary dataset
library("readxl")

lung_data <- read_excel("dataset/dataset_3.xlsx")
lung_da <- data.frame(lung_data)
lung_convert_data <- apply(lung_da, 2, as.numeric)
lung_remove_first_column <- lung_convert_data[,-1]
lung_remove_column_name <- unname(lung_remove_first_column)
lung_mutational_data <- lung_remove_column_name/colSums(lung_remove_column_name)
lung_mutational_data <- t(lung_mutational_data)
lung_mutational_data <- data.frame(lung_mutational_data)

############## Convert lung mutational data into .mat file for Matlab implementation
writeMat("lung_mutational_data.mat", lung_data = lung_mutational_data)

path_dat <- "../data/extdata_PintMF/"
file1 <- readRDS(file.path(path_dat, "OMICSSIMLA_v2.rds"))
file2 <- readRDS(file.path(path_dat , "OMICSSIMLA.rds"))
file <- c(file1, file2)
n_batch <- length(file)

print(str(file[[1]]))


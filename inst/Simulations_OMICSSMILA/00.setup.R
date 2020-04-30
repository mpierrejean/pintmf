path_dat <- "inst/extdata/"

file <- readRDS(file.path(path_dat, "OMICSSIMLA.rds"))

n_batch <- length(file)

print(str(file[[1]]))


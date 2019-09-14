# library imports
library(UpSetR)

# example of list input (list of named vectors)
listInput <- list(TF1 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF2 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF3 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF4 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF5 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF6 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF7 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF8 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF9 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)),
                  TF10 = ceiling(runif(ceiling(runif(1, 1, 100)), 1, 25)))


upset(data = fromList(listInput), nsets = 10, nintersects = 40,
      mb.ratio = c(0.5, 0.5), order.by = "freq", point.size = 4,
      empty.intersections = "on")

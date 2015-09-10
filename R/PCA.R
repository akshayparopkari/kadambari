# Abstract: Performs PCA and plots predicted results with labeling groups and
# Normal probability ellipsoids from relative abundance files.
#
# Date: 09/04/2015
#
# Author: Akshay Paropkari

require(MASS)
require(FactoMineR)
require(rgl)
require(bpca)
require(ggbiplot)
require(graphics)

# Read data from input file
rel.abd.data <- read.delim(file = 'insert your input file')
pca.data <- rel.abd.data[, 3:ncol(rel.abd.data)] # Category and OTU rel abd
pca.cat <- rel.abd.data[, 2] # Category column in input file

# Train PCA with 1/2 dataset
preg.pca <- prcomp(x = pca.data[1:floor(ncol(pca.data)/2)],
                   center = T,
                   scale. = T)

# Test predicting PCA results on full dataset
pca.predict <- predict(preg.pca, pca.data)

summary(preg.pca)

# Plot 2D PCA vector plot
g <- ggbiplot(preg.pca,
              scale = 1,
              var.scale = 1,
              obs.scale = 1,
              groups = pca.cat,
              ellipse = T,
              var.axes = T,
              varname.size = 4,
              varname.adjust = 2,
              cex = 5) +
              scale_color_discrete(name = '') +
              theme(legend.direction = 'horizontal',
                    legend.position = 'top')

print(g)
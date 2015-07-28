# Please use this script at your own risk!
#
# Abstract: Performs LDA, computes its accuracy and statistical significance
#           from relative abundance files.
#
# Date: 07/09/2015
#
# Author: Akshay Paropkari

# import LDA libraries
require(MASS)
require(scatterplot3d)
require(rgl)

# import data file
rel.abd.data <- read.delim(rel_abd_file)

# filter out input data, iff it contains surplus
lda.data <- rel.abd.data[, 2:220]

# calculate LDA and predict classes
r <- lda(Category ~ ., rel.abd.data[, 2:111])  # train with 50% of the dataset
prop <- (r$svd ^ 2 / sum(r$svd ^ 2)) * 100     # axis variance in %
r.pred <- predict(r, lda.data)                 # test LDA on input data
r.df <- data.frame(r.pred)                     # save results into dataframe

# calculate confusion matrix to measure the accuracy of LDA predictions
lda_accuracy <- table(r.pred$class, rel.abd.data$Category)

# significance test on all linear discriminent
ld1.aov <- anova(object = lm(r.pred$x[, 1] ~ rel.abd.data$Category))
ld2.aov <- anova(object = lm(r.pred$x[, 2] ~ rel.abd.data$Category))
ld3.aov <- anova(object = lm(r.pred$x[, 3] ~ rel.abd.data$Category))

# manova on lda results
fit <- manova(cbind(r.pred$x[, 1], r.pred$x[, 2], r.pred$x[, 3]) ~ r.pred$class)
print(summary(object = fit, test = 'Wilks'))

# creating colors column for each category
r.df$color <- NULL
r.df$color[r.df$class == 'Cat1'] <- '#4dac26'
r.df$color[r.df$class == 'Cat2'] <- '#b8e186'
r.df$color[r.df$class == 'Cat3'] <- '#d01c8b'
r.df$color[r.df$class == 'Cat4'] <- '#f1b6da'

# plot LDA prediction results
par(mar = c(5, 5, 4, 2))
plot(r.pred$x[, 1:2],
     pch = 19,
     col = r.df$color,
     cex = 1.25,
     xlab = paste('LD1 (',prop[1],'%)'),  # add % variance due to LD1
     ylab = paste('LD2 (',prop[2],'%)'))  # add % variance due to LD2
legend('topright',
       pch = 19,
       col = c('#f1b6da', '#d01c8b', '#b8e186', '#4dac26'),
       legend = c('Cat4', 'Cat3', 'Cat2', 'Cat1'))
abline(h = 0, v = 0, col = 'blue')

# plot 3D LDA prediction results
scatterplot3d(x = r.pred$x[, 3],
              y = r.pred$x[, 1],
              z = r.pred$x[, 2],
              color = r.df$color,
              pch = 19,
              xlab = paste('LD3 (', prop[3], '%)'),  # add % variance due to LD1
              ylab = paste('LD1 (', prop[1], '%)'),  # add % variance due to LD2
              zlab = paste('LD2 (', prop[2], '%)'),  # add # variance due to LD3
              box = FALSE,
              cex.symbols = 1.75,
              angle = 20)
legend('topright',
       pch = 19,
       col = c('#f1b6da', '#d01c8b', '#b8e186', '#4dac26'),
       legend = c('Cat4', 'Cat3', 'Cat2', 'Cat1'))
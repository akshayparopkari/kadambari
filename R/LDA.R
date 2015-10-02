# Abstract: Performs LDA, computes its accuracy and statistical significance
#           from unifrac distance matrix files.
#
# Date: 10/02/2015
#
# Author: Akshay Paropkari

# import LDA libraries
require(MASS)
require(scatterplot3d)
require(rgl)

# import unifrac distance data file
unifrac.data <- read.table(file = unifrac_dm_file,
                           header = T,
                           sep = '\t',
                           fill = F,
                           strip.white = T)

# train with 50% of the dataset
r <- lda(Condition ~ .,
         data = unifrac.data[, 1:floor(x = ncol(x = unifrac.data)/2)])
prop <- (r$svd ^ 2 / sum(r$svd ^ 2)) * 100     # axis variance in %
r.pred <- predict(r, unifrac.data)             # test LDA on input data
r.df <- data.frame(r.pred)                     # save results into dataframe

# calculate confusion matrix to measure the accuracy of LDA predictions
lda_accuracy <- table(r.pred$class, unifrac.data$Condition)
print(lda_accuracy)

# significance test on all linear discriminent
ld1.aov <- anova(object = lm(r.pred$x[, 1] ~ unifrac.data$Condition))
ld2.aov <- anova(object = lm(r.pred$x[, 2] ~ unifrac.data$Condition))
ld3.aov <- anova(object = lm(r.pred$x[, 3] ~ unifrac.data$Condition))
print(ld1.aov)
print(ld2.aov)
print(ld3.aov)

# manova on lda results
fit <- manova(cbind(r.pred$x[, 1], r.pred$x[, 2], r.pred$x[, 3]) ~ r.pred$class)
print(summary(object = fit, test = 'Wilks'))

# creating colors column for each condition
r.df$color <- NULL
r.df$color[r.df$class == 'Cat1'] <- '#4dac26'
r.df$color[r.df$class == 'Cat2'] <- '#b8e186'
r.df$color[r.df$class == 'Cat3'] <- '#d01c8b'
r.df$color[r.df$class == 'Cat4'] <- '#f1b6da'

# plot 2D LDA prediction results
par(mar = c(5, 5, 4, 2))
plot(r.pred$x[, 1:2],
     pch = 21,
     bg = r.df$color,
     col = '#000000',
     cex = 2.5,
     xlab = paste('LD1 (', round(prop[1], 3),'%)'),  # % variance due to LD1
     ylab = paste('LD2 (', round(prop[2], 3),'%)'))  # % variance due to LD2
legend('topright',
       pch = 19,
       col = c('#f1b6da', '#d01c8b', '#b8e186', '#4dac26'),
       legend = c('Cat4', 'Cat3', 'Cat2', 'Cat1'),
       pt.cex = 2)
abline(h = 0, v = 0, col = 'blue')

 # plot 3D LDA prediction results
scatterplot3d(x = r.pred$x[, 1],
              y = r.pred$x[, 3],
              z = r.pred$x[, 2],
              pch = 21,
              bg = r.df$color,
              color = '000000',
              xlab = paste('LD1 (', round(prop[1], 3), '%)'),  # % variance due to LD1
              ylab = paste('LD3 (', round(prop[3], 3), '%)'),  # % variance due to LD3
              zlab = paste('LD2 (', round(prop[2], 3), '%)'),  # # variance due to LD2
              box = F,
              cex.symbols = 2,
              angle = 30)
legend('topright',
    pch = 19,
    col = c('#f1b6da', '#d01c8b', '#b8e186', '#4dac26'),
    legend = c('Cat4', 'Cat3', 'Cat2', 'Cat1'),
    pt.cex = 2)

# plot interactive 3D plot
plot3d(x = r.pred$x[, 1],
       y = r.pred$x[, 3],
       z = r.pred$x[, 2],
       col = r.df$color,
       pch = 21,
       xlab = paste('LD1 (', round(prop[1], 3), '%)'),  # add % variance due to LD1
       ylab = paste('LD3 (', round(prop[3], 3), '%)'),  # add % variance due to LD3
       zlab = paste('LD2 (', round(prop[2], 3), '%)'),  # add # variance due to LD2
       box = T,
       size = 10,
       axes = T)

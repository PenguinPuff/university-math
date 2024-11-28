---
title: "Exercise04_Clustering"
output:
  pdf_document: default
  html_document: default
date: "`r format(Sys.time(), '%d %B, %Y')`"
editor_options: 
  markdown: 
    wrap: sentence
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,warning = FALSE,message = FALSE)
```

 

## Lab 01: K-Means Clustering

*This exercise has been adapted from James, G., Witten, D., Hastie, T., & Tibshirani, R. (2017). An Introduction to Statistical Learning with Applications in R.*

The function `kmeans()` performs K-means clustering in `R`.\
We begin with a simple simulated example in which there truly are two clusters in the data: the first 25 observations have a mean shift relative to the next 25 observations.

```{r clustering simulated example}
set.seed(2)
x <- matrix(rnorm(50 * 2), ncol = 2)
x[1:25, 1] <- x[1:25, 1] + 3
x[1:25, 2] <- x[1:25, 2] - 4
```

 

We now perform K-means clustering with K=2.

```{r kmeans clustering}
km.out <- kmeans(x, 2, nstart = 20)
```

 

The cluster assignments of the 50 observations are contained in `km.out$cluster`.

```{r cluster assignments}
km.out$cluster
```

 

The K-means clustering perfectly separated the observations into two clusters even though we did not supply any group information to `kmeans()`.
We can plot the data, with each observation colored according to its cluster assignment.

```{r chunk31}
plot(x, col = (km.out$cluster + 1),
    main = "K-Means Clustering Results with K = 2",
    xlab = "", ylab = "", pch = 20, cex = 2)
```

When performing K-means clustering, it is important to set a random seed using the `set.seed()` function.
This way, the initial cluster assignments can be replicated, and the output will be fully reproducible.

 

### Determining K - The Silhouette Coefficient

Especially for more complex data sets, it is not intuitively clear, what number of clusters should be input into the clustering command.
There are two methods to determine the value of `k`: the elbow method and the silhouette coefficient.
In the following, we focus on the latter.

To show their effectiveness, we load a new dataset that is available within the R distribution.
We select the `mtcars` dataset, which was extracted from the 1974 Motor Trend US magazine, and 11 aspects of automobile design and performance for 32 automobiles.

We need to scale the dataset prior to conducting the clustering.

```{r scaling}
scaled_data <- as.data.frame(scale(mtcars))
```

 

Silhouette analysis allows you to calculate how similar each observation is with the cluster it is assigned relative to other clusters.
This metric ranges from -1 (observation in the wrong cluster) to 1 (observation well matched to assigned cluster) for each observation.
We can determine the number of clusters k using the average silhouette width.
We pick the `k` which maximizes that score.

To run the silhouette coefficient, we need to import the `cluster`, `tidyverse`, and `purrr` libraries.

```{r silhouette}
library(cluster) 
library(tidyverse) 
library(purrr) 

# Use map_dbl to run many models with varying value of k
sil_width <- map_dbl(2:6,  function(k){
  model <- pam(x = scaled_data, k = k)
  model$silinfo$avg.width
})
# Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = 2:6,
  sil_width = sil_width
)
# Plot the relationship between k and sil_width
ggplot(sil_df, aes(x = k, y = sil_width)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 2:6)
```

In this example, based on the slihouette coefficient, we select `k = 4`.

 

## Lab 02: Hierarchical Clustering

The `hclust()` function implements hierarchical clustering in `R`.
In the following example, we use the `mtcars` data from the previous lab to plot the hierarchical clustering dendrogram with Euclidean distance as the dissimilarity measure.

We start by calculating the distance matrix, which is a 32 x 32 matrix for the given dataset.

```{r distance matrix}
dist_mat <- dist(scaled_data, method = 'euclidean')
```

 

Now, we can plot the corresponding hierarchical clustering dendrogram.

```{r hierarchical clustering}
plot(hclust(dist(scaled_data), method = "complete"),
    main = "Hierarchical Clustering with Scaled Features")
```

 

If you visually want to see the clusters on the dendrogram you can use the `abline()` function to draw the cut line and superimpose rectangular compartments for each cluster on the tree with the `rect.hclust()` function:

```{r print clusters}
plot(hclust(dist(scaled_data), method = "complete"),main = "Hierarchical Clustering with Scaled Features - Clusters")
rect.hclust(hclust(dist(scaled_data), method = "complete") , k = 4, border = 2:6)
abline(h = 4, col = 'red')
```

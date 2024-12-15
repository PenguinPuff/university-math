x1 = c(1, 1, 0, 5, 6, 4)
x2 = c(4, 3, 4, 1, 2, 0)
x = cbind(x1, x2)
plot(x[,1], x[,2],
xlab = expression(X[1]),
ylab = expression(X[2]),
bty = "l",
     pch = 16)
text(x[,1]+0.15, x[,2], 1:6)


set.seed(2100)
labels = sample(2, nrow(x), replace = T)
x = cbind(x, labels)
x

plot(x[,1], x[,2],
xlab = expression(X[1]),
ylab = expression(X[2]),
bty = "l",
col = labels+2,
pch = 16)
text(x[,1]+0.15, x[,2], 1:6)
legend("topright", legend = c("Cluster 1", "Cluster 2"),
col = c(1,2)+2, pch = 16, bty = "n")

# Cluster one
centroid_one = colMeans(x[x[,3] == 1,, drop = FALSE])[1:2]
centroid_one

# Cluster two
centroid_two = colMeans(x[x[,3] == 2,, drop = FALSE])[1:2]
centroid_two

plot(x[,1], x[,2],
xlab = expression(X[1]),
ylab = expression(X[2]),
col = labels+2,
bty = "l",
pch = 16)
text(x[,1]+0.15, x[,2], 1:6)
legend("topright", legend = c("Cluster 1", "Cluster 2"),
col = c(1,2)+2, pch = 16, bty = "n")
points(centroid_one[1], centroid_one[2], pch = 3, col = 3, cex = 2)
points(centroid_two[1], centroid_two[2], pch = 3, col = 4, cex = 2)



labels = c(2, 2, 2, 1, 1, 1)
x[,3] = labels
plot(x[, 1], x[, 2],
xlab = expression(X[1]),
ylab = expression(X[2]),
bty = "l",
col = labels + 2,
pch = 16)
text(x[,1]+0.15, x[,2], 1:6)
legend("topright", legend = c("Cluster 1", "Cluster 2"),
col = c(1,2)+2, pch = 16, bty = "n")
points(centroid_one[1], centroid_one[2], pch = 3, col = 3, cex = 2)
points(centroid_two[1], centroid_two[2], pch = 3, col = 4, cex = 2)




# Cluster one
centroid_one = colMeans(x[x[,3] == 1,, drop = FALSE])[1:2]
centroid_one


# Cluster two
centroid_two = colMeans(x[x[,3] == 2,, drop = FALSE])[1:2]
centroid_two


plot(x[,1], x[,2],
xlab = expression(X[1]),
ylab = expression(X[2]),
col = labels+2,
bty = "l",
pch = 16)
text(x[,1]+0.15, x[,2], 1:6)
legend("topright", legend = c("Cluster 1", "Cluster 2"),
col = c(1,2)+2, pch = 16, bty = "n")
points(centroid_one[1], centroid_one[2], pch = 3, col = 3, cex = 2)
points(centroid_two[1], centroid_two[2], pch = 3, col = 4, cex = 2)





# Scale the data
data_scaled = scale(USArrests)
# Set seed for reproducibility
set.seed(2100)
# Create the complete linkage cluster based on scaled data.
# 'dist()' computes the Euclidean distance between the observations.
complete_cluster_scaled = hclust(dist(data_scaled), method="complete")
# Plot the complete linkage cluster with scaled data
plot(complete_cluster_scaled)


# Set seed for reproducibility
set.seed(2021)
# Number of classes
n_classes = 3
# Number of obs in each class
n_obs = 20
# Number of variables
n_var = 50
# Generate the data
# Can increase the mean to separate the clusters of PC more.
data =
matrix(sapply(1:3,
function(x) rnorm(n_obs*n_var, mean = 12*sqrt(x))),
ncol=n_var) #, byrow = TRUE) yield completely different values.
dim(data)



pr_results = prcomp(data)
plot(pr_results$x[, 1:2],
col = true_labels,
bty = "l",
xlab = expression(Z[1]),
ylab = expression(Z[2]),
pch = 16)

K = 3
km_results <- kmeans(data, K)
table(true_labels, km_results$cluster)


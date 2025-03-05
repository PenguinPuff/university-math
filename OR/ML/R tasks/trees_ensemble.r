---
title: "Exercise03_Trees&Ensembles"
output:
  html_document: default
  pdf_document: default
date: "`r format(Sys.time(), '%d %B, %Y')`"
editor_options: 
  markdown: 
    wrap: sentence
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,warning = FALSE,message = FALSE)
```

## Lab 01: Classification Trees

*This exercise has been adapted from James, G., Witten, D., Hastie, T., & Tibshirani, R. (2017). An Introduction to Statistical Learning with Applications in R.* We use their library `ISLR2` for relevant data sets.

### Dataset & Libraries

The `tree` library is used to construct classification and regression trees.

```{r tree library}
library(tree)
```

We first use classification trees to analyze the `Carseats` data set, which contains information about car seat sales in 400 different stores.
The data set contains categorical predictors such as `Shelveloc` indicating the quality of the shelving location of the seat at each location as well as ordinal such as the `Income` level around the corresponding store.

```{r library}
library(ISLR2)
attach(Carseats)
names(Carseats)
summary(Carseats)
```

 

In these data, `Sales` is a continuous variable, which we want to change to a binary variable.
We use the `ifelse()` function to create a variable, called `High`, which takes on a value of `Yes` if there are at least $8$ `Sales` and takes on a value of `No` otherwise.
The threshold of high sales has been selected randomly based on the mean value for sales per store.

```{r sales}
High <- factor(ifelse(Sales <= 8, "No", "Yes"))
```

Finally, we use the `data.frame()` function to merge `High` with the rest of the `Carseats` data.

```{r carseats update}
Carseats <- data.frame(Carseats, High)
```

 

### Fitting the Classification Tree

We now use the `tree()` function to fit a **classification tree** in order to predict `High` using all variables but `Sales`.
The syntax of the `tree()` function is quite similar to that of the `lm()` function.

```{r classification tree}
tree.carseats <- tree(High ~ . - Sales, Carseats)
```

The `summary()` function lists the variables that are used as internal nodes in the tree, the number of terminal nodes, and the (training) error rate.

```{r summary of classification tree}
summary(tree.carseats)
```

The output also shows us the training error rate of $9\%$ as well as the deviance of $46\%$.
It is important to understand how the respective libraries define characteristics such as deviance.
For this example, the deviance is closely related to the definition of entropy, but not exactly the same (for more information on the `tree` library, please refer to its [documentation](https://cran.r-project.org/web/packages/tree/tree.pdf))

 

### Plotting

We can use the `plot()` function to display the tree structure, and the `text()` function to display the node labels.
The argument `pretty = 0` instructs `R` to include the category names for any qualitative predictors, rather than simply displaying a letter for each category.

```{r plot tree, fig.height=12, fig.width=18}
plot(tree.carseats)
text(tree.carseats, pretty = 0)
```

The most important indicator of `Sales` appears to be shelving location, since the first branch differentiates `Good` locations from `Bad` and `Medium` locations.

 

### Evaluating

In order to properly evaluate the performance of a classification tree, we split the observations into a training set and a test set, build the tree using the training set, and evaluate its performance on the test data.
The `predict()` function can be used for this purpose.
In the case of a classification tree, the argument `type = "class"` instructs R to return the actual class prediction.

```{r performance}
set.seed(2)
train <- sample(1:nrow(Carseats), 200)
carseats.test <- Carseats[-train, ]
High.test <- High[-train]
tree.carseats <- tree(High ~ . - Sales, Carseats, subset = train)
tree.pred <- predict(tree.carseats, carseats.test, type = "class")
table(tree.pred, High.test)
(104 + 50) / 200
```

We see that we are able to correctly predict for around $77\,\%$ of the locations in the test data set (classifcation accuracy).

## Lab 02: Ensemble Methods

### Bagging & Random Forests

we apply bagging and random forests to the `Carseats` data, using the `randomForest` package.
The exact results obtained in this section may depend on the version of `R` and the version of the `randomForest` package installed on your computer.

Recall that bagging is simply a special case of a random forest with $m=p$.
Therefore, the `randomForest()` function can be used to perform both random forests and bagging.
We perform bagging as follows indicating with `mtry=10` that we want to consider all ten predictors for each split of the tree:

```{r bagging}
library(randomForest)
set.seed(1)
bag.carseats <- randomForest(High ~ . - Sales, Carseats, subset = train, mtry = 10, importance = TRUE)
bag.carseats
```

 

Growing a **random forest** proceeds in exactly the same way, except that we use a smaller value of the `mtry` argument.

```{r random forest}
set.seed(1)
rf.carseats <- randomForest(High ~ . - Sales, Carseats, subset = train, mtry = 3, importance = TRUE)
rf.carseats
```

The reported **Out-of-bag (OOB) error estimate** is the mean prediction error on each training sample using only trees that did not have this training individual in their bootstrap sample.

 

Using the `importance()` function, we can view the **importance of each variable**.
Plots of these importance measures can be produced using the `varImpPlot()` function.

```{r importance}
importance(rf.carseats)
```

Two measures of variable importance are reported.
The first is based upon the mean decrease of accuracy in predictions on the OOB samples when a given variable is permuted.
The second is a measure of the total decrease in node impurity that results from splits over that variable, averaged over all trees, measured by the deviance.

```{r importance plot}
varImpPlot(rf.carseats)
```

### Boosting

Here we use the `gbm` package, and within it the `gbm()` function, to fit boosted regression trees to the `Boston` data set.
We run `gbm()` with the option `distribution = "bernoulli"` since this is a binary classification problem.
The argument `n.trees = 5000` indicates that we want $5000$ trees, and the option `interaction.depth = 4` limits the depth of each tree.

```{r boosting}
library(gbm)
set.seed(1)
boost.carseats <- gbm(High ~ . - Sales, data=Carseats[train, ], distribution = "gaussian", n.trees = 5000, interaction.depth = 4)
summary(boost.carseats)
```

The `summary()` function produces a relative influence plot and also outputs the relative influence statistics.
We see that `Price` and `ShelveLoc` are the most important variables.

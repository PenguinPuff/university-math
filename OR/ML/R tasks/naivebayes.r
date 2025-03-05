---
title: "Exercise02_NaiveBayes"
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

## Lab 01: Naive Bayes for Classification of Stock Market Data

### Data Set

*This exercise has been adapted from James, G., Witten, D., Hastie, T., & Tibshirani, R. (2017). An Introduction to Statistical Learning with Applications in R.* We use their library `ISLR2` for relevant data sets in the following.

For the implementation of a Naive Bayes Classifier, we utilize the `Smarket` data set, which contains percentage returns for the S&P 500 stock index over 1,250 days from 2001 to 2005.

```{r read smarket data}
library(ISLR2)
names(Smarket)
summary(Smarket)
```

 

For each date, we have recorded the percentage returns for each of the five previous trading days, `Lag1` through `Lag5`.
We have also recorded `Volume` (the number of shares traded on the previous day, in billions), `Today` (the percentage return on the date in question) and `Direction` (whether the market was Up or Down on this date).

 

To access the objects in the data set easier, we can `attach` it to the R search path.
This allows us to access the objects by simply stating their names.

```{r attach}
attach(Smarket)
```

 

To generate a **good training - testing data split**, we take the first four years as our training data.
The elements of the `train` vector that correspond to observations that occurred before 2005 are set to `TRUE`, whereas those that correspond to observations in 2005 are set to `FALSE`.
The object `train` is a Boolean vector, since its elements are `TRUE` and `FALSE`.
Boolean vectors can be used to obtain a subset of the rows or columns of a matrix.
For instance, the command `Smarket[train, ]` would pick out a submatrix of the stock market data set, corresponding only to the dates before 2005, since those are the ones for which the elements of `train` are `TRUE`.
Additionally, we can use the `!` symbol to reverse all of the elements of a Boolean vector.
Thus, `!train` is a Boolean vector whose elements corresponding to observations in 2005 are set to `TRUE` while all prior are set to `FALSE`.

```{r test train split}
train <- (Year < 2005)
Smarket.2005 <- Smarket[!train, ]
Direction.2005 <- Direction[!train]
```

 

### Naive Bayes in R

Naive Bayes is implemented in `R` using the `naiveBayes()` function, which is part of the `e1071` library.
By default, this implementation of the naive Bayes classifier models each quantitative feature using a Gaussian distribution.
For this example, we assume that we can classify the direction of the stock market based on `Lag1` and `Lag2` alone.

```{r naivebayes}
library(e1071)
nb.fit <- naiveBayes(Direction ~ Lag1 + Lag2, data = Smarket,
     subset = train)
nb.fit
```

 

The output contains the estimated mean and standard deviation for each variable in each class.
For example, the mean for `lag1` is $0.0428$ for `Direction=Down`, and the standard deviation is $1.23$.

 

The `predict()` function is straightforward.

```{r predictI}
nb.class <- predict(nb.fit, Smarket.2005)
```

 

Now we compare the results of the prediction with the actual movements of the market over that time period, i.e., in the year 2005.

```{r predictII}
table(nb.class, Direction.2005)
```

 

We can also determine the rate of correctly assigning the classification in the test set as well as its opposite value, the test set error rate.

```{r predictIII}
mean(nb.class == Direction.2005)
mean(nb.class != Direction.2005)
```

 

Naive Bayes in this setting outputs accurate predictions over $59\%$ of the time.
The `predict()` function can also generate estimates of the probability that each observation belongs to a particular class.

```{r predictIV}
nb.preds <- predict(nb.fit, Smarket.2005, type = "raw")
nb.preds[1:5, ]
```

Again, if no data set is supplied to the `predict()` function, then the probabilities are computed for the training data that was used to fit the model.

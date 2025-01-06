# A script to produce a set of regularized logistic regression models
#
# Figure 8D
#

library(glmnet)
library(pROC)
library(ggplot2)
library(ggpubr)
library(this.path)
library(tidyverse)
library(fitdistrplus)


rm(list=ls())

# Import expression sets containing BTM module scores
exp_FC <- readRDS(paste0(this.dir(), "/exp_FC.RDS"))

# Regularized logistic regression models on the training data and evaluation on the validation set --------

logpred <- function(X_train, X_validate, y_train, y_validate, percentile, sep_test, test.seq, names, FC_exp, folds) {
  
  # percentile (character, such as "20%") is a value of the percentile range of residuals 
  # that we want to include in the training data
  
  # sep_test (logical, character) is an indicator whether we want to subset the validation and test data 
  # using the same percentile or use all samples for making predictions
  
  # test.seq (integer vector) is the index of datasets in the exp_FC list that we assign to be 
  # used as test data 
  
  # names is the named vector of simplified names for objects in the exp_FC list
  
  # Set seed to ensure reproducibility
  set.seed(7)
  
  subs <- sep_test[[2]] # For the correct file path
  
  # Range of alpha values
  alpha_vals <- seq(0, 1, 0.01)
  
  # Subset the samples in the training data
  q <- quantile(y_train, prob = seq(0, 1, 0.1))
  idx <- y_train < q[percentile[1]] | y_train > q[percentile[2]]
  y_train <- y_train[idx]
  y_train <- ifelse(y_train > 0, 1, 0) # Convert to binary outcome
  X_train <- X_train[, idx]
  
  # Transpose predictor variables and normalize by column
  X_train <- scale(t(X_train))
  X_validate <- scale(t(X_validate))
  
  # Elastic net regularization
  cv_models <- list()
  for (alpha in alpha_vals) {
    cat("Alpha = ", alpha, "\n")
    cv_models[[paste0("alpha_", alpha)]] <- cv.glmnet(X_train, y_train, type.measure = "auc", alpha = alpha, family = "binomial", nfolds = folds)
  }
  
  # Evaluate performance on validation dataset and store AUC
  auc_vals <- vector()
  for (model in cv_models) {
    lambda_best <- model$lambda.min
    preds <- predict(model, newx = X_validate, s = lambda_best, type = "response")
    auc_vals <- c(auc_vals, auc(roc(y_validate, as.vector(preds))))
  }
  
  # Best alpha-lambda combinations for top 10 models
  top10_indices <- order(auc_vals, decreasing = TRUE)[1:10]
  top10_alphas <- alpha_vals[top10_indices]
  top10_lambdas <- sapply(top10_indices, function(index) cv_models[[index]]$lambda.min)
  
  # Re-training on the combined dataset and evaluation on test datasets: 
  
  # Combine training and validation datasets
  X_combined <- rbind(X_train, X_validate)
  y_combined <- c(y_train, y_validate)
  
  # Normalize combined data
  X_combined <- scale(X_combined)
  
  # Re-train models
  final_models <- list()
  feature_weights <- list()
  auc_results <- matrix(0, nrow=10, ncol=(length(test.seq)+2))  # columns for training, validation, and each test dataset
  
  for (i in 1:10) {
    alpha <- top10_alphas[i]
    lambda <- top10_lambdas[i]
    
    # Train the model on combined data
    model <- glmnet(X_combined, y_combined, alpha = alpha, lambda = lambda, family = "binomial")
    final_models[[i]] <- model
    
    # Store features and weights
    feature_weights[[i]] <- coef(model)
    
    # AUC for training and validation
    preds_train <- predict(model, newx = X_train, type = "response")
    preds_validate <- predict(model, newx = X_validate, type = "response")
    auc_results[i, 1] <- auc(roc(y_train, as.vector(preds_train)))
    auc_results[i, 2] <- auc(roc(y_validate, as.vector(preds_validate)))
    
    # AUC for test datasets
    for (j in test.seq) {
      X_test <- scale(t(exprs(exp_FC[[j]])))  # transpose and normalize
      y_test <- exp_FC[[j]]$ab.resid.binary
      preds_test <- predict(model, newx = X_test, type = "response")
      auc_results[i, j-1] <- auc(roc(y_test, as.vector(preds_test)))
    }
  }
  
  # Assign the final models list to a variable specific to that cutoff value
  assign(paste("final.models.", percentile[3]), final_models)
  assign(paste("feature.weights.", percentile[3]), feature_weights)
  assign(paste("auc.resuts.", percentile[3]), auc_results)
  
  system(paste0("mkdir ", this.dir(), "/", folds, ".fold"))
  system(paste0("mkdir ", this.dir(), "/", folds, ".fold/", subs))
  system(paste0("mkdir ", this.dir(), "/", folds, ".fold/", subs, "/Cutoff.", percentile[3]))
  
  saveRDS(final_models, file = paste0(this.dir(), "/", folds, ".fold/", subs, "/Cutoff.", percentile[3], "/final.models.RDS"))
  saveRDS(feature_weights, file = paste0(this.dir(), "/", folds, ".fold/", subs, "/Cutoff.", percentile[3], "/feature.weights.RDS"))
  saveRDS(auc_results, file = paste0(this.dir(), "/", folds, ".fold/", subs, "/Cutoff.", percentile[3], "/AUC.results.RDS"))
  
  for (i in 1:10) {
    pdf(file = paste0(this.dir(), "/", folds, ".fold/", subs, "/Cutoff.", percentile[3], "/Model.", i, ".cutoff=", percentile[3], ".pdf"), width = 11, height = 8)
    # Plotting for training set
    preds <- predict(final_models[[i]], newx = X_train, type = "response")
    roc_obj_train <- roc(y_train, as.vector(preds))
    plot.roc(roc_obj_train, main=paste("Cutoff = ", percentile[1], " Model", i))
    
    # Add validation and test sets to the ROC plot
    preds <- predict(final_models[[i]], newx = X_validate, type = "response")
    roc_obj_val <- roc(y_validate, as.vector(preds))
    lines.roc(roc_obj_val, col="blue")  # validation in blue
    
    legend_labels <- c(
      paste("Training AUC =", round(auc(roc_obj_train), 2)),
      paste("Validation AUC =", round(auc(roc_obj_val), 2))
    )
    
    for (j in test.seq) {
      X_test <- scale(t(exprs(exp_FC[[j]])))  # normalize and transpose
      y_test <- exp_FC[[j]]$ab.resid.binary
      if (sep_test[[1]] == TRUE) {
        y_test <- exp_FC[[j]]$ab_resid
        q <- quantile(y_test, prob = seq(0, 1, 0.1), na.rm = TRUE)
        idx <- y_test < q[percentile[1]] | y_test > q[percentile[2]]
        y_test <- exp_FC[[j]]$ab.resid.binary[idx]
        X_test <- X_test[idx, ]
      }
      
      if (all(y_test == 0) | all(y_test == 1)) {
        next
      }
      
      name <- as.character(names[names(exp_FC)[j]])
      preds <- predict(final_models[[i]], newx = X_test, type = "response")
      roc_obj <- roc(y_test, as.vector(preds))
      lines.roc(roc_obj, col=j)  # using different colors for each test set
      
      legend_labels <- c(legend_labels, 
                         paste(name, " AUC = ", round(auc(roc_obj), 2))
      )
    }
    legend("bottomright", legend=legend_labels, col=c("Black", "Blue", test.seq), lwd=2, cex = 0.5)
    dev.off()
    
    # Generate bar plot for model coefficients
    x <- as.matrix(coef(final_models[[i]]))
    df = data.frame('Coefficients' = x[, 1], 'BTM' = row.names(x))
    df <- subset(df, df$Coefficients != 0)
    df$BTM = factor(df$BTM, levels = df$BTM[order(df$Coefficients)])
    df$color <- "steelblue"
    df$color[df$Coefficients > 0] <- "orangered"
    ggplot(df, aes(x = BTM, y = Coefficients, fill = color)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      xlab("") +
      ylab("Model Coefficient") +
      theme_minimal(base_size=20) + 
      theme(legend.position = "none") + 
      labs(title = paste0("Model ",  i, ", cutoff = ", percentile[3]))
    ggsave(file = paste0(this.dir(), "/", folds, ".fold/", subs, "/Cutoff.", percentile[3], "/Cutoff = ", percentile[3], ".Model.", i, ".model.weights.pdf"), width = 11, height = 8)
  }
  
} # End or function

names <- c("Vax010", 
           "Immunity",
           "CHI", 
           "Pfizer", 
           "RTS,S", 
           "MPSV4")

names(names) <- names(exp_FC)

# Extract predictors and response for training and validation sets
X_train <- merge(exprs(exp_FC[[1]]), exprs(exp_FC[[2]]), by = 0) # Using both Immunity and Vax010 as a single training set. This gives us 79 subjects total
y_train <- c(exp_FC[[1]]$ab_resid, exp_FC[[2]]$ab_resid)
X_validate <- exprs(exp_FC[[3]])  # Using CHI as validation, 19 subjects
y_validate <- exp_FC[[3]]$ab.resid.binary

row.names(X_train) <- X_train[, 1]
X_train <- X_train[, -1]

sep_test = list(FALSE, "Without.subsetting")
# sep_test = list(TRUE, "With.subsetting")


test.seq = c(4:6)

cutoffs <- list(c("10%", "90%", 10), 
                c("20%", "80%", 20), 
                c("30%", "70%", 30), 
                c("40%", "60%", 40), 
                c("50%", "50%", 50)
)

folds <- 10
# folds <- 5

for (cuts in 1:length(cutoffs)) {
  percentile <- cutoffs[[cuts]]
  logpred(X_train, X_validate, y_train, y_validate, percentile, sep_test, test.seq, names, FC_exp, folds)
}

#
# Model #6 at 40% cutoff is shown in Figure 7D
#


# Access the significance of AUC values by permutation analysis -----------

# Will use one model for this, model #6 at 40% cutoff
models <- readRDS(paste0(this.dir(), "/10.fold/Without.subsetting/Cutoff.40/final.models.RDS"))

aucs <- data.frame(matrix(nrow = 1000, ncol = 6))
colnames(aucs) <- names

for (i in 1:1000) {
  for (j in 1:6) {
    X_test <- scale(t(exprs(exp_FC[[j]])))  # normalize and transpose
    y_test <- exp_FC[[j]]$ab.resid.binary
    
    # Permute ab residuals randomly 
    y_test <- y_test[sample.int(length(y_test), size = length(y_test))]
    
    # Predict on the scrambled dataset and calculate AUC
    preds <- predict(models[[6]], newx = X_test, type = "response")
    roc_obj <- roc(y_test, as.vector(preds))
    aucs[i, j] <- auc(roc_obj)
  }
}

auc.probs.model.6 <- data.frame("Dataset" <- colnames(aucs), 
                                "Actual.AUC" = c(0.76,
                                                 0.76, 
                                                 0.94,
                                                 0.73, 
                                                 0.86, 
                                                 0.81), 
                                "p-value" = numeric(6))

for (i in 1:6) {
  fit <- fitdistr(aucs[, i], "normal")
  auc.probs.model.6[i, 3] <- pnorm(auc.probs.model.6[i, 2], 
                                   mean = fit$estimate[1], 
                                   sd = fit$estimate[2], 
                                   lower.tail = FALSE)
}

write.table(auc.probs.model.6, file = paste0(this.dir(), "/Model.6.predictions.AUC.p.values.txt"), 
            row.names = FALSE, 
            sep = "\t")
# Note: due to the stochastic nature of the permutation test, resulting p-values are 
# close estimates of the actual significance, but may (and will) vary slightly 
# from one run to the next. 



# Load necessary libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOSemSim)
library(circlize)
library(GOplot)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(xgboost)
library(pROC)

# Define file paths
file_path <- "...Own your path..."
sample_description_path <- "...Own your path..."

# Read raw counts data
raw_counts <- as.matrix(read.table(file_path, header=TRUE, row.names=1, sep="\t"))

# Read sample description file
sample_description <- read.csv(sample_description_path, header=TRUE, stringsAsFactors=FALSE)

# Filter for IL6 and PBS samples
filtered_samples <- sample_description[sample_description$Condition %in% c("IL6", "PBS"), ]
raw_counts_filtered <- raw_counts[, filtered_samples$Sample.ID]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered, colData = filtered_samples, design = ~ Condition)

# Perform differential expression analysis
dds <- DESeq(dds)
contrast <- c("Condition", "IL6", "PBS")
res <- results(dds, contrast = contrast)

# Filter significant DE genes (adjusted p-value < 0.05)
sig_res <- res[which(res$padj < 0.05), ]
significant_genes <- rownames(sig_res)

# Prepare data for machine learning
gene_expression <- as.data.frame(t(assay(dds)[significant_genes, ]))
gene_expression$Condition <- filtered_samples$Condition
gene_expression$Condition <- as.factor(gene_expression$Condition)

# Set seed for reproducibility
set.seed(123)

# Split data into training (70%) and testing (30%)
train_index <- createDataPartition(gene_expression$Condition, p = 0.7, list = FALSE)
train_data <- gene_expression[train_index, ]
test_data <- gene_expression[-train_index, ]

# Define trainControl for 5-fold cross-validation
train_control <- trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)

# Train models with cross-validation on training set
lasso_model <- train(Condition ~ ., data = train_data, method = "glmnet", family = "binomial", trControl = train_control, metric = "ROC", tuneLength = 10)
rf_model <- train(Condition ~ ., data = train_data, method = "rf", trControl = train_control, metric = "ROC", tuneLength = 10)
svm_model <- train(Condition ~ ., data = train_data, method = "svmRadial", trControl = train_control, metric = "ROC", tuneLength = 10)
xgb_model <- train(Condition ~ ., data = train_data, method = "xgbTree", trControl = train_control, metric = "ROC", tuneLength = 10)
lr_model <- train(Condition ~ ., data = train_data, method = "glm", family = binomial, trControl = train_control, metric = "ROC")

# Evaluate models on test set
lasso_pred <- predict(lasso_model, newdata = test_data)
rf_pred <- predict(rf_model, newdata = test_data)
svm_pred <- predict(svm_model, newdata = test_data)
xgb_pred <- predict(xgb_model, newdata = test_data)
lr_pred <- predict(lr_model, newdata = test_data)

# Print confusion matrices for each model
print("Confusion Matrix for LASSO:")
print(confusionMatrix(lasso_pred, test_data$Condition))

print("Confusion Matrix for Random Forest:")
print(confusionMatrix(rf_pred, test_data$Condition))

print("Confusion Matrix for SVM:")
print(confusionMatrix(svm_pred, test_data$Condition))

print("Confusion Matrix for XGBoost:")
print(confusionMatrix(xgb_pred, test_data$Condition))

print("Confusion Matrix for Logistic Regression:")
print(confusionMatrix(lr_pred, test_data$Condition))

# Collect ROC data for all models
roc_data <- data.frame()

for (model in list(lasso_model, rf_model, svm_model, xgb_model, lr_model)) {
  model_name <- model$method
  roc_curve <- roc(response = model$pred$obs, predictor = model$pred$IL6)
  roc_data <- rbind(roc_data, data.frame(model = model_name, specificity = rev(roc_curve$specificities), sensitivity = rev(roc_curve$sensitivities)))
}

# Plot ROC curves
ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = model)) + 
  geom_line(size = 1) + 
  geom_abline(linetype = "dashed") + 
  labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC Curves for Different Models") +
  theme_minimal() + 
  theme(legend.title = element_blank())

# Identify important genes from each model
lasso_imp <- varImp(lasso_model, scale = FALSE)
rf_imp <- varImp(rf_model, scale = FALSE)
svm_imp <- varImp(svm_model, scale = FALSE)
xgb_imp <- varImp(xgb_model, scale = FALSE)
lr_imp <- varImp(lr_model, scale = FALSE)

# Print top important genes
print("Top important genes identified by LASSO:")
print(lasso_imp)
print("Top important genes identified by Random Forest:")
print(rf_imp)
print("Top important genes identified by SVM:")
print(svm_imp)
print("Top important genes identified by XGBoost:")
print(xgb_imp)
print("Top important genes identified by Logistic Regression:")
print(lr_imp)

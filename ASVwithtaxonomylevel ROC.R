# 필요한 라이브러리 로드
library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)

# ASV level 데이터 로드
asv_data <- read.table("clr_transformed_ASV_table.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
asv_KV_KC <- asv_data %>% select(contains("KV"), contains("KC"))
group_asv <- factor(c(rep("KV", sum(grepl("KV", colnames(asv_KV_KC)))), 
                      rep("KC", sum(grepl("KC", colnames(asv_KV_KC))))), 
                    levels = c("KV", "KC"))

# Taxonomy level 데이터 로드
taxonomy_data <- read.table("microbiometable_taxafix_clr_transformed.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
taxonomy_KV_KC <- taxonomy_data %>% select(contains("KV"), contains("KC"))
group_taxonomy <- factor(c(rep("KV", sum(grepl("KV", colnames(taxonomy_KV_KC)))), 
                           rep("KC", sum(grepl("KC", colnames(taxonomy_KV_KC))))), 
                         levels = c("KV", "KC"))

# ASV level에서 LOOCV로 ROC 계산
asv_probs <- as.data.frame(matrix(nrow = length(group_asv), ncol = 2))
colnames(asv_probs) <- c("KV", "KC")
for (i in 1:length(group_asv)) {
  train_idx <- setdiff(1:length(group_asv), i)
  test_idx <- i
  rf_asv <- randomForest(t(asv_KV_KC[, train_idx]), group_asv[train_idx], ntree = 500)
  asv_probs[test_idx, ] <- predict(rf_asv, t(asv_KV_KC[, test_idx]), type = "prob")
}

# Taxonomy level에서 LOOCV로 ROC 계산
taxonomy_probs <- as.data.frame(matrix(nrow = length(group_taxonomy), ncol = 2))
colnames(taxonomy_probs) <- c("KV", "KC")
for (i in 1:length(group_taxonomy)) {
  train_idx <- setdiff(1:length(group_taxonomy), i)
  test_idx <- i
  rf_taxonomy <- randomForest(t(taxonomy_KV_KC[, train_idx]), group_taxonomy[train_idx], ntree = 500)
  taxonomy_probs[test_idx, ] <- predict(rf_taxonomy, t(taxonomy_KV_KC[, test_idx]), type = "prob")
}

# ASV level에서 multiclass ROC 계산
roc_asv <- multiclass.roc(group_asv, asv_probs)

# Taxonomy level에서 multiclass ROC 계산
roc_taxonomy <- multiclass.roc(group_taxonomy, taxonomy_probs)

# ROC curve plotting
roc_plot <- ggroc(roc_asv, aes(color = "ASV Level")) +
  geom_line(size = 1) +
  ggroc(roc_taxonomy, aes(color = "Taxonomy Level")) +
  geom_line(size = 1) +
  labs(title = "Multiclass ROC Curves for ASV Level and Taxonomy Level", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
  scale_color_manual(values = c("ASV Level" = "blue", "Taxonomy Level" = "red")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(roc_plot)

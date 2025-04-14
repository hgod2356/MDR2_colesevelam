# Load necessary libraries
library(compositions)

# Load the ASV table
asv_data <- read.table("ASV_tablewithtaxonomy_gutmicrobiota.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Separate the Taxonomy column from numeric data
taxonomy <- asv_data$Taxonomy
asv_numeric <- asv_data[, -ncol(asv_data)]

# Replace zeros with a small value to handle clr transformation
asv_numeric[asv_numeric == 0] <- 0.5

# Apply clr transformation
clr_transformed <- clr(asv_numeric)

# Combine clr-transformed data with Taxonomy
clr_data <- cbind(as.data.frame(clr_transformed), Taxonomy = taxonomy)

# Save or view the transformed data
write.table(clr_data, "clr_transformed_ASV_table.txt", sep = "\t", row.names = TRUE, quote = FALSE)


# 필요한 라이브러리 로드
library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)

# 데이터 로드
bacteria_data <- read.table("clr_transformed_ASV_table.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# KV와 KC가 포함된 샘플만 선택
bacteria_KV_KC <- bacteria_data %>% select(contains("KV"), contains("KC"))

# 그룹 벡터 생성 (KV를 control, KC를 case로 지정)
group_KV_KC <- factor(c(rep("KV", sum(grepl("KV", colnames(bacteria_KV_KC)))), 
                        rep("KC", sum(grepl("KC", colnames(bacteria_KV_KC))))),
                      levels = c("KV", "KC"))  # KV를 control, KC를 case로 지정

# 예측 확률을 저장할 벡터 생성
bacteria_predictions <- numeric(length(group_KV_KC))

# LOOCV 적용
for (i in 1:length(group_KV_KC)) {
  # 훈련 및 테스트 세트 분리
  train_idx <- setdiff(1:length(group_KV_KC), i)
  test_idx <- i
  
  # 박테리아 데이터에 대한 Random Forest 모델 생성 및 예측 확률 추출
  rf_bacteria <- randomForest(t(bacteria_KV_KC[, train_idx]), group_KV_KC[train_idx], ntree = 500)
  bacteria_predictions[test_idx] <- predict(rf_bacteria, t(bacteria_KV_KC[, test_idx]), type = "prob")[, 2]
}

# ROC 커브 생성 및 AUC 계산 (KV vs KC)
roc_bacteria <- roc(group_KV_KC, bacteria_predictions, levels = c("KV", "KC"))
plot(roc_bacteria, main = "ROC Curve for Bacteria Data (KV as Control, KC as Case)", col = "blue", legacy.axes = TRUE)
print(paste("AUC for Bacteria Data:", auc(roc_bacteria)))

# ggroc()를 사용하여 ROC 커브 그리기
roc_plot <- ggroc(roc_bacteria, color = "#0072b2", size = 1, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#ee6677") +  # 기준선 추가
  labs(title = "ROC Curve for RandomForest", 
       x = "1-Specificity", 
       y = "Sensitivity") +
  annotate("text", x = 0.75, y = 0.25, 
           label = paste("AUC =", round(auc(roc_bacteria), 3)), 
           color = "black", size = 4) +  # AUC 값 표시
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

print(roc_plot)

## 테두리를 추가한 plot 코드 + x,y축 이름 변경
roc_plot <- ggroc(roc_bacteria, color = "#0072b2", size = 1, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#ee6677") +  # 기준선 추가
  labs(title = "ROC Curve for RandomForest", 
       x = "False Positive Rate", 
       y = "True Positive Rate") +
  annotate("text", x = 0.8, y = 0.15, 
           label = paste("AUC =", round(auc(roc_bacteria), 3)), 
           color = "black", size = 4, hjust = 0.5) +  # AUC 값 표시
  annotate("text", x = 0.8, y = 0.05, 
           label = "Control: KO-VEH\nCase: KO-COL", 
           color = "black", size = 4, hjust = 0.5) +  # Control/Case 정보 추가
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "gray", fill = NA, size = 1)  # 테두리 추가
  )

print(roc_plot)

# 예측 확률 및 중요도 누적용 데이터 구조 초기화
bacteria_predictions <- numeric(length(group_KV_KC))
importance_matrix <- matrix(0, nrow = nrow(bacteria_KV_KC), ncol = 1)  # feature 수에 맞게 초기화
rownames(importance_matrix) <- rownames(bacteria_KV_KC)

# LOOCV 적용
for (i in 1:length(group_KV_KC)) {
  # 훈련 및 테스트 세트 분리
  train_idx <- setdiff(1:length(group_KV_KC), i)
  test_idx <- i
  
  # Random Forest 모델 생성 및 중요도 추출
  rf_bacteria <- randomForest(t(bacteria_KV_KC[, train_idx]), group_KV_KC[train_idx], ntree = 500, importance = TRUE)
  bacteria_predictions[test_idx] <- predict(rf_bacteria, t(bacteria_KV_KC[, test_idx]), type = "prob")[, 2]
  
  # 중요도 값을 누적
  current_importance <- as.vector(importance(rf_bacteria, type = 1)) # type = 1은 GINI 중요도, type = 2는 Mean Decrease Accuracy
  importance_matrix <- importance_matrix + matrix(current_importance, nrow = nrow(bacteria_KV_KC), ncol = 1)
}

# 중요도 평균 계산
mean_importance <- importance_matrix / length(group_KV_KC)

# 중요도 높은 순서대로 정렬
sorted_importance <- sort(mean_importance[, 1], decreasing = TRUE)
print("Top Important Features:")
print(sorted_importance[1:20])  # 상위 20개 feature 출력

# 중요도 전체를 데이터 프레임으로 변환하고 0 이상인 feature만 필터링
importance_df <- data.frame(
  Feature = names(sorted_importance),
  Importance = sorted_importance
)
importance_df_filtered <- importance_df %>% filter(Importance > 0)

# 중요도 바 차트 시각화
ggplot(importance_df_filtered, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 가로로 바 차트 표시
  labs(title = "Feature Importance", x = "Feature", y = "Mean Decrease Gini") +
  theme_minimal()


###
# 필요한 라이브러리 로드
library(dplyr)
library(ggplot2)

# ASV 테이블 로드 (원본 데이터에서 Taxonomy 열 추출)
asv_data <- read.table("clr_transformed_ASV_table.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# ASV 데이터에서 Taxonomy 열 추출
taxonomy <- asv_data$Taxonomy
importance_df_filtered$Taxonomy <- taxonomy[match(importance_df_filtered$Feature, rownames(asv_data))]

# 같은 Taxonomy로 중요도 합산
importance_taxonomy <- importance_df_filtered %>%
  group_by(Taxonomy) %>%
  summarize(Total_Importance = sum(Importance, na.rm = TRUE)) %>%
  arrange(desc(Total_Importance))

# 중요도 시각화 (합산된 Taxonomy 기준)
ggplot(importance_taxonomy, aes(x = reorder(Taxonomy, Total_Importance), y = Total_Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Feature Importance by Taxonomy", x = "Taxonomy", y = "Total Mean Decrease Gini") +
  theme_minimal()

# 중요도 시각화 (합산된 Taxonomy 기준) - Akkermansia muciniphila 강조
ggplot(importance_taxonomy, aes(x = reorder(Taxonomy, Total_Importance), y = Total_Importance)) +
  geom_bar(stat = "identity", aes(fill = Taxonomy == "d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Akkermansiaceae; g__Akkermansia; s__Akkermansia_muciniphila")) +  # 조건부 색상 설정
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "skyblue"), 
                    labels = c("Others", "Akkermansia")) +  # 색상 맵핑
  labs(title = "Feature Importance by Taxonomy", x = "Taxonomy", y = "Total Mean Decrease Gini") +
  theme_minimal() +
  theme(legend.position = "none")


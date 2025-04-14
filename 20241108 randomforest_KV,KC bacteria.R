setwd("D:/Experiment (실험)/11.MDR2 KO/2.Colesevelam/Random forest")
# 필요한 라이브러리 로드
library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)

# 데이터 로드
bacteria_data <- read.table("microbiometable_taxafix_clr_transformed.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

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

#true positive rate = sensitivity and false positive rate = 1 - specificity.

# 테두리를 추가한 plot 코드
roc_plot <- ggroc(roc_bacteria, color = "#0072b2", size = 1, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#ee6677") +  # 기준선 추가
  labs(title = "ROC Curve for RandomForest", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
  annotate("text", x = 0.8, y = 0.1, 
           label = paste("AUC =", round(auc(roc_bacteria), 3)), 
           color = "black", size = 4) +  # AUC 값 표시
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.border = element_rect(color = "gray", fill = NA, size = 1)  # 테두리 추가
  )

print(roc_plot)

## 테두리를 추가한 plot 코드, labeling 추가
roc_plot <- ggroc(roc_bacteria, color = "#0072b2", size = 1, legacy.axes = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#ee6677") +  # 기준선 추가
  labs(title = "ROC Curve for RandomForest", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
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


#############important feature
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

# 중요도 바 차트 시각화 (Akkermansia 강조)
ggplot(importance_df_filtered, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", aes(fill = ifelse(grepl("Akkermansia", Feature), "red", "darkgrey"))) +  # Akkermansia만 빨간색, 나머지 파란색
  coord_flip() +  # 가로로 바 차트 표시
  labs(title = "Feature Importance", x = "Feature", y = "Mean Decrease Gini") +
  scale_fill_identity() +  # 색상 조정
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# 중요도 바 차트 시각화 (Akkermansia 강조, Y축 이름 기울임꼴 처리)
ggplot(importance_df_filtered, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", aes(fill = ifelse(grepl("Akkermansia", Feature), "red", "darkgrey"))) +  # Akkermansia만 빨간색, 나머지 파란색
  coord_flip() +  # 가로로 바 차트 표시
  labs(title = "Feature Importance", x = "Feature", y = "Mean Decrease Gini") +
  scale_fill_identity() +  # 색상 조정
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),  # 크기 조정
    axis.text.y = element_text(face = "italic")  # Y축 레이블을 기울임꼴로 설정
  )




#############KV와 KC 많은 애들을 따로 표기
# 데이터 로드 및 필요한 그룹 선택
bacteria_data <- read.table("microbiometable_taxafix_clr_transformed.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
bacteria_KV_KC <- bacteria_data %>% select(contains("KV"), contains("KC"))

# 그룹 벡터 생성 (KV를 control, KC를 case로 지정)
group_KV_KC <- factor(c(rep("KV", sum(grepl("KV", colnames(bacteria_KV_KC)))), 
                        rep("KC", sum(grepl("KC", colnames(bacteria_KV_KC))))),
                      levels = c("KV", "KC"))

# Random Forest 모델로 중요 feature 계산
rf_bacteria <- randomForest(t(bacteria_KV_KC), group_KV_KC, ntree = 500, importance = TRUE)
feature_importance <- importance(rf_bacteria, type = 1)

# 중요 feature 추출 및 정렬
sorted_importance <- sort(feature_importance[, 1], decreasing = TRUE)
important_features <- names(sorted_importance)

# KV와 KC에서 각 중요한 feature의 평균값 계산
bacteria_data_KV <- bacteria_KV_KC %>% select(contains("KV"))
bacteria_data_KC <- bacteria_KV_KC %>% select(contains("KC"))

feature_means <- data.frame(
  Feature = important_features,
  KV_Mean = rowMeans(bacteria_data_KV[important_features, ]),
  KC_Mean = rowMeans(bacteria_data_KC[important_features, ])
)

# 두 그룹의 평균 차이를 계산하고, 차이가 큰 순서대로 정렬
feature_means <- feature_means %>%
  mutate(Difference = KC_Mean - KV_Mean) %>%
  arrange(desc(abs(Difference)))

# 각 그룹에 더 많이 나타나는 feature에 따라 색상 지정
feature_means$Group <- ifelse(feature_means$Difference > 0, "KC", "KV")

# 중요도 바 차트 시각화
ggplot(feature_means, aes(x = reorder(Feature, Difference), y = Difference, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 가로로 바 차트 표시
  scale_fill_manual(values = c("KV" = "#228833", "KC" = "#E7AC02")) +  # 색상 지정
  labs(title = "Feature Differences", x = "Feature", y = "Difference (KC Mean - KV Mean)") +
  theme_minimal()

# Group 열의 값 변경
feature_means$Group <- recode(feature_means$Group, "KV" = "KO-VEH", "KC" = "KO-COL")

# 중요도 바 차트 시각화 (KO-VEH이 위, KO-COL이 아래로 설정)
ggplot(feature_means, aes(x = reorder(Feature, Difference), y = Difference, fill = Group)) +
  geom_bar(stat = "identity", color = "gray") +  # 테두리 색상을 회색으로 설정
  coord_flip() +  # 가로로 바 차트 표시
  scale_fill_manual(values = c("KO-VEH" = "#228833", "KO-COL" = "#E7AC02"), 
                    breaks = c("KO-VEH", "KO-COL"),  # KO-VEH가 위에, KO-COL이 아래에 오도록 설정
                    labels = c("KO-VEH", "KO-COL")) +  # 범례에 표시될 그룹 이름 지정
  labs(title = "Feature Differences", x = "Feature", y = "Difference (KO-COL Mean - KO-VEH Mean)") +
  theme_minimal()

ggplot(feature_means, aes(x = reorder(Feature, Difference), y = Difference, fill = Group)) +
  geom_bar(stat = "identity", color = "gray") +  # 테두리 색상을 회색으로 설정
  coord_flip() +  # 가로로 바 차트 표시
  scale_fill_manual(values = c("KO-VEH" = "#228833", "KO-COL" = "#E7AC02"), 
                    breaks = c("KO-VEH", "KO-COL"),  # KO-VEH가 위에, KO-COL이 아래에 오도록 설정
                    labels = c("KO-VEH", "KO-COL")) +  # 범례에 표시될 그룹 이름 지정
  labs(title = "Feature Differences", x = "Feature", y = "Difference (KO-COL Mean - KO-VEH Mean)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic")  # Y축 레이블을 기울임꼴로 설정
  )

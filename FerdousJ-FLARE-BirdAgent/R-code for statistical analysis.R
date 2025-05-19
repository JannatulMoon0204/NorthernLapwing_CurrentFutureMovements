#####Jannatul 22/01/25#####
##Edited 19/01/2025
options(scipen = 999)
setwd("D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New")
library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(reshape2)

#####Load the data#####
current<-"D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Current.xlsx" ##real time scenario
hypo1<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Hypo_41to11.xlsx" ##Hypothetical 1
hypo2<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Hypo_41to12.xlsx" ##Hypothetical 2
hypo3<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Hypo_41to26.xlsx" ##Hypothetical 3
sens_current<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Sensitivity_current.xlsx" ##Real-time sensitivity scores
sens_h1<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Sensitivity_Hypo_41to11.xlsx" ## Hypothetical 1 sensitivity scores
sens_h2<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Sensitivity_Hypo_41to12.xlsx" ## Hypothetical 2 sensitivity scores
sens_h3<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Sensitivity_Hypo_41to26.xlsx" ## Hypothetical 3 sensitivity scores
merged_data<- "D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Habitat_Frequency_Table.csv" ##Merged model results
all<-"D:/3rd semester/TUD job applicatons/UFZ/Agri birds/Model/Shapes/analysis/New results/New/Merged_Habitat_Frequency_and_Sensitivity_Data.csv" ##Merged sensitivity scores
######Data tables#####
data1 <- read_excel(current) ###current frequency
data2 <- read_excel(hypo1) ###hypothetical frequency 41-11
data3 <- read_excel(hypo2) ###hypothetical frequency 41-12
data4 <- read_excel(hypo3) ###hypothetical frequency 41-26
data5 <- read_excel(sens_current)###sensitivity current
data6 <- read_excel(sens_h1)###sensitivity 41-11
data7 <- read_excel(sens_h2)###sensitivity 41-12
data8 <- read_excel(sens_h3)###sensitivity 41-26
data9<- read.csv(merged_data)
all<-read.csv(all)
view(data9)


#####Group by Habitat Class and calculate the frequency for each#####
prepared_data1 <- data1 %>%
  group_by(`Habitat Class`, `Bird ID`) %>% 
  summarise(Frequency = n(), .groups = "drop")
prepared_data1
View(prepared_data1)
prepared_data2 <- data2 %>%
  group_by(`Habitat Class`, `Bird ID`) %>% 
  summarise(Frequency = n(), .groups = "drop")
prepared_data2
View(prepared_data2)
prepared_data3 <- data3 %>%
  group_by(`Habitat Class`, `Bird ID`) %>% 
  summarise(Frequency = n(), .groups = "drop")
prepared_data3
prepared_data4 <- data4 %>%
  group_by(`Habitat Class`, `Bird ID`) %>% 
  summarise(Frequency = n(), .groups = "drop")
prepared_data4
##Set manual color codes for each habitat type
manual_color<-c(
  "0" = "#1F78B4",
  "1" = "#B2DF8A",
  "2" = "#FDBF6F",
  "3" = "#FF7F00",
  "4" = "#CAB2D6",
  "5" = "#6A3D9A",
  "6" = "#FFFF99",
  "8" = "#A6CEE3",
  "9" = "#33A02C",
  "10" = "#FB9A99",
  "11" = "#1B9E77",
  "12" = "#33A02C",
  "16" = "#377EB8",
  "17" = "#4DAF4A",
  "18" = "#E41A1C",
  "19" = "#984EA3",
  "20" = "#A65628",
  "22" = "#999999",
  "23" = "#66C2A5",
  "24" = "#D95F02",
  "25" = "#7570B3",
  "26" = "#E7298A",
  "27" = "#B15928",
  "28" = "#FFFFFF",
  "29" = "#FFFF33",
  "30" = "#A6D854",
  "31" = "#FFD92F",
  "32" = "#E5C494",
  "33" = "#B3B3B3",
  "34" = "#8DD3C7",
  "35" = "#FFFFB3",
  "36" = "#BEBADA",
  "37" = "#FB8072",
  "38" = "#80B1D3",
  "39" = "#FDB462",
  "40" = "#B3DE69",
  "41" = "#FCCDE5",
  "42" = "#D9D9D9",
  "43" = "#BC80BD",
  "45" = "#CCEBC5",
  "46" = "#FFED6F",
  "47" = "#D95F02",
  "61" = "#E6AB02",
  "62" = "#A6761D",
  "251" = "#666666",
  "252" = "#999999",
  "253" = "#CCCCCC",
  "321" = "#8DD3C7",
  "322" = "#FFFFB3",
  "323" = "#BEBADA",
  "331" = "#FB8072",
  "332" = "#80B1D3",
  "333" = "#FDB462")
##Set manual labels for habitat codes
habitat_labels <- c("0" = "Ocean", "1" = "Agricultural grass","2"="Maize","3"="Potatoes",
                    "4"="Beets","5"="Grains","6"="Other crops","8"="Greenhouse (horti)",
                    "9"="Orchards","10"="Flower bulbs","11"="Deciduous forest","12"="Coniferous forest",
                    "16"="Fresh water","17"="Salt water","18"="Construction PB","19"="Construction SB",
                    "20"="Forest PB","22"="Forest SB","23"="Grass PB","24"="Bare ground",
                    "25"="Roads and railways","26"="Construction in RA","27"="Other LU in RA","28"="Grass SB",
                    "29"="Solar parks","30"="Salt marshes","31"="Open sand coastal","32"="Dunes LV",
                    "33"="Dunes HV","34"="Dune heath","35"="River sand","36"="Heather", "37"="Grassed (M) heathland",
                    "38"="Grassed (H) heathland","39"="Raised bog","40"="Forest RBA","41"="Other Swamp vegetation",
                    "42"="Reed vegetation","43"="Swamp forest","45"="Natural agri grassland",
                    "46"="Coastal grass","47"="Other grass","61"="Tree nurseries","62"="Fruit farms",
                    "251"="Infrastructure","252"="Semi-paved roads","253"="Narrow roads",
                    "321"="Shrub vegetation RBA(L)","322"="Shrub in swamp(L)",
                    "323"="Other shrub veg(L)","331"="Shrub vegetation RBA(H)",
                    "332"="Shrub in swamp(H)","333"="Other shrub veg(H)")
###PB= Primary-built up area, SB= Secondary built-up area, LV=Low vegetation
###RA= Rural area, LU=Land use, HV= High vegetation,Grassed (H/M)=Heavily/Moderately
###RBA=Raised bog area, (L/H)=Low/High
#####Figure14:  Habitat usage frequencies across different scenarios#####
##Prepare data in long format 
freq_long <- melt(all,
                  id.vars = "Habitat_type",
                  measure.vars = c("Current", "Hypo_41to11", "Hypo_41to12", "Hypo_41to26"),
                  variable.name = "Scenario",
                  value.name = "Frequency")
view(freq_long)
freq_long <- left_join(freq_long, all[, c("Habitat_type", "Sensitivity_Current")], by = "Habitat_type")

## Ordering the bars according to the current frequency of use in descending format
habitat_order <- freq_long %>%
  filter(Scenario == "Current") %>%
  arrange(desc(Frequency)) %>%  
  pull(Habitat_type)

## Applying the order to the dataset
freq_long <- freq_long %>%
  mutate(Habitat_type = factor(Habitat_type, levels = habitat_order))
scenario_colors <- c(
  "Current" = "#4E79A7",      
  "Hypo_41to11" = "#F28E2B", 
  "Hypo_41to12" = "#59A14F",  
  "Hypo_41to26" = "#E15759"   
)
##Manual label for scenarios
scenario_labels<- c("Current"="Real time scenario",
                    "Hypo_41to11"="Hypothetical: 33% expansion\nof deciduous forests",
                    "Hypo_41to12"="Hypothetical: 33% expansion\nof coniferous forests",
                    "Hypo_41to26"="Hypothetical: 39% expansion\nof construction in rural areas")
freq_long
##Plot1
Frequency<-ggplot(freq_long, aes(x = factor(Habitat_type), y = log1p(Frequency), fill = Scenario)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Habitat Usage Frequencies Across Scenarios",
       x = "Habitat Class", y = "Habitat Usage Frequency (log1p)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values = scenario_colors,labels=scenario_labels)+
  scale_x_discrete(labels=habitat_labels)+ theme(
    legend.position = c(0.91, 0.89),  # Adjust position as needed
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
Frequency
##save in png format with 300dpi resolution
ggsave("habitat usage_ordered_by_freq.png", width = 12, height = 6, dpi = 300)


#####Chi-sq for habitat frequency####
# Removing habitat class column and convert to matrix
freq_matrix <- as.matrix(data9[,-1])
rownames(freq_matrix) <- data9$Habitat_type
chi_test <- chisq.test(freq_matrix)
chi_test
##Error message due to lower number (<5) of data points in multiple classes
##Hence, applying simulated p-value with monte carlo randomization factor 10000
chi_result<-chisq.test(freq_matrix, simulate.p.value = TRUE, B = 10000)
chi_result

# Convert standardized residuals into a data frame
resid_df <- as.data.frame(as.table(chi_result$stdres))

#####Figure15:Chi-square residuals showing which habitats had affected the results most#####
residual<-ggplot(resid_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Chi-square Test: Standardized Residuals",
       x = "Scenario", y = "Habitat Type", fill = "Residual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1,size = 14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16))+
  scale_y_discrete(labels=habitat_labels)+
  scale_x_discrete(labels=scenario_labels)
residual
ggsave("chi_square_heatmap_corrected.png", width = 12, height = 11, dpi = 300)


#####Sensitivity score analysis#####
######Wilcoxon rank-sum#######
# Add a column to each dataframe indicating the scenario
data5$Scenario <- "Current"
data6$Scenario <- "Hypo_41to11"
data7$Scenario <- "Hypo_41to12"
data8$Scenario <- "Hypo_41to26"
combined_df <- rbind(data5, data6, data7, data8)
combined_df
view(combined_df)
wilcox.test(`Sensitivity Score` ~ Scenario, data = subset(combined_df, Scenario %in% c("Current", "Hypo_41to11")))
wilcox.test(`Sensitivity Score`~ Scenario, data = subset(combined_df, Scenario %in% c("Current", "Hypo_41to12")))
wilcox.test(`Sensitivity Score` ~ Scenario, data = subset(combined_df, Scenario %in% c("Current", "Hypo_41to26")))

##the comparison list to use for visualizing tthe p-values of targeted pairs
comparisons <- list(
  c("Current", "Hypo_41to11"),
  c("Current", "Hypo_41to12"),
  c("Current", "Hypo_41to26")
)
#####Figure16:Comparison of sensitivity scores over different scenarios#####
ggplot(combined_df, aes(x = Scenario, y = log(`Sensitivity Score`), fill = Scenario)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = comparisons,
                     label = "p.signif", size = 5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Sensitivity Scores",
    y = "Sensitivity Score (log values)",
    x = "Scenario"
  ) +
  theme(legend.position = "none")+
  scale_x_discrete(labels=scenario_labels)+
  scale_fill_manual(values = scenario_colors)
ggsave("wilcox violin_corrected.png", width = 12, height = 6, dpi = 300)
#####Figure17:Heatmap showing the changes in sensitivity scores over different scenarios#####
combined_df$`Habitat Class` <- as.factor(combined_df$`Habitat Class`)
combined_df$label <- round(combined_df$`Sensitivity Score`, 2)
ggplot(combined_df, aes(x = Scenario, y = `Habitat Class`, fill = `Sensitivity Score`)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3, color = "black") +
  scale_fill_gradient2(low = "#003399",mid="#e6ffff",high = "red", midpoint = 1.5) +
  labs(title = "Sensitivity Score Heatmap by Habitat and Scenario",
       x = "Scenario", y = "Habitat Class") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  scale_y_discrete(labels=habitat_labels)+scale_x_discrete(labels=scenario_labels)

ggsave("sensitivity heatmap_corrected.png", width = 12, height = 9, dpi = 300)


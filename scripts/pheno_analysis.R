## Phenotype Analysis ##
## Yvonne DEC 2022 ##

## Gets rid of previously stored variables ##
rm(list=ls())


############################
#### Load the libraries ####
############################

library(openxlsx)
library(dplyr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(vcd)
library(tidyr)



####################
#### Get set up ####
####################

## Set the directory ##
## Should be where you have the input data stored and where you want the output data ##

file.dir <- paste0("/Users/carolinarios/Documents/ELA/ELAPCB126.SCO")
setwd(file.dir)


## Label this with how you want to identify the data - typically Chemical and Pop ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
project <- "PCB126.SCO"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


## Make vector of ELA numbers that we don't want to include in analysis based on high background mortality, etc. ## ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#high_bkgd_elas <- "13012"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


## Set an output prefix with the date so we can track when output files were generated ##
output.prefix <- toupper(format(Sys.Date(),"%Y%b%d"))

## Define the output file directories ##
out.dir <- paste0(file.dir, "/Phenotype_Output_Final_", output.prefix)
stats.dir <- paste0(out.dir, "/stats/")
data.dir <- paste0(out.dir, "/data/")
plots.dir <- paste0(out.dir, "/plots/")

## Create the output file directories ##
# out.dir must be made manually #
dir.create(file.path(out.dir))
dir.create(stats.dir)
dir.create(paste0(stats.dir, "Perct Incidence/"))
dir.create(paste0(stats.dir, "Pheno Rating/"))
dir.create(data.dir)
dir.create(plots.dir)



##########################
#### Read in the data ####
##########################

data <- read.csv("Phenotype_data_all.csv")



###########################
#### Annotate the data ####
###########################

## Annotate the file with Concentration metadata ##
data <- data %>%
  mutate(NUM =
           case_when((EMB < 200) ~ "100",
                     (EMB > 200 & EMB < 300) ~ "200",
                     (EMB > 300 & EMB < 400) ~ "300",
                     (EMB > 400 & EMB < 500) ~ "400",
                     (EMB > 500 & EMB < 600) ~ "500",
                     (EMB > 600 & EMB < 700) ~ "600",
                     (EMB > 700 & EMB < 800) ~ "700",
                     (EMB > 800 & EMB < 900) ~ "800",)) %>%
  relocate(NUM, .after = EMB)


## Read in the map ##
map <- read.xlsx(paste0(file.dir, "/ELA_Map.xlsx"))
map$Group <- as.character(map$Group)
map$NUM <- as.character(map$NUM)
map$Chemical

data_mapped <- merge(map, data, by = c("Group", "NUM", "Chemical"), all.x = TRUE)

unique(data_mapped$NUM)
unique(data_mapped$Chemical)
unique(data_mapped$Group)


## Look at Standard.phenotype.score column ##
#############################################

## Remove rows where there is no pheno data ##
## Manual checking showed standard.phenotype.score is a good metric for whether fish were evaluated ##
data_filtered <- data_mapped %>%
  filter(Standard.phenotype.score != "") %>%
  filter(Standard.phenotype.score != " ")

## check number of rows removed ##
nrow(data_mapped)
nrow(data_filtered)

## Check for unique entries in case fish were not scored but marked as dead, NS, etc. ##
unique(data_filtered$Standard.phenotype.score)

## Make sure that any non-numbers (not 0, 1, 2, 3, 4, etc.) in the 'unique' line above is included in the below code ##
## We don't want to include these fish in pheno analysis because they were not evaluated for pheno ##
data_filtered2 <- data_filtered %>%
  filter(Standard.phenotype.score != "-") %>%
  filter(Standard.phenotype.score != "dead") %>%
  filter(Standard.phenotype.score != "DEAD") %>%
  filter(Standard.phenotype.score != "SACD") %>%
  filter(Standard.phenotype.score != "sac'd") %>%
  filter(Standard.phenotype.score != "SAC") %>%
  filter(Standard.phenotype.score != "MASSIVE MISDEVO") %>%
  filter(Standard.phenotype.score != "MASSICE MISDEVO") %>%
  filter(Standard.phenotype.score != "NS MASSIVE MISDEVO") %>%
  filter(Standard.phenotype.score != "misdevo") %>%
  filter(Standard.phenotype.score != "misdev") %>%
  filter(Standard.phenotype.score != "POPPED") %>%
  filter(Standard.phenotype.score != "NS") %>%
  filter(Standard.phenotype.score != "Na") %>%
  filter(Standard.phenotype.score != "/") %>%
  filter(Standard.phenotype.score != "N/A")

## Check number of rows removed ##
nrow(data_filtered2)

## Make sure we only have standard numbers (0, 1, 2, 3, 4) ##
unique(data_filtered2$Standard.phenotype.score)

unique(data_filtered2$NUM)               
unique(data_filtered2$Conc_uM)  


## Look at Heart.specific.score and other columns ##
####################################################

## Check whether there are weird, non-numeric heart-specific scores ##
unique(data_filtered2$Heart.specific.score)

## Make empty cells ("", aka NA after converted to numeric) into zeros ##
## Given the previous checks, we can consider all remaining empty cells as actual data ##
## First make a vector with all the endpoint column names (copy and paste from colnames output) ##
colnames(data_filtered2)
endpoints <- c('Pericard..edema', 'Severe.pericard..edema',
               'Heart.Abnormality..', 'Tail.hemorrhage',
               'Head.hem.', 'Body.head.size',
               'Other', 'Standard.phenotype.score',
               'Elongated.SV', 'Offset.Chambers',
               'Other.1', 'Heart.specific.score')

data_filtered3 <- data_filtered2 %>%
  mutate_at(c(endpoints), as.numeric) %>%
  mutate_at(c(endpoints), ~ tidyr::replace_na(.,0))

unique(data_filtered3$Heart.specific.score)

## Manually look at the comments to see if we need to delete anything ##
comments <- data_filtered3 %>%
  filter(Comments!= "")

# ## Remove fish as needed... ##
# fish_to_remove <- data_filtered3 %>%
#   filter((Group == XXX & EMBRYO.. == XXX) | (Group == XXX & EMBRYO.. == XXX))
# nrow(fish_to_remove)
# 
# write.csv(fish_to_remove, paste0(data.dir, project, "_Fish_removed_after_manual_QC_of_comments_", output.prefix, ".csv"), row.names=F)
# 
# nrow(data_filtered3)
# data_filtered3 <- anti_join(data_filtered3, fish_to_remove)
# nrow(data_filtered3)


#############################################
## Get counts of phenotype per dose level  ##
#############################################

#Subset data for just SCO PCB126 

data_filtered_SCO.PCB126 <- data_filtered3 %>%
  filter(Chemical == "PCB126" & POP.x == "SCO")

data_ctl <- subset(data_filtered_SCO.PCB126, Conc_nM == 0)
data_0.006nM <- subset(data_filtered_SCO.PCB126, Conc_nM == 0.006)
data_0.06nM <- subset(data_filtered_SCO.PCB126, Conc_nM == 0.06)
data_0.6nM <- subset(data_filtered_SCO.PCB126, Conc_nM == 0.6)
data_6nM <- subset(data_filtered_SCO.PCB126, Conc_nM == 6)
data_60nM <- subset(data_filtered_SCO.PCB126, Conc_nM == 60)
data_600nM <- subset(data_filtered_SCO.PCB126, Conc_nM == 600)

## Get percent counts/n by concentration for each standard phenotype

PE <- c(sum(data_ctl$Pericard..edema)/length(data_ctl$Pericard..edema),
            sum(data_0.006nM$Pericard..edema)/length(data_0.006nM$Pericard..edema), 
            sum(data_0.06nM$Pericard..edema)/length(data_0.06nM$Pericard..edema), 
            sum(data_6nM$Pericard..edema)/length(data_6nM$Pericard..edema), 
            sum(data_60nM$Pericard..edema)/length(data_60nM$Pericard..edema), 
            sum(data_600nM$Pericard..edema)/length(data_600nM$Pericard..edema))

Sev_PE <- c(sum(data_ctl$Severe.pericard..edema)/length(data_ctl$Severe.pericard..edema),
        sum(data_0.006nM$Severe.pericard..edema)/length(data_0.006nM$Severe.pericard..edema), 
        sum(data_0.06nM$Severe.pericard..edema)/length(data_0.06nM$Severe.pericard..edema), 
        sum(data_6nM$Severe.pericard..edema)/length(data_6nM$Severe.pericard..edema), 
        sum(data_60nM$Severe.pericard..edema)/length(data_60nM$Severe.pericard..edema), 
        sum(data_600nM$Severe.pericard..edema)/length(data_600nM$Severe.pericard..edema))

H_Abnorm <- c(sum(data_ctl$Heart.Abnormality..)/length(data_ctl$Heart.Abnormality..),
        sum(data_0.006nM$Heart.Abnormality..)/length(data_0.006nM$Heart.Abnormality..), 
        sum(data_0.06nM$Heart.Abnormality..)/length(data_0.06nM$Heart.Abnormality..), 
        sum(data_6nM$Heart.Abnormality..)/length(data_6nM$Heart.Abnormality..), 
        sum(data_60nM$Heart.Abnormality..)/length(data_60nM$Heart.Abnormality..), 
        sum(data_600nM$Heart.Abnormality..)/length(data_600nM$Heart.Abnormality..))

T_Hem <- c(sum(data_ctl$Tail.hemorrhage)/length(data_ctl$Tail.hemorrhage),
              sum(data_0.006nM$Tail.hemorrhage)/length(data_0.006nM$Tail.hemorrhage), 
              sum(data_0.06nM$Tail.hemorrhage)/length(data_0.06nM$Tail.hemorrhage), 
              sum(data_6nM$Tail.hemorrhage)/length(data_6nM$Tail.hemorrhage), 
              sum(data_60nM$Tail.hemorrhage)/length(data_60nM$Tail.hemorrhage), 
              sum(data_600nM$Tail.hemorrhage)/length(data_600nM$Tail.hemorrhage))

H_hem <- c(sum(data_ctl$Head.hem.)/length(data_ctl$ead.hem.),
           sum(data_0.006nM$ead.hem.)/length(data_0.006nM$Head.hem.), 
           sum(data_0.06nM$Head.hem.)/length(data_0.06nM$Head.hem.), 
           sum(data_6nM$Head.hem.)/length(data_6nM$Head.hem.), 
           sum(data_60nM$Head.hem.)/length(data_60nM$Head.hem.), 
           sum(data_600nM$Head.hem.)/length(data_600nM$Head.hem.))

Size <- c(sum(data_ctl$Body.head.size)/length(data_ctl$Body.head.size),
           sum(data_0.006nM$Body.head.size)/length(data_0.006nM$Body.head.size), 
           sum(data_0.06nM$Body.head.size)/length(data_0.06nM$Body.head.size), 
           sum(data_6nM$Body.head.size)/length(data_6nM$Body.head.size), 
           sum(data_60nM$Body.head.size)/length(data_60nM$Body.head.size), 
           sum(data_600nM$Body.head.size)/length(data_600nM$Body.head.size))

#Get average (counts/n) for heart specific scores

E_SV <- c(sum(data_ctl$Elongated.SV)/length(data_ctl$Elongated.SV),
          sum(data_0.006nM$Elongated.SV)/length(data_0.006nM$Elongated.SV), 
          sum(data_0.06nM$Elongated.SV)/length(data_0.06nM$Elongated.SV), 
          sum(data_6nM$Elongated.SV)/length(data_6nM$Elongated.SV), 
          sum(data_60nM$Elongated.SV)/length(data_60nM$Elongated.SV), 
          sum(data_600nM$Elongated.SV)/length(data_600nM$Elongated.SV))

Off_Ch <- c(sum(data_ctl$Offset.Chambers)/length(data_ctl$Offset.Chambers),
          sum(data_0.006nM$Offset.Chambers)/length(data_0.006nM$Offset.Chambers), 
          sum(data_0.06nM$Offset.Chambers)/length(data_0.06nM$Offset.Chambers), 
          sum(data_6nM$Offset.Chambers)/length(data_6nM$Offset.Chambers), 
          sum(data_60nM$Offset.Chambers)/length(data_60nM$Offset.Chambers), 
          sum(data_600nM$Offset.Chambers)/length(data_600nM$Offset.Chambers))

Other <- c(sum(data_ctl$Other)/length(data_ctl$Other),
            sum(data_0.006nM$Other)/length(data_0.006nM$Other), 
            sum(data_0.06nM$Other)/length(data_0.06nM$Other), 
            sum(data_6nM$Other)/length(data_6nM$Other), 
            sum(data_60nM$Other)/length(data_60nM$Other), 
            sum(data_600nM$Other)/length(data_600nM$Other))

library(ggplot2)
#data_filtered_trans_SCO.PCB126 <- data_filtered_SCO.PCB126
#data_filtered_trans_SCO.PCB126$POP.x <- as.factor(data_filtered_trans_SCO.PCB126$POP.x)
#data_filtered_trans_SCO.PCB126$Conc_nM <- as.factor(data_filtered_trans_SCO.PCB126$Conc_nM)



#data_filtered_trans_PCB126.3 <- data_filtered_trans_PCB126 %>% 
  #group_by(Conc_nM, POP.x) %>% 
  #summarize_if(is.numeric, list(mean = mean, sd = sd))

#data_filtered_trans_PCB126.2 <- data_filtered_trans_PCB126 %>% 
  #group_by(Conc_nM, POP.x) %>% 
  #summarise(PE = sum(Pericard..edema), n = n())

#ggplot(data_filtered_trans_PCB126,aes(x=Conc_nM, fill=POP.x)) + 
  #geom_boxplot(aes(y=Pericard..edema),color = 'blue',position = position_identity()) +
  #expand_limits(y=0) # This expands the y-axis to include 0 so we can better evaluate the importance of the effect!

#ggplot(data_filtered_trans_PCB126, aes(fill=PE, y= PE, x=Conc_nM) +
  #geom_bar(stat="identity", color="black", position=position_dodge())+
  #geom_text(aes(label=PE), color="black", 
            #size=3.5, position = position_dodge(width = 1), 
            #vjust = -.25)+
  #theme_bw()+
  #labs(x = "Population and Concentration",
       #y = "Incidences of Pericardial Edema",
       #)
#)



# Create bar plot for individual traits

## Create barplots of PE incidence by population in response to PCB exposure  
#PE <- ggplot(data_filtered_trans_PCB126, aes(x = Conc_nM, fill = factor(Pericard..edema))) +
  #geom_bar(position = "stack", stat = "count") +
  #labs(title = "Pericardial Edema in Response to PCB126",
       #x = "PCB126 (nM)",
       #y = "Count") +
  #scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  #facet_grid(~POP.x) +
  #theme_minimal()
#PE
#ggsave("PE_PCB126.png", plot = PE, width = 12, height = 6)

#Sev_PE <- ggplot(data_filtered_trans_PCB126, aes(x = Conc_nM, fill = factor(Severe.pericard..edema))) +
  #geom_bar(position = "stack", stat = "count") +
  #labs(title = "Severe Pericardial Edema in Response to PCB126",
   #    x = "PCB126 (nM)",
    #   y = "Count") +
  #scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
  #facet_grid(~POP.x) +
  #theme_minimal()
#Sev_PE
#ggsave("Sev_PE_PCB126.png", plot = PE, width = 12, height = 6)

#H_Abnorm <- ggplot(data_filtered_trans_PCB126, aes(x = Conc_nM, fill = factor(Heart.Abnormality..))) +
 # geom_bar(position = "stack", stat = "count") +
  #labs(title = "Heart Abnormality in Response to PCB126",
   #    x = "PCB126 (nM)",
    #   y = "Count") +
#  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
#  facet_grid(~POP.x) +
#  theme_minimal()
#H_Abnorm
#ggsave("H_Abnorm_PCB126.png", plot = PE, width = 12, height = 6)

#T_Hem <- ggplot(data_filtered_trans_PCB126, aes(x = Conc_nM, fill = factor(Tail.hemorrhage))) +
#  geom_bar(position = "stack", stat = "count") +
#  labs(title = "Tail Hemorrhage in Response to PCB126",
#       x = "PCB126 (nM)",
#      y = "Count") +
#  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
#  facet_grid(~POP.x) +
#  theme_minimal()
#T_Hem
#ggsave("T_Hem_PCB126.png", plot = PE, width = 12, height = 6)

#T_Hem <- ggplot(data_filtered_trans_PCB126, aes(x = Conc_nM, fill = factor(Tail.hemorrhage))) +
#  geom_bar(position = "stack", stat = "count") +
#  labs(title = "Tail Hemorrhage in Response to PCB126",
#       x = "PCB126 (nM)",
#       y = "Count") +
#  scale_fill_manual(values = c("0" = "blue", "1" = "red")) +
#  facet_grid(~POP.x) +
#  theme_minimal()
#T_Hem
#ggsave("T_Hem_PCB126.png", plot = PE, width = 12, height = 6)




#################################################
## Generate Barplots for PCB traits and scores ##
#################################################

library(openxlsx)
library(dplyr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(vcd)
library(tidyr)
library(forcats)
library(RColorBrewer)

data_filtered_PCB126 <- data_filtered3 %>%
  filter(Chemical == "PCB126")

data_filtered_trans_PCB126 <- data_filtered_PCB126

data_filtered_trans_PCB126$POP.x <- as.factor(data_filtered_trans_PCB126$POP.x)

data_filtered_trans_PCB126$Conc_nM <- as.factor(data_filtered_trans_PCB126$Conc_nM)

data_filtered_trans_PCB126$EMB <- as.factor(data_filtered_trans_PCB126$EMB)

data_filtered_trans_PCB126$Plate <- as.factor(data_filtered_trans_PCB126$Plate)

data_filtered_trans_PCB126$Conc_uM <- as.factor(data_filtered_trans_PCB126$Conc_uM)

data_filtered_trans_PCB126_long <- data_filtered_trans_PCB126 %>%
  gather(key = "trait", value = "value", -Group, -NUM, -Standard.phenotype.score, -Chemical, -Conc_nM,
         -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Comments, -POP.y, -Heart.specific.score, -Archived, -Heart.vid, -Initials)

# Reorder factor levels for traits
data_filtered_trans_PCB126_long$trait <- fct_reorder(data_filtered_trans_PCB126_long$trait, data_filtered_trans_PCB126_long$value)

# Ensure proper ordering of factor levels for 'population'
data_filtered_trans_PCB126_long$POP.x <- factor(data_filtered_trans_PCB126_long$POP.x, 
                                                levels = unique(data_filtered_trans_PCB126_long$POP.x))

color_palette_cnt <- brewer.pal(n = nlevels(factor(data_filtered_trans_PCB126_long$value)), name = "Blues")

ALL_PCB126 <- ggplot(data_filtered_trans_PCB126_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Traits Across PCB126 Concentrations",
       x = "[PCB126] nM",
       y = "Count") +
  facet_wrap(~POP.x + trait, scales = "free_y", ncol=10) +
  scale_fill_manual(values = color_palette_cnt)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
ALL_PCB126
ggsave("ALL_PCB126.png", plot = ALL_PCB126, width = 24, height = 8)

## Create boxplot for PCB Scores

data_filtered_trans_PCB126_score_long <- data_filtered_trans_PCB126 %>%
  gather(key = "score", value = "value", -Group, -NUM, -Chemical, -Conc_nM, -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Pericard..edema, -Severe.pericard..edema, -Heart.Abnormality.., -Tail.hemorrhage, -Head.hem., -Body.head.size,
         -Other, -Elongated.SV, -Offset.Chambers, -Other.1, -Comments, -POP.y, -Archived, -Heart.vid, -Initials)

# Reorder factor levels for traits
data_filtered_trans_PCB126_score_long$score <- fct_reorder(data_filtered_trans_PCB126_score_long$score, data_filtered_trans_PCB126_score_long$value)

# Ensure proper ordering of factor levels for 'population'
data_filtered_trans_PCB126_score_long$POP.x <- factor(data_filtered_trans_PCB126_score_long$POP.x, 
                                                levels = unique(data_filtered_trans_PCB126_score_long$POP.x))

color_palette_scores <- brewer.pal(n = nlevels(factor(data_filtered_trans_PCB126_score_long$value)), name = "YlOrRd")

Score_PCB126 <- ggplot(data_filtered_trans_PCB126_score_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Scores Across PCB126 Concentrations",
       x = "[PCB126] nM",
       y = "Count") +
  facet_wrap(~POP.x + score, scales = "free_y", ncol=2) +
  scale_fill_manual(values = color_palette_scores)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
Score_PCB126
ggsave("Score_PCB126.png", plot = Score_PCB126, width = 18, height = 8)

#############################
## Create bar plot for BkF ##
#############################

data_filtered_BkF <- data_filtered3 %>%
  filter(Chemical == "BkF")

data_filtered_trans_BkF <- data_filtered_BkF

data_filtered_trans_BkF$POP.x <- as.factor(data_filtered_trans_BkF$POP.x)

data_filtered_trans_BkF$Conc_nM <- as.factor(data_filtered_trans_BkF$Conc_nM)

data_filtered_trans_BkF$EMB <- as.factor(data_filtered_trans_BkF$EMB)

data_filtered_trans_BkF$Plate <- as.factor(data_filtered_trans_BkF$Plate)

data_filtered_trans_BkF$Conc_uM <- as.factor(data_filtered_trans_BkF$Conc_uM)

data_filtered_trans_BkF_long <- data_filtered_trans_BkF %>%
  gather(key = "trait", value = "value", -Group, -NUM, -Standard.phenotype.score, -Chemical, -Conc_nM,
         -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Comments, -POP.y, -Heart.specific.score, -Archived, -Heart.vid, -Initials)
# Reorder factor levels for traits
data_filtered_trans_BkF_long$trait <- fct_reorder(data_filtered_trans_BkF_long$trait, data_filtered_trans_BkF_long$value)

# Ensure proper ordering of factor levels for 'population'
data_filtered_trans_BkF_long$POP.x <- factor(data_filtered_trans_BkF_long$POP.x, 
                                                levels = unique(data_filtered_trans_BkF_long$POP.x))
# Generate color palette
color_palette <- brewer.pal(n = nlevels(factor(data_filtered_trans_BkF_long$value)), name = "Blues")

# Generate boxplot
ALL_BkF <- ggplot(data_filtered_trans_BkF_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Traits Across BkF Concentrations",
       x = "[BkF] uM",
       y = "Count") +
  facet_wrap(~POP.x + trait, scales = "free_y", ncol=10) +
  scale_fill_manual(values = color_palette_cnt)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
ALL_BkF
ggsave("ALL_BkF.png", plot = ALL_BkF, width = 24, height = 8)

## Create barplot for scores

data_filtered_trans_BkF_score_long <- data_filtered_trans_BkF %>%
  gather(key = "score", value = "value", -Group, -NUM, -Chemical, -Conc_nM, -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Pericard..edema, -Severe.pericard..edema, -Heart.Abnormality.., -Tail.hemorrhage, -Head.hem., -Body.head.size,
         -Other, -Elongated.SV, -Offset.Chambers, -Other.1, -Comments, -POP.y, -Archived, -Heart.vid, -Initials)
# Reorder factor levels for traits
data_filtered_trans_BkF_score_long$score <- fct_reorder(data_filtered_trans_BkF_score_long$score, data_filtered_trans_BkF_score_long$score)

# Ensure proper ordering of factor levels for 'population'
data_filtered_trans_BkF_score_long$POP.x <- factor(data_filtered_trans_BkF_score_long$POP.x, 
                                             levels = unique(data_filtered_trans_BkF_score_long$POP.x))

# Generate boxplot
Score_BkF <- ggplot(data_filtered_trans_BkF_score_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Scores Across BkF Concentrations",
       x = "[BkF] nM",
       y = "Score") +
  facet_wrap(~POP.x + score, scales = "free_y", ncol=2) +
  scale_fill_manual(values = color_palette_scores)  +
  theme_minimal()
Score_BkF
ggsave("Score_BkF.png", plot = Score_BkF, width = 18, height = 8)


#############################################
## Identify Unique Phenotypic Combinations ##
#############################################

# ID unique phenotype combinations in PCB exposed fish

data_filtered_trans_PCB126_concat <- data_filtered_trans_PCB126 %>%
  mutate(combination = paste(POP.x, Chemical, Conc_nM, Pericard..edema, Severe.pericard..edema,
         Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
         Elongated.SV, Offset.Chambers, Other.1, sep = "_"))

unique_combinations_PCB <- unique(data_filtered_trans_PCB126_concat$combination)

combination_counts_PCB <- data_filtered_trans_PCB126_concat %>%
  count(combination, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

print(unique_combinations_PCB)

print(combination_counts_PCB)

## Create bar plots of counts of unique combinations by population

# Unique combination of SCO PCB126 exposure

color_palette_conc <- brewer.pal(n = nlevels(factor(data_filtered_trans_PCB126$Conc_uM)), name = "YlOrRd")

data_filtered_trans_SCO.PCB126_concat <- data_filtered_trans_PCB126_concat %>%
  filter(POP.x == "SCO")

combination_counts_SCO.PCB <- data_filtered_trans_SCO.PCB126_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_SCO.PCB126 <- ggplot(combination_counts_SCO.PCB, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for PCB126 exposed SCO Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Unique_SCO.PCB126.png", plot = Unique_SCO.PCB126, width = 18, height = 8)

# Unique combination of NBH PCB126 exposure

data_filtered_trans_NBH.PCB126_concat <- data_filtered_trans_PCB126_concat %>%
  filter(POP.x == "NBH")

combination_counts_NBH.PCB <- data_filtered_trans_NBH.PCB126_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_NBH.PCB126 <- ggplot(combination_counts_NBH.PCB, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for PCB126 exposed NBH Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Unique_NBH.PCB126
ggsave("Unique_NBH.PCB126.png", plot = Unique_NBH.PCB126, width = 18, height = 8)

# Unique combination of KC PCB126 exposure

data_filtered_trans_KC.PCB126_concat <- data_filtered_trans_PCB126_concat %>%
  filter(POP.x == "KC")

combination_counts_KC.PCB <- data_filtered_trans_KC.PCB126_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_KC.PCB126 <- ggplot(combination_counts_KC.PCB, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for PCB126 exposed KC Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Unique_KC.PCB126
ggsave("Unique_KC.PCB126.png", plot = Unique_KC.PCB126, width = 18, height = 8)

# Unique combination of ER PCB126 exposure

data_filtered_trans_ER.PCB126_concat <- data_filtered_trans_PCB126_concat %>%
  filter(POP.x == "ER")

combination_counts_ER.PCB <- data_filtered_trans_ER.PCB126_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_ER.PCB126 <- ggplot(combination_counts_ER.PCB, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for PCB126 exposed KC Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Unique_ER.PCB126
ggsave("Unique_ER.PCB126.png", plot = Unique_ER.PCB126, width = 18, height = 8)

# ID unique phenotype combinations in BkF exposed fish

data_filtered_trans_BkF_concat <- data_filtered_trans_BkF %>%
  mutate(combination = paste(POP.x, Chemical, Conc_nM, Pericard..edema, Severe.pericard..edema,
                             Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
                             Elongated.SV, Offset.Chambers, Other.1, sep = "_"))


## Create bar plots of counts of unique combinations by population

### Unique combination of SCO BkF exposure

color_palette_conc <- brewer.pal(n = nlevels(factor(data_filtered_trans_BkF$Conc_uM)), name = "YlOrRd")

data_filtered_trans_SCO.BkF_concat <- data_filtered_trans_BkF_concat %>%
  filter(POP.x == "SCO")

combination_counts_SCO.BkF <- data_filtered_trans_SCO.BkF_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_SCO.BkF <- ggplot(combination_counts_SCO.BkF, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for BkF exposed SCO Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Unique_SCO.BkF.png", plot = Unique_SCO.BkF, width = 18, height = 8)

### Unique combination of NBH BkF exposure

data_filtered_trans_NBH.BkF_concat <- data_filtered_trans_BkF_concat %>%
  filter(POP.x == "NBH")

combination_counts_NBH.BkF <- data_filtered_trans_NBH.BkF_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_NBH.BkF <- ggplot(combination_counts_NBH.BkF, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for BkF exposed NBH Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Unique_NBH.BkF
ggsave("Unique_NBH.BkF.png", plot = Unique_NBH.BkF, width = 18, height = 8)

### Unique combination of KC BkF exposure

data_filtered_trans_KC.BkF_concat <- data_filtered_trans_BkF_concat %>%
  filter(POP.x == "KC")

combination_counts_KC.BkF <- data_filtered_trans_KC.BkF_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_KC.BkF <- ggplot(combination_counts_KC.BkF, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for BkF exposed KC Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Unique_KC.BkF
ggsave("Unique_KC.BkF.png", plot = Unique_KC.BkF, width = 18, height = 8)

### Unique combination of ER BkF exposure

data_filtered_trans_ER.BkF_concat <- data_filtered_trans_BkF_concat %>%
  filter(POP.x == "ER")

combination_counts_ER.BkF <- data_filtered_trans_ER.BkF_concat %>%
  count(combination, Conc_nM, Pericard..edema, Severe.pericard..edema,
        Heart.Abnormality.., Tail.hemorrhage, Head.hem., Body.head.size, Other, 
        Elongated.SV, Offset.Chambers, Other.1,)

Unique_ER.BkF <- ggplot(combination_counts_ER.BkF, aes(x = combination, y = n, fill=factor(Conc_nM))) +
  geom_bar(stat = "identity") +
  labs(title = "Counts of Unique Combinations for BkF exposed KC Embryos",
       x = "Combination",
       y = "Count") +
  scale_fill_manual(values = color_palette_conc) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Unique_ER.BkF
ggsave("Unique_ER.BkF.png", plot = Unique_ER.BkF, width = 18, height = 8)

####################################################
## Random subsetting of samples based on group ID ##
####################################################

subset_data <- data_filtered3 %>% group_by(Group, NUM) %>% slice_sample(n=5)
nrow(subset_data)

##########################################################
## Assess subsetting using barplots by phenotype counts ##
##########################################################

subset_data_PCB126 <- subset_data %>%
  filter(Chemical == "PCB126")
write.csv(subset_data_PCB126, "/Users/carolinarios/Documents/ELA/subset_data_PCB126.csv")

subset_data_trans_PCB126 <- subset_data_PCB126

subset_data_trans_PCB126$POP.x <- as.factor(subset_data_trans_PCB126$POP.x)

subset_data_trans_PCB126$Conc_nM <- as.factor(subset_data_trans_PCB126$Conc_nM)

subset_data_trans_PCB126$EMB <- as.factor(subset_data_trans_PCB126$EMB)

subset_data_trans_PCB126$Plate <- as.factor(subset_data_trans_PCB126$Plate)

subset_data_trans_PCB126$Conc_uM <- as.factor(subset_data_trans_PCB126$Conc_uM)

subset_data_trans_PCB126_long <- subset_data_trans_PCB126 %>%
  gather(key = "trait", value = "value", -Group, -NUM, -Standard.phenotype.score, -Chemical, -Conc_nM,
         -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Comments, -POP.y, -Heart.specific.score, -Archived, -Heart.vid, -Initials)

# Reorder factor levels for traits
subset_data_trans_PCB126_long$trait <- fct_reorder(subset_data_trans_PCB126_long$trait, subset_data_trans_PCB126_long$value)

# Ensure proper ordering of factor levels for 'population'
subset_data_trans_PCB126_long$POP.x <- factor(subset_data_trans_PCB126_long$POP.x, 
                                                levels = unique(subset_data_trans_PCB126_long$POP.x))

color_palette_subcnt <- brewer.pal(n = nlevels(factor(subset_data_trans_PCB126_long$value)), name = "Greens")

ALL_PCB126_subset <- ggplot(subset_data_trans_PCB126_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Traits Across PCB126 Concentrations in Subset Individuals",
       x = "[PCB126] nM",
       y = "Count") +
  facet_wrap(~POP.x + trait, scales = "free_y", ncol=10) +
  scale_fill_manual(values = color_palette_subcnt)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
ALL_PCB126_subset
ggsave("ALL_PCB126_subset.png", plot = ALL_PCB126_subset, width = 24, height = 8)

## Create boxplot for PCB Scores

subset_data_trans_PCB126_score_long <- subset_data_trans_PCB126 %>%
  gather(key = "score", value = "value", -Group, -NUM, -Chemical, -Conc_nM, -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Pericard..edema, -Severe.pericard..edema, -Heart.Abnormality.., -Tail.hemorrhage, -Head.hem., -Body.head.size,
         -Other, -Elongated.SV, -Offset.Chambers, -Other.1, -Comments, -POP.y, -Archived, -Heart.vid, -Initials)

# Reorder factor levels for traits
subset_data_trans_PCB126_score_long$score <- fct_reorder(subset_data_trans_PCB126_score_long$score, subset_data_trans_PCB126_score_long$value)

# Ensure proper ordering of factor levels for 'population'
subset_data_trans_PCB126_score_long$POP.x <- factor(subset_data_trans_PCB126_score_long$POP.x, 
                                                      levels = unique(subset_data_trans_PCB126_score_long$POP.x))

color_palette_subsetscores <- brewer.pal(n = nlevels(factor(subset_data_trans_PCB126_score_long$value)), name = "YlGnBu")

Score_PCB126_subset <- ggplot(subset_data_trans_PCB126_score_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Scores Across PCB126 Concentrations in Subset Individuals",
       x = "[PCB126] nM",
       y = "Count") +
  facet_wrap(~POP.x + score, scales = "free_y", ncol=2) +
  scale_fill_manual(values = color_palette_scores)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
Score_PCB126_subset
ggsave("Score_PCB126_subset.png", plot = Score_PCB126_subset, width = 18, height = 8)

#############################
## Create bar plot for BkF ##
#############################

subset_data_BkF <- subset_data %>%
  filter(Chemical == "BkF")
write.csv(subset_data_BkF, "/Users/carolinarios/Documents/ELA/subset_data_BkF.csv")

subset_data_trans_BkF <- subset_data_BkF

subset_data_trans_BkF$POP.x <- as.factor(subset_data_trans_BkF$POP.x)

subset_data_trans_BkF$Conc_nM <- as.factor(subset_data_trans_BkF$Conc_nM)

subset_data_trans_BkF$EMB <- as.factor(subset_data_trans_BkF$EMB)

subset_data_trans_BkF$Plate <- as.factor(subset_data_trans_BkF$Plate)

subset_data_trans_BkF$Conc_uM <- as.factor(subset_data_trans_BkF$Conc_uM)

subset_data_trans_BkF_long <- subset_data_trans_BkF %>%
  gather(key = "trait", value = "value", -Group, -NUM, -Standard.phenotype.score, -Chemical, -Conc_nM,
         -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Comments, -POP.y, -Heart.specific.score, -Archived, -Heart.vid, -Initials)
# Reorder factor levels for traits
subset_data_trans_BkF_long$trait <- fct_reorder(subset_data_trans_BkF_long$trait, subset_data_trans_BkF_long$value)

# Ensure proper ordering of factor levels for 'population'
subset_data_trans_BkF_long$POP.x <- factor(subset_data_trans_BkF_long$POP.x, 
                                             levels = unique(subset_data_trans_BkF_long$POP.x))
# Generate color palette
color_palette <- brewer.pal(n = nlevels(factor(subset_data_trans_BkF_long$value)), name = "Greens")

# Generate boxplot
ALL_subset_BkF <- ggplot(subset_data_trans_BkF_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Traits Across BkF Concentrations in Subset Individuals",
       x = "[BkF] uM",
       y = "Count") +
  facet_wrap(~POP.x + trait, scales = "free_y", ncol=10) +
  scale_fill_manual(values = color_palette_cnt)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
ALL_subset_BkF
ggsave("ALL_subset_BkF.png", plot = ALL_subset_BkF, width = 24, height = 8)

## Create barplot for scores

subset_data_trans_BkF_score_long <- subset_data_trans_BkF %>%
  gather(key = "score", value = "value", -Group, -NUM, -Chemical, -Conc_nM, -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Pericard..edema, -Severe.pericard..edema, -Heart.Abnormality.., -Tail.hemorrhage, -Head.hem., -Body.head.size,
         -Other, -Elongated.SV, -Offset.Chambers, -Other.1, -Comments, -POP.y, -Archived, -Heart.vid, -Initials)
# Reorder factor levels for traits
subset_data_trans_BkF_score_long$score <- fct_reorder(subset_data_trans_BkF_score_long$score, subset_data_trans_BkF_score_long$value)

# Ensure proper ordering of factor levels for 'population'
subset_data_trans_BkF_score_long$POP.x <- factor(subset_data_trans_BkF_score_long$POP.x, 
                                                   levels = unique(subset_data_trans_BkF_score_long$POP.x))

# Generate boxplot
Score_subset_BkF <- ggplot(subset_data_trans_BkF_score_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Scores Across BkF Concentrations in Subset Individuals",
       x = "[BkF] nM",
       y = "Score") +
  facet_wrap(~POP.x + score, scales = "free_y", ncol=2) +
  scale_fill_manual(values = color_palette_scores)  +
  theme_minimal()
Score_subset_BkF
ggsave("Score_subset_BkF.png", plot = Score_subset_BkF, width = 18, height = 8)


###################
##resubset scores##
###################



################################################
##Boxplots (scores only) resubset individuals)## 
################################################

#Barplots for BkF

resubset_data_BkF <- read.csv("/Users/carolinarios/Documents/ELA/resubset_data_BkF.csv")

resubset_data_trans_BkF <- resubset_data_BkF

resubset_data_trans_BkF$POP.x <- as.factor(resubset_data_trans_BkF$POP.x)

resubset_data_trans_BkF$Conc_nM <- as.factor(resubset_data_trans_BkF$Conc_nM)

resubset_data_trans_BkF$EMB <- as.factor(resubset_data_trans_BkF$EMB)

resubset_data_trans_BkF$Plate <- as.factor(resubset_data_trans_BkF$Plate)

resubset_data_trans_BkF$Conc_uM <- as.factor(resubset_data_trans_BkF$Conc_uM)

resubset_data_trans_BkF_long <- resubset_data_trans_BkF %>%
  gather(key = "trait", value = "value", -Group, -NUM, -Standard.phenotype.score, -Chemical, -Conc_nM,
         -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Comments, -POP.y, -Heart.specific.score, -Archived, -Heart.vid, -Initials)
# Reorder factor levels for traits
resubset_data_trans_BkF_long$trait <- fct_reorder(resubset_data_trans_BkF_long$trait, resubset_data_trans_BkF_long$value)

# Ensure proper ordering of factor levels for 'population'
resubset_data_trans_BkF_long$POP.x <- factor(resubset_data_trans_BkF_long$POP.x, 
                                           levels = unique(resubset_data_trans_BkF_long$POP.x))
# Generate color palette
color_palette <- brewer.pal(n = nlevels(factor(resubset_data_trans_BkF_long$value)), name = "Greens")

## Create barplot for scores

resubset_data_trans_BkF_score_long <- resubset_data_trans_BkF %>%
  gather(key = "score", value = "value", -Group, -NUM, -Chemical, -Conc_nM, -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Pericard..edema, -Severe.pericard..edema, -Heart.Abnormality.., -Tail.hemorrhage, -Head.hem., -Body.head.size,
         -Other, -Elongated.SV, -Offset.Chambers, -Other.1, -Comments, -POP.y, -Archived, -Heart.vid, -Initials)
# Reorder factor levels for traits
resubset_data_trans_BkF_score_long$score <- fct_reorder(resubset_data_trans_BkF_score_long$score, resubset_data_trans_BkF_score_long$value)

# Ensure proper ordering of factor levels for 'population'
resubset_data_trans_BkF_score_long$POP.x <- factor(resubset_data_trans_BkF_score_long$POP.x, 
                                                 levels = unique(resubset_data_trans_BkF_score_long$POP.x))

# Generate boxplot
Score_resubset_BkF <- ggplot(resubset_data_trans_BkF_score_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Scores Across BkF Concentrations in Resubset Individuals",
       x = "[BkF] nM",
       y = "Score") +
  facet_wrap(~POP.x + score, scales = "free_y", ncol=2) +
  scale_fill_manual(values = color_palette_scores)  +
  theme_minimal()
Score_resubset_BkF
ggsave("Score_resubset_BkF.png", plot = Score_resubset_BkF, width = 18, height = 8)

#Barplot for PCB 126 resubset

resubset_data_PCB126 <- read.csv("/Users/carolinarios/Documents/ELA/resubset_data_PCB126.csv")

resubset_data_trans_PCB126 <- resubset_data_PCB126

resubset_data_trans_PCB126$POP.x <- as.factor(resubset_data_trans_PCB126$POP.x)

resubset_data_trans_PCB126$Conc_nM <- as.factor(resubset_data_trans_PCB126$Conc_nM)

resubset_data_trans_PCB126$EMB <- as.factor(resubset_data_trans_PCB126$EMB)

resubset_data_trans_PCB126$Plate <- as.factor(resubset_data_trans_PCB126$Plate)

resubset_data_trans_PCB126$Conc_uM <- as.factor(resubset_data_trans_PCB126$Conc_uM)

resubset_data_trans_PCB126_long <- resubset_data_trans_PCB126 %>%
  gather(key = "trait", value = "value", -Group, -NUM, -Standard.phenotype.score, -Chemical, -Conc_nM,
         -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Comments, -POP.y, -Heart.specific.score, -Archived, -Heart.vid, -Initials)

# Reorder factor levels for traits
resubset_data_trans_PCB126_long$trait <- fct_reorder(resubset_data_trans_PCB126_long$trait, resubset_data_trans_PCB126_long$value)

# Ensure proper ordering of factor levels for 'population'
resubset_data_trans_PCB126_long$POP.x <- factor(resubset_data_trans_PCB126_long$POP.x, 
                                              levels = unique(resubset_data_trans_PCB126_long$POP.x))

color_palette_subcnt <- brewer.pal(n = nlevels(factor(resubset_data_trans_PCB126_long$value)), name = "Greens")

ALL_PCB126_resubset <- ggplot(resubset_data_trans_PCB126_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Traits Across PCB126 Concentrations in Resubset Individuals",
       x = "[PCB126] nM",
       y = "Count") +
  facet_wrap(~POP.x + trait, scales = "free_y", ncol=10) +
  scale_fill_manual(values = color_palette_subcnt)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
ALL_PCB126_resubset
ggsave("ALL_PCB126_resubset.png", plot = ALL_PCB126_resubset, width = 24, height = 8)

## Create boxplot for PCB Scores

resubset_data_trans_PCB126_score_long <- resubset_data_trans_PCB126 %>%
  gather(key = "score", value = "value", -Group, -NUM, -Chemical, -Conc_nM, -Conc_uM, -POP.x, -POP_NUM, -Plate, -EMB,
         -Pericard..edema, -Severe.pericard..edema, -Heart.Abnormality.., -Tail.hemorrhage, -Head.hem., -Body.head.size,
         -Other, -Elongated.SV, -Offset.Chambers, -Other.1, -Comments, -POP.y, -Archived, -Heart.vid, -Initials)

# Reorder factor levels for traits
resubset_data_trans_PCB126_score_long$score <- fct_reorder(resubset_data_trans_PCB126_score_long$score, resubset_data_trans_PCB126_score_long$value)

# Ensure proper ordering of factor levels for 'population'
resubset_data_trans_PCB126_score_long$POP.x <- factor(resubset_data_trans_PCB126_score_long$POP.x, 
                                                    levels = unique(resubset_data_trans_PCB126_score_long$POP.x))

color_palette_subsetscores <- brewer.pal(n = nlevels(factor(resubset_data_trans_PCB126_score_long$value)), name = "YlGnBu")

Score_PCB126_resubset <- ggplot(resubset_data_trans_PCB126_score_long, aes(x = Conc_nM, fill = factor(value))) +
  geom_bar(position = "stack", stat = "count") +
  labs(title = "Distribution of Scores Across PCB126 Concentrations in Resubset Individuals",
       x = "[PCB126] nM",
       y = "Count") +
  facet_wrap(~POP.x + score, scales = "free_y", ncol=2) +
  scale_fill_manual(values = color_palette_scores)  +
  theme(strip.text.x = element_text(size =40)) +
  theme_minimal()
Score_PCB126_resubset
ggsave("Score_PCB126_resubset.png", plot = Score_PCB126_resubset, width = 18, height = 8)

subset_data_PCB126


# Assuming you have a data frame named 'your_data' with multiple traits and a grouping variable 'condition'
# Install and load the dplyr package if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# Example: Iterating over doses and populations
doses <- unique(subset_data$Conc_uM)
populations <- unique(subset_data$POP.x)
chemicals <- unique(subset_data$Chemical)

# Identify binary and three-level traits
binary_traits <- c("Pericard..edema", "Severe.pericard..edema", "Heart.Abnormality..", "Tail.hemorrhage", 
                   "Head.hem.", "Body.head.size ")
three_level_traits <- c("Elongated.SV", "Offset.Chambers")

for (dose in doses) {
  for (population in populations) {
    for (chemical in chemicals)
      cat("Dose:", dose, "Population:", population, "Chemical:", chemical, "\n")
    
      for (trait in c(binary_traits, three_level_traits)) {
        cat("  Trait:", trait, "\n")
      
      # Extract data before and after subsetting
        data_before <- subset(data_filtered3, Conc_uM == dose & POP.x == population, Chemical == chemical)
        data_after <- subset(subset_data, Conc_uM == dose & POP.x == population, Chemical == chemical)
      
      # Perform the appropriate statistical test
        if (trait %in% binary_traits) {
          test_result <- t.test(data_before[[trait]], data_after[[trait]], paired = TRUE)
        } else {
          test_result <- wilcox.test(data_before[[trait]], data_after[[trait]], paired = TRUE)
        }
      
      # Print the test result
        print(test_result)
        cat("\n")
    }
  }
}

# Perform paired t-tests for each trait
t_test_results <- lapply(phenotypes, function(phenotype) {
  t_test_result <- t.test(data_filtered3[[phenotype]], subset_data[[phenotype]], paired = TRUE)
  return(t_test_result)
})

# Print the results
for (i in seq_along(phenotype)) {
  cat("Phenotype:", phenotypes[i], "\n")
  print(t_test_results[[i]])
  cat("\n")
}

##########################################################
##random subset of PCB126 embryo numbers for extractions##
##########################################################

resubset_extractions_PCB126_round1 <- resubset_data_PCB126 %>% group_by(Group, NUM) %>% slice_sample(n=2)
nrow(resubset_extractions_PCB126_round1)

resubset_extractions_PCB126_remaining <- anti_join(resubset_data_PCB126, 
                                                   resubset_extractions_PCB126_round1)

resubset_extractions_PCB126_round2 <- resubset_extractions_PCB126_remaining %>% group_by(Group, NUM) %>% slice_sample(n=2)
nrow(resubset_extractions_PCB126_round2)

resubset_extractions_PCB126_round3 <- anti_join(resubset_extractions_PCB126_remaining, 
                                                resubset_extractions_PCB126_round2)

nrow(resubset_extractions_PCB126_round3)

resubset_extractions_PCB126_round1 <- resubset_extractions_PCB126_round1 %>%
  mutate(EMBID = Group+EMB)

resubset_extractions_PCB126_round2 <- resubset_extractions_PCB126_round2 %>%
  mutate(EMBID = Group+EMB)

resubset_extractions_PCB126_round3 <- resubset_extractions_PCB126_round3 %>%
  mutate(EMBID = Group+EMB)


resubset_extractions_BkF_round1 <- resubset_data_BkF %>% group_by(Group, NUM) %>% slice_sample(n=3)
nrow(resubset_extractions_BkF_day1)

resubset_extractions_BkF_round2<-anti_join(resubset_data_BkF, resubset_extractions_BkF_day1)
nrow(resubset_extractions_BkF_round2)

resubset_extractions_BkF_round1 <- resubset_extractions_BkF_round1 %>%
  mutate(EMBID = Group+EMB)
resubset_extractions_BkF_round2 <- resubset_extractions_BkF_round2 %>%
  mutate(EMBID = Group+EMB)


#look at ID3 depenbdent survival
#directories
base.dir<- "c010-datasets/Internal/2018-Ali/LiteratureSearch"
analysis.dir <- file.path(base.dir, "analysis")
data.dir <- file.path(base.dir, "data")

#libraries 
library(ggpubr)

#load data
survival <- read.table(file.path(data.dir, "Overall_Survival_metabric.txt"),sep="\t", header=TRUE)
survival$Sample.Id <- survival$Case.ID
expression <- read.table(file.path(data.dir, "ID3_expresssion_metabric.txt"),sep="\t", header=TRUE)
#get quantiles
lowCutoff <- median(expression$ID3..mRNA.expression..microarray. - sd(expression$ID3..mRNA.expression..microarray.))
highCutoff <- median(expression$ID3..mRNA.expression..microarray. + sd(expression$ID3..mRNA.expression..microarray.))

#first_qu <- quantile(expression$ID3..mRNA.expression..microarray.)[2]
#third_qu <- quantile(expression$ID3..mRNA.expression..microarray.)[4]

#subset
expression_rad_sub_fistqu <- expression[expression$ID3..mRNA.expression..microarray. <= lowCutoff,]
dim(expression_rad_sub_fistqu)
expression_rad_sub_thirdqu <- expression[expression$ID3..mRNA.expression..microarray. >= highCutoff,]
dim(expression_rad_sub_thirdqu)

#expression_rad_sub_fistqu <- expression[expression$ID3..mRNA.expression..microarray. <= first_qu,]
#dim(expression_rad_sub_fistqu)
#expression_rad_sub_thirdqu <- expression[expression$ID3..mRNA.expression..microarray. >= third_qu,]
#dim(expression_rad_sub_thirdqu)

#get survival data for those patients
library(dplyr)
expression_rad_sub_fistqu<- left_join(expression_rad_sub_fistqu, survival, by="Sample.Id")
summary(expression_rad_sub_fistqu$Survival.Rate)
expression_rad_sub_thirdqu<- left_join(expression_rad_sub_thirdqu, survival, by="Sample.Id")
summary(expression_rad_sub_thirdqu$Survival.Rate)

#plot caplan meier survival curves
library(survminer)
library(survival)

expression_rad_sub_thirdqu$ID3_expression <- "high"
expression_rad_sub_fistqu$ID3_expression <- "low"
plot <- rbind(expression_rad_sub_thirdqu,expression_rad_sub_fistqu )

plot$status2 <- ifelse(plot$Status=="deceased", 1,0)

fit <- survfit(Surv(Time..months.,status2) ~ ID3_expression, data = plot)
 
pdf(file.path(analysis.dir, "allpatients_survival_rad_median+sd.pdf"))
ggsurvplot(fit, xlab="Time [month]", palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE)
dev.off()



#NO radiotherapy

expression_rad_sub <- expression_brca_free[expression_brca_free$Radio.Therapy=="NO",]


#get quantiles
first_qu <- quantile(expression_rad_sub$ID3..mRNA.expression..microarray.)[2]
third_qu <- quantile(expression_rad_sub$ID3..mRNA.expression..microarray.)[4]
#subset
expression_rad_sub_fistqu <- expression_rad_sub[expression_rad_sub$ID3..mRNA.expression..microarray. <= first_qu,]
dim(expression_rad_sub_fistqu)
expression_rad_sub_thirdqu <- expression_rad_sub[expression_rad_sub$ID3..mRNA.expression..microarray. >= third_qu,]
dim(expression_rad_sub_thirdqu)

#get survival data for those patients
library(dplyr)
expression_rad_sub_fistqu<- left_join(expression_rad_sub_fistqu, survival, by="Case.ID")
summary(expression_rad_sub_fistqu$Survival.Rate)
expression_rad_sub_thirdqu<- left_join(expression_rad_sub_thirdqu, survival, by="Case.ID")
summary(expression_rad_sub_thirdqu$Survival.Rate)

#plot caplan meier survival curves
library(survminer)
library(survival)

expression_rad_sub_thirdqu$ID3_expression <- "high"
expression_rad_sub_fistqu$ID3_expression <- "low"
plot <- rbind(expression_rad_sub_thirdqu,expression_rad_sub_fistqu )

plot$status2 <- ifelse(plot$Status=="deceased", 1,0)

fit <- survfit(Surv(Time..months.,status2) ~ ID3_expression, data = plot)
 
pdf(file.path(analysis.dir, "noRadiotherapy_survival_rad_25.pdf"))
ggsurvplot(fit, xlab="Time [month]", palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE)
dev.off()


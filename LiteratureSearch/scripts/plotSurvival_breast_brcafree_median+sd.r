#look at ID3 depenbdent survival
#directories
base.dir<- "c010-datasets/Internal/2018-Ali/LiteratureSearch"
analysis.dir <- file.path(base.dir, "analysis")
data.dir <- file.path(base.dir, "data")

#libraries 
library(ggpubr)

#load data
survival <- read.table(file.path(data.dir, "Overall_Survival_metabric.txt"),sep="\t", header=TRUE)
expression <- read.table(file.path(data.dir, "ID3_expresssion_metabric.txt"),sep="\t", header=TRUE)
brca_mut <-  read.table(file.path(data.dir, "BRCA1_BRCA2_patient_data.txt"),sep="\t", header=TRUE)

#define patients with brca1 alteration
brca_mut_df<- sapply(brca_mut, function(x)as.character(x))
brca_mut_df[brca_mut_df==""]<- NA
brca_mut_df[brca_mut_df=="NA"]<- NA
rownames(brca_mut_df)<- brca_mut$track_name
brca1 <- brca_mut_df[,grep("BRCA1", colnames(brca_mut_df))]
brca1_na<- is.na(brca1)
brca1_na<-rowSums(brca1_na)
brca1_mut<- names(brca1_na[brca1_na!=5])
length(brca1_mut)

brca2 <- brca_mut_df[,grep("BRCA2", colnames(brca_mut_df))]
brca2_na<- is.na(brca2)
brca2_na<-rowSums(brca2_na)
brca2_mut<- names(brca2_na[brca2_na!=5])
length(brca2_mut)

brca_mut <- unique(c(brca1_mut,brca2_mut ))
length(brca_mut)


#rename expression colnames
colnames(expression)<- c("Case.ID", "Radio.Therapy", "ID3..mRNA.expression..microarray.", "Mutations" )
#subset by brca mutation
dim(expression)
expression_brca_free <- expression[!expression$Case.ID %in% brca_mut,  ]
dim(expression_brca_free)

#plot percentage of brca mut in total
brca_table <- as.data.frame(table(expression$Case.ID %in% brca_mut))
pdf(file.path(analysis.dir, "BRCA_alternated_Patients.pdf"), height=3.5, width=3.5)
ggbarplot(brca_table, x="Var1", y="Freq", fill="Var1", col="grey", palette="jco", #add="jitter",
xlab="BRCA alterations", ylab="Number of Patients")+rremove("legend")# +geom_hline(yintercept=c(first_qu, third_qu), linetype = 2)
dev.off()
#subset radio therapy
expression_rad_sub <- expression_brca_free[expression_brca_free$Radio.Therapy=="YES",]
dim(expression_rad_sub)
expression_rad_sub2 <- expression[expression$Radio.Therapy=="YES",]
dim(expression_rad_sub2)
#get distribution of ID3 expression
summary(expression_rad_sub$ID3..mRNA.expression..microarray.)

#plot distribution of id3 expression
pdf(file.path(analysis.dir, "BRCAfreesub_Metabric_ID3_Radio_Distr2.pdf"), height=3.5, width=3.5)
ggboxplot(expression_rad_sub, x="Radio.Therapy", y="ID3..mRNA.expression..microarray.", fill="Radio.Therapy", col="grey", palette="jco", #add="jitter",
xlab="Radio Therapy", ylab="ID3 expression")+rremove("legend")# +geom_hline(yintercept=c(first_qu, third_qu), linetype = 2)
dev.off()
#plot number of cases
library(reshape2)
cases <- melt(as.data.frame(table(expression_rad_sub$Radio.Therapy)))
pdf(file.path(analysis.dir, "BRCAfreesub_Metabric_BreastCancerCases.pdf"), height=3.5, width=3.5)
ggbarplot(cases, x="Var1", y="value", fill="Var1", col="grey", palette="jco", #add="jitter",
xlab="Radio Therapy", ylab="Number of Cases")+rremove("legend")
dev.off()

#subset
lowCutoff <- median(expression_rad_sub$ID3..mRNA.expression..microarray. - sd(expression_rad_sub$ID3..mRNA.expression..microarray.))
highCutoff <- median(expression_rad_sub$ID3..mRNA.expression..microarray. + sd(expression_rad_sub$ID3..mRNA.expression..microarray.))

expression_rad_sub_fistqu <- expression_rad_sub[expression_rad_sub$ID3..mRNA.expression..microarray. <= lowCutoff,]
dim(expression_rad_sub_fistqu)
expression_rad_sub_thirdqu <- expression_rad_sub[expression_rad_sub$ID3..mRNA.expression..microarray. >= highCutoff,]
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
 
pdf(file.path(analysis.dir, "BRCAfreesub_survival_rad_median+sd.pdf"))
ggsurvplot(fit, xlab="Time [month]", palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE)
dev.off()
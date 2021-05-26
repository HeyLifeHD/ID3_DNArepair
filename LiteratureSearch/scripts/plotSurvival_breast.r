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
#rename expression colnames
colnames(expression)<- c("Case.ID", "Radio.Therapy", "ID3..mRNA.expression..microarray.", "Mutations" )
#subset radio therapy
dim(expression)
expression_rad_sub <- expression[expression$Radio.Therapy=="YES",]
dim(expression_rad_sub)

#get distribution of ID3 expression
summary(expression_rad_sub$ID3..mRNA.expression..microarray.)
#plot distribution of id3 expression
first_qu <- quantile(expression_rad_sub$ID3..mRNA.expression..microarray.)[2]
third_qu <- quantile(expression_rad_sub$ID3..mRNA.expression..microarray.)[4]

pdf(file.path(analysis.dir, "Metabric_ID3_Radio_Distr2.pdf"), height=3.5, width=3.5)
ggboxplot(expression, x="Radio.Therapy", y="ID3..mRNA.expression..microarray.", fill="Radio.Therapy", col="grey", palette="jco", #add="jitter",
xlab="Radio Therapy", ylab="ID3 expression")+rremove("legend")# +geom_hline(yintercept=c(first_qu, third_qu), linetype = 2)
dev.off()
#plot number of cases
library(reshape2)
cases <- melt(as.data.frame(table(expression$Radio.Therapy)))
pdf(file.path(analysis.dir, "Metabric_BreastCancerCases.pdf"), height=3.5, width=3.5)
ggbarplot(cases, x="Var1", y="value", fill="Var1", col="grey", palette="jco", #add="jitter",
xlab="Radio Therapy", ylab="Number of Cases")+rremove("legend")
dev.off()

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

expression_rad_sub_thirdqu$ID3_expression <- "high"
expression_rad_sub_fistqu$ID3_expression <- "low"
plot <- rbind(expression_rad_sub_thirdqu,expression_rad_sub_fistqu )

plot$status2 <- ifelse(plot$Status=="deceased", 1,0)

fit <- survfit(Surv(Time..months.,status2) ~ ID3_expression, data = plot)
 
pdf(file.path(analysis.dir, "survival_rad_25.pdf"))
ggsurvplot(fit, xlab="Time [month]", palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE)
dev.off()


#get 10

#get quantiles
first_qu <- quantile(expression_rad_sub$ID3..mRNA.expression..microarray., probs=c(.05))
third_qu <- quantile(expression_rad_sub$ID3..mRNA.expression..microarray., probs=c(.95))
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

expression_rad_sub_thirdqu$ID3_expression <- "high"
expression_rad_sub_fistqu$ID3_expression <- "low"
plot <- rbind(expression_rad_sub_thirdqu,expression_rad_sub_fistqu )

plot$status2 <- ifelse(plot$Status=="deceased", 1,0)

fit <- survfit(Surv(Time..months.,status2) ~ ID3_expression, data = plot)
 
pdf(file.path(analysis.dir, "survival_rad_5.pdf"))
ggsurvplot(fit, xlab="Time [month]", palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE)
dev.off()



library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/")
test(mycgds)
setVerbose(mycgds, TRUE)

# Get list of cancer studies at server
getCancerStudies(mycgds)[,c(1,2)]

getGeneticProfiles(mycgds,'brca_metabric')[,c(1:2)]
getCaseLists(mycgds,'brca_metabric')[,c(1:2)]
getProfileData(mycgds, "NF1", c("gbm_tcga_gistic","gbm_tcga_mrna"), "gbm_tcga_all")[c(1:5),]
getProfileData(mycgds, c("MDM2","MDM4"), "gbm_tcga_mrna", "gbm_tcga_all")[c(25:30),]
 getClinicalData(mycgds, "brca_metabric")
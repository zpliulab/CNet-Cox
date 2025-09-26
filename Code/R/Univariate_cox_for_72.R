library(survival)
library(plyr)

rm(list = ls())

## set pathway
path <- '/Users/lilingyu/E/PhD/R/'
# path <- '/home/lly/R/'


# Creat Files ------------------------------------------------------------------
setwd(paste(path, 'CNetCox/Data/Result', sep=''))


setwd("D:\\E\\博士\\R_程序\\BRCA\\Data\\TCGA_NEW\\Result")

data = read.table("TCGA_feature_data.txt", header = T, check.names = FALSE)
Data <- data.frame(t(data))
# View(Data[,1:10])

str(Data)

Data$status<-factor(Data$status)
summary(Data$status)


## -----------------------
y <- Surv(time=Data$time, event=Data$status==0)  #0Ϊ????

Uni_cox_model<- function(x){
  # x <- variable_names
  FML <- as.formula(paste0 ("y~",x))
  cox <- coxph(FML,data=Data)
  cox1 <-summary(cox)
  coef <- cox1$coefficients[,1]
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <- round(cox1$conf.int[,3],2)
  CI95 <- round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' =x,
                             names <-rownames(cox1$conf.int),
                             'coef' = coef,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'Uni_P' = PValue
  )
  return(Uni_cox_model)}  

names(Data)

variable_names <- colnames(Data)[c(3:74)]  

Uni_cox <- lapply(variable_names, Uni_cox_model)
Uni_cox <- ldply(Uni_cox, data.frame)

Uni_cox$HR.CI95 <- paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")")
Uni_cox <- Uni_cox[,-4:-6] 

View(Uni_cox)

which(Uni_cox[,4] <= 0.05)
Uni_cox[which(Uni_cox[,4] <= 0.05),c(1,3)]


## ----------------------------

Uni_cox$Characteristics[Uni_cox$Uni_P<0.05]

mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$Uni_P<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=Data)
cox4 <- summary(mul_cox) 
coef <- cox4$coefficients[,c(1,5)]
mul_coef <- coef[which(coef[,2] <= 0.05),1]
View(mul_coef)
# write.csv(mul_coef, file = "univariate_cox_coef.csv")

mul_corf <- cox4$coefficients[,1] 
mul_HR <- round(cox4$coefficients[,2],2) 
mul_PValue<- round(cox4$coefficients[,5],4) 
mul_CI1<-round(cox4$conf.int[,3],2)
mul_CI2<-round(cox4$conf.int[,4],2)

## HR(95%CI)+P
mul_HR.CI95 <- paste(mul_HR,"(",mul_CI1,'-',mul_CI2,")")
mul_cox1 <- data.frame("coef"=mul_corf, "mul_HR.CI95"=mul_HR.CI95,"P"=mul_PValue)

## HR+95%CI+P
#mul_CI<-paste(mul_CI1,'-',mul_CI2)
#mul_cox1<- data.frame("HR"=mul_HR,"mul_CI"=mul_CI, "P"=mul_PValue)


## ------------------------------

Uni_cox <- Uni_cox[,-1]

colnames(Uni_cox)[1] <- 'Characteristics'

mul_cox1 <- cbind(rownames(mul_cox1), mul_cox1, row.names=NULL); names(mul_cox1 )[1]<-"Characteristics"

table2 <- merge.data.frame(Uni_cox, mul_cox1, by="Characteristics", all = T, sort = T)

table3 <- table2[,c(1,2,4,3,5,6,7)]
colnames(table3) <- c("Feature", "Coef_Uni", "HR(95%CI)_Uni", "P_Uni",
                     "Coef_Mul", "HR(95%CI)_Mul", "P_Mul")
View(table3)

# write.csv(table3, file = "univariate_cox.csv",row.names = F)


## latex ----------------------------------------------------------------
library(stargazer)
stargazer(table3, summary=FALSE, rownames=FALSE) 

stargazer(table3[,], summary=FALSE, rownames=FALSE) 

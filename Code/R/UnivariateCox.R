
library(survival)
library(plyr)

rm(list = ls())


# 导入数据
setwd("D:\\E\\博士\\R_程序\\BRCA\\Data\\TCGA_NEW\\result")

data = read.table("TCGA_feature_data.txt", header = T, check.names = FALSE)
Data <- data.frame(t(data))
# View(Data[,1:10])

#4.查看数据数据性质
str(Data)

#5.查看结局，0=复发，1未复发
Data$status<-factor(Data$status)
summary(Data$status)



################### 一、批量单因素回归 -----------------------

#1.构建模型的y
y <- Surv(time=Data$time, event=Data$status==0)  #0为复发

#2.批量单因素回归模型建立：Uni_cox_model
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

#3.将想要进行的单因素回归变量输入模型

#3-(1)查看变量的名字和序号
names(Data)

#3-(2)输入变量序号
variable_names <- colnames(Data)[c(3:74)] #例：这里选择了3-10号变量

#4.输出结果
Uni_cox <- lapply(variable_names, Uni_cox_model)
Uni_cox <- ldply(Uni_cox, data.frame)

#5.优化表格，这里举例HR+95% CI+P 风格
Uni_cox$HR.CI95 <- paste0(Uni_cox$HR," (",Uni_cox$CI5,'-',Uni_cox$CI95,")")
Uni_cox <- Uni_cox[,-4:-6] #HR (95% CI)+P风格

# Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
# Uni_cox<-Uni_cox[,-4:-6] #方法2：HR+95% CI+P 风格

#查看单因素cox表格
View(Uni_cox)

which(Uni_cox[,4] <= 0.05)
Uni_cox[which(Uni_cox[,4] <= 0.05),c(1,3)]


##################### 二、单因素回归p<0.05做多因素回归 ----------------------------

#1.提取单因素p<0.05变量
Uni_cox$Characteristics[Uni_cox$Uni_P<0.05]

#2.多因素模型建立
mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$Uni_P<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=Data)
cox4 <- summary(mul_cox) 
coef <- cox4$coefficients[,c(1,5)]
mul_coef <- coef[which(coef[,2] <= 0.05),1]
View(mul_coef)
# write.csv(mul_coef, file = "univariate_cox_coef.csv")

#3.提取多因素回归的信息
mul_corf <- cox4$coefficients[,1] 
mul_HR <- round(cox4$coefficients[,2],2) 
mul_PValue<- round(cox4$coefficients[,5],4) 
mul_CI1<-round(cox4$conf.int[,3],2)
mul_CI2<-round(cox4$conf.int[,4],2)

#4.多因素结果优化并成表：mul_cox1
##4-1：HR(95%CI)+P风格
mul_HR.CI95 <- paste(mul_HR,"(",mul_CI1,'-',mul_CI2,")")
mul_cox1 <- data.frame("coef"=mul_corf, "mul_HR.CI95"=mul_HR.CI95,"P"=mul_PValue)

##4-2：HR+95%CI+P风格
#mul_CI<-paste(mul_CI1,'-',mul_CI2)
#mul_cox1<- data.frame("HR"=mul_HR,"mul_CI"=mul_CI, "P"=mul_PValue)


#####################  三、单因素多因素整合成一个表  ------------------------------
#1.删除单因素表Uni_cox的第一列
Uni_cox <- Uni_cox[,-1]

#2.删除第一列后的单因素表Uni_cox的第一列命名为Characteristics
colnames(Uni_cox)[1] <- 'Characteristics'

#3.多因素表mul_cox1的行名放入单元格中，也命名为Characteristics
mul_cox1 <- cbind(rownames(mul_cox1), mul_cox1, row.names=NULL); names(mul_cox1 )[1]<-"Characteristics"

#4.Uni_cox表和mul_cox1表以Characteristics列为标准合并
table2 <- merge.data.frame(Uni_cox, mul_cox1, by="Characteristics", all = T, sort = T)

#5.查看最终表格
table3 <- table2[,c(1,2,4,3,5,6,7)]
colnames(table3) <- c("Feature", "Coef_Uni", "HR(95%CI)_Uni", "P_Uni",
                     "Coef_Mul", "HR(95%CI)_Mul", "P_Mul")
View(table3)

# write.csv(table3, file = "univariate_cox.csv",row.names = F)


# 转出为latex ----------------------------------------------------------------
library(stargazer)
stargazer(table3, summary=FALSE, rownames=FALSE) 

stargazer(table3[,], summary=FALSE, rownames=FALSE) 

setwd("Desktop/Final/parasite/") 
dir()
data<- read.table("opg-hiding.txt", header = T)
#数据独立性、正态性、方差齐性检验
datax<-table(data$Treat,data$opg)
chisq.test(datax)#卡方独立性检验
shapiro.test(data$opg)#夏皮罗—威尔克正态检验
bartlett.test(data$opg~data$Treat,data)#Bartlett方差齐性检验

#P<0.05数据不满足正态性
wilcox.test(data$opg~data$Treat,paired=FALSE)

#kruskal.test(data$opg~data$Treat,data)#非参数检验

#t.test(data$opg~data$Treat)


data2<- read.table("opg-Confilict phase.txt", header = T)
#数据独立性、正态性、方差齐性检验
datax<-table(data2$Treat,data2$opg)
chisq.test(datax)#卡方独立性检验
shapiro.test(data2$opg)#夏皮罗—威尔克正态检验
bartlett.test(data2$opg~data2$Treat,data2)#Bartlett方差齐性检验

#P<0.05数据不满足正态性
wilcox.test(data2$opg~data2$Treat,paired=FALSE)

#kruskal.test(data2$opg~data2$Treat,data2)#非参数检验

#t.test(data2$opg~data2$Treat)

data3<- read.table("opg-Disintegrated phase.txt", header = T)
#数据独立性、正态性、方差齐性检验
datax<-table(data3$Treat,data3$opg)
chisq.test(datax)#卡方独立性检验
shapiro.test(data3$opg)#夏皮罗—威尔克正态检验
bartlett.test(data3$opg~data3$Treat,data3)#Bartlett方差齐性检验

#P<0.05数据不满足正态性
wilcox.test(data3$opg~data3$Treat,paired=FALSE)

#kruskal.test(data3$opg~data3$Treat,data3)#非参数检验

#t.test(data3$opg~data3$Treat)

data4<- read.table("opg-Independent phase.txt", header = T)
#数据独立性、正态性、方差齐性检验
datax<-table(data4$Treat,data4$opg)
chisq.test(datax)#卡方独立性检验
shapiro.test(data4$opg)#夏皮罗—威尔克正态检验
bartlett.test(data4$opg~data3$Treat,data4)#Bartlett方差齐性检验

#P<0.05数据不满足正态性
wilcox.test(data4$opg~data4$Treat,paired=FALSE)

#kruskal.test(data4$opg~data3$Treat,data4)#非参数检验

#t.test(data4$opg~data4$Treat)

data5<- read.table("Control.txt", header = T)
attach(data5)
wilcox.test(opg~Treat,paired=FALSE)



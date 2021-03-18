#### Distance Sampling - Parrots - MArch 2021 ####

# Load packages ####
library(Distance)

# Load dataset ####
df_parrots <- read.csv2("~/Documents/GitHub/DS_parrots_2021/df_parrots.csv",stringsAsFactors=TRUE)

head(df_parrots)
str(df_parrots)

table(df_parrots$Ecoregion,df_parrots$Pais)

# Separate dataframe into a list of dataframes, one for each species/ecoregion/country combination. ####

list1<-list()
for(i in 1:length(levels(df_parrots$Pais))){
  dfcountry<-subset(df_parrots,Pais==levels(df_parrots$Pais)[i])
  dfcountry<-droplevels(dfcountry)
  for(z in 1:length(levels(dfcountry$Ecoregion))){
    dfecoregion <- subset(dfcountry,Ecoregion==levels(dfcountry$Ecoregion)[z])
    dfecoregion<-droplevels(dfecoregion)
    for(k in 1:length(levels(dfecoregion$Species))){
      dfspecies<-subset(dfecoregion,Species==levels(dfecoregion$Species)[k])
      dfspecies<-droplevels(dfspecies)
      list1[[length(list1) + 1]]<-dfspecies
    }
  }
}
rm(list=setdiff(ls(), "list1"))

# There are 212 list elements. Need to fit DS model to each

# Setup dataframe to store detection probabilities and other data

output<-matrix(NA,nrow=212,ncol=13)
output<-data.frame(output)
colnames(output)<-c("Country",
                    "Ecoregion",
                    "Species",
                    "Km",
                    "Count",
                    "Encounters_(<=500m)",
                    "Group_sizeXdistance_cor_test_p",
                    "DS_cutpoints",
                    "Detection_probability",
                    "Detection_function",
                    "GoF_DS_Xsq",
                    "N_(count/p)",
                    "D_(N/Km)")
output
str(output)

for(i in 1:length(list1)){
  output$Country[i]<-as.character(levels(list1[[i]]$Pais))
  output$Ecoregion[i]<-levels(list1[[i]]$Ecoregion)
  output$Species[i]<-levels(list1[[i]]$Species)
  output$Km[i]<-list1[[i]]$Km.total[1]
  output$Count[i]<-sum(list1[[i]]$size)
}
output
write.csv(output,"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/output_table.csv",row.names = F)

rm(list=setdiff(ls(), "list1"))



# 1 ####
index<-1

df.all<-list1[[index]]
df.all
df<-subset(df.all,distance<501)
df
hist(df$distance,xlab="Distance (m)")

plot(df$distance~df$size) # explore if detection distance is correlated with group size
cor.test(df$distance,df$size,method="spearman")

cutpoints <- c(0,100,200,300,400,500)
hist(df$distance,xlab="Distance (m)",breaks=cutpoints)

DS1 <- ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints)
(gof1 <- gof_ds(DS1,plot=T))
plot(DS1,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x1<-summary(DS1))

DS2 <- ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints) # Half-normal key function with cosine(2) adjustments selected.
(gof2 <- gof_ds(DS2,plot=T))
plot(DS2,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x2<-summary(DS2))

DS3 <- ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints) # Half-normal key function selected.

DS4 <- ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints) # Half-normal key function selected.

AIC(DS1,DS2)

output <- read.csv("~/Documents/GitHub/DS_parrots_2021/output_table.csv",stringsAsFactors = F,check.names = F)
output[index,]$Detection_probability <-round(x2$ds$average.p,digits=2)
output[index,]$Detection_function <-x2$ddf$name.message
output[index,]$GoF_DS_Xsq <- round(gof2$chisquare$chi1$p,digits=2)
output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$p,digits=2)
output[index,]$`D_(N/Km)`<-round(output[index,]$N/output[index,]$Km,digits=2)
output[index,]$`Encounters_(<=500m)` <-nrow(df)
output[index,]
write.csv(output,"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/output_table.csv",row.names = F)

pdf(paste0("/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/figures/",paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "),".pdf"))
plot(DS2,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
dev.off()

rm(list=setdiff(ls(), "list1"))

# 2 ####
index<-2

df.all<-list1[[index]]
df.all
df<-subset(df.all,distance<501)
df
hist(df$distance,xlab="Distance (m)")

plot(df$distance~df$size) # explore if detection distance is correlated with group size
cor.test(df$distance,df$size,method="spearman")

cutpoints <- c(0,50,100,150,200,250,300,350,400,450,500)
hist(df$distance,xlab="Distance (m)",breaks=cutpoints)

DS1 <- ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints)
(gof1 <- gof_ds(DS1,plot=T))
plot(DS1,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x1<-summary(DS1))

DS2 <- ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints) # Half-normal key function selected.

DS3 <- ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints) # Half-normal key function selected.

DS4 <- ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints) # Half-normal key function selected.

output <- read.csv("~/Documents/GitHub/DS_parrots_2021/output_table.csv",stringsAsFactors = F,check.names = F)
output[index,]$p <-round(x1$ds$average.p,digits=2)
output[index,]$Detection_function <-x1$ddf$name.message
output[index,]$GoF_DS_Xsq <- round(gof1$chisquare$chi1$p,digits=2)
output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$p,digits=2)
output[index,]$`D_(N/Km)`<-round(output[index,]$N/output[index,]$Km,digits=2)
output[index,]$`Encounters_(<=500m)` <-nrow(df)
output[index,]
write.csv(output,"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/output_table.csv",row.names = F)

pdf(paste0("/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/figures/",paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "),".pdf"))
plot(DS1,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
dev.off()

rm(list=setdiff(ls(), "list1"))

# 3 ####
index<-3

df.all<-list1[[index]]
df.all
df<-subset(df.all,distance<501)
df
hist(df$distance,xlab="Distance (m)")

plot(df$distance~df$size) # explore if detection distance is correlated with group size
cor.test(df$distance,df$size,method="spearman")

cutpoints <- seq(0,500,by=20)
hist(df$distance,xlab="Distance (m)",breaks=cutpoints)

DS1 <- ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints)
(gof1 <- gof_ds(DS1,plot=T))
plot(DS1,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x1<-summary(DS1))

DS2 <- ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints) # Half-normal key function with cosine(2) adjustments selected.
(gof2 <- gof_ds(DS2,plot=T))
plot(DS2,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x2<-summary(DS2))

DS3 <- ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints) # Half-normal key function with Hermite(4) adjustments
(gof3 <- gof_ds(DS3,plot=T))
plot(DS3,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x3<-summary(DS3))

DS4 <- ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints) # Half-normal key function selected.

AIC(DS1,DS2,DS3) #DS2 and DS3 and essentially the same detection curve. Selecting DS2

output <- read.csv("~/Documents/GitHub/DS_parrots_2021/output_table.csv",stringsAsFactors = F,check.names = F)
output[index,]$p <-round(x2$ds$average.p,digits=2)
output[index,]$Detection_function <-x2$ddf$name.message
output[index,]$GoF_DS_Xsq <- round(gof2$chisquare$chi1$p,digits=2)
output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$p,digits=2)
output[index,]$`D_(N/Km)`<-round(output[index,]$N/output[index,]$Km,digits=2)
output[index,]$`Encounters_(<=500m)` <-nrow(df)
output[index,]
write.csv(output,"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/output_table.csv",row.names = F)

pdf(paste0("/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/figures/",paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "),".pdf"))
plot(DS2,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
dev.off()

rm(list=setdiff(ls(), "list1"))

# 4 ####
index<-4

df.all<-list1[[index]]
df.all
df<-subset(df.all,distance<501)
df
hist(df$distance,xlab="Distance (m)")

plot(df$distance~df$size) # explore if detection distance is correlated with group size
cor.test(df$distance,df$size,method="spearman")

cutpoints <- seq(0,500,by=100)
hist(df$distance,xlab="Distance (m)",breaks=cutpoints)

DS1 <- ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints)
(gof1 <- gof_ds(DS1,plot=T))
plot(DS1,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
(x1<-summary(DS1))

DS2 <- ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints) # Half-normal key function selected.

DS3 <- ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints) # Half-normal key function selected.

DS4 <- ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints) # Half-normal key function selected.

AIC(DS1,DS2)

output <- read.csv("~/Documents/GitHub/DS_parrots_2021/output_table.csv",stringsAsFactors = F,check.names = F)
output[index,]$p <-round(x1$ds$average.p,digits=2)
output[index,]$Detection_function <-x1$ddf$name.message
output[index,]$GoF_DS_Xsq <- round(gof1$chisquare$chi1$p,digits=2)
output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$p,digits=2)
output[index,]$`D_(N/Km)`<-round(output[index,]$N/output[index,]$Km,digits=2)
output[index,]$`Encounters_(<=500m)` <-nrow(df)
output[index,]
write.csv(output,"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/output_table.csv",row.names = F)

pdf(paste0("/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/figures/",paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "),".pdf"))
plot(DS2,breaks=cutpoints,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
dev.off()

rm(list=setdiff(ls(), "list1"))

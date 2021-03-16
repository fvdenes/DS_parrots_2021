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


# There are 212 list elements. Need to fit DS model to each

# Setup dataframe to store detection probabilities and other data

output<-matrix(NA,nrow=212,ncol=8)
output<-data.frame(output)
colnames(output)<-c("Country", "Ecoregion","Species","Km","count","p","N", "D (ind/Km)")
output
str(output)

for(i in 1:length(list1)){
  output$Country[i]<-as.character(levels(list1[[i]]$Pais))
  output$Ecoregion[i]<-levels(list1[[i]]$Ecoregion)
  output$Species[i]<-levels(list1[[i]]$Species)
  output$Km[i]<-list1[[i]]$Km.total[1]
  output$count[i]<-sum(list1[[i]]$size)
}
output

# 1 ####
index<-1
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
plot(DS1,breaks=cutpoints)
(x<-summary(DS1))

output[index,]$p<-x$ds$average.p
output[index,]$N<-output[index,]$count/output[index,]$p
output[index,]$`D (ind/Km)`<-output[index,]$N/output[index,]$Km
output[index,]



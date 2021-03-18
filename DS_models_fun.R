# A function to fit DS models, evaluate GoF, print plots and save results for each of the spp/ecoregion/country combinations

output.folder<-"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/"
df.all<-list1[[1]]

DSmodels<-function(df.all,print.plots=T,save.output=T,output.folder,max.dist=500,cutpoints=NULL){
  df<-subset(df.all,distance<max.dist+1)

  cor_test<-cor.test(df$distance,df$size,method="spearman")

  if(is.null(cutpoints)==T){
    cutpoints<-list(
      seq(0,max.dist,by=100),
      seq(0,max.dist,by=50),
      seq(0,max.dist,by=25)
    )
  }


  if(cor_test$p.value>0.05){

    listDSmodels<-list(NA)
    for(i in 1:length(cutpoints)){
      DS1 <- ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints[[i]])
      DS2 <- ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints[[i]])
      DS3 <- ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints[[i]])
      DS4 <- ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints[[i]])


      listDSmodels[[i]]<-c(list(DS1),list(DS2,DS3,DS4)[which(c(DS2$ddf$name.message,DS3$ddf$name.message,DS4$ddf$name.message)!=DS1$ddf$name.message)])
    }

    listDSmodels2  <- unlist(listDSmodels,recursive=F)

    GoF<-lapply(listDSmodels2,gof_ds)
    GoF_p<-NA
    for(i in 1:length(GoF)){
      GoF_p[i]<-GoF[[i]]$chisquare$chi1$p
    }

    SelMod<-listDSmodels2[[which.max(GoF_p)]]
    summarySelMod<-summary(SelMod)

  }



  output <- read.csv(paste0(output.folder,"output_table.csv"),stringsAsFactors = F,check.names = F)

  combination<-paste0(df[1,1],df[1,2],df[1,4])
  output.combinations<-apply(output,1,function(x)paste0(x[1],x[2],x[3]))
  index<-which(combination%in%output.combinations)

  output[index,]$Group_sizeXdistance_cor_test_p <- round(cor_test$p.value,digits=2)
  output[index,]$DS_cutpoints<-paste0(SelMod$ddf$meta.data$breaks[2]-SelMod$ddf$meta.data$breaks[1],"m")
  output[index,]$Detection_probability <-round(summarySelMod$ds$average.p,digits=2)
  output[index,]$Detection_function <-summarySelMod$ddf$name.message
  output[index,]$GoF_DS_Xsq <- round(max(GoF_p),digits=2)
  output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$Detection_probability,digits=2)
  output[index,]$`D_(N/Km)`<-round(output[index,]$`N_(count/p)`/output[index,]$Km,digits=2)
  output[index,]$`Encounters_(<=500m)` <-nrow(df)
  output[index,]
  write.csv(output,paste0(output.folder,"output_table.csv"),row.names = F)

  pdf(paste0(output.folder,"figures/",paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "),".pdf"))
  hist(df$distance,xlab="Distance (m)",main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
  plot(xlab="Distance (m)",ylab="Group size",df$distance~df$size,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
  plot(SelMod,breaks=SelMod$ddf$meta.data$breaks,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
  dev.off()

}

listDSmodels2[[6]]$ddf$meta.data$breaks

# 1 ####
index<-1

df.all<-list1[[index]]
df.all
df<-subset(df.all,distance<501)
df
hist(df$distance,xlab="Distance (m)")

plot(xlab="Distance (m)",ylab="Group size",df$distance~df$size,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - ")) # explore if detection distance is correlated with group size
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

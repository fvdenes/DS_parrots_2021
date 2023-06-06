# A function to fit DS models, evaluate GoF, print plots and save results for each of the spp/ecoregion/country combinations
library(Distance)
library(pbapply)

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

  if(max(df$distance)<100){
    cutpoints<-cutpoints[2:3]
  }

  if(max(df$distance)<50){
    cutpoints<-cutpoints[2]
  }

  if(cor_test$p.value>0.05){

    listDSmodels<-list(NA)
    for(i in 1:length(cutpoints)){
      DS1 <- try(suppressMessages(ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
      DS2 <- try(suppressMessages(ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
      DS3 <- try(suppressMessages(ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
      DS4 <- try(suppressMessages(ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))

      listDSmodels[[i]]<-try(c(list(DS1),list(DS2,DS3,DS4)[which(c(DS2$ddf$name.message,DS3$ddf$name.message,DS4$ddf$name.message)!=DS1$ddf$name.message)]))
    }

    listDSmodels2 <- unlist(listDSmodels,recursive=F)
    listDSmodels2 <- listDSmodels2[which(lapply(listDSmodels2,class)=="dsmodel")]

    GoF<-lapply(listDSmodels2,gof_ds)
    GoF_p<-NA
    for(i in 1:length(GoF)){
      GoF_p[i]<-GoF[[i]]$chisquare$chi1$p
    }

    SelMod<-listDSmodels2[[which.max(GoF_p)]]
    summarySelMod<-summary(SelMod)

  }

  if(cor_test$p.value<=0.05){

    listDSmodels<-list(NA)
    for(i in 1:length(cutpoints)){
      DS1 <- try(suppressMessages(ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
      DS2 <- try(suppressMessages(ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
      DS3 <- try(suppressMessages(ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
      DS4 <- try(suppressMessages(ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))
      DS5 <- try(suppressMessages(ds(data=df,key="hn",adjustment=NULL, formula = ~size, cutpoints = cutpoints[[i]],quiet=T)))

      listDSmodels[[i]]<- try(c(list(DS1),list(DS2,DS3,DS4)[which(c(DS2$ddf$name.message,DS3$ddf$name.message,DS4$ddf$name.message)!=DS1$ddf$name.message)],list(DS5)))
    }

    listDSmodels2 <- unlist(listDSmodels,recursive=F)
    listDSmodels2 <- listDSmodels2[which(lapply(listDSmodels2,class)=="dsmodel")]

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
  index<-which(output.combinations%in%combination)

  output[index,]$Group_sizeXdistance_cor_test_p <- round(cor_test$p.value,digits=2)
  output[index,]$DS_cutpoints<-paste0(SelMod$ddf$meta.data$breaks[2]-SelMod$ddf$meta.data$breaks[1],"m")
  output[index,]$Detection_probability <-round(summarySelMod$ds$average.p,digits=2)

  if(any(is.na(GoF_p)==T)) {
    output[index,]$GoF_DS_Xsq <- max(GoF_p[-which(is.na(GoF_p))])
  } else output[index,]$GoF_DS_Xsq <- max(GoF_p)

  output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$Detection_probability,digits=2)
  output[index,]$`D_(N/Km)`<-round(output[index,]$`N_(count/p)`/output[index,]$Km,digits=2)
  output[index,]$`Encounters_(<=500m)` <-nrow(df)

  if(length(row.names(summarySelMod$ds$coeff$key.scale))>1){
    if(row.names(summarySelMod$ds$coeff$key.scale)[2]=="size"){
      output[index,]$Detection_function <-paste0(summarySelMod$ddf$name.message," with group size covariate")
      output[index,]$Group_size_beta <- round(summarySelMod$ddf$par[2], digits=2)
    }
    if(row.names(summarySelMod$ds$coeff$key.scale)[2]!="size"){
      output[index,]$Detection_function <-summarySelMod$ddf$name.message
      output[index,]$Group_size_beta <- "-"
    }
  }

  if(length(row.names(summarySelMod$ds$coeff$key.scale))==1){
    output[index,]$Detection_function <-summarySelMod$ddf$name.message
    output[index,]$Group_size_beta <- "-"
  }

  write.csv(output,paste0(output.folder,"output_table.csv"),row.names = F)

  pdf(paste0(output.folder,"figures/",paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "),".pdf"))
  hist(df$distance,xlab="Distance (m)",main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
  plot(xlab="Distance (m)",ylab="Group size",df$size~df$distance,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
  plot(SelMod,breaks=SelMod$ddf$meta.data$breaks,main=paste(df$Species[1],df$Ecoregion[1],df$Pais[1],sep = " - "))
  dev.off()

}

# Run the function
pblapply(list1,DSmodels,output.folder="/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/")


# Troubleshooting a few cases

#8

df.all<-list1[[5]]

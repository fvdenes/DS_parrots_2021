#### Distance Sampling - Parrots - June 2021 ####

# Load packages ####
library(Distance)
library(pbapply)

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

output<-matrix(NA,nrow=212,ncol=14)
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
                    "D_(N/Km)",
                    "Group_size_beta")
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
write.csv(output,"/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/review/output_table.csv",row.names = F)

rm(list=setdiff(ls(), "list1"))





# A function to fit DS models, evaluate GoF, print plots and save results for each of the spp/ecoregion/country combinations ####
library(Distance)
library(pbapply)

DSmodels<-function(df.all,print.plots=T,save.output=T,output.folder,max.dist=500,cutpoints=NULL){
  df<-subset(df.all,distance<max.dist+1)

  output <- read.csv(paste0(output.folder,"output_table.csv"),stringsAsFactors = F,check.names = F)
  combination<-paste0(df[1,1],df[1,2],df[1,4])
  print(combination)


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

      DS5 <- try(suppressMessages(ds(data=df,key="hr",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
      #DS6 <- try(suppressMessages(ds(data=df,key="hr",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
      #DS7 <- try(suppressMessages(ds(data=df,key="hr",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
      #DS8 <- try(suppressMessages(ds(data=df,key="hr",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))

      DS9 <- try(suppressMessages(ds(data=df,key="unif",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
      #DS10 <- try(suppressMessages(ds(data=df,key="unif",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
      #DS11 <- try(suppressMessages(ds(data=df,key="unif",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))


      #modlist<-list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11),class)=="dsmodel")]
      try(modlist<-list(DS1,DS2,DS3,DS4,DS5,DS9)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS9),class)=="dsmodel")])

      listDSmodels[[i]]<-try(modlist[[which.min(lapply(modlist,function(x) AIC(x)$AIC))]])

      if(class(listDSmodels[[i]])=="try-error"){
        DS1 <- try(suppressMessages(ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
        DS2 <- try(suppressMessages(ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
        DS3 <- try(suppressMessages(ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
        DS4 <- try(suppressMessages(ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))

        DS5 <- try(suppressMessages(ds(data=df,key="hr",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
        #DS6 <- try(suppressMessages(ds(data=df,key="hr",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
        #DS7 <- try(suppressMessages(ds(data=df,key="hr",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
        #DS8 <- try(suppressMessages(ds(data=df,key="hr",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))

        DS9 <- try(suppressMessages(ds(data=df,key="unif",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
        #DS10 <- try(suppressMessages(ds(data=df,key="unif",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
        #DS11 <- try(suppressMessages(ds(data=df,key="unif",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))


        #modlist<-list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11),class)=="dsmodel")]
        try(modlist<-list(DS1,DS2,DS3,DS4,DS5,DS9)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS9),class)=="dsmodel")])

        listDSmodels[[i]]<-try(modlist[[which.min(lapply(modlist,function(x) AIC(x)$AIC))]])
      }
    }


    listDSmodels2 <- listDSmodels

    GoF<-lapply(listDSmodels2,gof_ds)
    GoF_chi<-NA
    GoF_p<-NA
    for(i in 1:length(GoF)){
      GoF_chi[i]<-GoF[[i]]$chisquare$chi1$chisq
      GoF_p[i]<-GoF[[i]]$chisquare$chi1$p
    }

    SelMod<-listDSmodels2[[which.max(GoF_p)]]

    if(any(is.na(GoF_chi)==T)) {
      if(length(which(is.na(GoF_chi)))==1){
        SelMod2 <- try(listDSmodels2[[which(GoF_p<max(GoF_p)&GoF_p>min(GoF_p))]])
      }
      if(length(which(is.na(GoF_chi)))==2){
        SelMod2 <- NA
      }
    } else{
      SelMod2<-try(listDSmodels2[[which(GoF_p<max(GoF_p)&GoF_p>min(GoF_p))]])
      SelMod3<-try(listDSmodels2[[which.min(GoF_p)]])}

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

      DS6 <- try(suppressMessages(ds(data=df,key="hr",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
      #DS7 <- try(suppressMessages(ds(data=df,key="hr",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
      #DS8 <- try(suppressMessages(ds(data=df,key="hr",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
      #DS9 <- try(suppressMessages(ds(data=df,key="hr",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))
      DS10 <- try(suppressMessages(ds(data=df,key="hr",adjustment=NULL, formula = ~size, cutpoints = cutpoints[[i]],quiet=T)))

      DS11 <- try(suppressMessages(ds(data=df,key="unif",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
      #DS12 <- try(suppressMessages(ds(data=df,key="unif",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
      #DS13 <- try(suppressMessages(ds(data=df,key="unif",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))

      #modlist<-list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11,DS12,DS13)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11,DS12,DS13),class)=="dsmodel")]

      modlist<-list(DS1,DS2,DS3,DS4,DS5,DS6,DS10,DS11)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS6,DS10,DS11),class)=="dsmodel")]

      listDSmodels[[i]]<-try(modlist[[which.min(lapply(modlist,function(x) AIC(x)$AIC))]])

      if(class(listDSmodels[[i]])=="try-error"){
        DS1 <- try(suppressMessages(ds(data=df,key="hn",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
        DS2 <- try(suppressMessages(ds(data=df,key="hn",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
        DS3 <- try(suppressMessages(ds(data=df,key="hn",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
        DS4 <- try(suppressMessages(ds(data=df,key="hn",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))
        DS5 <- try(suppressMessages(ds(data=df,key="hn",adjustment=NULL, formula = ~size, cutpoints = cutpoints[[i]],quiet=T)))

        DS6 <- try(suppressMessages(ds(data=df,key="hr",adjustment=NULL,cutpoints = cutpoints[[i]],quiet=T)))
        #DS7 <- try(suppressMessages(ds(data=df,key="hr",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
        #DS8 <- try(suppressMessages(ds(data=df,key="hr",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
        #DS9 <- try(suppressMessages(ds(data=df,key="hr",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))
        DS10 <- try(suppressMessages(ds(data=df,key="hr",adjustment=NULL, formula = ~size, cutpoints = cutpoints[[i]],quiet=T)))

        DS11 <- try(suppressMessages(ds(data=df,key="unif",adjustment="cos",cutpoints = cutpoints[[i]],quiet=T)))
        #DS12 <- try(suppressMessages(ds(data=df,key="unif",adjustment="herm",cutpoints = cutpoints[[i]],quiet=T)))
        #DS13 <- try(suppressMessages(ds(data=df,key="unif",adjustment="poly",cutpoints = cutpoints[[i]],quiet=T)))

        #modlist<-list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11,DS12,DS13)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS6,DS7,DS8,DS9,DS10,DS11,DS12,DS13),class)=="dsmodel")]

        modlist<-list(DS1,DS2,DS3,DS4,DS5,DS6,DS10,DS11)[which(lapply(list(DS1,DS2,DS3,DS4,DS5,DS6,DS10,DS11),class)=="dsmodel")]

        listDSmodels[[i]]<-try(modlist[[which.min(lapply(modlist,function(x) AIC(x)$AIC))]])
      }
    }

    listDSmodels2 <- listDSmodels

    GoF<-lapply(listDSmodels2,gof_ds)
    GoF_chi<-NA
    GoF_p<-NA

     for(i in 1:length(GoF)){
      GoF_chi[i]<-GoF[[i]]$chisquare$chi1$chisq
      GoF_p[i]<-GoF[[i]]$chisquare$chi1$p
    }

    SelMod<-listDSmodels2[[which.max(GoF_p)]]

    if(any(is.na(GoF_chi)==T)) {
      if(length(which(is.na(GoF_chi)))==1){
       SelMod2 <- try(listDSmodels2[[which(GoF_p<max(GoF_p)&GoF_p>min(GoF_p))]])
      }
      if(length(which(is.na(GoF_chi)))==2){
        SelMod2 <- NA
      }
    } else{
    SelMod2<-try(listDSmodels2[[which(GoF_p<max(GoF_p)&GoF_p>min(GoF_p))]])
    SelMod3<-try(listDSmodels2[[which.min(GoF_p)]])}

    summarySelMod<-summary(SelMod)

  }



  output.combinations<-apply(output,1,function(x)paste0(x[1],x[2],x[3]))
  index<-which(output.combinations%in%combination)

  output[index,]$Group_sizeXdistance_cor_test_p <- round(cor_test$p.value,digits=2)
  output[index,]$DS_cutpoints<-paste0(SelMod$ddf$meta.data$breaks[2]-SelMod$ddf$meta.data$breaks[1],"m")
  output[index,]$Detection_probability <-round(summarySelMod$ds$average.p,digits=2)

  output[index,]$GoF_DS_Xsq <- GoF_chi[which.max(GoF_p)]

  output[index,]$`N_(count/p)`  <- round(output[index,]$Count/output[index,]$Detection_probability,digits=2)
  output[index,]$`D_(N/Km)`<-round(output[index,]$`N_(count/p)`/output[index,]$Km,digits=2)
  output[index,]$`Encounters_(<=500m)` <-nrow(df)

  if(summarySelMod$ds$key=="unif"){
    output[index,]$Detection_function <-summarySelMod$ddf$name.message
    output[index,]$Group_size_beta <- "-"
  }

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
  mtext(paste("Spearman's rank correlation: p =",round(cor_test$p.value,2) ),cex=0.7,side=3,line=-2,adj=0.95)

  try(chi1<-round(GoF_chi[which.max(GoF_p)],4))
  try(p1<-round(max(GoF_p),4))
  try(plot(SelMod,breaks=SelMod$ddf$meta.data$breaks,main=bquote(paste(.(paste0(df$Species[1]," - ",df$Ecoregion[1]," - ",df$Pais[1]," ",paste0("(","N = ",nrow(df),", ","key function = ",SelMod$call$key,")")))))))
  try(mtext(bquote(paste("GoF: ",chi^2 == .(chi1),", ","p = ",.(p1))),cex=0.7,side=3,line=-2,adj=0.95))

  try(chi2<-round(GoF_chi[which(GoF_p<max(GoF_p)&GoF_p>min(GoF_p))],4))
  try(p2<- round(GoF_p[which(GoF_p<max(GoF_p)&GoF_p>min(GoF_p))],4))
  try(plot(SelMod2,breaks=SelMod2$ddf$meta.data$breaks,main=bquote(paste(.(paste0((df$Species[1])," - ",df$Ecoregion[1]," - ",df$Pais[1]," ",paste0("(","N = ",nrow(df),", ","key function = ",SelMod2$call$key,")")))))))
  try(mtext(bquote(paste("GoF: ",chi^2 == .(chi2),", ","p = ",.(p2))),cex=0.7,side=3,line=-2,adj=0.95))

  try(chi3<-round(GoF_chi[which.min(GoF_p)],4))
  try(p3<- round(min(GoF_p),4))
  try(plot(SelMod3,breaks=SelMod3$ddf$meta.data$breaks,main=bquote(paste(.(paste0(df$Species[1]," - ",df$Ecoregion[1]," - ",df$Pais[1]," ",paste0("(","N = ",nrow(df),", ","key function = ",SelMod3$call$key,")")))))))
  try(mtext(bquote(paste("GoF: ",chi^2 == .(chi3),", ","p = ",.(p3))),cex=0.7,side=3,line=-2,adj=0.95))

  dev.off()

}

# Run the function
pblapply(list1,DSmodels,output.folder="/Users/franciscodenes/Documents/GitHub/DS_parrots_2021/review/")





#troubleshoot
df.all<-list1[[21]]
max.dist=500
cutpoints=NULL

# assessing frequency of key function in selected models for each combination:
output_table <- read.csv("~/Documents/GitHub/DS_parrots_2021/review/output_table.csv")
c<-table(output_table$Detection_function,output_table$Encounters_...500m.)
d<-as.data.frame(c)
e<-d[order(d$Var1),]
e_hn<-e[e$Var1=="half-normal key function",]
e_hz<-e[e$Var1=="hazard-rate key function",]

chisq.test( table(output_table$Detection_function,output_table$Encounters_...500m.))


plot(e_hz$Freq~e_hn$Var2,xlab="Encounters",ylab="Frequency as best model",main="Hazard rate key function")
plot(e_hn$Freq~e_hn$Var2,xlab="Encounters",ylab="Frequency as best model",main="Half-normal key function")

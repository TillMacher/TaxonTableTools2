# ==============================================================================
# Title: Empty R Script
# Description: Empty R Script
# Authors:
#   Pierre Bady <pierre.bady@free.fr>
#   Maxime Logez <maxime.logez@cemagref.fr>
# Date: 09/02/08 10:36:20
# Revision: 03/12/09 14:39:21
# Version: 1.0
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2008-2009  Pierre BADY
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# ==============================================================================
require(MASS)
#------------------------------------------------
# additional functions
# 09/02/08 10:10:38 by pbady
#------------------------------------------------
# 11/25/08 16:25:42 by mlogez
CheckSpNames <- function(x){
  # fonction qui a pour but d'homog?n?iser les noms des esp?ces
  # elle supprime les blancs des d?buts et fin de cha?ne
  # elle remplace les points entre deux cha?nes de caract?res par des blancs
  # elle remplace les blancs multiples au sein de la cha?ne de caract?re par
  # un seul blanc
  # elle passe la premi?re lettre de la cha?ne de cacract?re en majuscule
  # et toutes les autres lettres en minuscule
  # et elle passe la premi?re lettre de la deuxi?me esp?ce en majuscule
  # dans le cas des hybrides si les deux noms d'esp?ces sont s?par?s par un
  # 'x' pr?c?d? et suivi d'un espace ou d'un .
  # exemple : c("  Cottus.GOBIO "," cOttus   Gobio", "    cOttus     Gobio.x  SaLmo trutta Fario    ")
  # devient : "Cottus gobio" "Cottus gobio" "Cottus gobio x Salmo trutta fario"
  if (!is.character(x))
    x <- as.character(x)
  # supression des blancs de d?but et fins
  x <- gsub("(^\\s+|\\s+$)","",x,perl=TRUE)
  # remplace les points par des blancs, remplace les blancs multiples au sein de la cha?ne de caract?re par un seul blanc
  x <- gsub("(?:(?<=\\w))(\\.+|\\s+)(?=\\w+)"," ",x,perl=TRUE)
  #  passe toutes les lettres en minuscules
  x <- tolower(x)
  # remplace la premi?re lettre de l'esp?ce par sa majuscule
  # prend en charge les hybride a condition qu'il soit ?crit de cette mani?re : Esp1 x Esp2
  x <- gsub("(?<=^| x )(\\w{1})","\\U\\1",x,perl=TRUE)
  x
}
# 09/02/08 10:10:38 by pbady
unzload <- function(dirzip,filename){
  con <- unz(dirzip,filename,open="r")
  load(con,envir = .GlobalEnv)
  close(con)
}
ftable2df <- function(x,...){
  if(inherits(file,"ftable"))
    stop("'x' must be an \"ftable\" object")
  at1 <- attributes(x)
  nl <- at1$dim[1]
  nc <- at1$dim[2]
  mat1 <- matrix(x,nl,nc)
  niv <- unlist(lapply(at1$row.vars,function(x) length(x)))
  cump1 <- as.vector(cumprod(niv))
  cump2 <- as.vector(cumprod(niv[length(niv):1]))
  k <- c(length(niv):2)-1
  fcol <- NULL
  for(i in length(niv):1) {
    if(i < length(niv)) {
      w <- as.vector(sapply(unlist(at1$row.vars[i]),function(x) rep(x,cump2[k[i]])))
    }
    else w <- unlist(at1$row.vars[i])
    if(i == 1)
      vec <- w
    else  vec <- rep(w,cump1[i-1])
    fcol <- cbind(vec,fcol)
  }
  row.names(fcol)<- 1:nl
  res <- cbind(as.data.frame(fcol),as.data.frame(mat1))
  names(res) <- c(names(niv),(at1$col.vars[[1]]))
  return(res)
}
xtabs2df <- function (x){
    if (!inherits(x, "xtabs"))
        stop("non convenient argument")
    x <- as.data.frame(unclass(x))
    return(x)
}
try.test <- function(x) !inherits(x, "try-error")
rescale.names <- function(tab,ref,method ="row.row"){
  if(!inherits(tab,"data.frame"))
    stop("non convenient argument")
  if(!inherits(ref,"data.frame"))
    stop("non convenient argument")
  if(method =="row.row"){
        vec <- intersect(row.names(ref),row.names(tab))
        return(tab[vec,])
  }else if(method =="row.col"){
        vec <- intersect(names(ref),row.names(tab))
        return(tab[vec,])
  }else if(method =="col.col"){
        vec <- intersect(names(ref),names(tab))
        return(tab[,vec])
  }else if(method =="col.row"){
        vec <- intersect(row.names(ref),names(tab))
        return(tab[,vec])
  }else stop("non convenient method!")
}
scale01 <- function (x){
  xx <- as.matrix(x)
  x01 <- as.numeric(xx > 0)
  dim(x01) <- dim(xx)
  x01 <- as.data.frame(x01)
  names(x01) <- names(x)
  row.names(x01) <- row.names(x)
  return(x01)
}
### see the library ade4 (Chessel et al.)
acm.util <- function(x) {
  rnames <- names(x)
  cl <- as.factor(x)
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (unclass(cl) - 1)] <- 1
  dimnames(x) <- list(rnames(x),levels(cl))
  return(x)
}
### see the library ade4 (Chessel et al.)
acm.disjonctif <-function (df){
    acm.util <- function(i) {
        cl <- df[, i]
        cha <- names(df)[i]
        n <- length(cl)
        cl <- as.factor(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1:n) + n * (unclass(cl) - 1)] <- 1
        dimnames(x) <- list(row.names(df), paste(cha, levels(cl),
            sep = "."))
        return(x)
    }
    G <- lapply(1:ncol(df), acm.util)
    G <- data.frame(G, check.names = FALSE)
    return(G)
}
#------- functions "Replace" -------------
# argument class: data.frame
Replace <- function (x,old,new){
    UseMethod("Replace")
}
Replace.data.frame <- function(x,old,new){
  if(!inherits(x,"data.frame"))
    stop("non convenient argument")
  f1 <- function(x,old,new){
   if(is.na(old)){
    if(is.factor(x)){
      x <- as.character(x)
      x[is.na(x)] <- rep(new,length(x[is.na(x)]))
      x <- as.factor(x)
    }else
      x[is.na(x)] <- rep(new,length(x[is.na(x)]))
   return(x)
   }else {
    if(is.factor(x)){
      x <- as.character(x)
      x[x==old & !is.na(x)] <- rep(new,length(x[x==old & !is.na(x)]))
      x <- as.factor(x)
    }else
      x[x==old & !is.na(x)] <- rep(new,length(x[x==old & !is.na(x)]))
   return(x)
   }
  }
  tab <- data.frame(lapply(x,function(y) f1(y,old,new)))
  names(tab) <- names(x)
  row.names(tab) <- row.names(x)
  return(tab)
}
# argument class: factor
Replace.factor <- function(x,old, news){
  if(!inherits(x,"factor"))
    stop("non convenient argument")
  x <- as.character(x)
  if(is.na(old)){
    ted <- is.na(x)
    nx <- length(x[ted])
    x[ted] <- rep(news,nx)
    x <- as.factor(x)
  }else{
    ted <- x==old & !is.na(x)
    nx <- length(x[ted])
    x[ted] <- rep(news,nx)
    x <- as.factor(x)
  }
  return(x)
}
# argument class: numeric or character
Replace.default <- function(x,old, news){
  if(!any(!inherits(x,"numeric"),!inherits(x,"character")))
      stop("non convenient argument")
  if(is.na(old)){
    ted <- is.na(x)
    nx <- length(x[ted])
    x[ted] <- rep(news,nx)
  }else{
    ted <- x==old & !is.na(x)
    nx <- length(x[ted])
    x[ted] <- rep(news,nx)
  }
  return(x)
}
#------------------------------------------------
# computation of observed values
# 09/02/08 10:06:43 by pbady
#------------------------------------------------
# check ok, ? revoir pour optimiser le code
# 09/02/08 18:16:08
MetricCompute <- function(faun,guild,labels,type,Tol=1e-7){
  if(!inherits(faun,"data.frame"))
    stop("non convenient argument")
  if(!inherits(guild,"data.frame"))
    stop("non convenient argument")
  f2 <- function(x,tol) ifelse(abs(x) < tol & !is.na(x),0,x)
# definition of ST-species
  salmospecies <- c("Alburnoides.bipunctatus",
    "Cobitis.calderoni",
    "Coregonus.lavaretus",
    "Cottus.gobio",
    "Cottus.poecilopus",
    "Eudontomyzon.mariae",
    "Hucho.hucho",
    "Lampetra.planeri",
    "Phoxinus.phoxinus",
    "Salmo.salar",
    "Salmo.trutta.fario",
    "Salmo.trutta.lacustris",
    "Salmo.trutta.macrostigma",
    "Salmo.trutta.trutta",
    "Salmo.trutta.marmoratus",
    "Salvelinus.fontinalis",
    "Salvelinus.namaycush",
    "Salvelinus.umbla",
    "Thymallus.thymallus")
# computations
  faun <- rescale.names(faun,guild,"col.row")
  guild <- rescale.names(guild,faun,"row.col")
  guild <- acm.disjonctif(guild)
  if(!missing(labels))
    guild <- guild[,labels]
  faun01 <- scale01(faun)
  den <- apply(faun,1,sum)
  ric <- apply(faun01,1,sum)
  gfaunA <- as.data.frame(as.matrix(faun)%*%as.matrix(guild))
  gfaunR <- as.data.frame(as.matrix(faun01)%*%as.matrix(guild))
  names(gfaunA) <- paste("dens",names(gfaunA),sep=".")
  names(gfaunR) <- paste("ric",names(gfaunR),sep=".")
  w <- cbind(gfaunA,gfaunR)
  if(!missing(type))
    w <- w[,paste(type,labels,sep=".")]
  w <- as.data.frame(apply(w,2,function(x,y=Tol) f2(x,tol=y)))
  psalmo <- apply(faun[,names(faun)%in%salmospecies],1,sum)
  psalmo <- psalmo/den
  attr(w,"Richness") <- ric
  attr(w,"Density") <- den
  attr(w,"Psalmo") <- psalmo
  attr(w,"Nbmetric") <- ncol(guild)
  return(w)
}
#------------------------------------------------
# preparation of the faunistic table
# 11/25/08 14:22:42 by pbady
#------------------------------------------------
# the names of variables must be verified !!
# 11/25/08 14:27:06
prepTab <- function(df,param=c("code","species","number")){
  if(!inherits(df,'data.frame'))
    stop("non convenient argument !")
  df <- df[,param]
# need to modify the formula for more genericity !
# 11/26/08 14:57:31
  names(df) <- c("code","species","number")
  df <- xtabs2df(xtabs(number ~ code + species, data=df))
  return(df)
}
ValidatedLength <- function(target,complement,obs,Tol = 1e-7){
  valid1 <- abs((target+complement)-obs[row.names(complement),names(complement)]) < Tol  # &  obs[,names(complement)]) > 0
# pour le sup > 0 voir avec max et didier
# 11/28/08 13:44:16
  target <- as.data.frame(sapply(1:ncol(target),function(x) ifelse(valid1[,x],target[,x],NA)))
  names(target) <- names(complement)
  return(target)
}
#------------------------------------------------
# preparation of the environmental variables
# 11/25/08 14:22:42 by pbady
#------------------------------------------------
##### nom de variable ? v?rifier
prepEnv <- function(df,param=c("code","Actual.river.slope","temp.jul","temp.jan","Floodplain.site","Water.source.type",
  "Geomorph.river.type","Distance.from.source","AREA.ctch","Natural.sediment","Fished.area","Method",
  "Sampling.location"),richness,abundance,psalmo,geomodels,nameOrder,...)
{
  if(!inherits(df,'data.frame'))
    stop("non convenient argument !")
  df <- df[,param]
  df <-  df[!duplicated(df$code),]
  row.names(df) <- df$code
  if(!missing(nameOrder))
    df <- df[nameOrder,]
  # slope
  df$lslope <- ifelse(df$Actual.river.slope < 0.002,log(0.001),log(df$Actual.river.slope))
  # Tdif
  df$Tjul <- df$temp.jul
  df$Tjan <- df$temp.jan
  df$Tdif <- df$Tjul-df$Tjan
  # natsed
  df$natsed <- as.factor(c("large","medium",rep("small",3))[match(df$Natural.sediment,
    c("Boulder/Rock","Gravel/Pebble/Cobble","Organic","Sand","Silt"))])
  # fldpl
  df$fldpl <- df$Floodplain.site
  # watersource
  df$watersource <- as.factor(c("Non-Pluvial","Non-Pluvial","Non-Pluvial",
    "Pluvial")[match(df$Water.source.type,c("Glacial","Groundwater","Nival","Pluvial"))])
  # geomorph
  df$geomorph <- as.factor(c("braided","meand","meand","constraint","sinuous")[match(df$Geomorph.river.type,
    c("Braided","Meand regular","Meand tortous","Naturally constraint no mob","Sinuous"))])
  # method   "Boat"   "Mixed"  "Wading"
  df$method <- as.factor(c("boat","noboat","noboat")[match(df$Method,c("Boat","Mixed","Wading"))])
  # dist
  df$ldist <- ifelse(df$Distance.from.source < 0.5,log(0.05),log(df$Distance.from.source))
  # drainage area
  df$lDR <- ifelse(df$AREA.ctch < 0.5,log(0.05),log(df$AREA.ctch))
  # df$lDR <- ifelse(df$Size.of.catchment < 0.5,log(0.05),log(df$Size.of.catchment))
  syngeomorph1 <- predict(geomodels[[1]],newdata=df)
  syngeomorph2 <- predict(geomodels[[2]],newdata=df)
  df <- data.frame(df,syngeomorph1=syngeomorph1,syngeomorph2=syngeomorph2,richesse=richness,captures=abundance,psalmo=psalmo)
  class(df) <- c("env","data.frame")
  return(df)
}
#------------------------------------------------
# Regionalisation
# 11/26/08 10:57:11
# we need two variables: EFT (Andy's typologie,
# see FAME project) and EcoRegion (ECOREG).
# D.PONT, M.LOGEZ & P.BADY
# 11/26/08 11:41:04
#------------------------------------------------
prepRegEcoOLD <- function(df,...){
  # ECOREG1
  ECOREG1<- as.character(df$ECOREG)
  ECOREG1 <- c("Alp","Bal","Bor","C.h","C.p","E.p","Eng","Fen","Hun","Ibe","Ita",
    "Pon","Pyr","Car","W.h","W.p")[match(ECOREG1,sort(unique(ECOREG1)))]
  # correction for mediterranean sites
  ECOREG1<-ifelse(df$MEDIT==1 & df$Latitude<45 & ECOREG1=="Ita","Med",ECOREG1)
  ECOREG1<-ifelse(df$MEDIT==1 & df$Latitude<45 & ECOREG1=="W.p","Med",ECOREG1)
  ECOREG1<-ifelse(df$MEDIT==1 & df$Latitude<45 & ECOREG1=="Ibe","Med",ECOREG1)
  # ECOREG2
  ECOREG2<-as.character(ECOREG1)
  ECOREG2<-ifelse(ECOREG1=="Alp","Alp",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Pyr","Alp",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="E.p","Est",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Pon","Est",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Hun","Est",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Fen","Nor",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Bor","Nor",ECOREG2)
  # EFT.new
  EFT.new<-rep("Salmonid",length(ECOREG2))
  EFT.new<-ifelse(ECOREG2=="Est","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.5","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.6","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.14","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.15","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$MEDIT==1,"Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.NA",NA,EFT.new)
  # EFT <- ifelse(EFT.new!="Trout" & !is.na(EFT.new),"NoTrout",EFT.new)
  # results
# modification of the varable which indicates the functional type (EFT.new)
  EFT2 <- EFT.new
  EFT2 <- ifelse(ECOREG2=="Alp","Salmonid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Bal","Cyprinid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Est","Cyprinid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Med","Cyprinid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Nor","Salmonid",EFT2)
  EFT2 <- ifelse(is.na(EFT.new),NA,EFT2)
  EFT2 <- ifelse(EFT.new=="NA",NA,EFT2)
# data preparation
  res <- data.frame(EFT.Type=df$EFT.Type,EFT.new,ECOREG=df$ECOREG,ECOREG1,ECOREG2,EFT2)
  row.names(res) <- row.names(df)
  res
}
prepRegEco <- function(df,...){
  # ECOREG1
  ECOREG1<- as.character(df$ECOREG)
  ECOREG1 <- c("Alp","Bal","Bor","C.h","C.p","E.p","Eng","Fen","Hun","Ibe","Ita",
    "Pon","Pyr","Car","W.h","W.p")[match(ECOREG1,c("Alps","Baltic province",
 "Borealic uplands","Central highlands","Central plains","Eastern plains",
 "England","Fenno-scandian shield", "Hungarian lowlands","Ibero-Macaronesian region",
 "Italy and Corsica","Pontic province","Pyrenees","The Carpathiens","Western highlands",
 "Western plains"))]
   # correction for mediterranean sites
  ECOREG1<-ifelse(df$MEDIT==1 & df$Latitude<45 & ECOREG1=="Ita","Med",ECOREG1)
  ECOREG1<-ifelse(df$MEDIT==1 & df$Latitude<45 & ECOREG1=="W.p","Med",ECOREG1)
  ECOREG1<-ifelse(df$MEDIT==1 & df$Latitude<45 & ECOREG1=="Ibe","Med",ECOREG1)
  # ECOREG2
  ECOREG2<-as.character(ECOREG1)
  ECOREG2<-ifelse(ECOREG1=="Alp","Alp",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Pyr","Alp",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="E.p","Est",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Pon","Est",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Hun","Est",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Fen","Nor",ECOREG2)
  ECOREG2<-ifelse(ECOREG1=="Bor","Nor",ECOREG2)
  # EFT.new
  EFT.new<-rep("Salmonid",length(ECOREG2))
  EFT.new<-ifelse(ECOREG2=="Est","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.5","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.6","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.14","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.15","Cyprinid",EFT.new)
  EFT.new<-ifelse(df$MEDIT==1,"Cyprinid",EFT.new)
  EFT.new<-ifelse(df$EFT.Type=="T.NA",NA,EFT.new)
  # EFT <- ifelse(EFT.new!="Trout" & !is.na(EFT.new),"NoTrout",EFT.new)
  # results
# modification of the varable which indicates the functional type (EFT.new)
  EFT2 <- EFT.new
  EFT2 <- ifelse(ECOREG2=="Alp","Salmonid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Bal","Cyprinid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Est","Cyprinid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Med","Cyprinid",EFT2)
  EFT2 <- ifelse(ECOREG2=="Nor","Salmonid",EFT2)
  EFT2 <- ifelse(is.na(EFT.new),NA,EFT2)
  EFT2 <- ifelse(EFT.new=="NA",NA,EFT2)
# modification temporaire 03/05/09 16:24:17
# test version Didier
# ECOREG2<-ifelse(ECOREG2=="Ita" & EFT2=="Cyprinid","W.p",ECOREG2)
# ECOREG2<-ifelse(ECOREG2=="W.h" & EFT2=="Cyprinid","W.p",ECOREG2)
# ECOREG2<-ifelse(ECOREG2=="C.h" & EFT2=="Cyprinid","C.p",ECOREG2)
# data preparation
  res <- data.frame(EFT.Type=df$EFT.Type,EFT.new,ECOREG=df$ECOREG,ECOREG1,ECOREG2,EFT2)
  row.names(res) <- row.names(df)
  res
}
#------------------------------------------------
# computation of expected values
# 09/02/08 10:06:43 by pbady
#------------------------------------------------
# check ok 09/02/08 17:19:13
ExpectedMetric  <- function(env,models){
  if(!inherits(env,"env"))
    stop("object 'env' expected !")
  p <- length(models)
  n <- nrow(env)
  res <- matrix(NA,nc=p,nr=n)
  for(i in 1:p){
    res[,i] <- predict(models[[i]],newdata=env,type="response")
  }
  res <- as.data.frame(res)
  names(res) <- names(models)
  row.names(res) <- row.names(env)
  return(res)
}
#------------------------------------------------
# computation of standardized distances
# between expected and observed values
# 12/03/08 16:26:09 by pbady and mlogez
# pour v?rifier
# 12/03/08 17:38:17 by pbady
# dans efitools-0.1
#------------------------------------------------
DistMetricOLD <- function(obs,fitted,fac,center,scale){
  if(!inherits(obs,"data.frame"))
    stop("object 'data.frame' expected !")
  if(!inherits(fitted,"data.frame"))
    stop("object 'data.frame' expected !")
  obs <- log(obs+1)
  fitted <- log(fitted+1)
  indrow <- match(fac,row.names(center))
  indcol <- match(names(obs),names(center))
  Mcenter <- center[indrow,indcol]
  res <- obs-fitted
  res <- as.matrix(res-Mcenter)
  res <- sweep(res, 2, scale, "/", check.margin = FALSE)
  res <- as.data.frame(res)
  names(res) <- names(obs)
  row.names(res) <- row.names(obs)
  class(res) <- c("distmet",class(res))
  return(res)
}
# modification after modification
# of the guild table
# 03/12/09 14:15:44
DistMetric <- function(obs,fitted,fac,center,scale,labels=NULL){
  if(!inherits(obs,"data.frame"))
    stop("object 'data.frame' expected !")
  if(!inherits(fitted,"data.frame"))
    stop("object 'data.frame' expected !")
  if(!is.null(labels)){
    fitted <- fitted[,labels]
    obs <- obs[,labels]
    scale <-  scale[labels]
    center <-  center[,labels]
    }
  obs <- log(obs+1)
  fitted <- log(fitted+1)
  indrow <- match(fac,row.names(center))
  indcol <- match(names(obs),names(center))
  Mcenter <- center[indrow,indcol]
  res <- obs-fitted
  res <- as.matrix(res-Mcenter)
  res <- sweep(res, 2, scale, "/", check.margin = FALSE)
  res <- as.data.frame(res)
  names(res) <- names(obs)
  row.names(res) <- row.names(obs)
  class(res) <- c("distmet",class(res))
  return(res)
}
#------------------------------------------------
# computation of individual score
# 12/03/08 18:35:12 by pbady
# ? modifier !!! 01/07/09 15:49:00
#------------------------------------------------
IDScoring <- function(distmet,stdpack,labels){
  if(!inherits(distmet,"distmet"))
    stop("non convienent argument!")
  if(!inherits(stdpack,"stdpack"))
    stop("non convienent argument!")
  limits <- stdpack$limits[labels]
  object <- distmet[,labels]
  if(is.null(ncol(object))){
    res <- punif(object,limits[1,labels],limits[2,labels])
    names(res) <- row.names(distmet)
  }else{
    res <- data.frame(sapply(1:ncol(object),function(x) punif(object[,x],limits[1,x],limits[2,x])))
    names(res) <- labels
    row.names(res) <- row.names(distmet)
    }
  qclass <- stdpack$qclass
  resclass <- ClassifIDScoring(res,qclass[,labels],labels=labels)
  res <- cbind(as.data.frame(res),resclass)
  names(res) <- c(labels,paste(labels,"class",sep="."))
  class(res) <- c("idscore",class(res))
  return(res)
}
####### 01/07/09 15:44:25
ClassifIDScoring <- function(object,qclassobject,labels=names(object),...){
  if(inherits(object,"data.frame")){
    if(length(labels) >1){
      wclass <- as.data.frame(sapply(labels,
        function(x) 6-as.numeric(cut(object[,x],qclassobject[,x],right=FALSE))))
      row.names(wclass) <- row.names(object)
      names(wclass) <- labels
    }else if(length(labels) == 1){
      wclass <- 6-as.numeric(cut(object[,labels],qclassobject[,labels],right=FALSE))
      names(wclass) <- row.names(object)
    }else stop("non convenient dimension !")
  }else if(inherits(object,"numeric")){
      wclass <- 6-as.numeric(cut(object,qclassobject,right=FALSE))
      names(wclass) <- names(object)
  }else stop("non convenient argument !")
  return(wclass)
}
#----------------------------------------------------
# aggregation based on the mean of individual score
# 12/10/08 17:32:10 by pbady
# the both Index are always computed if it's possible
# 01/30/09 10:59:32 by pbady
#----------------------------------------------------
IDaggregation <- function(object,qclassMean,...){
  if(!inherits(object,"data.frame"))
    stop("non convienent argument!")
  object$EFT3 <- recodeRiverZone(object$EFT2,object$psalmo,object$ECOREG2)
  object$comments.river.zone <- commentRiverZone(object$EFT3,object$psalmo)
  object$comments.fish.index <- commentFishIndex(object$EFT3,object$method)
  object$comments.sampling.effort <- commentSamplingEffort(object$captures)
  object$comments.sampling.location <- commentSamplingLocation(object$Sampling.location)
# indices for Salm and Cypr zones
  object$IndexSalm <- ifelse(!is.na(object$EFT3),
    (object$dens.WQO2.O2INTOL + object$dens.HTOL.HINTOL150)/2,NA)
  object$IndexCypr <- ifelse(!is.na(object$EFT3),
    (object$ric.HabSp.RHPAR + object$dens.Repro.LITH)/2,NA)
# classes for individual index
  object$IndexSalm.class <- ClassifIDScoring(object$IndexSalm,qclassMean[,"Salmonid"])
  object$IndexCyprboat.class <- ClassifIDScoring(object$IndexCypr,qclassMean[,"Cyprinid.boat"])
  object$IndexCyprnoboat.class <- ClassifIDScoring(object$IndexCypr,qclassMean[,"Cyprinid.noboat"])
# global index
  object$IndexMean<- rep(NA,nrow(object))
  object$IndexMean <- ifelse(object$EFT3=="Salmonid" & !is.na(object$EFT3),object$IndexSalm,object$IndexMean)
  object$IndexMean <- ifelse(object$EFT3=="Cyprinid" & !is.na(object$EFT3),object$IndexCypr,object$IndexMean)
# global index en classes
  object$IndexMean.class <- rep(NA,nrow(object))
  object$IndexMean.class <- ifelse(object$EFT3=="Salmonid" & !is.na(object$EFT3) & !is.na(object$method),
    object$IndexSalm.class,object$IndexMean.class)
  object$IndexMean.class <- ifelse(object$EFT3=="Cyprinid" & object$method=="boat"& !is.na(object$EFT3) & !is.na(object$method),
    object$IndexCyprboat.class,object$IndexMean.class)
  object$IndexMean.class <- ifelse(object$EFT3=="Cyprinid" & object$method=="noboat" & !is.na(object$EFT3) & !is.na(object$method),
    object$IndexCyprnoboat.class,object$IndexMean.class)
# observed and fitted data in /m2
  obspred1 <- c("dens.HTOL.HINTOL150.obs","dens.WQO2.O2INTOL.obs","dens.Repro.LITH.obs",
    "dens.HTOL.HINTOL150.pred","dens.WQO2.O2INTOL.pred","dens.Repro.LITH.pred")
  obspred <- c("dens.HTOL.HINTOL150.obs","dens.WQO2.O2INTOL.obs","ric.HabSp.RHPAR.obs","dens.Repro.LITH.obs",
    "dens.HTOL.HINTOL150.pred","dens.WQO2.O2INTOL.pred","ric.HabSp.RHPAR.pred","dens.Repro.LITH.pred")
  for(i in obspred1)
    object[,i] <- 100*object[,i]/object$Fished.area
# outputnames preparation
  outputnames <- c("Sample.code","Day","Month","Year","Longitude","Latitude","Fished.area",obspred,
    "dens.HTOL.HINTOL150","dens.WQO2.O2INTOL","ric.HabSp.RHPAR","dens.Repro.LITH",
    "dens.HTOL.HINTOL150.class","dens.WQO2.O2INTOL.class","ric.HabSp.RHPAR.class","dens.Repro.LITH.class",
    "Method","richesse","captures","comments.sampling.effort","ECOREG2","EFT3","psalmo","comments.river.zone",
    "IndexSalm","IndexCypr","IndexMean","IndexMean.class","comments.fish.index",
    "Sampling.location","comments.sampling.location")
  object <- object[,outputnames]
  return(object)
}
##
RecodeNameOutput <- function(object,...){
  if(!inherits(object,"data.frame"))
    stop("non convienent argument!")
  obspred <- c("dens.HTOL.HINTOL150.obs","dens.WQO2.O2INTOL.obs","ric.HabSp.RHPAR.obs","dens.Repro.LITH.obs",
    "dens.HTOL.HINTOL150.pred","dens.WQO2.O2INTOL.pred","ric.HabSp.RHPAR.pred","dens.Repro.LITH.pred")

  outputnames <- c("Sample.code","Day","Month","Year","Longitude","Latitude",
    obspred,"dens.HTOL.HINTOL150","dens.WQO2.O2INTOL","ric.HabSp.RHPAR","dens.Repro.LITH",
    "Method","Sampling.location","comments.sampling.location","captures","comments.sampling.effort",
    "ECOREG2","EFT3","psalmo","comments.river.zone","IndexSalm","IndexCypr","IndexMean",
    "IndexMean.class","comments.fish.index")
  object <- object[,outputnames]

  names(object)[7:32] <- c("Obs.dens.HINTOL.inf.150","Obs.dens.O2INTOL","Obs.ric.RH.PAR","Obs.dens.LITH",
    "Exp.dens.HINTOL.inf150","Exp.dens.O2INTOL","Exp.ric.RH.PAR","Exp.dens.LITH",
    "Ids.dens.HINTOL.inf.150","Ids.dens.O2INTOL","Ids.ric.RH.PAR","Ids.dens.LITH",
    "Method","Sampling.location","comments.sampling.location","Captures","comments.sampling.effort",
    "Ecoregion","EFT.river.typology","ST-Species","Comments.river.zone","Aggregated.score.Salmonid.zone",
    "Aggregated.score.Cyprinid.zone","FishIndex","FishIndex.class","Comments.fish.index")
  return(object)
}
RecodeNameOutputDiad <- function(object,...){
  if(!inherits(object,"data.frame"))
    stop("non convienent argument!")
  obspred <- c("dens.HTOL.HINTOL150.obs","dens.WQO2.O2INTOL.obs","ric.HabSp.RHPAR.obs","dens.Repro.LITH.obs",
    "dens.HTOL.HINTOL150.pred","dens.WQO2.O2INTOL.pred","ric.HabSp.RHPAR.pred","dens.Repro.LITH.pred")

  outputnames <- c("Sample.code","Day","Month","Year","Longitude","Latitude",
    obspred,"dens.HTOL.HINTOL150","dens.WQO2.O2INTOL","ric.HabSp.RHPAR","dens.Repro.LITH",
    "Method","Sampling.location","comments.sampling.location","captures","comments.sampling.effort",
    "ECOREG2","EFT3","psalmo","comments.river.zone","IndexSalm","IndexCypr","IndexMean",
    "IndexMean.class","comments.fish.index")
  res <- object[,outputnames]

  names(res)[7:32] <- c("Obs.dens.HINTOL.inf.150","Obs.dens.O2INTOL","Obs.ric.RH.PAR","Obs.dens.LITH",
    "Exp.dens.HINTOL.inf150","Exp.dens.O2INTOL","Exp.ric.RH.PAR","Exp.dens.LITH",
    "Ids.dens.HINTOL.inf.150","Ids.dens.O2INTOL","Ids.ric.RH.PAR","Ids.dens.LITH",
    "Method","Sampling.location","comments.sampling.location","Captures","comments.sampling.effort",
    "Ecoregion","EFT.river.typology","ST-Species","Comments.river.zone","Aggregated.score.Salmonid.zone",
    "Aggregated.score.Cyprinid.zone","FishIndex","FishIndex.class","Comments.fish.index")
  res$Hist.ric.diadromous  <- object$nactual
  res$Present.ric.diadromous  <- object$nhistoric
  res$Ids.ric.diadromous  <- object$migmet
  return(res)
}
#-------------------------------------------
# commentaires
# 03/02/09 15:21:10 by pbady
#-------------------------------------------
recodeRiverZone <- function(zone,STspecies,ecoreg){
# correction related to % of salmonid
  csalmo <- cut(STspecies,c(-0.1,0.20,0.50,0.80,1.01))
  w <- rep(NA,length(zone))
  w <- ifelse(zone=="Salmonid" & !is.na(zone),c("Salmonid","Salmonid","Salmonid",
    "Salmonid")[match(as.character(csalmo),c("(-0.1,0.2]","(0.2,0.5]","(0.5,0.8]","(0.8,1.01]"))],w)
  w <- ifelse(zone=="Cyprinid" & !is.na(zone),c("Cyprinid","Salmonid","Salmonid",
    "Salmonid")[match(as.character(csalmo),c("(-0.1,0.2]","(0.2,0.5]","(0.5,0.8]","(0.8,1.01]"))],w)
# correction related to ecoregion
  w<-ifelse(ecoreg=="Alp" & !is.na(ecoreg),"Salmonid",w)
  w<-ifelse(ecoreg=="Nor" & !is.na(ecoreg),"Salmonid",w)
  w<-ifelse(ecoreg=="Est" & !is.na(ecoreg),"Cyprinid",w)
  w<-ifelse(ecoreg=="Bal" & !is.na(ecoreg),"Cyprinid",w)
  w<-ifelse(ecoreg=="Med" & !is.na(ecoreg),"Cyprinid",w)
  w
}
recodeRiverZone2 <- function(zone,STspecies,ecoreg){
# correction related to % of salmonid
  csalmo <- cut(STspecies,c(-0.1,0.20,0.50,0.80,1.01))
  w <- rep(NA,length(zone))
  w <- ifelse(zone=="Salmonid" & !is.na(zone),"Salmonid",w)
  w <- ifelse(zone=="Cyprinid" & !is.na(zone) & csalmo=="(0.8,1.01]","Salmonid",w)
  w <- ifelse(zone=="Cyprinid" & !is.na(zone) & csalmo=="(0.5,0.8]","Salmonid",w)
  w <- ifelse(zone=="Cyprinid" & !is.na(zone) & csalmo=="(0.2,0.5]","Salmonid",w)
  w <- ifelse(zone=="Cyprinid" & !is.na(zone) & csalmo=="(-0.1,0.2]","Cyprinid",w)
# correction related to ecoregion
  w<-ifelse(ecoreg=="Alp" & !is.na(ecoreg),"Salmonid",w)
  w<-ifelse(ecoreg=="Nor" & !is.na(ecoreg),"Salmonid",w)
  w<-ifelse(ecoreg=="Est" & !is.na(ecoreg),"Cyprinid",w)
  w<-ifelse(ecoreg=="Bal" & !is.na(ecoreg),"Cyprinid",w)
  w<-ifelse(ecoreg=="Med" & !is.na(ecoreg),"Cyprinid",w)
  w
}
commentRiverZone <- function(zone,STspecies){
  csalmo <- cut(STspecies,c(-0.1,0.20,0.50,0.80,1.01))
  w <- rep(NA,length(zone))
  w <- ifelse(zone=="Salmonid" & !is.na(zone),
    c(rep("To.be.checked.by.user",3),"Initial.classification.correct")[match(as.character(csalmo),
    c("(-0.1,0.2]","(0.2,0.5]","(0.5,0.8]","(0.8,1.01]"))],w)
  w <- ifelse(zone=="Cyprinid" & !is.na(zone),
    w <- c("Initial.classification.correct",
    rep("Initial.missclassification.possible.Modified river zone to.be.checked.by.user",3))[match(as.character(csalmo),
    c("(-0.1,0.2]","(0.2,0.5]","(0.5,0.8]","(0.8,1.01]"))],w)
  return(w)
}
commentFishIndex <- function(zone,method){
  w <- rep(NA,length(zone))
  w <- c("Sampled.by.wading/wading-boating.Fish Index adequate",
  "Sampled.by.boating.Fish Index.not.adequate.Only preliminary.corrected.result")[match(as.character(method),c("noboat","boat"))]
  return(w)
}
commentSamplingEffort <- function(NBind){
  w <- rep(NA,length(NBind))
  w <- ifelse(NBind < 30 & !is.na(NBind),"The number of fish caught is low. Results have to be used carefully!",w)
  w <- ifelse(NBind >= 30 & !is.na(NBind),"Nothing to report",w)
  return(w)
}
commentSamplingLocation <- function(SampLoc){
  w <- rep("No information !",length(SampLoc))
  w <- ifelse(SampLoc=="Backwaters" & !is.na(SampLoc),"results have to be used carefully!",w)
  w <- ifelse(SampLoc=="Mixed" & !is.na(SampLoc),"results have to be used carefully!",w)
  w <- ifelse(SampLoc=="Main channel" & !is.na(SampLoc),"Nothing to report",w)
  return(w)
}
#------------------------------------------
# metric based on diadromus species
# 09/08/08 16:38:15 by pbady
#------------------------------------------
getSimilarity <- function(df){
  if(!inherits(df,"data.frame"))
    stop("object 'data.frame' expected !")
  diadromousNames <- c("Acipenser.naccarii",
    "Acipenser.gueldenstaedti",
    "Acipenser.nudiventris",
    "Acipenser.stellatus",
    "Acipenser.sturio/A..oxyrinchus",
    "Alosa.alosa",
    "Alosa.fallax",
    "Alosa.immaculata",
    "Anguilla.anguilla",
    "Huso.huso",
    "Lampetra.fluviatilis",
    "Petromyzon.marinus",
    "Osmerus.eperlanus",
    "Salmo.salar",
    "Salmo.trutta.trutta",
    "Platichthys.flesus",
    "Coregonus.spp..diadr..form")
  actualnames <- paste(diadromousNames,"present",sep=".")
  historicalnames <- paste(diadromousNames,"historic",sep=".")
  nactual <- unlist(apply(df[,is.element(names(df),c(actualnames))],1,sum,na.rm=TRUE))
  nhistorical <- unlist(apply(df[,is.element(names(df),c(actualnames,historicalnames))],1,sum,na.rm=TRUE))
  migmet <- nactual/nhistorical
  res <- as.data.frame(cbind(nactual,nhistorical,migmet))
  names(res) <- c("nactual","nhistoric","migmet")
  row.names(res) <- row.names(df)
  return(res)
}
#------------------------------------------
# table merge
# 09/08/08 16:38:15 by pbady
#------------------------------------------
# change the names of the function addTab to addTable
addTab <- function(tab1,tab2){
  res <- merge(x=tab1,y=tab2, by="row.names",all=TRUE,sort = FALSE)
  if("Row.names"%in%names(res)){
    row.names(res) <- res$Row.names
    res <- res[,(!names(res)%in% "Row.names")]
    }
  return(res)
}
######################################## END ##################################























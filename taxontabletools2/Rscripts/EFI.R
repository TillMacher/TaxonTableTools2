#!/usr/bin/env Rscript

# ==============================================================================
# Title: main.R
# Description: main function of the computational R engine of EFI+ software
# Authors:
#   Pierre Bady <pierre.bady@free.fr>
#   Maxime Logez <maxime.logez@cemagref.fr>
#   Didier Pont <didier.pont@cemagref.fr>
# Date: 11/25/08 16:21:04
# Last modified: 03/12/09 14:39:21
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

# command to run the R script:
# R --vanilla --slave --args "input.txt" "output.txt" < main.R > info.out

# command to run the R script with diadromous species datafiles: Not Yet Implemented!
# R --vanilla --slave --args "input.txt" "output.txt" "diadromous.txt" < main.R > info.out

# print information
cat("==============================================================\n")
cat("Information Report\n")
cat("Date: ",as.character(Sys.time()),"\n")
cat("==============================================================\n\n")
cat("Authors:\n")
cat("Pierre Bady <pierre.bady@free.fr>\n")
cat("Maxime Logez <maxime.logez@cemagref.fr>\n")
cat("Didier Pont <didier.pont@cemagref.fr>\n\n")
cat("Created: 11/25/08 16:21:04\n")
cat("Last modified: 03/12/09 14:39:21\n")
cat("Version: 1.0\n")
cat("License: GPL version 2 or newer\n")
cat("Copyright (C) 2008-2009  Pierre BADY\n\n")
cat("This program is free software; you can redistribute it and/or\n")
cat("modify it under the terms of the GNU General Public License\n")
cat("as published by the Free Software Foundation; either version 2\n")
cat("of the License, or (at your option) any later version.\n\n")
cat("This program is distributed in the hope that it will be useful,\n")
cat("but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
cat("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
cat("GNU General Public License for more details.\n")
cat("\n==============================================================\n\n")
#--------------------------------------------------
# code Error
# 12/12/08 10:52:18 by pbady
#--------------------------------------------------
# Error -1: non convenient command!
# Error 0: non convenient arguments!
# Error 1: non convenient input file!
# Error 2: efiminus.R is absent!
# Error 3: internal.zip is absent!
# Error 4: non convenient input file!
# Error 5: non convenient diadromous input file!
# Error 6: non convenient output file!
#--------------------------------------------------

# first time point
t0 <- proc.time()
#----------------------------------------------------------
# definition of directory
#----------------------------------------------------------
# directory <- getwd()
# directory for AXP0328 (Cemagref desktop computer)
# directory <- "E:/DepotSVN/SVNefiminus/trunk/efitools-0.2"
# directory for BADOX (personal laptop computer)
# directory <- "C:/Documents and Settings/bad/Mes documents/data/SVNdepot/SVNefiminus/trunk/efitools-1.0"
# directory for triton (personal desktop computer)
# directory <- "H:/DepotSVN/SVNefiminus/trunk/efitools-0.2"

# directory for workstation antony
directory <- "/Users/tillmacher/Desktop/_dev/FEI+"
#--------------------------------------------------
# source loading (internal inputs and functions)
#--------------------------------------------------
if(file.exists(paste(directory,"/internal/efiminus.R",sep=""))){
  source(paste(directory,"/internal/efiminus.R",sep=""))
}else{
  print("Error 2: efiminus.R is absent!")
  stop("Error 2: efiminus.R is absent!")
}
if(file.exists(paste(directory,"/internal/internal.zip",sep=""))){
  unzload(paste(directory,"/internal/internal.zip",sep=""),"geomodels.rda")
  unzload(paste(directory,"/internal/internal.zip",sep=""),"models.rda")
  unzload(paste(directory,"/internal/internal.zip",sep=""),"lengthmodels.rda")
  unzload(paste(directory,"/internal/internal.zip",sep=""),"guildNew.rda")
  unzload(paste(directory,"/internal/internal.zip",sep=""),"stdpackNoTrout.rda")
  unzload(paste(directory,"/internal/internal.zip",sep=""),"stdpackTrout.rda")
  unzload(paste(directory,"/internal/internal.zip",sep=""),"stdpacklength.rda")
    unzload(paste(directory,"/internal/internal.zip",sep=""),"qclassMean.rda")
}else{
  print("Error 3: internal.zip is absent!")
  stop("Error 3: internal.zip is absent!")
}

#----------------------------------------------------------
# Handle command-line arguments
#----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values (for safety)
inputname <- NULL
outputname <- NULL

# Parse arguments
for (i in seq_along(args)) {
  if (args[i] == "-i" && i < length(args)) {
    inputname <- args[i + 1]
  }
  if (args[i] == "-o" && i < length(args)) {
    outputname <- args[i + 1]
  }
}

# Check if required arguments are provided
if (is.null(inputname) || is.null(outputname)) {
  stop("Usage: Rscript script.R -i <input_file> -o <output_file>")
}

#----------------------------------------------------------
# Check and install required packages
#----------------------------------------------------------
# Function to check and install missing packages
install_if_missing <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
}

# Set a CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# List of required packages
packages <- c("readxl", "writexl")

# Install missing packages
sapply(packages, install_if_missing)

input <- try(read_excel(inputname), silent = TRUE)
if (inherits(input, "try-error")) {
  print("Error 1: non convenient input file!")
  stop("Error 1: non convenient input file!")
}

# transformation of Missing Data in R format (NA)
input <- Replace(input,-999,NA)
input <- Replace(input,"NoData",NA)

# row code preparation
input$code <- paste(input$Sample.code,input$Day,input$Month,input$Year,sep="_")

# Print unique species names before applying the function
print(unique(input$Species))

# Identify species that become empty after applying CheckSpNames
result <- sapply(input$Species, CheckSpNames)
empty_indices <- which(result == "" | is.na(result))

# Print the problematic species
if (length(empty_indices) > 0) {
  print("Problematic species found:")
  print(input$Species[empty_indices])
} else {
  print("No problematic species found.")
}

# Species names checking
input$Species <- CheckSpNames(input$Species)
input$species <- gsub(" ",".",input$Species)

# minimal needed variables
input <- try(input[,c("Sample.code","Day","Month","Year","Longitude","Latitude",
  "Actual.river.slope","Temp.jul","Temp.jan","Floodplain.site","Water.source.type",
  "Geomorph.river.type","Distance.from.source","Area.ctch","Natural.sediment",
  "Ecoreg","Eft.type","Medit","Species","Total.number.run1",
  "Number.length.below.150","Number.length.over.150","Fished.area","Method","Sampling.location",
  "code","species")])

if(!try.test(input)){
  print("Error 4: non convenient input file!")
  stop("Error 4: non convenient input file!")
}

names(input) <- c("Sample.code","Day","Month","Year","Longitude","Latitude",
 "Actual.river.slope","temp.jul","temp.jan","Floodplain.site",
  "Water.source.type","Geomorph.river.type","Distance.from.source",
  "AREA.ctch","Natural.sediment","ECOREG","EFT.Type","MEDIT","Species",
  "Total.number.run1","Number.length.below.150","Number.length.over.150",
  "Fished.area","Method","Sampling.location","code","species")

#--------------------------------------------------
# faunistic table pr?paration
# give a table [code x species]
#--------------------------------------------------
faun <- prepTab(input,param=c("code","species","Total.number.run1"))

#--------------------------------------------------
# faunistic table pr?paration under constraint
# of the fish length
# give a table [code x species]
#--------------------------------------------------
length1 <- prepTab(input,param=c("code","species","Number.length.below.150"))
length2 <- prepTab(input,param=c("code","species","Number.length.over.150"))
#----------------------------------------------------------
# metric computation
# last modified: 03/12/09 14:28:05
#----------------------------------------------------------
labels1 <- c("WQO2.O2INTOL","HTOL.HINTOL","HabSp.RHPAR","Repro.LITH")
type1 <- c("dens","dens","ric","dens")
obs1 <- MetricCompute(faun,guild,labels=labels1,type=type1)

# estimation of the richness and abundance (captures)
abundance1 <- attributes(obs1)$Density
richness1 <- attributes(obs1)$Richness
psalmo1 <- attributes(obs1)$Psalmo

# for the metric based on length fish, we only use type="dens"
labels2 <- c("WQO2.O2INTOL","HTOL.HINTOL")
type2 <- rep("dens",2)
obslength1 <- MetricCompute(length1,guild,labels=labels2,type=type2)
obslength2 <- MetricCompute(length2,guild,labels=labels2,type=type2)

#--------------------------------------------------
# checking the length table
# obslength1 + obslength2 must equal to obs1
# 11/25/08 18:32:41
#--------------------------------------------------
obslength1 <- ValidatedLength(target=obslength1,complement=obslength2,obs=obs1)
print("Observed value computation: ok!")
#--------------------------------------------------
# Environmental table preparation
# 11/25/08 18:32:41
#--------------------------------------------------
varenv <- c("code","Actual.river.slope","temp.jul","temp.jan","Floodplain.site","Water.source.type",
  "Geomorph.river.type","Distance.from.source","AREA.ctch","Natural.sediment",
  "ECOREG","EFT.Type","Latitude","Longitude","MEDIT","Fished.area","Method","Sampling.location",
  "Sample.code","Day","Month","Year","Longitude","Latitude")
env1 <- prepEnv(input,param=varenv,richness=richness1,abundance=abundance1,psalmo=psalmo1,
  geomodels=geomodels,nameOrder=row.names(obs1))

#--------------------------------------------------
# computation of expected values
# 11/25/08 18:32:41
#--------------------------------------------------
fit1 <- ExpectedMetric(env1[row.names(obs1),],models)
fitlength <- ExpectedMetric(env1[row.names(obslength1),],lengthmodels)
print("Expected value computation: ok!")
#--------------------------------------------------
# regionalisation for the two metric types
# 12/05/08 15:33:08
#--------------------------------------------------
region1 <- prepRegEco(env1)
region1 <- rescale.names(region1,obs1)
regionlength <- rescale.names(region1,obslength1)

#--------------------------------------------------
# metric based on the length
# 12/05/08 17:56:41 by pbady
#--------------------------------------------------
# avec Trout
distLength <- DistMetric(obslength1,fitlength,regionlength$ECOREG2,stdpacklength$center,stdpacklength$scale,labels=c("dens.WQO2.O2INTOL","dens.HTOL.HINTOL"))
idLength <- IDScoring(distLength,stdpacklength,labels=c("dens.HTOL.HINTOL"))

# avec Trout
distTrout <- DistMetric(obs1,fit1,region1$ECOREG2,stdpackTrout$center,stdpackTrout$scale,labels=c("dens.WQO2.O2INTOL","dens.HTOL.HINTOL"))
idTrout<- IDScoring(distTrout,stdpackTrout,labels=c("dens.WQO2.O2INTOL"))

# avec No Trout
distNoTrout <- DistMetric(obs1,fit1,region1$ECOREG2,stdpackNoTrout$center,stdpackNoTrout$scale,labels=c("ric.HabSp.RHPAR","dens.Repro.LITH"))
idNoTrout <- IDScoring(distNoTrout,stdpackNoTrout,labels=c("ric.HabSp.RHPAR","dens.Repro.LITH"))

# print information
print("Metric computation: ok!")

# aggregation des tableaux
output <- addTab(idTrout,idNoTrout)
#names(output)[1] <-"dens.WQO2.O2INTOL"
output <- addTab(idLength,output)
names(output)[c(1,2)] <-c("dens.HTOL.HINTOL150","dens.HTOL.HINTOL150.class")
output <- addTab(output,region1)
output <- addTab(output,env1[,c("Sample.code","Day","Month","Year","Longitude",
  "Latitude","Fished.area","richesse","captures","psalmo","Method","method","Sampling.location")])

# observed values
obs1 <- obs1[,c("dens.WQO2.O2INTOL","ric.HabSp.RHPAR","dens.Repro.LITH")]
names(obs1) <- paste(names(obs1),"obs",sep=".")
output <- addTab(output,obs1)

# expected values
fit1 <- fit1[,c("dens.WQO2.O2INTOL","ric.HabSp.RHPAR","dens.Repro.LITH")]
names(fit1) <- paste(names(fit1),"pred",sep=".")
output <- addTab(output,fit1)

# observed values based on fish length
w <- obslength1[,c("dens.HTOL.HINTOL")]
names(w) <- row.names(obslength1)
output <- addTab(w,output)
names(output)[1] <- "dens.HTOL.HINTOL150.obs"

# expected values based on fish length
w <- fitlength[,c("dens.HTOL.HINTOL")]
names(w) <- row.names(fitlength)
output <- addTab(w,output)
names(output)[1] <- "dens.HTOL.HINTOL150.pred"

# print information
print("Table aggregation: ok!")
#--------------------------------------------------
# aggregation
# 12/05/08 17:56:41 by pbady
#--------------------------------------------------
output <- IDaggregation(output,qclassMean)

# print information
print("Index computation: ok!")

#--------------------------------------------------
# metricbased on diadromous species
# 12/01/08 18:22:43
# in construction
#--------------------------------------------------
# diadromous1 <- try(read.table(diadromnames,header=T,sep=";",comment.char=""))
# if(!try.test(diadromous1)){
#  print("Error 5: non convenient diadromous input file!")
#  stop("Error 5: non convenient diadromous input file!")
#}
# diadmet1 <- getSimilarity(diadromous1)
# output <- addTab(output,diadmet1)

#--------------------------------------------------
# recoding the names of the object 'output'
# 12/12/08 10:51:33 by pbady
#--------------------------------------------------
output<- RecodeNameOutput(output)
print("Names recoding: ok!")

#--------------------------------------------------
# write output file
# 12/12/08 10:51:33 by pbady
#--------------------------------------------------
# Save output as .xlsx
outtest <- try(writexl::write_xlsx(output, path = paste0(outputname)))
# Error handling
if (inherits(outtest, "try-error")) {
  print("Error 6: Non-convenient output file!")
  stop("Error 6: Non-convenient output file!")
}
#--------------------------------------------------
# End of the programs
# 12/18/08 14:56:21 by pbady
#---------------------------------------------
# second time point
t1 <- proc.time()
print(paste("Time: ",t1[3]-t0[3]," s",sep=""))
print(paste("Sample number: ",nrow(output),sep=""))











setwd("P:/Lukas/GitHubRepositories/DLCAnalyzer-master/data/EPM/")

library(sp)
library(imputeTS)
library(ggplot2)
library(ggmap)
library(data.table)

speed_cutoff <- 5             #cutoff in cm/s ,above which an animal is considered moving
likelyhood_cutoff <- 0.95       #DLC likelyhood below which points are removed and insted imputed
scale_existence <- 1.8
arena_size <- 65.5*65.5             #size of the y-maze in cm x cm
integration_window <- 5         #number of frames (+-) over which certain scorings are averages
fps <- 25                       #frames per second of the video
input_folder <- "Output_DLC/"      #folder in which all CSV input files are
output_folder <- "Results/"
maxspeed <- 20                  #maximum pixels a point can move in 1 frame before being removed
cut_start <- 100                #frames cut off at the start of the video
cut_end <- 100                  #frames cut off at the end of the video

#================ MAIN LOOP =========================

dir.create(output_folder)
files <- list.files(input_folder)
parameters <- data.frame(speed_cutoff,scale_existence,likelyhood_cutoff,arena_size,fps,maxspeed)

Full_Report <- NULL
PlotsBC <- vector("list")
PlotsNose <- vector("list")
SkeletonData <- NULL

for(j in files){
    #Read data and remove header
  data_header <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, nrows = 2)
  header_1 <- data.frame(sapply(data_header[1,], as.character), stringsAsFactors=FALSE)
  header_2 <- data.frame(sapply(data_header[2,], as.character), stringsAsFactors=FALSE)
  print(paste("reading file ",j, sep = ""))
  data <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, skip = 2)
  names(data) <- c("frame",paste(header_1[-1,1],header_2[-1,1],sep = "_"))

  Results <- parseEPMdata(data,parameters)
  
  SkeletonData_animal <- cbind(File = j, frame = data$frame, create_skeleton_v1(Results$long_data))
  
  

  dataBC <- Results$long_data[Results$long_data$Feature == "bodycentre",]
  dataNose <- Results$long_data[Results$long_data$Feature == "nose",]
  dataHC <- Results$long_data[Results$long_data$Feature == "headcentre",]
  dataNeck <- Results$long_data[Results$long_data$Feature == "neck",]
  dataBC <- dataBC[cut_start:(nrow(dataBC)-cut_end),]
  dataNose <- dataNose[cut_start:(nrow(dataNose)-cut_end),]
  dataHC <- dataHC[cut_start:(nrow(dataHC)-cut_end),]
  dataNeck <- dataNeck[cut_start:(nrow(dataNeck)-cut_end),]
  
  
  Report <- data.frame(file = j, 
                       raw.distance = sum(dataBC$speed, na.rm = T) * Results$px_to_cm,
                       distance.moving = sum(dataBC[dataBC$is.moving ,"speed"], na.rm = T) * Results$px_to_cm,
                       distance.moving.nose = sum(dataNose[dataNose$is.moving ,"speed"], na.rm = T) * Results$px_to_cm,
                       raw.speed = mean(dataBC$speed, na.rm = T) * Results$px_to_cm * fps,
                       speed.moving = mean(dataBC[dataBC$is.moving ,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       speed.in.center = mean(dataBC[dataBC$is.moving & dataBC$is.in.center ,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       speed.in.open = mean(dataBC[dataBC$is.moving & (dataBC$is.in.open.A1 | dataBC$is.in.open.A2) ,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       speed.in.closed = mean(dataBC[dataBC$is.moving & (dataBC$is.in.closed.A1 | dataBC$is.in.closed.A1),"speed"], na.rm = T) * Results$px_to_cm * fps,
                       time.moving = sum(dataBC$is.moving, na.rm = T) / fps,
                       time.stationary = sum(!dataBC$is.moving, na.rm = T) / fps,
                       time.in.center = sum(dataBC$is.in.center) / fps,
                       time.in.open.A1 = sum(dataBC$is.in.open.A1) / fps,
                       time.in.open.A2 = sum(dataBC$is.in.open.A2) / fps,
                       time.in.closed.A1 = sum(dataBC$is.in.closed.A1) / fps,
                       time.in.closed.A2 = sum(dataBC$is.in.closed.A2) / fps,
                       time.in.open = sum(dataBC$is.in.open.A1 | dataBC$is.in.open.A2) / fps,
                       time.in.closed = sum(dataBC$is.in.closed.A1 | dataBC$is.in.closed.A2) / fps,
                       transitions.center = calculate_transitions(dataBC$is.in.center,integration_window),
                       transitions.open = calculate_transitions((dataBC$is.in.open.A1 | dataBC$is.in.open.A2),integration_window),
                       transitions.closed = calculate_transitions((dataBC$is.in.closed.A1 | dataBC$is.in.closed.A2),integration_window),
                       nose.dip = calculate_transitions((!dataHC$is.in.arena & !dataNeck$is.in.closed.A1 & !dataNeck$is.in.closed.A2), integration_window) / 2
                       )
  
  Full_Report <- rbind(Full_Report, Report)
  SkeletonData <- rbind(SkeletonData,SkeletonData_animal)
  
  #create 2d plots of features
  p1 <- plot_density_path(dataBC,10) + 
    geom_path(data = rbind(Results$arena,Results$arena[1,]),aes(x,y), color = "darkolivegreen4")
  PlotsBC[[j]] <- p1
  
  p2 <- plot_density_path(dataNose,10) + 
    geom_path(data = rbind(Results$arena,Results$arena[1,]),aes(x,y), color = "darkolivegreen4")
  PlotsNose[[j]] <- p2 
} 


write.table(Full_Report,paste(output_folder,"Results.csv", sep = ""), sep = ";", row.names = F)
write.table(SkeletonData,paste(output_folder,"SkeletonData.csv", sep = ""), sep = ";", row.names = F)

pdf(paste(output_folder,"BodyCenterDensity.pdf",sep = ""), height = 10, width = 12)
for(j in PlotsBC){
  print(j)
}
dev.off() 

pdf(paste(output_folder,"NoseDensity.pdf",sep = ""), height = 10, width = 12)
for(j in PlotsNose){
  print(j)
}
dev.off() 

#======================================= REQUIRED LIBRARIES (install first) ==============================

library(sp)
library(imputeTS)
library(ggplot2)
library(ggmap)
library(data.table)
library(keras)

#======================================= PARAMETERS ==============================

setwd("G:/DLC/DLC paper 2019 Datapipeline final/OFT") #sets working directory
speed_cutoff <- 5             #cutoff in cm/s ,above which an animal is considered moving
scale_center <- 0.5             #of overall area (side ratio, not volumetric) that is considered the center
scale_periphery <- 0.8          #of overall area (side ratio, not volumetric) that is considered the periphery
scale_existence <- 1.3          #of overall area (side ratio, not volumetric) that is considered the known universe. points outside will be removed
scale_corners <- 0.2            #of overall area (side ratio, not volumetric) that is considered corners
likelyhood_cutoff <- 0.95       #DLC likelyhood below which points are removed and insted imputed
arena_size <- 42*42             #size of the arena in cm x cm, required to calculate the pixle to cm ratio
integration_window <- 5         #number of frames (+-) over which certain scorings are averaged
fps <- 25                       #frames per second of the video
input_folder <- "Output_DLC/"      #folder in which all CSV input files are
output_folder <- "Results/"     #Folder where Results (table,pdfs) will be saved
maxspeed <- 20                  #maximum pixels a point can move in 1 frame before being removed and interpolated
cut_start <- 100                #frames cut off at the start of the video, helps to remove start artefacts
cut_end <- 100                  #frames cut off at the end of the video, helps to remove end artefacts
Network_path <- "Networks/Oliver_OFT_model.hdf5"        #select hdf5 file that contains the network for feature selection
Network_parameters_path <- "Networks/Oliver_OFT_para.RDS"  #select RDS file that contains parameters for the same network

#======================================= MAIN LOOP =========================

dir.create(output_folder)
files <- list.files(input_folder)
parameters <- data.frame(scale_corners,speed_cutoff,scale_center,scale_periphery,scale_existence,likelyhood_cutoff,arena_size,fps,maxspeed)

#Initialize
Full_Report <- NULL
PlotsBC <- vector("list")
PlotsNose <- vector("list")
PlotsBehavior <- vector("list")
model <- load_model_hdf5(Network_path)
model_para <- readRDS(Network_parameters_path)

#process all files
for(j in files){
  
  #Read data and remove header
  data_header <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, nrows = 2)
  header_1 <- data.frame(sapply(data_header[1,], as.character), stringsAsFactors=FALSE)
  header_2 <- data.frame(sapply(data_header[2,], as.character), stringsAsFactors=FALSE)
  print(paste("reading file ",j, sep = ""))
  data <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, skip = 2)
  names(data) <- c("frame",paste(header_1[-1,1],header_2[-1,1],sep = "_"))

  #parse data and transform in long data format
  print(paste("processing file ",j, sep = ""))
  Results <- parseOFTdata(data,parameters)
  
  #use long data to create a Skeleton representation over time, then Zscore normalize this representation
  SkeletonData <- create_skeleton_v1(Results$long_data)
  SkeletonData_Norm <- ZScore_normalize_skeleton_v1(SkeletonData)
  
  #Expand Skeletion data over integration window to generate test data compatible with pretrained model, then use model to classify behaviors
  x_test <- create_test_set_OFT(SkeletonData_Norm, model_para$integration_window)
  behavior <- model %>% predict_classes(as.matrix(x_test$x))
  behavior <- data.frame(type = model_para$Feature_names[behavior + 1], 
                     time = (model_para$integration_window):(length(behavior)+model_para$integration_window - 1) / parameters$fps)
  
  dataBC <- Results$long_data[Results$long_data$Feature == "bodycentre",]
  dataNose <- Results$long_data[Results$long_data$Feature == "nose",]
  dataBC <- dataBC[cut_start:(nrow(dataBC)-cut_end),]
  dataNose <- dataNose[cut_start:(nrow(dataNose)-cut_end),]
  
  #Create final report
  Report <- data.frame(file = j, 
                       raw.distance = sum(dataBC$speed, na.rm = T) * Results$px_to_cm,
                       distance.moving = sum(dataBC[dataBC$is.moving ,"speed"], na.rm = T) * Results$px_to_cm,
                       distance.moving.nose = sum(dataNose[dataNose$is.moving ,"speed"], na.rm = T) * Results$px_to_cm,
                       raw.speed = mean(dataBC$speed, na.rm = T) * Results$px_to_cm * fps,
                       speed.moving = mean(dataBC[dataBC$is.moving ,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       speed.in.center = mean(dataBC[dataBC$is.moving & dataBC$is.in.center ,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       speed.in.periphery = mean(dataBC[dataBC$is.moving & dataBC$is.in.periphery ,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       speed.in.corners = mean(dataBC[dataBC$is.moving & (dataBC$is.in.C1 | dataBC$is.in.C2 | dataBC$is.in.C3 | dataBC$is.in.C4),"speed"], na.rm = T) * Results$px_to_cm * fps,
                       time.moving = sum(dataBC$is.moving, na.rm = T) / fps,
                       time.stationary = sum(!dataBC$is.moving, na.rm = T) / fps,
                       time.in.center = sum(dataBC$is.in.center) / fps,
                       time.in.periphery = sum(dataBC$is.in.periphery) / fps,
                       time.in.periphery.nose = sum(dataNose$is.in.periphery) / fps,
                       time.in.C1 = sum(dataBC$is.in.C1) / fps,
                       time.in.C2 = sum(dataBC$is.in.C2) / fps,
                       time.in.C3 = sum(dataBC$is.in.C3) / fps,
                       time.in.C4 = sum(dataBC$is.in.C4) / fps,
                       time.in.corners = sum(dataBC$is.in.C1 | dataBC$is.in.C2 | dataBC$is.in.C3 | dataBC$is.in.C4) / fps,
                       transitions.center = calculate_transitions(dataBC$is.in.center,integration_window),
                       transitions.periphery = calculate_transitions(dataBC$is.in.periphery,integration_window),
                       transitions.corners = calculate_transitions((dataBC$is.in.C1 | dataBC$is.in.C2 | dataBC$is.in.C3 | dataBC$is.in.C4),integration_window),
                       Supported.Rears.count = sum(abs(integratevector(avgbool(behavior$type == "Supported",10))))/2,
                       Supported.Rears.time = sum((avgbool(behavior$type == "Supported",10))) / parameters$fps,
                       Unsupported.Rears.count = sum(abs(integratevector(avgbool(behavior$type == "Unsupported",10))))/2,
                       Unsupported.Rears.time = sum((avgbool(behavior$type == "Unsupported",10))) / parameters$fps
                       )
  Full_Report <- rbind(Full_Report, Report)
  
  
  #create 2d plots of features
  p1 <- plot_density_path(dataBC,4) + 
    geom_path(data = rbind(Results$arena,Results$arena[1,]),aes(x,y), color = "darkolivegreen4") +
    geom_path(data = rbind(Results$arena.periphery,Results$arena.periphery[1,]),aes(x,y), color = "darkolivegreen4") +
    geom_path(data = rbind(Results$arena.center,Results$arena.center[1,]),aes(x,y), color = "darkolivegreen4")
  PlotsBC[[j]] <- p1
  
  p2 <- plot_density_path(dataNose,2) + 
    geom_path(data = rbind(Results$arena,Results$arena[1,]),aes(x,y), color = "darkolivegreen4") +
    geom_path(data = rbind(Results$arena.periphery,Results$arena.periphery[1,]),aes(x,y), color = "darkolivegreen4") +
    geom_path(data = rbind(Results$arena.center,Results$arena.center[1,]),aes(x,y), color = "darkolivegreen4")
  PlotsNose[[j]] <- p2 
  
  #create behavior plots
  p3 <- ggplot() + geom_path(data=behavior,aes(time,as.numeric(avgbool(type == "Supported",10))), color = "red") + 
    geom_path(data=behavior,aes(time, 1.1 + as.numeric(avgbool(type == "Unsupported",10))), color = "blue") + theme_bw() + ylab("active behavior (supported rear = red, unsupported = blue)") + xlab("time / s") + ggtitle(j)
  
  PlotsBehavior[[j]] <- p3
} 


#======================================= EXPORT ================================
write.table(Full_Report,paste(output_folder,"Results.csv", sep = ""), sep = ";", row.names = F)

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

pdf(paste(output_folder,"Behaviors.pdf",sep = ""), height = 6, width = 12)
for(j in PlotsBehavior){
  print(j)
}
dev.off() 


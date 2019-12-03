setwd("P:/Oliver/DLC PAPER/PIPELINE FINAL/FST")

library(sp)
library(imputeTS)
library(ggplot2)
library(ggmap)
library(data.table)


fps <- 25                       #frames per second of the video
input_folder <- "Output_DLC/"      #folder in which all CSV input files are
output_folder <- "Results/"
maxspeed <- 20                  #maximum pixels a point can move in 1 frame before being removed
likelyhood_cutoff <- 0.95
diameter_toptobottom <- 20
swim_cutoff <- 15              #numeric value that is used to determine if an animals is floating or not. rates the average change of bodyarea (in px^2)/frame over time
speed_cutoff <- 5
cut_start <- 300
cut_end <- 100

#================ MAIN LOOP =========================

dir.create(output_folder)
files <- list.files(input_folder)
parameters <- data.frame(fps,maxspeed,speed_cutoff,likelyhood_cutoff,diameter_toptobottom)

Full_Report <- NULL
PlotsBC <- vector("list")
PlotsNose <- vector("list")

for(j in files){
    #Read data and remove header
  data_header <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, nrows = 2)
  header_1 <- data.frame(sapply(data_header[1,], as.character), stringsAsFactors=FALSE)
  header_2 <- data.frame(sapply(data_header[2,], as.character), stringsAsFactors=FALSE)
  print(paste("reading file ",j, sep = ""))
  data <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, skip = 2)
  names(data) <- c("frame",paste(header_1[-1,1],header_2[-1,1],sep = "_"))

  Results <- parseFSTdata(data,parameters)
  
  SkeletonData <- create_skeleton_FST_v1(Results$long_data)
  
  dataBC <- Results$long_data[Results$long_data$Feature == "bodycentre",]
  dataNose <- Results$long_data[Results$long_data$Feature == "nose",]
  dataBC$is.floating <- (avgmean(abs(integratevector(get_polygon_areas(Results$long_data,c("headcentre","bcr","tailbase","bcl")))),20) < swim_cutoff)
  
  dataBC <- dataBC[cut_start:(nrow(dataBC)-cut_end),]
  dataNose <- dataNose[cut_start:(nrow(dataNose)-cut_end),]

  Report <- data.frame(file = j, 
                       raw.distance = sum(dataBC$speed, na.rm = T) * Results$px_to_cm,
                       raw.speed = mean(dataBC$speed, na.rm = T)  * Results$px_to_cm * fps,
                       distance.swimming = sum(dataBC[!dataBC$is.floating,"speed"], na.rm = T) * Results$px_to_cm,
                       speed.swimming = mean(dataBC[!dataBC$is.floating,"speed"], na.rm = T) * Results$px_to_cm * fps,
                       time.floating = sum(dataBC$is.floating) / parameters$fps,
                       percentage.floating = sum(dataBC$is.floating) / nrow(dataBC) * 100
                       )
  
  Full_Report <- rbind(Full_Report, Report)
  
  
  #create 2d plots of features
  p1 <- plot_density_path(dataBC,4)
  PlotsBC[[j]] <- p1
  
  p2 <- plot_density_path(dataNose,2)
  PlotsNose[[j]] <- p2 
} 


write.table(Full_Report,paste(output_folder,"Results_explorative14.csv", sep = ""), sep = ";", row.names = F)

pdf(paste(output_folder,"BodyCenterDensity_final.pdf",sep = ""), height = 10, width = 12)
for(j in PlotsBC){
  print(j)
}
dev.off() 

pdf(paste(output_folder,"NoseDensity_final.pdf",sep = ""), height = 10, width = 12)
for(j in PlotsNose){
  print(j)
}
dev.off() 

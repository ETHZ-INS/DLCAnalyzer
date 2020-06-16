library(sp)
library(imputeTS)
library(ggplot2)
library(ggmap)
library(data.table)
library(stringi)
library(keras)
library(abind)


#set wroking directory
setwd("DLC paper 2019 Datapipeline final/OFT/")

speed_cutoff <- 5             #cutoff in cm/s ,above which an animal is considered moving
scale_center <- 0.5             #of overall area (side ratio, not volumetric) that is considered the center
scale_periphery <- 0.8          #of overall area (side ratio, not volumetric) that is considered the periphery
scale_existence <- 1.3          #of overall area (side ratio, not volumetric) that is considered the known universe. points outside will be removed
scale_corners <- 0.2
likelyhood_cutoff <- 0.95       #DLC likelyhood below which points are removed and insted imputed
arena_size <- 42*42             #size of the arena in cm x cm
integration_window <- 15         #number of frames (+-) over which certain scorings are averages
fps <- 25                       #frames per second of the video
input_folder <- "Output_DLC/"      #folder in which all CSV input files are
output_folder <- "Results/"
maxspeed <- 20                  #maximum pixels a point can move in 1 frame before being removed
cut_start <- 100                #frames cut off at the start of the video
cut_end <- 100                  #frames cut off at the end of the video

#================ MAIN LOOP =========================

lab_data <- read.table("Labels/AllLabDataOFT_final.csv", sep = ";", header = T)
metadata <- read.table("metadata/metadata_OFT.csv", sep = ";", header = T)

dir.create(output_folder)
files <- list.files(input_folder)
parameters <- data.frame(scale_corners,speed_cutoff,scale_center,scale_periphery,scale_existence,likelyhood_cutoff,arena_size,fps,maxspeed)

SkeletonData <- list()
SkeletonData_Norm <- list()
for(j in unique(lab_data$DLCFile)){
  #Read data and remove header
  data_header <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, nrows = 2)
  header_1 <- data.frame(sapply(data_header[1,], as.character), stringsAsFactors=FALSE)
  header_2 <- data.frame(sapply(data_header[2,], as.character), stringsAsFactors=FALSE)
  print(paste("reading file ",j, sep = ""))
  data <- read.table(paste(input_folder,j,sep = ""), sep =",", header = T, skip = 2)
  names(data) <- c("frame",paste(header_1[-1,1],header_2[-1,1],sep = "_"))
  
  Results <- parseOFTdata(data,parameters)
  
  #Create Skeleton data and Zscore normalized Skeleton data
  SkeletonData[[j]] <- data.frame(create_skeleton_v1(Results$long_data))
  SkeletonData_Norm[[j]] <- ZScore_normalize_skeleton_v1(SkeletonData[[j]])
} 

#Select experimenter Behavior and Files for training to estimate network performance of experimenters
Experimenter <- c("Jin","Furkan","Oliver")
Behaviors <- c("Supported", "Unsupported","Grooming")


Files <- unique(lab_data$ID)
AllRes <- list()
ps <- list()
ResDFs <- list()
for(activeExperimenter in Experimenter){
for(activeFile in Files){
  print(activeExperimenter)
  print(paste("testing",activeFile,sep = " "))
  Files_train <- Files[Files != activeFile]
  Files_test <- activeFile 

lab_data_train <- lab_data[lab_data$Experimenter %in% activeExperimenter & lab_data$type %in% Behaviors & lab_data$ID %in% Files, ]
lab_data_train <- na.omit(lab_data_train)


x_train <- NULL
y_train <- NULL
for(i in Files_train){
  print(paste("adding train file",i,sep = " "))
  ldat <- lab_data_train[lab_data_train$ID == i,]
  Results <- create_train_set_OFT(SkeletonData_Norm[[paste(metadata[metadata$IDNumber == i, "DLCFile"])]], ldat, parameters,integration_window)
  x_train <- rbind(x_train,Results$x)
  y_train <- append(y_train,Results$y)
}

x_train <- as.matrix(x_train)
N_input <- ncol(x_train)
N_features <- length(unique(y_train))
Feature_names <- levels(as.factor(y_train))
y_train_cat <- to_categorical(-1 + as.integer(as.factor(y_train)))

new_order <- sample(1:nrow(y_train_cat))
x_train <- x_train[new_order,]
y_train_cat <- y_train_cat[new_order,]


model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 256, activation = 'relu', input_shape = c(N_input),kernel_regularizer = regularizer_l2(l = 0)) %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units = 128, activation = 'relu',kernel_regularizer = regularizer_l2(l = 0)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = N_features, activation = 'softmax')

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

history <- model %>% fit(
  x_train, y_train_cat, 
  epochs = 10, batch_size = 32, 
  validation_split = 0
)


x_test <- NULL
y_test <- NULL
lab_data_test<- lab_data[lab_data$Experimenter %in% activeExperimenter & lab_data$type %in% Behaviors & lab_data$ID %in% Files_test, ]
for(i in Files_test){
  print(paste("evaluating",i,sep = " "))
  ldat <- lab_data_test[lab_data_test$ID == i,]
  Results <- create_train_set_OFT(SkeletonData_Norm[[paste(metadata[metadata$IDNumber == i, "DLCFile"])]], ldat, parameters,integration_window)
  x_test <- rbind(x_test,Results$x)
  y_test <- append(y_test,Results$y)
}

x_test <- as.matrix(x_test)
y_test_cat <- to_categorical(-1 + as.integer(as.factor(y_test)))
AllRes[[paste(activeFile,activeExperimenter,sep="_")]] <- model %>% evaluate(x_test, y_test_cat)


res <- model %>% predict_classes(x_test)

res2 <- data.frame(type = Feature_names[res + 1], 
                   time = (integration_window):(length(res)+integration_window - 1) / parameters$fps,
                   original.scored = y_test)

p1 <- ggplot() + geom_path(data=res2,aes(time,as.numeric(avgbool(type == "Supported",10))), color = "red") + 
  geom_path(data=res2,aes(time,2.1 + as.numeric(avgbool(type == "Unsupported",10))), color = "blue") + 
  geom_path(data=res2,aes(time,4.2 + as.numeric(avgbool(type == "Grooming",10))), color = "green") + 
  geom_path(data=res2,aes(time,1.05 + as.numeric(original.scored == "Supported")), color = "red4") + 
  geom_path(data=res2,aes(time,3.15  + as.numeric(original.scored == "Unsupported")), color = "blue4") + 
  geom_path(data=res2,aes(time,5.25 + as.numeric(original.scored == "Grooming")), color = "green4") +
  theme_bw()

ps[[paste(activeFile,activeExperimenter,sep="_")]] <- p1
ResDFs[[paste(activeFile,activeExperimenter,sep="_")]] <- res2

}
}

saveRDS(ResDFs,"MLResults/ResDFs_final_reshuff.RDS")
saveRDS(AllRes,"MLResults/AllRes_final:reshuff.RDS")

pdf("MLResults/Results_final_reshuff.pdf", width = 15, height = 8)
for(i in 1:length(ps)){
  p1 <- ps[i]
  print( p1[[names(p1)]] + ggtitle(names(p1)))
}
dev.off()

Stats <- NULL
for(i in 1:length(ResDFs)){
  Stats <- rbind(Stats, data.frame(File = names(ResDFs)[i],
                                   Supported_DLC = sum(abs(integratevector(avgbool(ResDFs[[names(ResDFs)[i]]]$type == "Supported",10))))/2,
                                   Supported_Exp = sum(abs(integratevector(avgbool(ResDFs[[names(ResDFs)[i]]]$original.scored == "Supported",10))))/2,
                                   Unsupported_DLC = sum(abs(integratevector(avgbool(ResDFs[[names(ResDFs)[i]]]$type == "Unsupported",10))))/2,
                                   Unsupported_Exp = sum(abs(integratevector(avgbool(ResDFs[[names(ResDFs)[i]]]$original.scored == "Unsupported",10))))/2,
                                   Grooming_DLC = sum(abs(integratevector(avgbool(ResDFs[[names(ResDFs)[i]]]$type == "Grooming",10))))/2,
                                   Grooming_Exp = sum(abs(integratevector(avgbool(ResDFs[[names(ResDFs)[i]]]$original.scored == "Grooming",10))))/2,
                                   Supported_DLC_time = sum((avgbool(ResDFs[[names(ResDFs)[i]]]$type == "Supported",10))) / parameters$fps,
                                   Supported_Exp_time = sum((avgbool(ResDFs[[names(ResDFs)[i]]]$original.scored == "Supported",10)))/ parameters$fps,
                                   Unsupported_DLC_time = sum((avgbool(ResDFs[[names(ResDFs)[i]]]$type == "Unsupported",10))) / parameters$fps,
                                   Unsupported_Exp_time = sum((avgbool(ResDFs[[names(ResDFs)[i]]]$original.scored == "Unsupported",10))) / parameters$fps,
                                   Grooming_DLC_time = sum((avgbool(ResDFs[[names(ResDFs)[i]]]$type == "Grooming",10))) / parameters$fps,
                                   Grooming_Exp_time = sum((avgbool(ResDFs[[names(ResDFs)[i]]]$original.scored == "Grooming",10))) / parameters$fps
                                   
  ))
}

write.table(Stats,"MLResults/Results_final_reshuff.csv", sep = ";", row.names = F)


# Train final networks for individual experimenters
Experimenter <- c("Jin","Furkan","Oliver")
Behaviors <- c("Supported", "Unsupported","Grooming")
Files <- unique(lab_data$ID)

Networks <- list()

for(activeExperimenter in Experimenter){
  lab_data_train <- lab_data[lab_data$Experimenter %in% activeExperimenter & lab_data$type %in% Behaviors & lab_data$ID %in% Files, ]
  lab_data_train <- na.omit(lab_data_train)
  
  x_train <- NULL
  y_train <- NULL
  for(i in Files){
    print(paste("adding train file",i,sep = " "))
    ldat <- lab_data_train[lab_data_train$ID == i,]
    Results <- create_train_set_OFT(SkeletonData_Norm[[paste(metadata[metadata$IDNumber == i, "DLCFile"])]], ldat, parameters,integration_window)
    x_train <- rbind(x_train,Results$x)
    y_train <- append(y_train,Results$y)
  }
  
  x_train <- as.matrix(x_train)
  N_input <- ncol(x_train)
  N_features <- length(unique(y_train))
  Feature_names <- levels(as.factor(y_train))
  y_train_cat <- to_categorical(-1 + as.integer(as.factor(y_train)))
  
  new_order <- sample(1:nrow(y_train_cat))
  x_train <- x_train[new_order,]
  y_train_cat <- y_train_cat[new_order,]
  
  
  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = 256, activation = 'relu', input_shape = c(N_input),kernel_regularizer = regularizer_l2(l = 0)) %>% 
    layer_dropout(rate = 0.4) %>% 
    layer_dense(units = 128, activation = 'relu',kernel_regularizer = regularizer_l2(l = 0)) %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = N_features, activation = 'softmax')
  
  model  %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  history <- model %>% fit(
    x_train, y_train_cat, 
    epochs = 10, batch_size = 32, 
    validation_split = 0
  )
  
  #Store networks and parameters
  Res <- list()
  Res$Network <- model
  Res$N_input <- N_input
  Res$N_features <-N_features
  Res$Feature_names <- Feature_names
  Res$integration_window <- integration_window
  Networks[[paste(activeExperimenter)]] <- Res
}

saveRDS(Networks$Jin,"Networks/Jin_OFT_reshuff_para.RDS")
saveRDS(Networks$Oliver,"Networks/Oliver_OFT_reshuff_para.RDS")
saveRDS(Networks$Furkan,"Networks/Furkan_OFT_reshuff_para.RDS")

Networks$Jin$Network %>% save_model_hdf5("Networks/Jin_OFT_reshuff_model.hdf5")
Networks$Oliver$Network %>% save_model_hdf5("Networks/Oliver_OFT_reshuff_model.hdf5")
Networks$Furkan$Network %>% save_model_hdf5("Networks/Furkan_OFT_reshuff_model.hdf5")

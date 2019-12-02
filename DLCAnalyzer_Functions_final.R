
scale_arena <- function(arena, scale){
  data.frame(x = c(arena$x[1] + (arena$x[3]-arena$x[1]) * (1-scale)/2,
                   arena$x[2] + (arena$x[4]-arena$x[2]) * (1-scale)/2,
                   arena$x[3] + (arena$x[1]-arena$x[3]) * (1-scale)/2,
                   arena$x[4] + (arena$x[2]-arena$x[4]) * (1-scale)/2),
             y = c(arena$y[1] + (arena$y[3]-arena$y[1]) * (1-scale)/2,
                   arena$y[2] + (arena$y[4]-arena$y[2]) * (1-scale)/2,
                   arena$y[3] + (arena$y[1]-arena$y[3]) * (1-scale)/2,
                   arena$y[4] + (arena$y[2]-arena$y[4]) * (1-scale)/2))
}

scale_from_center <- function(arena, center, scale){
  arena$x <- center$x + (arena$x - center$x) * scale
  arena$y <- center$y + (arena$y - center$y) * scale
  return(arena)
}


scale_arena_corners <- function(arena, scale_corners){

  scale_one_corner <- function(x,scale_corners){
    c(x[1], x[1] + scale_corners * (x[2] - x[1]),x[1] + scale_corners * (x[3] - x[1]), x[1] + scale_corners * (x[4] - x[1]))
  }
  
  out <- list()
  out$C1 <- data.frame(x = scale_one_corner(arena[c(1,2,3,4),"x"],scale_corners), y = scale_one_corner(arena[c(1,2,3,4),"y"],scale_corners))
  out$C2 <- data.frame(x = scale_one_corner(arena[c(2,3,4,1),"x"],scale_corners), y = scale_one_corner(arena[c(2,3,4,1),"y"],scale_corners))
  out$C3 <- data.frame(x = scale_one_corner(arena[c(3,4,1,2),"x"],scale_corners), y = scale_one_corner(arena[c(3,4,1,2),"y"],scale_corners))
  out$C4 <- data.frame(x = scale_one_corner(arena[c(4,1,2,3),"x"],scale_corners), y = scale_one_corner(arena[c(4,1,2,3),"y"],scale_corners))
  return(out)
}

distance_xy <- function(x,feat1, feat2){
  df1 <- x[x$Feature == feat1,]
  df2 <- x[x$Feature == feat2,]
  sqrt((df1$x- df2$x)^2 + (df1$y - df2$y)^2)
}

angle_features <- function(x,feat1,feat2,feat3,feat4){
  a1 <- x[x$Feature == feat1,]
  a2 <- x[x$Feature == feat2,]
  b1 <- x[x$Feature == feat3,]
  b2 <- x[x$Feature == feat4,]
  ax <- a1$x - a2$x
  ay <- a1$y - a2$y
  bx <- b1$x - b2$x
  by <- b1$y - b2$y
  
  res <- rep(0,nrow(a1))
  for(i in 1:nrow(a1)){
    res[i] <- acos(sum(c(ax[i],ay[i])*c(bx[i],by[i])) / ( sqrt(sum(c(ax[i],ay[i]) * c(ax[i],ay[i]))) * sqrt(sum(c(bx[i],by[i]) * c(bx[i],by[i]))) ) )
  }
  return(res)
}

calculate_transitions <- function(x,integration_window){
  x <- avgbool(x,integration_window)
  sum(append(0, (x[2:length(x)]!= x[1:length(x)-1])))
}

avgbool <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- mean(x[max(0,i-window):min(length(x), i + window)], na.rm = T) > 0.5
  }
  return(res)
}

avgmean <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- mean(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

SDPlot <- function(x,nSD){
  ifelse((mean(x) - x) / sd(x) > -nSD, x, mean(x) + nSD * sd(x))
}

get_polygon_areas <-function(data,Points){
  x <- NULL
  y <- NULL
  for(i in Points){
    x <- cbind(x,data[data$Feature == i,"x"] )
    y <- cbind(y,data[data$Feature == i,"y"] )
  }
  area_polygons(x,y)
}


area_polygons <- function(x,y){
  area <- 0
  for(i in 1:dim(x)[2]){
    if(i != dim(x)[2]){
      area <- area + x[,i]*y[,i+1] - y[,i]*x[,i+1]
    }else{
      area <- area + x[,i]*y[,1] - y[,i]*x[,1]
    }
  }
  return(abs(area/2))
}

integratevector <- function(x){
  append(0, x[2:length(x)] - x[1:(length(x)-1)])
}

NormalizeZscore <- function(x){
  apply(x, 2, FUN = function(x){(x - mean(x)) / (sd(x))})
}

NormalizeZscore_median <- function(x){
  apply(x, 2, FUN = function(x){(x - median(x)) / (sd(x))})
}


ZScore_normalize_skeleton_v1 <- function(x){
  cbind(NormalizeZscore(x[,-((ncol(x)-1):ncol(x))]),x[,((ncol(x)-1):ncol(x))])
}

ZScore_normalize_skeleton_EPM <- function(x){
  cbind(NormalizeZscore(x[,-((ncol(x)-7):ncol(x))]),x[,((ncol(x)-7):ncol(x))])
}


#This function parses data from Csv files from DLC, interpolates missing values and calculates all measures that are important to generate the final report
parseOFTdata <- function(data,parameters){
  
  out <- list()
  
  #Read corners and define areas 
  corner_tl <- data.frame(x = median(data$tl_x),y = median(data$tl_y))
  corner_bl <- data.frame(x = median(data$bl_x),y = median(data$bl_y))
  corner_tr <- data.frame(x = median(data$tr_x),y = median(data$tr_y))
  corner_br <- data.frame(x = median(data$br_x),y = median(data$br_y))
  arena <- data.frame(x = c(corner_tl$x, corner_tr$x, corner_br$x, corner_bl$x),y = c(corner_tl$y, corner_tr$y, corner_br$y, corner_bl$y))
  arena_center <- scale_arena(arena,parameters$scale_center) 
  arena_corners <- scale_arena_corners(arena, parameters$scale_corners)
  arena_periphery <- scale_arena(arena,parameters$scale_periphery)
  arena_existence <- scale_arena(arena,parameters$scale_existence)
  px_to_cm <- sqrt(parameters$arena_size / area_polygons(t(arena$x),t(arena$y)))

  #create long data
  long_data <- NULL
  for(i in 1:((ncol(data) - 1) / 3)){
    long_data <- rbind(long_data,data.frame(frame = data[,1], x = data[,3*i-1], y = data[,3*i], Likelyhood = data[,3*i+1] , Feature = strsplit(paste(names(data)[3*i-1]),"_")[[1]][1]))
  }
  
  #Remove points with bad likelyhood or points that are outside of existence area
  long_data[long_data$Likelyhood < parameters$likelyhood_cutoff,"x"] <- NA
  long_data[long_data$Likelyhood < parameters$likelyhood_cutoff,"y"] <- NA
  long_data$is.existing <- point.in.polygon(long_data$x,long_data$y,arena_existence$x,arena_existence$y)
  long_data[!long_data$is.existing,"x"] <- NA
  long_data[!long_data$is.existing,"y"] <- NA
 
  
  Results <- NULL
  #get impute data within feature, get velocity and acceleration for each feature. define if it is moving at this frame or not
  for(i in levels(long_data$Feature)){
    datafeat <- long_data[long_data$Feature == i,]
    datafeat <- na.interpolation(datafeat)
    datafeat$delta_x <- integratevector(datafeat$x)
    datafeat$delta_y <- integratevector(datafeat$y)
    datafeat$speed <- sqrt(datafeat$delta_x ^2 + datafeat$delta_y^2)
    datafeat[datafeat$speed > parameters$maxspeed ,c("x","y","delta_x","delta_y","speed")] <- NA
    datafeat <- na.interpolation(datafeat)
    datafeat$acceleration <- integratevector(datafeat$speed)
    datafeat$is.moving <- (datafeat$speed * parameters$fps * px_to_cm) > parameters$speed_cutoff
    Results <- rbind(Results,datafeat)
  }
  
  #define areas and create booleans for the is in area checks
  Results$is.in.center <- point.in.polygon(Results$x,Results$y,arena_center$x,arena_center$y)
  Results$is.in.arena <- point.in.polygon(Results$x,Results$y,arena$x,arena$y)
  Results$is.in.periphery <- 1-point.in.polygon(Results$x,Results$y,arena_periphery$x,arena_periphery$y)
  Results$is.in.C1 <- point.in.polygon(Results$x,Results$y,arena_corners$C1$x,arena_corners$C1$y)
  Results$is.in.C2 <- point.in.polygon(Results$x,Results$y,arena_corners$C2$x,arena_corners$C2$y)
  Results$is.in.C3 <- point.in.polygon(Results$x,Results$y,arena_corners$C3$x,arena_corners$C3$y)
  Results$is.in.C4 <- point.in.polygon(Results$x,Results$y,arena_corners$C4$x,arena_corners$C4$y)
  Results$RealTime <- Results$frame / parameters$fps
  
  out$long_data <- Results
  out$px_to_cm <- px_to_cm
  out$arena <- arena
  out$arena.existence <- arena_existence
  out$arena.center <- arena_center
  out$arena.periphery <- arena_periphery
  out$arena.corners <- arena_corners
  return(out)
}

parseFSTdata <- function(data,parameters){
  
  out <- list()
  
  #create long data
  long_data <- NULL
  for(i in 1:((ncol(data) - 1) / 3)){
    long_data <- rbind(long_data,data.frame(frame = data[,1], x = data[,3*i-1], y = data[,3*i], Likelyhood = data[,3*i+1] , Feature = strsplit(paste(names(data)[3*i-1]),"_")[[1]][1]))
  }
  
  #Remove points with bad likelyhood or points that are outside of existence area
  long_data[long_data$Likelyhood < parameters$likelyhood_cutoff,"x"] <- NA
  long_data[long_data$Likelyhood < parameters$likelyhood_cutoff,"y"] <- NA
  
  top <- data.frame(x = median(data$t_x),y = median(data$t_y))
  bottom <- data.frame(x = median(data$b_x),y = median(data$b_y))
  px_to_cm <- parameters$diameter_toptobottom / sqrt((top$x- bottom$x)^2 + (top$y - bottom$y)^2)
  
  Results <- NULL
  #get impute data within feature, get velocity and acceleration for each feature. define if it is moving at this frame or not
  for(i in levels(long_data$Feature)){
    datafeat <- long_data[long_data$Feature == i,]
    datafeat <- na.interpolation(datafeat)
    datafeat$delta_x <- integratevector(datafeat$x)
    datafeat$delta_y <- integratevector(datafeat$y)
    datafeat$speed <- sqrt(datafeat$delta_x ^2 + datafeat$delta_y^2)
    datafeat[datafeat$speed > parameters$maxspeed ,c("x","y","delta_x","delta_y","speed")] <- NA
    datafeat <- na.interpolation(datafeat)
    datafeat$acceleration <- integratevector(datafeat$speed)
    datafeat$is.moving <- (datafeat$speed * parameters$fps) > parameters$speed_cutoff
    Results <- rbind(Results,datafeat)
  }
  
  Results$RealTime <- Results$frame / parameters$fps
  
  out$long_data <- Results
  out$px_to_cm <- px_to_cm
  return(out)
}


parseEPMdata <- function(data,parameters){
  
  out <- list()
  
  #Read corners and define areas 
  corner_tl <- data.frame(x = median(data$tl_x),y = median(data$tl_y))
  corner_tr <- data.frame(x = median(data$tr_x),y = median(data$tr_y))
  
  corner_bl <- data.frame(x = median(data$bl_x),y = median(data$bl_y))
  corner_br <- data.frame(x = median(data$br_x),y = median(data$br_y))
  
  corner_lt <- data.frame(x = median(data$lt_x),y = median(data$lt_y))
  corner_lb <- data.frame(x = median(data$lb_x),y = median(data$lb_y)) 
  
  corner_rt <- data.frame(x = median(data$rt_x),y = median(data$rt_y))
  corner_rb <- data.frame(x = median(data$rb_x),y = median(data$rb_y)) 
  
  corner_ctl <- data.frame(x = median(data$ctl_x),y = median(data$ctl_y))
  corner_ctr <- data.frame(x = median(data$ctr_x),y = median(data$ctr_y))
  corner_cbl <- data.frame(x = median(data$cbl_x),y = median(data$cbl_y))
  corner_cbr <- data.frame(x = median(data$cbr_x),y = median(data$cbr_y))
  
  arena <- data.frame(x = c(corner_tl$x, 
                            corner_tr$x, 
                            corner_ctr$x,
                            corner_rt$x,
                            corner_rb$x,
                            corner_cbr$x,
                            corner_br$x, 
                            corner_bl$x,
                            corner_cbl$x,
                            corner_lb$x,
                            corner_lt$x,
                            corner_ctl$x
                            ),
                      y = c(corner_tl$y, 
                            corner_tr$y, 
                            corner_ctr$y,
                            corner_rt$y,
                            corner_rb$y,
                            corner_cbr$y,
                            corner_br$y, 
                            corner_bl$y,
                            corner_cbl$y,
                            corner_lb$y,
                            corner_lt$y,
                            corner_ctl$y
                      ))
  center <- data.frame(x = mean(c(corner_ctr$x,corner_ctl$x,corner_cbr$x,corner_cbl$x)),y = mean(c(corner_ctr$y,corner_ctl$y,corner_cbr$y,corner_cbl$y)))
  arena_center <- data.frame(x = c(corner_ctl$x,
                                   corner_ctr$x,
                                   corner_cbr$x,
                                   corner_cbl$x),
                             y = c(corner_ctl$y,
                                   corner_ctr$y,
                                   corner_cbr$y,
                                   corner_cbl$y))
  closed_arms <- list()
  closed_arms$A1 <- data.frame(x = c(corner_tl$x,
                                   corner_tr$x,
                                   corner_ctr$x,
                                   corner_ctl$x),
                             y = c(corner_tl$y,
                                   corner_tr$y,
                                   corner_ctr$y,
                                   corner_ctl$y))
  closed_arms$A2 <- data.frame(x = c(corner_bl$x,
                                   corner_br$x,
                                   corner_cbr$x,
                                   corner_cbl$x),
                             y = c(corner_bl$y,
                                   corner_br$y,
                                   corner_cbr$y,
                                   corner_cbl$y))
  open_arms <- list()
  open_arms$A1 <- data.frame(x = c(corner_rb$x,
                                     corner_rt$x,
                                     corner_ctr$x,
                                     corner_cbr$x),
                               y = c(corner_rb$y,
                                     corner_rt$y,
                                     corner_ctr$y,
                                     corner_cbr$y))
  open_arms$A2 <- data.frame(x = c(corner_lb$x,
                                   corner_lt$x,
                                   corner_ctl$x,
                                   corner_cbl$x),
                             y = c(corner_lb$y,
                                   corner_lt$y,
                                   corner_ctl$y,
                                   corner_cbl$y))
  
  arena_existence <- scale_from_center(arena, center,parameters$scale_existence)
  arena_cubic <- data.frame(x = c(mean(arena[4:5,"x"]),mean(arena[4:5,"x"]),mean(arena[10:11,"x"]),mean(arena[10:11,"x"])),
                            y = c(mean(arena[1:2,"y"]),mean(arena[7:8,"y"]),mean(arena[7:8,"y"]),mean(arena[1:2,"y"])))
  px_to_cm <- sqrt(parameters$arena_size / area_polygons(t(arena_cubic$x),t(arena_cubic$y)))
  
  #create long data
  long_data <- NULL
  for(i in 1:((ncol(data) - 1) / 3)){
    long_data <- rbind(long_data,data.frame(frame = data[,1], x = data[,3*i-1], y = data[,3*i], Likelyhood = data[,3*i+1] , Feature = strsplit(paste(names(data)[3*i-1]),"_")[[1]][1]))
  }
  
  #Remove points with bad likelyhood or points that are outside of existence area
  long_data[long_data$Likelyhood < parameters$likelyhood_cutoff,"x"] <- NA
  long_data[long_data$Likelyhood < parameters$likelyhood_cutoff,"y"] <- NA
  long_data$is.existing <- point.in.polygon(long_data$x,long_data$y,arena_existence$x,arena_existence$y)
  long_data[!long_data$is.existing,"x"] <- NA
  long_data[!long_data$is.existing,"y"] <- NA
  
  
  Results <- NULL
  #get impute data within feature, get velocity and acceleration for each feature. define if it is moving at this frame or not
  for(i in levels(long_data$Feature)){
    datafeat <- long_data[long_data$Feature == i,]
    datafeat <- na.interpolation(datafeat)
    datafeat$delta_x <- integratevector(datafeat$x)
    datafeat$delta_y <- integratevector(datafeat$y)
    datafeat$speed <- sqrt(datafeat$delta_x ^2 + datafeat$delta_y^2)
    datafeat[datafeat$speed > parameters$maxspeed ,c("x","y","delta_x","delta_y","speed")] <- NA
    datafeat <- na.interpolation(datafeat)
    datafeat$acceleration <- integratevector(datafeat$speed)
    datafeat$is.moving <- (datafeat$speed * parameters$fps * px_to_cm) > parameters$speed_cutoff
    Results <- rbind(Results,datafeat)
  }
  
  #define areas and create booleans for the is in area checks
  Results$is.in.center <- point.in.polygon(Results$x,Results$y,arena_center$x,arena_center$y)
  Results$is.in.arena <- point.in.polygon(Results$x,Results$y,arena$x,arena$y)
  Results$is.in.open.A1 <- point.in.polygon(Results$x,Results$y,open_arms$A1$x,open_arms$A1$y)
  Results$is.in.open.A2 <- point.in.polygon(Results$x,Results$y,open_arms$A2$x,open_arms$A2$y)
  Results$is.in.closed.A1 <- point.in.polygon(Results$x,Results$y,closed_arms$A1$x,closed_arms$A1$y)
  Results$is.in.closed.A2 <- point.in.polygon(Results$x,Results$y,closed_arms$A2$x,closed_arms$A2$y)
  Results$RealTime <- Results$frame / parameters$fps
  
  out$long_data <- Results
  out$px_to_cm <- px_to_cm
  out$arena <- arena
  out$arena.existence <- arena_existence
  out$arena.center <- arena_center
  out$open.arms <- open_arms
  out$closed.arms <- closed_arms
  return(out)
}

#This function takes a long data frame created with ParseOFTdata and creates a skeleton description from it
create_skeleton_v1 <- function(x){
  MLdata <- data.frame(S1 = distance_xy(x,"nose","headcentre"))
  MLdata$S2 <- distance_xy(x,"headcentre","neck")
  MLdata$S3 <- distance_xy(x,"neck","bodycentre")
  MLdata$S4 <- distance_xy(x,"bodycentre","bcr")
  MLdata$S5 <- distance_xy(x,"bodycentre","bcl")
  MLdata$S6 <- distance_xy(x,"bodycentre","tailbase")
  MLdata$S7 <- distance_xy(x,"tailbase","hipr")
  MLdata$S8 <- distance_xy(x,"tailbase","hipl")
  MLdata$S9 <- distance_xy(x,"tailbase","tailcentre")
  MLdata$S10 <- distance_xy(x,"tailcentre","tailtip")
  MLdata$A1 <- angle_features(x,"tailbase","tailcentre","tailcentre","tailtip")
  MLdata$A2 <- angle_features(x,"hipr","tailbase","tailbase","hipl")
  MLdata$A3 <- angle_features(x,"tailbase","bodycentre","bodycentre","neck")
  MLdata$A4 <- angle_features(x,"bcr","bodycentre","bodycentre","bcl")
  MLdata$A5 <- angle_features(x,"bodycentre","neck","neck","headcentre")
  MLdata$A6 <- angle_features(x,"tailbase","bodycentre","neck","headcentre")
  MLdata$Ar1 <- get_polygon_areas(x,c("tailbase","hipr","hipl"))
  MLdata$Ar2 <- get_polygon_areas(x,c("hipr","hipl","bcl","bcr"))
  MLdata$Ar3 <- get_polygon_areas(x,c("bcr","earr","earl","bcl"))
  MLdata$Ar4 <- get_polygon_areas(x,c("earr","nose","earl"))
  MLdata$P1 <- as.numeric(x[x$Feature == "nose","is.in.arena"])
  MLdata$P2 <- as.numeric(x[x$Feature == "headcentre","is.in.arena"])
  MLdata <- as.data.frame(MLdata)
  
  return(MLdata)
}

create_skeleton_EPM <- function(x){
  MLdata <- data.frame(S1 = distance_xy(x,"nose","headcentre"))
  MLdata$S2 <- distance_xy(x,"headcentre","neck")
  MLdata$S3 <- distance_xy(x,"neck","bodycentre")
  MLdata$S4 <- distance_xy(x,"bodycentre","bcr")
  MLdata$S5 <- distance_xy(x,"bodycentre","bcl")
  MLdata$S6 <- distance_xy(x,"bodycentre","tailbase")
  MLdata$S7 <- distance_xy(x,"tailbase","hipr")
  MLdata$S8 <- distance_xy(x,"tailbase","hipl")
  MLdata$S9 <- distance_xy(x,"tailbase","tailcentre")
  MLdata$S10 <- distance_xy(x,"tailcentre","tailtip")
  MLdata$Sp1 <- x[x$Feature == "nose","speed"]
  MLdata$Sp2 <- x[x$Feature == "headcentre","speed"]
  MLdata$Sp3 <- x[x$Feature == "tailbase","speed"]
  MLdata$Sp4 <- x[x$Feature == "bodycentre","speed"]
  MLdata$A1 <- angle_features(x,"tailbase","tailcentre","tailcentre","tailtip")
  MLdata$A2 <- angle_features(x,"hipr","tailbase","tailbase","hipl")
  MLdata$A3 <- angle_features(x,"tailbase","bodycentre","bodycentre","neck")
  MLdata$A4 <- angle_features(x,"bcr","bodycentre","bodycentre","bcl")
  MLdata$A5 <- angle_features(x,"bodycentre","neck","neck","headcentre")
  MLdata$A6 <- angle_features(x,"tailbase","bodycentre","neck","headcentre")
  MLdata$Ar1 <- get_polygon_areas(x,c("tailbase","hipr","hipl"))
  MLdata$Ar2 <- get_polygon_areas(x,c("hipr","hipl","bcl","bcr"))
  MLdata$Ar3 <- get_polygon_areas(x,c("bcr","earr","earl","bcl"))
  MLdata$Ar4 <- get_polygon_areas(x,c("earr","nose","earl"))
  MLdata$P1 <- as.numeric(x[x$Feature == "nose","is.in.arena"])
  MLdata$P2 <- as.numeric(x[x$Feature == "headcentre","is.in.arena"])
  MLdata <- as.data.frame(MLdata)
  
  return(MLdata)
}


create_skeleton_FST_v1 <- function(x){
  MLdata <- data.frame(S1 = distance_xy(x,"nose","headcentre"))
  MLdata$S2 <- distance_xy(x,"headcentre","neck")
  MLdata$S3 <- distance_xy(x,"neck","bodycentre")
  MLdata$S4 <- distance_xy(x,"bodycentre","bcr")
  MLdata$S5 <- distance_xy(x,"bodycentre","bcl")
  MLdata$S6 <- distance_xy(x,"bodycentre","tailbase")
  MLdata$S7 <- distance_xy(x,"tailbase","hipr")
  MLdata$S8 <- distance_xy(x,"tailbase","hipl")
  MLdata$S9 <- distance_xy(x,"tailbase","tailcentre")
  MLdata$S10 <- distance_xy(x,"tailcentre","tailtip")
  MLdata$A1 <- angle_features(x,"tailbase","tailcentre","tailcentre","tailtip")
  MLdata$A2 <- angle_features(x,"hipr","tailbase","tailbase","hipl")
  MLdata$A3 <- angle_features(x,"tailbase","bodycentre","bodycentre","neck")
  MLdata$A4 <- angle_features(x,"bcr","bodycentre","bodycentre","bcl")
  MLdata$A5 <- angle_features(x,"bodycentre","neck","neck","headcentre")
  MLdata$A6 <- angle_features(x,"tailbase","bodycentre","neck","headcentre")
  MLdata$Ar1 <- get_polygon_areas(x,c("tailbase","hipr","hipl"))
  MLdata$Ar2 <- get_polygon_areas(x,c("hipr","hipl","bcl","bcr"))
  MLdata$Ar3 <- get_polygon_areas(x,c("bcr","earr","earl","bcl"))
  MLdata$Ar4 <- get_polygon_areas(x,c("earr","nose","earl"))
  MLdata <- as.data.frame(MLdata)
  return(MLdata)
}

equalizeSets <- function(x,y){
  N_obs <- NULL
  for(i in 1:ncol(y)){
    N_obs <- append(N_obs,sum(y[,i]))
  }
  
  keep <- rep(FALSE,nrow(y))
  
  for(i in 1:ncol(y)){
    keep[sample(which(y[,i] == 1))[1:min(N_obs)]] <- TRUE
  }
  out <- list()
  out$x <- x[keep,]
  out$y <- y[keep,]
  return(out)
}

create_train_set_OFT <- function(x, lab, parameters, integration_window){
  lab$from <- as.integer(lab$from * parameters$fps)
  lab$to <- as.integer(lab$to * parameters$fps)
  x_window <- x[1:(nrow(x) - 2*integration_window),]
  for(i in (-integration_window + 1):integration_window){
    x_window <- cbind(x_window, x[(integration_window + i + 1):(nrow(x) - integration_window + i),])
  }
  y <- rep("None",nrow(x_window))
  
  for(i in 1:nrow(lab)){
    y[(lab[i,"from"] - integration_window):(lab[i,"to"] - integration_window)] <- paste(lab[i,"type"])
  }
  y <- y[1:nrow(x_window)]
  
  out <- list()
  out$x <- x_window
  out$y <- y
  return(out)
}

create_test_set_OFT <- function(x, integration_window){
  x_window <- x[1:(nrow(x) - 2*integration_window),]
  for(i in (-integration_window + 1):integration_window){
    x_window <- cbind(x_window, x[(integration_window + i + 1):(nrow(x) - integration_window + i),])
  }

  out <- list()
  out$x <- x_window
  return(out)
}

plot_density_path <- function(data_plot,SDcutoff){
  xbreaks <- seq(floor(min(data_plot$x)), ceiling(max(data_plot$x)), by = 1)
  ybreaks <- seq(floor(min(data_plot$y)), ceiling(max(data_plot$y)), by = 1)
  data_plot$latbin <- xbreaks[cut(data_plot$x, breaks = xbreaks, labels=F)]
  data_plot$longbin <- ybreaks[cut(data_plot$y, breaks = ybreaks, labels=F)]
  ggplot(data = data_plot, aes(x,y)) + 
    stat_density_2d(data = data_plot, aes(latbin,longbin, fill=..density..), geom = "raster", contour = FALSE) + 
    scale_fill_gradient(name = "Time Density", low = "blue", high = "yellow") +
    geom_path(data = data_plot, aes(x,y, color = SDPlot((speed * Results$px_to_cm * fps),SDcutoff))) + 
    theme_bw() + 
    ggtitle(j) + 
    scale_color_gradient2(name = "speed (cm/s)", high = "white", low="black", mid = "black")
}

plotexperimenters <- function(x,b){
  AllPlotData <- NULL
  for(i in levels(x$Experimenter)){
    subdata <- x[x$Experimenter == i & x$type == b,]
    plot_dat <- data.frame(time = seq(from = 0, to = (max(x$to) + 1), by = 0.1))
    plot_dat$active.behavior <- 0
    plot_dat$Experimenter <- i
    for(j in 1:nrow(subdata)){
      plot_dat[between(plot_dat$time,lower = subdata[j,"from"], upper = subdata[j,"to"]), "active.behavior"] <- 1
    }
    AllPlotData <- rbind(AllPlotData,plot_dat)
  }
  return(AllPlotData)
}
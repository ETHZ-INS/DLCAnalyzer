#' Reads DLC Tracking data from a csv file and returns a TrackingData object
#' 
#' @param file path to a DLC tracking .csv file
#' @param fps frames per second of the recording. required to enable time resolved metrics
#' @return An object of type TrackingData
#' @examples
#' ReadDLCDataFromCSV("DLCData/Data.csv", fps = 25)
ReadDLCDataFromCSV <- function(file,fps = 1){
  out <- list()
  out$data <- list()
  data.header <- read.table(file, sep =",", header = T, nrows = 1)
  data.header <- data.frame(sapply(data.header, as.character), stringsAsFactors=FALSE)
  raw.data <- read.table(file, sep =",", header = T, skip = 2)
  for(i in seq(from = 2, to = nrow(data.header), by = 3)){
    out$data[[paste(data.header[i,])]] <- data.frame(frame = raw.data$coords, x = raw.data[,i], y = raw.data[,(i+1)], likelihood = raw.data[,(i+2)])
  }
  
  if(fps == 1){
    warning("no fps set. setting fps to 1. keep in mind that time based analyses are resolved in frames / second")
  }
  out$frames <- raw.data$coords
  out$fps <- fps
  out$seconds <- out$frames / fps
  
  out$median.data <- NULL
  for(i in names(out$data)){
    out$median.data <- rbind(out$median.data, data.frame(PointName = i, x = median(out$data[[i]]$x), y = median(out$data[[i]]$y)))
  }
  rownames(out$median.data) <- out$median.data$PointName
  
  out$point.info <- data.frame(PointName = names(out$data), PointType = "NotDefined")
  out$distance.units <- "pixel"
  out$labels <- list()
  out$filename <- last(strsplit(file,split = "/")[[1]])
  out$object.type = "TrackingData"

  return(out)
}

#' Adds a dataframe with point info to the tracking data and checks its validity
#' 
#' @param t an object of type TrackingData
#' @param pointinfo A data frame with additional point info. Requires variables: 'PointName' and `PointTyp'
#' @return An object of type TrackingData
#' @examples
#' AddPointInfo(Data, PointInfoDataFrame)
#' AddPointInfo(Data, data.frame(PointName = c("nose","tail","top","bottom), PointTyp = c("Mouse","Mouse","Maze","Maze)))
AddPointInfo <- function(t,pointinfo){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  out <- t
  if(is.null(pointinfo$PointName)){
    warning("point info does not contain required variable PointName. could not add point info")
    return(t)
  }
  if(is.null(pointinfo$PointType)){
    warning("point info does not contain required variable PointType could not add point info")
    return(t)
  }
  
  if(length(setdiff(names(t$data),pointinfo$PointName)) != 0){
    warning(paste("point info missing for following points:", setdiff(names(t$data),pointinfo$PointName), " "))
  }
  out$point.info <- pointinfo
  
  return(out)
}
 
IsTrackingData <- function(t){
  if(!is.null(t$object.type)){
    if(t$object.type == "TrackingData"){
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Cuts an object of type TrackingData into a shorter object of the same type
#' 
#' @param start if not null the first start frames will be removed
#' @param end if not null the last end frames will be removed
#' @param keep.frames a numeric vector. if not null, frames that intersect with this vector will be kept
#' @param remove.frames a numeric vector. if not null, frames that intersect with this vector will be removed
#' @return An object of type TrackingData
#' @examples
#' CutTrackingData(Data, start = 100, end = 100)
#' CutTrackingData(Data, keep.frames = c(10,11,12,13))
#' CutTrackingData(Data, remove.frames = c(21,22,23,24))
CutTrackingData <- function(t,start = NULL, end = NULL, remove.frames = NULL, keep.frames = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  keep <- t$frames
  if(!is.null(start)){
    keep <- keep[-(1:start)]
  }
  if(!is.null(end)){
    keep <- keep[-((length(keep) - end):length(keep))]
  }
  if(!is.null(remove.frames)){
    keep <- setdiff(keep, remove.frames)
  }
  if(!is.null(keep.frames)){
    keep <- intersect(keep, keep.frames)
  }
  t$frames <- keep
  if(!is.null(t$seconds)){
    t$seconds <- t$seconds[keep+1]
  }
  if(length(names(t$labels)) > 0){
    for(i in (names(t$labels))){
      t$labels[[i]] <- t$labels[[i]][keep+1]
    }
  }
  if(!is.null(t$features)){
    t$features <- t$features[keep + 1,]
  }
    
  for(i in 1:length(t$data)){
    t$data[[i]] <- t$data[[i]][t$data[[i]]$frame %in% keep,]
  }
  return(t)
}

#' Cleanes up an object of type Tracking data. missing data is replaced by data interpolation
#' 
#' @param t an object of type TrackingData
#' @param likelihoodcutoff points below this likelihoodcutoff from DLC will be interpolated
#' @param existence.pol points outside of the polygon existence.pol will be interpolated
#' @param maxdelta points that move more than maxdelta in one frame (cm or px, depending on calibrated data or not) will be interpolated
#' @return An object of type TrackingData
#' @examples
#' CleanTrackingData(t)
#' CleanTrackingData(t, likelihoodcutoff = 0.9)
#' CleanTrackingData(t, existence.pol = data.frame(x = c(0,100,100,0), y = c(100,100,0,0)))
#' CleanTrackingData(t, likelihoodcutoff = 1, maxdelta = 5)
CleanTrackingData <- function(t, likelihoodcutoff = 0.95, existence.pol = NULL, maxdelta = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  print(paste("interpolating points with likelihood < ", likelihoodcutoff, sep = ""))
  if(!is.null(existence.pol)){
    print("interpolating points which are outside of the existence area")
  }
  if(!is.null(maxdelta)){
  print(paste("interpolating points with a maximum delta of", maxdelta, t$distance.units,"per frame", sep = " "))
  }
  
  for(i in 1:length(t$data)){
    process <- t$data[[i]]$likelihood < likelihoodcutoff
    if(!is.null(existence.pol)){
      if(is.null(existence.pol$x) | is.null(existence.pol$y)){
        warning("warning. existence.pol is invalid. polygon data needs to include variable x and variable")
      }else if(length(existence.pol$x) != length(existence.pol$y)){
        warning("invalid polygon entered. polygon data needs to include variable x and variable y of equal length")
      }else{
      process <- process | !point.in.polygon(t$data[[i]]$x,t$data[[i]]$y,existence.pol$x, existence.pol$y)
      }
    }
    if(!is.null(maxdelta)){
      process <- process | (sqrt(integratevector(t$data[[i]]$x) ^2 + integratevector(t$data[[i]]$y)^2) > maxdelta)
    }
    t$data[[i]]$x[process] <- NA
    t$data[[i]]$y[process] <- NA
    t$data[[i]] <- na.interpolation(t$data[[i]])
  }
  return(t)
}

#' Adds an object median.data to a TrackingData object
#' 
#' @param t an object of type TrackingData
#' @param TypeName points of this TypeName are considered maze data. Requires existence of point information (See AddPointInfo() )
#' @return An object of type TrackingData
#' @examples
#' AddMazeData(t)
#' AddMazeData(t, TypeName = "SomeOtherMaze")
AddMazeData <- function(t,TypeName = "Maze"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!sum(t$point.info$PointType == TypeName)){
    warning(paste("Point type", TypeName,"does not exist. No maze data added", sep = " "))
    return(t)
  }
  t$median.data <- NULL
  for(i in t$point.info[t$point.info$PointType == TypeName,"PointName"]){
    t$median.data <- rbind(t$median.data, data.frame(PointName = i, x = median(t$data[[i]]$x), y = median(t$data[[i]]$y)))
  }
  rownames(t$median.data) <- t$median.data$PointName
  t$has.median.data <- TRUE
  return(t)
}

#' Calibrates an object of type TrackingData from pixel to metric space
#' 
#' @param t an object of type TrackingData to be calibrated
#' @param method use "distance" or "area". distance requires 2 points, area a polygon with > 2 points
#' @param in.metric the measured distance or area in metric units
#' @param points a vector of tracked points that are used for calibration
#' @return An object of type TrackingData
#' @examples
#' CalibrateTrackingData(t, method = "distance", in.metric = 40, points = c("top","bottom"))
#' CalibrateTrackingData(t, method = "area", in.metric = 1600, points = c("top.left","top.right","bottom.right","bottom.left"))
CalibrateTrackingData <- function(t, method, in.metric = NULL, points = NULL, ratio = NULL, new.units = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!method %in% c("distance","area","ratio")){
    warning("invalid method: valid methods are distance, area or ratio, can not calibrate")
    return(t)
  }
  if(method == "ratio"){
    if(is.numeric(ratio)){
    t$px.to.cm <- ratio
    }else{
      warning("method ratio needs a valid ratio to be entered. can not calibrate")
      return(t)
    }
  }else{
    if(is.numeric(in.metric) & is.null(points)){
      warning("method requires both points and a in.metric (numeric!) measurement")
      return(t)
    }
    if(sum(!(points %in% t$median.data$PointName))){
      warning("invalid points entered, can not calibrate")
      return(t)
    }

    if(method == "distance"){
      if(length(points) == 2){
        t$px.to.cm <- in.metric / Distance2d(t$median.data[points[1],],t$median.data[points[2],]) 
      }else{
        warning("invalid number of points: distance needs 2 points, can not calibrate")
        return(t)
      }
    }
    if(method == "area"){
      if(length(points) > 2){
        t$px.to.cm <- sqrt(in.metric / AreaPolygon2d(t$median.data[points,])) 
      }else{
        warning("invalid number of points: area need polygon of > 2 points, can not calibrate")
        return(t)
      }
    }
  }
  
  for(i in 1:length(t$data)){
    t$data[[i]]$x <- t$data[[i]]$x * t$px.to.cm
    t$data[[i]]$y <- t$data[[i]]$y * t$px.to.cm
  }
  if(!is.null(t$median.data)){
    t$median.data[,c("x","y")] <- t$median.data[,c("x","y")] * t$px.to.cm
  }
  if(!is.null(t$zones)){
    for(i in 1:length(t$zones)){
      t$zones[[i]]$x <- t$zones[[i]]$x * t$px.to.cm
      t$zones[[i]]$y <- t$zones[[i]]$y * t$px.to.cm
    }
  }
  if(!is.null(new.units)){
    t$distance.units <- new.units
  }else{
  t$distance.units <- "cm"
  }
  
  return(t)
} 

#' Scales a polygon by a factor
#' 
#' @param p a polygon (type list() or data.frame() with elements x and y are required)
#' @param factor scaling facor
#' @return a polygon
#' @examples
#' ScalePolygon(p = data.frame(x = c(0,0,1,1), y = c(0,1,1,0)), factor = 1.3)
ScalePolygon <-function(p, factor){
  if(is.null(p$x) | is.null(p$y)){
    stop("invalid input. Needs polygon of type list() or data.frame() with 2 variables, x = and y =")
  }
  center_x <- mean(p$x)
  center_y <- mean(p$y)
  return(data.frame(x = center_x + factor * (p$x - center_x), y = center_y + factor * (p$y - center_y)))
}

#' Recenters a polygon to a new center point
#' 
#' @param p a polygon (type list() or data.frame() with elements x and y are required)
#' @param new_center the new center point (type list() or data.frame() with elements x and y are required)
#' @return a polygon
#' @examples
#' RecenterPolygon(p = data.frame(x = c(0,0,1,1), y = c(0,1,1,0)), new_center = data.frame(x=1,y=1))
RecenterPolygon <-function(p, new_center){
  center_x <- mean(p$x)
  center_y <- mean(p$y)
  return(data.frame(x = p$x + new_center$x - center_x , y = p$y + new_center$y - center_y))
}

#' Adds Zones required for OFT analysis to an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param points a vector of point names describing the corners of the arena (!THE ORDER IS IMPORTANT)
#' @param scale_center a factor by which the arena is scaled to define the center
#' @param scale_periphery a factor by which the arena is scaled to define the periphery
#' @param scale_corners a factor by which the arena is scaled before recentering to each corner
#' @return a polygon
#' @examples
#' AddOFTZones(t, c("tl","tr","br","bl"), 0.5,0.4,0.8)
AddOFTZones <- function(t, points = c("tl","tr","br","bl"), scale_center = 0.5, scale_corners = 0.4, scale_periphery = 0.8){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(points,names(t$data) != 4))){
    warning("invalid number or type of points entered. exactly 4 existing points needed for OFT Zones")
    return(t)
  }
  zones <- list()
  zones.invert <- list()
  zones[["arena"]] <- t$median.data[points,c("x","y")]
  zones.invert[["arena"]] <- FALSE
  zones[["center"]] <- ScalePolygon(t$median.data[points,c("x","y")], scale_center)
  zones.invert[["center"]] <- FALSE
  zones[["periphery"]] <- ScalePolygon(t$median.data[points,c("x","y")], scale_periphery)
  zones.invert[["periphery"]] <- TRUE
  t$corner.names <- NULL
  for(i in points){
    zones[[paste("corner",i, sep = ".")]] <- RecenterPolygon(ScalePolygon(t$median.data[points,c("x","y")], scale_corners), t$median.data[i,c("x","y")])
    zones.invert[[paste("corner",i, sep = ".")]] <- FALSE
    t$corner.names <- append(t$corner.names,paste("corner",i, sep = "."))
  }
  t$zones <- zones
  t$zones.invert <- zones.invert
  return(t)
}

#' Performs an OFT analysis on an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param movement_cutoff a numeric value that denotes the cutoff point above which an animal is considered moving
#' @param integration_period a numeric value that denotes the duration over which metrics are smoothed.
#' @param points a string or vector of strings that denotes the name of points which will be analysed
#' @return a TrackingData object
#' @examples
#' OFTAnalysis(t, 5,5,"bodycentre")
OFTAnalysis <- function(t, movement_cutoff,integration_period, points){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined for OFT analysis. Returing simple analysis only")
  }
  
  t <- CalculateMovement(t,movement_cutoff,integration_period)
  
  t$Report <- list()
  if(!is.null(t$labels)){
    t$Report <- append(t$Report,LabelReport(t, integration_period))
  }
  
  for(k in points){
    dat <- t$data[[k]]
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(dat[,"speed"], na.rm = T)
    t$Report[[paste(k, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving,"speed"], na.rm = T)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(dat[,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "time.moving", sep = ".")]] <- sum(dat[,"is.moving"], na.rm = T) / t$fps
    t$Report[[paste(k, "total.time", sep = ".")]] <- length(dat[,"is.moving"]) / t$fps
    t$Report[[paste(k, "time.stationary", sep = ".")]] <- t$Report[[paste(k, "total.time", sep = ".")]] - t$Report[[paste(k, "time.moving", sep = ".")]]
    t$Report[[paste(k, "percentage.moving", sep = ".")]] <- t$Report[[paste(k, "time.moving", sep = ".")]] / t$Report[[paste(k, "total.time", sep = ".")]] * 100
    
    if(!is.null(t$zones)){
      t$Report <- append(t$Report, ZoneReport(t,k,"center", zone.name = paste(k,"center", sep = ".")))
      t$Report <- append(t$Report, ZoneReport(t,k,"periphery" , zone.name = paste(k,"periphery", sep = "."), invert = TRUE))
      t$Report <- append(t$Report, ZoneReport(t,k,t$corner.names, zone.name = paste(k,"corners", sep = ".")))
    }
    }
  return(t)
}

#' Performs an EPM analysis on an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param movement_cutoff a numeric value that denotes the cutoff point above which an animal is considered moving
#' @param integration_period a numeric value that denotes the duration over which metrics are smoothed.
#' @param points a string or vector of strings that denotes the name of points which will be analysed
#' @return a TrackingData object
#' @examples
#' EPMAnalysis(t, 5,5,"bodycentre")
EPMAnalysis <- function(t, movement_cutoff,integration_period, points,nosedips = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined for EPM analysis. Returing simple analysis only")
  }
  
  t <- CalculateMovement(t, movement_cutoff,integration_period)
  t$Report <- list()
  
  #THIS IS A VERA ARBITRARY SECTION FOR THIS TYPE OF ANALYSIS. IT WILL ONLY WORK WITH THE CORRECT POINTS AND ZONES
  if(nosedips){
    if((length(setdiff(c("headcentre","bodycentre","neck"),names(t$data))) != 0) | (length(setdiff(c("closed.top","closed.bottom","arena"),names(t$zones))) != 0)){
      warning("Not all points or zones needed for nosedip analysis. Requires points : headcentre,bodycentre,neck and zones closed.top, closed.bottom, arena")
    }else{
    t$labels$automatic.nosedip <- avgbool(!IsInZone(t,"headcentre","arena") & IsInZone(t,"bodycentre","arena") &!IsInZone(t,"neck",c("closed.top","closed.bottom")),integration_period)
    t$Report[["nose.dip"]] <- CalculateTransitions(t$labels$automatic.nosedip,integration_period) / 2
    t$labels$automatic.nosedip <- ifelse(t$labels$automatic.nosedip == 1,"Nosedip","None")
    }
  }
  
  for(k in points){
    dat <- t$data[[k]]
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(dat[,"speed"], na.rm = T)
    t$Report[[paste(k, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving,"speed"], na.rm = T)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(dat[,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "time.moving", sep = ".")]] <- sum(dat[,"is.moving"], na.rm = T) / t$fps
    t$Report[[paste(k, "total.time", sep = ".")]] <- length(dat[,"is.moving"]) / t$fps
    t$Report[[paste(k, "time.stationary", sep = ".")]] <- t$Report[[paste(k, "total.time", sep = ".")]] - t$Report[[paste(k, "time.moving", sep = ".")]]
    t$Report[[paste(k, "percentage.moving", sep = ".")]] <- t$Report[[paste(k, "time.moving", sep = ".")]] / t$Report[[paste(k, "total.time", sep = ".")]] * 100
    
    t$Report <- append(t$Report, ZoneReport(t,k,"center", zone.name = paste(k,"center", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,c("open.right","open.left"), zone.name = paste(k,"open", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,c("closed.top","closed.bottom"), zone.name = paste(k,"closed", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"closed.top",  zone.name = paste(k,"closed.top", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"closed.bottom",  zone.name = paste(k,"closed.bottom", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"open.left", zone.name = paste(k,"open.left", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"open.right", zone.name = paste(k,"open.right", sep = ".")))
  }
  return(t)
}

#' Checks if point p is in zone(s) z
#' 
#' @param t an object of type TrackingData
#' @param p string value of a point
#' @param z a string or vector of strings naming the zones to be checked
#' @param invert if TRUE instead it will be checked if the point is outside the zone
#' @return a boolean vector for the assessment at each frame
#' @examples
#' IsInZone(t, "bodycentre","center")
IsInZone <- function(t,p,z,invert = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(p,names(t$data))) != 1){
    warning(paste("Points not available in Tracking data:",paste(setdiff(points, names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  if(length(intersect(z,names(t$zones))) != length(z)){
    warning(paste("Zones in Tracking data:",paste(setdiff(z, names(t$zones)),collapse = " "), sep = " "))
    return(NULL)
  }
  
  
  zones <- t$zones[z]
  in.zone <- rep(FALSE,nrow(t$data[[p]]))
  for(i in zones){
    in.zone <- in.zone | (point.in.polygon(t$data[[p]]$x,t$data[[p]]$y,i$x,i$y) == 1)
  }
  if(invert){
    in.zone <- !in.zone
  }
  return(in.zone)
}

#' Creates a report for each behavior which includes total time of behavior and number of onset/offsets
#' 
#' @param t an object of type TrackingData
#' @param integration_period string value of a point
#' @return a list of metrics for each behavior
#' @examples
#' ClassificationReport(t, integration_period = 5)
LabelReport <- function(t, integration_period = 0){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(t$labels) == 0){
    warning("no label data present")
    return(NULL)
  }
  Report <- list()
  for(j in names(t$labels)){
  c <- SmoothLabel(t$labels[[j]],integration_period)
  c <- na.omit(c)
    for(i in unique(c)){
      Report[[paste(j,i,"time", sep = ".")]] <- sum(c == i) / t$fps
      Report[[paste(j,i,"count", sep = ".")]] <- sum(CalculateTransitions(c == i, 0)) / 2
    }
  }
  return(Report)
}

CalculateMovement <- function(t, movement_cutoff, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  for(i in 1:length(t$data)){
    t$data[[i]]$delta_x <- integratevector(t$data[[i]]$x)
    t$data[[i]]$delta_y <- integratevector(t$data[[i]]$y)
    t$data[[i]]$speed <- sqrt(t$data[[i]]$delta_x ^2 + t$data[[i]]$delta_y^2)
    t$data[[i]]$acceleration <- integratevector(t$data[[i]]$speed)
    t$data[[i]]$is.moving <- as.logical(avgbool(t$data[[i]]$speed > (movement_cutoff / t$fps), integration_period))
  }
  t$integration_period <- integration_period
  return(t)
}

CalculateAccelerations <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  for(i in 1:length(t$data)){
    t$data[[i]]$delta_x <- integratevector(t$data[[i]]$x)
    t$data[[i]]$delta_y <- integratevector(t$data[[i]]$y)
    t$data[[i]]$speed <- sqrt(t$data[[i]]$delta_x ^2 + t$data[[i]]$delta_y^2)
    t$data[[i]]$acceleration <- integratevector(t$data[[i]]$speed)
  }
  return(t)
}

ZoneReport <- function(t,point,zones, zone.name = NULL, invert = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  Report <- list()
  if(is.null(zone.name)){
    zone.name <- paste(zones, collapse = ".")
  }
  if(!point %in% names(t$data)){
    warning("Invalid point")
    return(NULL)
  }
  if(!sum(zones %in% names(t$zones))){
    warning("Invalid zone(s)")
    return(NULL)
  }
  
  dat <- t$data[[point]]
  in.zone <- rep(FALSE,nrow(dat))
  for(i in t$zones[zones]){
    in.zone <- in.zone | (point.in.polygon(dat$x,dat$y,i$x,i$y) == 1)
  }
  if(invert){
    in.zone <- !in.zone
  }
  Report[[paste(zone.name, "raw.distance", sep = ".")]] <- sum(dat[in.zone,"speed"], na.rm = T)
  Report[[paste(zone.name, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving & in.zone,"speed"], na.rm = T)
  Report[[paste(zone.name, "raw.speed", sep = ".")]] <- mean(dat[in.zone,"speed"], na.rm = T) * t$fps
  Report[[paste(zone.name, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving & in.zone,"speed"], na.rm = T) * t$fps
  Report[[paste(zone.name, "time.moving", sep = ".")]] <- sum(dat[in.zone,"is.moving"], na.rm = T) / t$fps
  Report[[paste(zone.name, "total.time", sep = ".")]] <- length(dat[in.zone,"is.moving"]) / t$fps
  Report[[paste(zone.name, "time.stationary", sep = ".")]] <- Report[[paste(zone.name, "total.time", sep = ".")]] - Report[[paste(zone.name, "time.moving", sep = ".")]]
  Report[[paste(zone.name, "percentage.moving", sep = ".")]] <- Report[[paste(zone.name, "time.moving", sep = ".")]] / Report[[paste(zone.name, "total.time", sep = ".")]] * 100
  Report[[paste(zone.name, "transitions", sep = ".")]] <- CalculateTransitions(in.zone, t$integration_period) 

  return(Report)
}

Distance2d <- function(a,b){
  sqrt((a$x - b$x)^2 + (a$y - b$y)^2)
}

MedianMouseArea <- function(t,points = c("nose","earr","bcr","hipr","tailbase","hipl","bcl","earl")){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  median(GetPolygonAreas(t,points), na.rm = T)
}

MedianMouseLength <- function(t, front  = "nose", back = "tailbase"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  median(GetDistances(t,front,back), na.rm = T)
}

AreaPolygon2d <- function(p){
  area <- 0
  for(i in 1:nrow(p)){
    if(i != nrow(p)){
      area <- area + p$x[i]*p$y[i+1] - p$y[i]*p$x[i+1]
    }else{
      area <- area + p$x[i]*p$y[1] - p$y[i]*p$x[1]
    }
  }
  return(abs(area/2))
}

GetPolygonAreas <-function(t,points){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(points,names(t$data))) != length(unique(points))){
    warning(paste("Points for distance measurement not available in Tracking data:",paste(setdiff(points, names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  x <- NULL
  y <- NULL
  for(i in points){
    x <- cbind(x,t$data[[i]]$x )
    y <- cbind(y,t$data[[i]]$y )
  }
  AreaPolygon3d(x,y)
}

GetDistances <- function(t,f,b){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(c(f,b),names(t$data))) != 2){
    warning(paste("Points for distance measurement not available in Tracking data:",paste(setdiff(c(f,b), names(t$data)),collapse = " "), sep = " "))
  }
  Distance2d(t$data[[f]],t$data[[b]])
}

AreaPolygon3d <- function(x,y){
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

AddLabelingData <- function(t, lab){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t$labels$manual <- rep("None",length(t$seconds))
  for(i in 1:nrow(lab)){
    t$labels$manual[t$seconds >= lab[i,"from"] & t$seconds <= lab[i,"to"]] <- as.character(lab[i,"type"])
  }
  return(t)
}
 
CreateSkeletonData_OFT <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(S1 = GetDistances(t,"nose","headcentre"))
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat$P1 <- as.integer(IsInZone(t,"nose","arena"))
  dat$P2 <- as.integer(IsInZone(t,"headcentre","arena"))
  dat <- as.data.frame(dat) 
  t$features <- dat
  
  return(t)
}


CreateSkeletonData_FST_v2 <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat <- as.data.frame(dat) 
  t$features <- abs(dat)
  
  return(t)
}

CreateAccelerationFeatures<- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat <- as.data.frame(dat) 
  t$features <- abs(dat)
  
  return(t)
}

ZscoreNormalizeSkeleton <- function(t, omit = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  change <- setdiff(names(t$features),omit)
  for(i in change){
    t$features[i] <- NormalizeZscore(t$features[i])
  }
  return(t)
}

PlotDensityPaths <- function(t,points,SDcutoff = 4, Title = "density path"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  out <- NULL
  for(i in points){
    data_plot <- t$data[[i]]
    xbreaks <- seq(floor(min(data_plot$x)), ceiling(max(data_plot$x)), by = 0.1)
    ybreaks <- seq(floor(min(data_plot$y)), ceiling(max(data_plot$y)), by = 0.1)
    data_plot$latbin <- xbreaks[cut(data_plot$x, breaks = xbreaks, labels=F)]
    data_plot$longbin <- ybreaks[cut(data_plot$y, breaks = ybreaks, labels=F)]
    out[[i]] <- ggplot(data = data_plot, aes(x,y)) + 
      stat_density_2d(data = data_plot, aes(latbin,longbin, fill=..density..), geom = "raster", contour = FALSE) + 
      scale_fill_gradient(name = "Time Density", low = "blue", high = "yellow") +
      geom_path(data = data_plot, aes(x,y, color = SDPlot((speed * t$fps),SDcutoff))) + 
      theme_bw() + 
      ggtitle(paste(i,Title, sep = " ")) + 
      scale_color_gradient2(name = paste("speed (",t$distance.units,"/s)",sep = ""), high = "white", low="black", mid = "black")
  }
  return(out)
}

SDPlot <- function(x,nSD){
  ifelse((mean(x) - x) / sd(x) > -nSD, x, mean(x) + nSD * sd(x))
}

SmoothLabel <- function(x, integration_period){
  types <-unique(x)
  mat <- NULL
  for (i in types){
    mat <- cbind(mat,periodsum(x == i, integration_period))
  }
  c <- apply(mat,1,FUN = which.max)
  return(types[c])
}

AddZonesToPlots <- function(p,z){
  for(i in 1:length(p)){
      for(j in z){
    p[[i]] <- p[[i]] + geom_path(data=j[c(1:nrow(j),1),],aes(x,y))
    }
  }
  return(p)
}

PlotZones <- function(t, zones = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("Zones not defined")
    return(NULL)
  }
  if(is.null(zones)){
    zones <- names(t$zones)
  }
  p <- ggplot()
  for(i in zones){
    dat <- t$zones[[i]]
    p <- p + geom_path(data=dat[c(1:nrow(dat),1),],aes(x,y))
  }
  return(p)
}

AddZones <- function(t,z){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    t$zones <- list()
    t$zones.invert <- list()
  }
  
  for(i in names(z)){
    t$zones[[i]] <- t$median.data[as.character(z[z[,i]!= "",i]),c("x","y")]
    t$zones.invert[[i]] <- FALSE
  }
  return(t)
}

AddBinData <- function(t, bindat = NULL, unit = "frame", binlength = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!is.null(bindat)){
    if(unit == "second"){
      bindat$from <- bindat$from * t$fps
      bindat$to <- bindat$to * t$fps
    }
    if(unit == "minute"){
      bindat$from <- bindat$from * t$fps * 60
      bindat$to <- bindat$to * t$fps * 60
    }
    t$bins <- bindat
  }else if(!is.null(binlength)){
    binlength <- binlength * ifelse(unit != "frame",t$fps,1) * ifelse(unit == "minute",60,1)
    t$bins <- NULL
    for(i in 1:ceiling(length(t$frames) / binlength)){
      t$bins <- rbind(t$bins, data.frame(bin = paste("bin",i,sep="."), 
                                         from = t$frames[(i-1)*binlength + 1], 
                                         to = t$frames[min(i*binlength,length(t$frames))]))
    }
  }else{
    warning("To add bin data either a the bindata or binlenght has to be added")
  }
  return(t)
}

BinAnalysis <- function(t,FUN, ...){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$bins)){
    warning("No bins defined, can not perform bin analysis")
    return(NULL)
  }
  
  Report <- NULL
  for(i in 1:nrow(t$bins)){
    tb <- CutTrackingData(t, keep.frames = t$bins[i,"from"]:t$bins[i,"to"])
    Report <- rbind(Report,data.frame(bin = t$bins[i,"bin"], FUN(tb,...)$Report))
  }
  return(Report)
}

FSTAnalysis <- function(t, cutoff_floating, integration_period = 0, points, Object){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  for(i in 1:length(t$data)){
    t$data[[i]]$delta_x <- integratevector(t$data[[i]]$x)
    t$data[[i]]$delta_y <- integratevector(t$data[[i]]$y)
    t$data[[i]]$speed <- sqrt(t$data[[i]]$delta_x ^2 + t$data[[i]]$delta_y^2)
    t$data[[i]]$acceleration <- integratevector(t$data[[i]]$speed)
  }
  
    temp <- t$data[as.character(t$point.info[t$point.info$PointType == Object,"PointName"])]
    acc <- NULL
    for(i in 1:length(temp)){
      acc <- cbind(acc, i = temp[[i]]$acceleration)
    }

    t$object <- list()
    t$object$movement <- abs(apply(acc,1,FUN = mean))
    t$object$is.floating <- avgmean(t$object$movement,integration_period) < cutoff_floating
    t$labels$cutoff.floating <- ifelse(t$object$is.floating == 1, "Floating","None")
  
    t$Report <- list()
    t$Report[["time.floating"]] <- sum(t$object$is.floating) / t$fps
    t$Report[["total.time"]] <- length(t$object$is.floating) / t$fps
    t$Report[["percentage.floating"]] <- sum(t$object$is.floating) / length(t$object$is.floating) * 100
    
    for(k in points){
      t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(t$data[[k]]$speed, na.rm = T)
      t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(t$data[[k]]$speed, na.rm = T) * t$fps
      t$Report[[paste(k, "distance.swiming", sep = ".")]] <- sum(t$data[[k]]$speed[!t$object$is.floating], na.rm = T)
      t$Report[[paste(k, "speed.swiming", sep = ".")]] <- mean(t$data[[k]]$speed[!t$object$is.floating], na.rm = T) * t$fps
      }
  
  return(t)
}  

HeadAngleAnalysis <- function(t, points = c("tailbase","neck","neck","nose"), angle_cutoff, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t$object <- list()
  t$object$head.angle.CW <- GetAngleClockwise(t,points[1],points[2],points[3],points[4])
  t$object$head.angle.total <- GetAngleTotal(t,points[1],points[2],points[3],points[4])
  t$object$head.angle.CW.degree <- (t$object$head.angle.CW + pi) * 180 / pi - 180
  t$object$head.angle.total.degree <- t$object$head.angle.total * 180 / pi
  t$object$head.tilted.CW <- t$object$head.angle.CW.degree > angle_cutoff
  t$object$head.tilted.CCW <- t$object$head.angle.CW.degree < -angle_cutoff
  
  t$Report <- list()
  t$Report[["average.head.angle.CW"]] <- mean(t$object$head.angle.CW.degree)
  t$Report[["time.head.tilted.CW"]] <- sum(avgbool(t$object$head.tilted.CW, integration_period)) / t$fps
  t$Report[["time.head.tilted.CCW"]] <- sum(avgbool(t$object$head.tilted.CCW, integration_period)) / t$fps
  
  return(t)
}

GetAngleClockwise <- function(t,p1,p2,p3,p4){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(sum(c(p1,p2,p3,p4) %in% names(t$data)) != 4){
    warning(paste("Points for angle measurement not available in Tracking data:",paste(setdiff(c(p1,p2,p3,p4), names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  
  a1 <- t$data[[p1]]
  a2 <- t$data[[p2]]
  b1 <- t$data[[p3]]
  b2 <- t$data[[p4]]
  ax <- a1$x - a2$x
  ay <- a1$y - a2$y
  bx <- b1$x - b2$x
  by <- b1$y - b2$y
  dot <- ax*bx + ay*by      # dot product between [x1, y1] and [x2, y2]
  det <- ax*by - ay*bx      # determinant
  res <- atan2(det, dot)  # atan2(y, x) or atan2(sin, cos)
  return(res)
}

GetAngleTotal <- function(t,p1,p2,p3,p4){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(sum(c(p1,p2,p3,p4) %in% names(t$data)) != 4){
    warning(paste("Points for angle measurement not available in Tracking data:",paste(setdiff(c(p1,p2,p3,p4), names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  
  a1 <- t$data[[p1]]
  a2 <- t$data[[p2]]
  b1 <- t$data[[p3]]
  b2 <- t$data[[p4]]
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

CalculateTransitions <- function(x,integration_period){
  x <- avgbool(x,integration_period)
  sum(append(0, (x[2:length(x)]!= x[1:length(x)-1])))
}

avgmean <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- mean(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

periodsum <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- sum(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

CreateTrainingSet <- function(t, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  x <- t$features
  x_window <- x[1:(nrow(x) - 2*integration_period),]
  
  if(integration_period > 0){
    for(i in (-integration_period + 1):integration_period){
      x_window <- cbind(x_window, x[(integration_period + i + 1):(nrow(x) - integration_period + i),])
    }
  }
  t$train_x <- as.matrix(x_window)
  t$train_y <- t$labels$manual[(integration_period + 1):(length(t$labels$manual) - integration_period)]
  t$ml_integration <- integration_period
  return(t)
}

CreateTestSet <- function(t, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  x <- t$features
  x_window <- x[1:(nrow(x) - 2*integration_period),]
  
  if(integration_period > 0){
    for(i in (-integration_period + 1):integration_period){
      x_window <- cbind(x_window, x[(integration_period + i + 1):(nrow(x) - integration_period + i),])
    }
  }
  t$train_x <- as.matrix(x_window)
  t$ml_integration <- integration_period
  return(t)
}

ClassifyBehaviors <- function(t,model, model_parameters){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t <- CreateTestSet(t, model_parameters$integration_period)
  t$labels$classifications <- model %>% predict_classes(t$train_x)
  t$labels$classifications <- c(rep(NA,model_parameters$integration_period), model_parameters$Feature_names[t$labels$classifications + 1], rep(NA,model_parameters$integration_period))
  return(t)
}

PlotLabels <- function(t, p.size = 2){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- NULL
  
  if(length(t$labels) == 0){
    warning("No labeling data available. can not plot")
    return(NULL)
  }
  
  for(i in names(t$labels)){
    dat <- rbind(dat,data.frame(seconds = t$seconds, behavior =  t$labels[[i]], type = i))
  }
  ggplot(data = na.omit(dat),aes(seconds, behavior, color = behavior)) + geom_point(size = p.size, shape = 124) + facet_grid(type~., scales = "free_y")
}

PlotZoneVisits <- function(t, points, zones = NULL, p.size = 2){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined")
    return(NULL)
  }
  if(is.null(zones)){
    zones <- names(t$zones)
  }
  if(length(setdiff(zones,names(t$zones))) > 0){
    warning("invalid zones")
    return(NULL)
  }
  if(length(setdiff(points,names(t$data))) > 0){
    warning("invalid points")
    return(NULL)
  }
  
  dat <- NULL
  for(j in points){
    for(i in zones){
      dat <- rbind(dat,data.frame(seconds = t$seconds, zone = ifelse(IsInZone(t,j,i,t$zones.invert[[i]]),i,NA), type = "automatic", points = j))
    }
  }
  
  ggplot(data = na.omit(dat),aes(seconds, zone, color = zone)) + geom_point(size = p.size, shape = 124) + facet_grid(points~.)
}

PlotPointData <- function(t, points = NULL, from = NULL, to = NULL, unit = "frame", type = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(points) & is.null(type)){
    points <- names(t$data)
  }
  
  if(unit == "second"){
    if(!is.null(from)){
      from <- t$frames[which(t$seconds >= from)[1]]
    }
    if(!is.null(to)){
      to <- t$frames[which(t$seconds >= to)[1]]
    }
  }
  if(is.null(from)){
    from = min(t$frames)
  }
  if(is.null(to)){
    to = max(t$frames)
  }
  
  if(!is.null(type)){
    points <- t$point.info[t$point.info$PointType == type,"PointName"]
  }
  
  range <- from:to
  
  p <- ggdraw()
  dim <- ceiling(sqrt(length(points)))
  nplot <- 0
  
  for(i in points){
    p <- p + draw_plot(ggplot(data = t$data[[i]][t$data[[i]]$frame %in% range,], aes(x,y, color = likelihood)) + geom_path() + ggtitle(i) + xlab(paste("x /",t$distance.units,sep = " ")) + ylab(paste("y /",t$distance.units,sep = " ")), 
                       x = (nplot %% dim / dim),
                       y = ((dim - 1)/ dim) - floor(nplot / dim) / dim, 
                       width = 1/dim, 
                       height = 1/dim)
    nplot <- nplot + 1
  }
  
  return(p)
}

RunPipeline <- function(files, path, FUN){
  out <- list()
  for(j in files){
  out[[paste(j)]] <- FUN(paste(path,j,sep = ""))
  }
  return(out)
}

CombineTrainingsData <- function(ts, shuffle =TRUE){
  if(IsTrackingData(ts)){
    temp <- ts
    ts <- list()
    ts[[paste(temp$filename)]] <- temp
  }
  
  train_x <- NULL
  train_y <- NULL
  for(i in names(ts)){
      if(is.null(ts[[i]]$train_x) | is.null(ts[[i]]$train_y)){
        warning(paste("File:",i,"is missing trainings data. data for this file not included"))
      }else{
        train_x <- rbind(train_x,ts[[i]]$train_x)
        train_y <- append(train_y,ts[[i]]$train_y)
      }
  }
  out <- PrepareMLData(train_x,train_y, shuffle)
  out$parameters$integration_period <- ts[[1]]$ml_integration
  return(out)
}

PrepareMLData <- function(x_train, y_train, shuffle = TRUE){
  if(nrow(x_train) != length(y_train)){
    warning("Unequal size of x_train and y_train")
    return(NULL)
  }
  out <- list()
  parameters <- list()
  x_train <- as.matrix(x_train)
  out$parameters$N_input <- ncol(x_train)
  out$parameters$N_features <- length(unique(y_train))
  out$parameters$Feature_names <- levels(as.factor(y_train))
  
  y_train_cat <- to_categorical(-1 + as.integer(as.factor(y_train)))
  
  if(shuffle){
  new_order <- sample(1:nrow(y_train_cat))
  x_train <- x_train[new_order,]
  y_train_cat <- y_train_cat[new_order,]
  }
  out$train_x <- x_train
  out$train_y <- y_train_cat
  
  return(out)
}

SmoothLabels <- function(t, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(t$labels) == 0){
    warning("No labels present. Returning original object")
      return(t)
  }
  for(i in names(t$labels)){
    t$labels[[i]] <- SmoothLabel(t$labels[[i]], integration_period)
  }
  return(t)
}

UnsupervisedClusteringKmeans <- function(ts, N_clusters = 20, Z_score_Normalize = TRUE){
  if(IsTrackingData(ts)){
    print("single file detected. Runing kmeans in single file mode")
    if(is.null(ts$train_x)){
      warning("No training or testing data present. Returning original data")
      return(ts)
    }
    if(Z_score_Normalize){
      test <- kmeans(NormalizeZscore(ts$train_x),centers = N_clusters)
    }
    else{
      test <- kmeans(ts$train_x,centers = N_clusters)
    }
    ts$labels$unsupervised <- c(rep(NA,j$ml_integration),as.character(test$cluster),rep(NA,j$ml_integration))
    return(ts)
  }
  
  allx <- NULL
  id <- NULL
  print("multiple files detected. Runing kmeans in multi file mode")
  for(j in names(ts)){
    if(is.null(ts[[j]]$train_x)){
      warning(paste("File:",j, "does not cotain any training data. Returning original data"))
      return(ts)
    }
    allx <- rbind(allx, ts[[j]]$train_x)
    id <- append(id, rep(paste(j),nrow(ts[[j]]$train_x)))
  }

  if(Z_score_Normalize){
    allx <- NormalizeZscore(allx)
  }
  test <- kmeans(allx,centers = N_clusters)
  
  for(j in names(ts)){
    print(j)
    clust <- test$cluster[id == j]
    ts[[j]]$labels$unsupervised <- c(rep(NA,ts[[j]]$ml_integration),as.character(clust),rep(NA,ts[[j]]$ml_integration))
  }
  return(ts)
}


PlotZoneSelection <- function(t,point,zones, invert = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!point %in% names(t$data)){
    warning("Invalid point")
    return(NULL)
  }
  if(!sum(zones %in% names(t$zones))){
    warning("Invalid zone(s)")
    return(NULL)
  }
  
  dat <- t$data[[point]]
  in.zone <- rep(FALSE,nrow(dat))
  for(i in t$zones[zones]){
    in.zone <- in.zone | (point.in.polygon(dat$x,dat$y,i$x,i$y) == 1)
  }
  if(invert){
    in.zone <- !in.zone
  }
  p <- ggplot(dat,aes(x,y, color = in.zone)) + geom_point()
  return(p)
}

MultiFileReport <- function(ts){
  if(IsTrackingData(ts)){
    stop("Expected a list() of TrackingData objects. You entered a single object")
  }
  out <- list()
  for(i in names(ts)){
    if(IsTrackingData(ts[[i]])){
      if(is.null(ts[[i]]$Report)){
        warning(paste("Object",i,"Does not contain any Report. omitting", sep = " "))
      }else{
        out <- rbindlist(list(out,append(c(file = ts[[i]]$filename), ts[[i]]$Report)),use.names = TRUE, fill = TRUE,idcol = F)
      }
    }else{
      warning(paste("List contains an element that is not of type TrackingData:",i,".No report produced for these", sep = " "))
    }
  }
  return(data.frame(out))
}

MultiFileBinanalysis <- function(ts, FUN, ...){
  if(IsTrackingData(ts)){
    stop("Expected a list() of TrackingData objects. You entered a single object")
  }
  out <- NULL
  for(i in names(ts)){
    if(IsTrackingData(ts[[i]])){
      binrep <- BinAnalysis(ts[[i]], FUN = FUN, ...)
      out <- rbind(out, data.frame(file = i, binrep))
    }else{
      warning(paste("Object",i,"is not of type Tracking data. omiting from Analysis", sep = " "))
    }
  }
  return(out)
}

PlotDensityPaths.Multi.PDF <- function(ts, points, filename = "DensityPathMulti",width = 10, height = 8, add_zones = FALSE, selected_zones = NULL, ...){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  out <- list()
  for(i in names(ts)){
    if(IsTrackingData(ts[[i]])){
        ps <- PlotDensityPaths(ts[[i]],points = points, Title = i, ...)
        
        if(add_zones){
          if(!is.null(ts[[i]]$zones)){
            if(is.null(selected_zones)){
              zones <- names(ts[[i]]$zones)
            }else{
              zones <- intersect(names(ts[[i]]$zones), selected_zones)
            }
            if(length(zones) == 0){
              warning(paste("in file:", i, "none of the indicated zones found in data", sep = " "))
            }
            else{
              ps <- AddZonesToPlots(ps,ts[[i]]$zones[zones])
            }
          }else{
            warning(paste("can not add zones to",i, ".Does not have zones defined", sep = " "))
          }
        }
        out <- append(out,ps)
    }else{
      warning(paste("List contains an element that is not of type TrackingData:",i,".No plot produced for these", sep = " "))
    }
  }
  pdf(paste(filename,".pdf",sep = ""), width = width, height = height)
  for(i in out){
    print(i)
  }
  dev.off()
  return(NULL)
}

PlotZoneVisits.Multi.PDF <- function(ts, points,filename = "ZoneVisitsMulti", width = 10, height = 8,...){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  
  pdf(paste(filename,".pdf",sep = ""), width = width, height = height)
  for(i in names(ts)){
    print(PlotZoneVisits(ts[[i]], points,...) + ggtitle(paste(i)))
  }
  dev.off()
  return(NULL)
}

PlotLabels.Multi.PDF <- function(ts, filename = "LabelsMulti", width = 10, height = 8,...){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  
  pdf(paste(filename,".pdf",sep = ""), width = width, height = height)
  for(i in names(ts)){
    print(PlotLabels(ts[[i]],...) + ggtitle(paste(i)))
  }
  dev.off()
  return(NULL)
}

OverviewPlot <- function(t, point){
  if(!IsTrackingData(t)){
    stop("Input needs to be a single object of type Tracking")
  }
  if(!point %in% names(t$data)){
    warning("point does not exist in Trackingobject")
  }
  
  rel_h <- c(0.2,1.5)
  title <- ggdraw() + draw_label(paste("Overview file",t$filename, sep = " "))
  p1 <- PlotDensityPaths(t,point)
  if(!is.null(t$zones)){
  p1 <- AddZonesToPlots(p1, t$zones)
  }
  p1 <- p1[[1]] + scale_y_reverse()
  
  if(length(t$labels) > 0){
    p2 <- PlotLabels(t) +  theme(legend.position = "none")
    rel_h <- append(rel_h, 0.5)
  }
  if(!is.null(t$zones)){
    p3 <- PlotZoneVisits(t,points = point) +  theme(legend.position = "none")
    rel_h <- append(rel_h, 0.5)
  }
  
  if((length(t$labels) > 0) & !is.null(t$zones)){
    return(plot_grid(title,p1,p2,p3,rel_heights = rel_h, ncol = 1))
  }else if(length(t$labels) > 0){
    return(plot_grid(title,p1,p2,rel_heights = rel_h, ncol = 1))
  }else if(!is.null(t$zones)){
    return(plot_grid(title,p1,p3,rel_heights = rel_h, ncol = 1))
  }
  else{
    return(plot_grid(title,p1,rel_heights = rel_h, ncol = 1))
  }
}

NormalizeZscore <- function(x){
  apply(x, 2, FUN = function(x){(x - mean(x)) / (sd(x))})
}

NormalizeZscore_median <- function(x){
  apply(x, 2, FUN = function(x){(x - median(x)) / (sd(x))})
}

integratevector <- function(x){
  if(length(x) < 2){
    stop("can  not integrate a vector of length < 2")
  }
  append(0, x[2:length(x)] - x[1:(length(x)-1)])
}

avgbool <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- mean(x[max(0,i-window):min(length(x), i + window)], na.rm = T) > 0.5
  }
  return(res)
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

EvaluateClassification <- function(ts, truth = "manual", compare = "classifications"){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  Report <- list()
  Report$files <- NULL
  for(i in names(ts)){
    if(is.null(ts[[i]]$labels)){
      stop(paste("file",i,"Does not have label data", sep = " "))
    }
    if(is.null(ts[[i]]$labels[[truth]])){
      stop(paste("file",i,"Does not have label type", truth, sep = " "))
    }
    if(is.null(ts[[i]]$labels[[compare]])){
      stop(paste("file",i,"Does not have label type", compare, sep = " "))
    }
    for(j in na.omit(unique(ts[[i]]$labels[[compare]]))){
      precision <- sum(ts[[i]]$labels[[compare]] == j & ts[[i]]$labels[[truth]] == j, na.rm = T) / sum(ts[[i]]$labels[[compare]] == j, na.rm = T)
      N_truth = sum(ts[[i]]$labels[[truth]] == j, na.rm = T)
      N_compare = sum(ts[[i]]$labels[[compare]] == j, na.rm = T)
      correct <- sum(ts[[i]]$labels[[compare]] == j & ts[[i]]$labels[[truth]] == j, na.rm = T)
      wrong <- sum(ts[[i]]$labels[[compare]] != j & ts[[i]]$labels[[truth]] == j, na.rm = T)
      recall <- sum(ts[[i]]$labels[[compare]] == j & ts[[i]]$labels[[truth]] == j, na.rm = T) / sum(ts[[i]]$labels[[truth]] == j, na.rm = T)
      Report$files <- rbind(Report$files,data.frame(file = i, 
                                                    label = j, 
                                                    accuracy = correct/(correct + wrong),
                                                    precision = correct / N_compare, 
                                                    recall = correct / N_truth,
                                                    correct = correct, 
                                                    wrong = wrong, 
                                                    N_truth = N_truth, 
                                                    N_compare = N_compare))
    }
    
  }
  
  Report$overall <- NULL
  for(i in unique(Report$files$label)){
    s <- apply(Report$files[Report$files$label == i,-c(1:5)],2,FUN = sum)
    Report$overall <- rbind(Report$overall, data.frame(label = i,
                                                       accuracy = s[1]/(s[1] +s[2]), 
                                                       precision =  s[1] / s[4], 
                                                       recall =   s[1] / s[3],
                                                       correct = s[1],
                                                       wrong = s[2],
                                                       N_truth = s[3],
                                                       N_compare = s[4]
                                                       ))
  }
  return(Report)
}

CorrelationPlotLabels <- function(ts, include = NULL){
  compare <- list()
  for(i in ts){
    compare <- rbindlist(list(compare,LabelReport(i)),use.names = TRUE, fill = TRUE,idcol = F)
  }
  compare <-as.data.frame(compare)
  
  if(is.null(include)){
    include <- names(compare)
  }
  
  corrplot(cor(as.matrix(na.replace(compare[,include]))),
           method="color",
           type = "upper",
           addCoef.col = "black")
}







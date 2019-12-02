setwd("G:/DLC/DLC paper 2019 Datapipeline final/OFT/")


library(ggplot2)
library(av)
library(cowplot)

render_dat <- readRDS("MLResults/ResDFs_final_reshuff.RDS")

render_dat <- render_dat$OFT_15_Furkan


p_text <- ggplot() + geom_text()
  
  
png("C:/Users/Deeplabcut/Documents/others/Lukas/pngFiles/input%03d.png", width = 1280, height = 360, res = 108)
for(i in seq(from = 2,to = nrow(render_dat), by = 1)){
  
  res2 <- render_dat[ifelse(i-1000 > 0, i-1000,1):(i),]
  
  p1 <- ggplot() + geom_path(data=res2,aes(time,as.numeric(avgbool(type == "Supported",10))), color = "red") + 
    geom_path(data=res2,aes(time,2.1 + as.numeric(avgbool(type == "Unsupported",10))), color = "blue") + 
    geom_path(data=res2,aes(time,1.05 + as.numeric(original.scored == "Supported")), color = "red4") + 
    geom_path(data=res2,aes(time,3.15  + as.numeric(original.scored == "Unsupported")), color = "blue4") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("time / s") + ylab("Behavior scored")
  
  
  
  print(i)
  print(p1)
}

dev.off()
png_files <- sprintf("C:/Users/Deeplabcut/Documents/others/Lukas/pngFiles/input%03d.png", 1:(i-1))
av::av_encode_video(png_files, 'OFT_15_Furkan.mp4', framerate = 25)

library(tidyverse)


samples <- read.csv("metadata/Danish_FullScale_WWTPs/Danish-WWTPs-lat-long.csv")

p <- ggplot() + borders("world", colour="black", fill="gray50") + coord_fixed(xlim=c(8, 13), ylim=c(54, 58))
p1 <- p + geom_point(data=samples, aes(x=lon, y=lat)) 
p1
p1 + theme_classic()

p1 + geom_text(data=samples, aes(x=lon, y=lat, label=plant), size = 3, col = "black", nudge_y=0.05)
                                                            
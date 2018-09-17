data <- readRDS("./rds/Freichel.rds")
own_ss <- readRDS("rds/own_sample_sheet.rds")

colnames(data)[2:97] <-own_ss$V3 # change to useful colnames in original dataset

#plot(d[,a==TRUE], d[,a=TRUE], log='xy')      #would plot one sample against another
#par(mfrow=c(length(a),length(a)))


#quality control via clusterization

hc <- hclust(dist(t(log(1+data[1:10,-1]))))

png('./Figures/hcluster.png',width = 40, height = 20, units = "cm", res=90)

plot(hc)

dev.off()



#dev.off()
#cor(log(1+d[-1]))
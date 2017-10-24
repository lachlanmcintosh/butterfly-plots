rm(list=ls())

library("plyr")
library("DNAcopy")
library("aroma.affymetrix")
library("PSCBS")
library("truncnorm")
library("Rsolnp")
library("shiny")
library("ggplot2")
library("gridExtra")
library("lattice")
library("reshape2")
library("RColorBrewer")
library("R.utils")
library(ggtern)
library(MASS)
library(Hmisc)
library(mgcv)
library(dpmixsim)
library(igraph)
#require("ggrepel")
# if (!require("shinyBS")) install.packages("shinyBS")
library("shinyBS")
source("GridITH_functions.R")

CN_track_description <- "A copy number track shows the major (in colour) and minor (in black) copy numbers on the vertical axis against location on the horizontal axis. This plot is interactive with the butterfly plot in that you can highlight points (by dragging the mouse) on one plot and they will also be highlighted on the other."
bar_plot_description <- ""
data_table_description <- ""
density_plots_description <- ""
density_change_plots_description <- ""
gam_plots_description <- ""
subclone_barchart_description <- ""
Mahanalobis_control_plotui_description <- ""
GridIT_dialPlotui_description <- ""
butterfly_plot_description <- "A butterfly plot shows the major (vertical axis) and minor (horizontal axis) copy numbers for each location in the genome. This plot is <b> interactive</b> with the copy number track plot in that you can highlight points (<b>by dragging the mouse</b>) on one plot and they will also be highlighted on the other."

options(shiny.maxRequestSize=1024^3) 

ninetyNinthPercentile <- function(mat,var,weights="length"){
  max(round(unlist(wtd.quantile(as.numeric(unlist(c(mat[,var])),na.rm=T),weights = as.numeric(unlist(c(mat[,weights])),na.rm=T),probs=0.9999)),2)+0.1,0)
}

firstPercentile <- function(mat,var,weights="length"){
  min(round(unlist(wtd.quantile(as.numeric(unlist(c(mat[,var])),na.rm=T),weights = as.numeric(unlist(c(mat[,weights])),na.rm=T),probs=0.0001)),2)-0.1,0)
}

roundedMax <- function(mat,iter){
  round(max(as.numeric(unlist(c(mat[,x(iter)],mat[y(iter)])),na.rm=T)),2)
}

# PRIVATE FUNCTIONS
readCBS = function(filename, printfolder = NULL, fixIMBA = FALSE,dropdown=TRUE){
  if(!dropdown){ 
    da = filename
  } else{
    da = loadObject(filename)
    # da <<- loadObject(filename) 
    mat = da$output
  }
  names(mat)[names(mat) == 'c1Mean'] = 'a1'
  names(mat)[names(mat) == 'c2Mean'] = 'a2'
  names(mat)[names(mat) == 'tcnStart'] = 'Start'
  names(mat)[names(mat) == 'tcnEnd'] = 'End'
  names(mat)[names(mat) == 'chromosome'] = 'chr'
  mat$length = mat$End - mat$Start
  mat = subset(mat, !is.na(length) & chr<23 & !is.na(a2) & !is.na(a1))
  mat$W = mat$length/sum(mat$length)
  mat$tot = pmax(mat$a1,0) + mat$a2
  mat$imba = 2*(mat$a2/(mat$tot)-.5)
  mat$chr <- as.factor(mat$chr)
  
  # chr_sizes <- sapply(unique(mat$chr),function(x) max(mat[which(mat$chr == x),"End"])) - 
  # sapply(unique(mat$chr),function(x) min(mat[which(mat$chr == x),"Start"]))
  # mat$temp <- sapply(mat$chr, function(x) sum(chr_sizes[1:x])) - 
  # sapply(mat$chr, function(x) chr_sizes[x])+ (mat$Start +mat$End)/2
  # mat$location <- mat$temp/max(mat$temp) *1.032
  
  # mat$height = rep(0.5, nrow(mat))
  
  #obj$matCBS = mat
  # if(fixIMBA) {
  #   #Find index of points to move and theta
  #   den = density(atan(mat$a1/mat$a2), bw = 0.05)
  #   #plot(den, main = 'Angle theta')
  #   der = (den$y[-1] - den$y[-length(den$y)])/(den$x[2]-den$x[1])
  #   pos = which(der >= 0)
  #   neg = which(der < 0)
  #   cut = neg[neg<pos[length(pos)]]
  #   cut = cut[length(cut)]
  #   #abline(v = den$x[cut], col = 'red')
  #   theta = den$x[pos[length(pos)]]
  #   #abline(v = theta, lty = 'dotted')
  #   moveix = which(atan(mat$a1/mat$a2) > den$x[cut])
  #   
  #   g = c(1, tan(theta))
  #   ng = sqrt(sum(g * g))
  #   g = g/ng
  #   p.prime = mat$a2 * g[1] + mat$a1 * g[2]
  #   x.prime = p.prime * cos(theta)
  #   y.prime = p.prime * sin(theta)
  #   #lims = c(0, max(mat$a1, mat$a2))
  #   #plot(mat$a2, mat$a1, xlim = lims, ylim = lims,
  #   #main = 'Red points move to blue', xlab = 'Major array CN',
  #   #ylab = 'Minor array CN', pch = 16, col = col.transp('grey'))
  #   #abline(0,1)
  #   #lines(c(0,100*g[1]), c(0, 100*g[2]), col = 'red')
  #   #points(mat$a2[moveix], mat$a1[moveix], pch = 16,
  #   #col = col.transp('red'))
  #   d = (x.prime - y.prime)/2
  #   #points((mat$a2-d)[moveix], (mat$a1+d)[moveix],
  #   #col = col.transp('blue'), pch = 16)
  #   #if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
  #   
  #   
  #   #Rearrange a1 and a2 according to size
  #   mat$a1[moveix] = pmin((mat$a2-d)[moveix], (mat$a1+d)[moveix])
  #   mat$a2[moveix] = pmax((mat$a2-d)[moveix], (mat$a1+d)[moveix])
  #   mat$fixedIMBA = TRUE
  #   mat$allelic_balance = FALSE
  #   mat$allelic_balance[moveix] = TRUE
  #   #par(mfrow = c(1,1))
  #   #obj
  # }
  mat$tempx <- mat$a_1 <- mat$a1
  mat$tempy <- mat$b_1 <- mat$a2
  mat$invalid_length <- FALSE
  mat$invalid_number <- FALSE
  da$output = mat
  
  return(da)
}
## end readCBS

z <- function(allele,iter,extras=NULL){
  return(paste(c(allele,as.character(iter),as.character(extras)),collapse="_"))
}
x <- function(iter,extras=NULL){
  return(z("a",iter,extras))
}
y <- function(iter,extras=NULL){
  return(z("b",iter,extras))
}
xnext <- function(iter,extras=NULL){
  return(z("a",iter+1,extras))
}
ynext <- function(iter,extras=NULL){
  return(z("b",iter+1,extras))
}

#PLOTTING FUNCTIONS
number_ticks <- function(n) {function(limits) pretty(limits, n)}

get_CN_track <- function(mat,xlimits,var1,var2,centromeres,num_rows){
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "Start", y = var1, xend = "End", yend = var1,col="chr"),size=2)
  
  if(!is.null(var2)) g <- g + geom_segment(aes_string(x = "Start", y = var2, xend = "End", yend = var2),size=1)
  
  g <- g + geom_vline(data=centromeres,aes(xintercept=x,col=chr))
  if(num_rows ==1){ 
    g <- g+ facet_grid(.~chr,scales="free_x",space="free")
  }
  else{ 
    g <- g+ facet_wrap(~chr,nrow=num_rows)
  }
  g <- g+
    scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+
    theme_bw()+
    theme(axis.text.x = element_blank(),legend.position="none",strip.background = element_blank())+ 
    xlab("50 MB ticks")+
    scale_y_continuous(breaks=number_ticks(20))+
    ylab("major/minor CN")+ggtitle("CN track") + coord_cartesian(ylim=c(xlimits[[1]],xlimits[[2]]))
  return(g)
}

summary_histogram <- function(mat,xlimits,var1,var2){
  g1 <- ggplot() + geom_histogram(aes(x=mat[,var1],weight=mat[,"length"],col="major"))+theme_bw()+coord_flip()+xlim(xlimits[[1]],xlimits[[2]])
  g2 <- ggplot() + geom_histogram(aes(x=mat[,var2],weight=mat[,"length"],col="minor"))+theme_bw()+coord_flip()+xlim(ylimits[[1]],ylimits[[2]])
  grid.arrange(g1,g2,nrow=1)
}

# get_box_plot <- function(mat,var1,var2,angle1,angle2,ratio,grid_spacing,xlimits){
#   T1 <- matrix(c(1,tan(angle1),tan(angle2),1),nrow=2,ncol=2)
#   
#   xs <- seq(0,xlimits/(1+tan(angle1)),by=grid_spacing)
#   ys <- seq(0,xlimits/(1+tan(angle1)),by=grid_spacing)
#   transformed_grid <- as.data.frame(as.matrix(expand.grid(x=xs, y=ys)) %*% t(T1),col.names = c("x","y")) # a grid of all lattice values
#   colnames(transformed_grid) <- c("x","y")
#   transformed_grid_right <- as.data.frame(as.matrix(expand.grid(x=xs+grid_spacing,y=ys)) %*% t(T1),col.names = c("x","y")) # the same grid starting one unit across
#   transformed_grid_above <- as.data.frame(as.matrix(expand.grid(x=xs,y=ys+grid_spacing)) %*% t(T1),col.names = c("x","y")) # the same grid starting one unit above
#   
#   original_grid <- as.data.frame(as.matrix(expand.grid(x=xs, y=ys)),col.names = c("x","y")) # a grid of all lattice values
#   colnames(original_grid) <- c("x","y")
#   original_grid_right <- as.data.frame(as.matrix(expand.grid(x=xs+grid_spacing,y=ys)),col.names = c("x","y")) # the same grid starting one unit across
#   original_grid_above <- as.data.frame(as.matrix(expand.grid(x=xs,y=ys+grid_spacing)),col.names = c("x","y")) # the same grid starting one unit above
#   
#   transformed_boxes <- rbind(cbind(transformed_grid,transformed_grid_right),cbind(transformed_grid,transformed_grid_above)) # a data frame of all vectors pointing to the right and pointing above
#   colnames(transformed_boxes) <- c("x","y","x2","y2")
#   
#   original_boxes <- rbind(cbind(original_grid,original_grid_right),cbind(original_grid,original_grid_above))
#   colnames(original_boxes) <- c("x","y","x2","y2")
#   
#   g1 <- ggplot()+xlim(0,xlimits)+ylim(0,xlimits)+#geom_point(aes(x=var1,y=var2),alpha=0.1)+
#     geom_point(aes(x=transformed_grid$x,y=transformed_grid$y,col="transformed"))+
#     geom_segment(aes(x=transformed_boxes$x,y=transformed_boxes$y,
#                      xend=transformed_boxes$x2,yend=transformed_boxes$y2,col="transformed"))+
#     geom_point(aes(x=original_grid$x,y=original_grid$y,col="original"))+
#     geom_segment(aes(x=original_boxes$x,y=original_boxes$y,xend=original_boxes$x2,yend=original_boxes$y2,col="original"))
#   
#   g1
# }

# get_density_plot <- function(mat,allele,miniter,maxiter,density_adjust,limits,type=NULL){ # can put type in here to be the changes.
#   if(allele == "x") dat <- data.frame(value =mat[,x(miniter)],iteration =rep(as.character(miniter),nrow(mat)))
#   if(allele == "y") dat <- data.frame(value =mat[,y(miniter)],iteration =rep(as.character(miniter),nrow(mat)))
# 
#   if(maxiter >1 ){
#     for(iter in (miniter+1):maxiter){
#       if(allele == "x") tem <- data.frame(value =mat[,x(iter)],iteration =rep(as.character(iter),nrow(mat)))
#       if(allele == "y") tem <- data.frame(value =mat[,y(iter)],iteration =rep(as.character(iter),nrow(mat)))
#       dat <- rbind(dat,tem)
#     }
#   }
#   ggplot(dat) + geom_density(aes(value,col=iteration),adjust=density_adjust)+ theme(legend.position="none")
# }

right_minus_left <- function(var){
  paste(var,"right","minus","left",sep="_")
}

right <- function(var){
  paste(var,"right",sep="_")
}

get_max_iter <- function(mat){
  values = colnames(mat)
  all = grep(pattern = "^a_[[:digit:]]$",x=values,perl=T)
  val <- max(as.numeric(sub(pattern = "^a_",replacement="",x=values[all],perl=T)))
  if(is.null(val)) return(1)
  return(val)
}
get_change_density_plot <- function(mat,allele,miniter,maxiter,density_adjust,limits,type=NULL){
  for(i in 1:maxiter){
    xnext = paste(x(i),"next",sep="_")
    xdiff = paste(x(i),"diff",sep="_")
    mat[,xnext] <- c(mat[2:nrow(mat),x(i)],mat[1,x(i)])
    mat[,xdiff] <- mat[,x(i)]-mat[,xnext]
    ynext = paste(y(i),"next",sep="_")
    ydiff = paste(y(i),"diff",sep="_")
    mat[,ynext] <- c(mat[2:nrow(mat),y(i)],mat[1,y(i)])
    mat[,ydiff] <- mat[,y(i)]-mat[,ynext]
  }
  
  g <- ggplot(mat)
  for(i in 1:maxiter){
    if(allele == "x") g <- g + geom_density(aes_string(paste(x(i),"diff",sep="_"),col=factor(i)),adjust = density_adjust)
    if(allele == "y") g <- g + geom_density(aes_string(paste(y(i),"diff",sep="_"),col=factor(i)),adjust = density_adjust)
  }
  
  # if(allele == "x") dat <- data.frame(value =mat[,x(miniter)],iteration =rep(as.character(miniter),nrow(mat)))
  # if(allele == "y") dat <- data.frame(value =mat[,y(miniter)],iteration =rep(as.character(miniter),nrow(mat)))
  # 
  # if(maxiter >1 ){
  #   for(iter in (miniter+1):maxiter){
  #     if(allele == "x") tem <- data.frame(value =mat[,x(iter)],iteration =rep(as.character(iter),nrow(mat)))
  #     if(allele == "y") tem <- data.frame(value =mat[,y(iter)],iteration =rep(as.character(iter),nrow(mat)))
  #     dat <- rbind(dat,tem)
  #   }
  # }
  # ggplot(dat) + geom_density(aes(value,col=iteration),adjust=density_adjust)
  if(allele =="y"){g <- g+coord_flip()}
  g+theme_minimal()+theme(legend.position="none")
}
get_density_plot <- function(mat,allele,miniter,maxiter,density_adjust,limits,type=NULL){
  g <- ggplot(mat)
  for(i in 1:maxiter){
    if(allele == "x") g <- g + geom_density(aes_string(x(i),col=factor(i)),adjust = density_adjust)
    if(allele == "y") g <- g + geom_density(aes_string(y(i),col=factor(i)),adjust = density_adjust)
  }
  if(allele =="y"){g <- g +coord_flip()}
  g+theme_minimal()+theme(legend.position="none")
}
# get_density_plot <- function(mat,allele,miniter,maxiter,density_adjust,limits,type=NULL){ # can put type in here to be the changes.
#   if(allele == "x") dat <- data.frame(value =mat[,x(miniter)],iteration =rep(as.character(miniter),nrow(mat)))
#   if(allele == "y") dat <- data.frame(value =mat[,y(miniter)],iteration =rep(as.character(miniter),nrow(mat)))
#   
#   if(maxiter >1 ){
#     for(iter in (miniter+1):maxiter){
#       if(allele == "x") tem <- data.frame(value =mat[,x(iter)],iteration =rep(as.character(iter),nrow(mat)))
#       if(allele == "y") tem <- data.frame(value =mat[,y(iter)],iteration =rep(as.character(iter),nrow(mat)))
#       dat <- rbind(dat,tem)
#     }
#   }
#   ggplot(dat) + geom_density(aes(value,col=iteration),adjust=density_adjust)
# }

get_density <- function(mat,iter,bandwidth,num_grid_points,min_density_height){
  # find the kernal density estimator
  twoddensity <- kde2d.weighted(x=mat[,x(iter)],y=mat[,y(iter)],n=num_grid_points,h=bandwidth,w=mat[,"tcnNbrOfHets"])
  
  # find the local maxima
  z <- as.matrix(twoddensity$z)
  local_maxima <- as.data.frame(which(z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[1:(nrow(z)-2),2:(ncol(z)-1)] &
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[3:(nrow(z)),2:(ncol(z)-1)] &
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[2:(nrow(z)-1),1:(ncol(z)-2)] &
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[2:(nrow(z)-1),3:(ncol(z))]&
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[3:(nrow(z)),3:(ncol(z))]&
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[1:(nrow(z)-2),1:(ncol(z)-2)]&
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[1:(nrow(z)-2),3:(ncol(z))]&
                                        z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[3:(nrow(z)),1:(ncol(z)-2)],
                                      arr.ind = T))
  # array indicies are off by 1,1 because which reports the bottom corner of the square that we are comparing.
  local_maxima[,1] <- local_maxima[,1] +1
  local_maxima[,2] <- local_maxima[,2] +1
  # local_maxima$z is the height of the density at that x,y coordinate.
  local_maxima$z <-  z[as.matrix(local_maxima)]
  
  # arbritrary thresholding on the local_maxima
  local_maxima <- local_maxima[which(local_maxima$z>min_density_height),]
  
  local_maxima$x <- twoddensity$x[local_maxima[,1]]
  local_maxima$y <- twoddensity$y[local_maxima[,2]]
  # order the points by the most dense
  local_maxima = local_maxima[order(-local_maxima$z),]
  
  mat$xchange <- mat[,x(iter)] - c(mat[2:nrow(mat),x(iter)],NA)
  mat$ychange <- mat[,y(iter)] - c(mat[2:nrow(mat),y(iter)],NA)
  
  xchangedens <- density(abs(mat$xchange),na.rm=T)
  ychangedens <- density(abs(mat$ychange),na.rm=T)
  
  xc <- xchangedens$y
  bestxspacinglist <- xchangedens$x[which(c(FALSE,xc[2:(length(xc)-1)] > xc[1:(length(xc)-2)] & xc[2:(length(xc)-1)] > xc[3:(length(xc))],FALSE))]
  # plot(xchangedens)
  
  if(length(bestxspacinglist) >=2){
    if(abs(bestxspacinglist[2]-bestxspacinglist[1])>0.2){
      bestxspacing = 1
    } else{
      bestxspacing = bestxspacinglist[2]-bestxspacinglist[1]
    }}
  else{
      bestxspacing = 1
    }
  
  yc <- ychangedens$y
  bestyspacinglist <- ychangedens$x[which(c(FALSE,yc[2:(length(yc)-1)] > yc[1:(length(yc)-2)] & yc[2:(length(yc)-1)] > yc[3:(length(yc))],FALSE))]
  if(length(bestyspacinglist) >=2){
    if(abs(bestyspacinglist[2]-bestyspacinglist[1])>0.2){
      bestyspacing = 1
    } else{
      bestyspacing = bestyspacinglist[2]-bestyspacinglist[1]
    }}
  else{
    bestyspacing = 1
  }

  dens <- density(abs(local_maxima$y/local_maxima$x),adjust=0.1)
  maxima <- which(c(F,c(dens$y[2:(length(dens$y)-1)]) > c(dens$y[1:(length(dens$y)-2)]) & c(dens$y[2:(length(dens$y)-1)]) > c(dens$y[3:(length(dens$y))]),F))
  
  diagonal_pts = local_maxima[which(abs(local_maxima$y/local_maxima$x -dens$x[maxima][1]) < 0.25),]
  if(nrow(diagonal_pts)>1){
    oneone = diagonal_pts[which(diagonal_pts$z == max(diagonal_pts$z)),]
  } else{
    oneone = local_maxima[1,]
  }
  if( ( (oneone$x-1)^2 + (oneone$y-1)^2 )^0.5 > 0.2){
    oneone$x = 1
    oneone$y = 1
  }
  
  local_maxima$xtrue <- round((local_maxima$x - oneone$x)/bestxspacing) + 1
  local_maxima$ytrue <- round((local_maxima$y - oneone$y)/bestyspacing) + 1
  local_maxima <- local_maxima[order(-local_maxima$z),]
  
  # ggplot()+geom_point(aes(mydata$a1,mydata$a2,size=mydata$length),alpha=0.25)+
  # geom_point(aes(x=local_maxima$x,y=local_maxima$y,size=local_maxima$z),col="red")+
  # coord_fixed()+geom_density2d(aes(x=mydata[,"a1"],y=mydata[,"a2"]))+ylim(0,1.5)
  
  local_maxima$name <- paste(as.character(local_maxima$xtrue),as.character(local_maxima$ytrue))
  
  local_maxima$order <- order(-local_maxima$z)
  local_maxima <- local_maxima[local_maxima$order,]
  local_maxima[which(duplicated(local_maxima$name)),c("xtrue","ytrue")] <- NA
  
  length(unique(local_maxima$xtrue))
  length(unique(local_maxima$ytrue))
  
  # ggplot()+geom_point(aes(mydata$a1,mydata$a2,size=mydata$length),alpha=0.25,col="white")+
  #   geom_point(aes(x=local_maxima$x,y=local_maxima$y,size=local_maxima$z),col="red")+
  #   coord_fixed()+geom_density2d(aes(x=mydata[,"a1"],y=mydata[,"a2"]))+ylim(0,1.5)+
  #   geom_text(aes(x=local_maxima$x,y=local_maxima$y,label=local_maxima$name),col="black",size=7)
  
  return(list(local_maxima = local_maxima,z=z))
}




use_labels_with_gam <- function(mat,iter,local_maxima,z){
  mat$x <- mat[,x(iter)]
  mat$y <- mat[,y(iter)]
  
  gx <- gam(xtrue ~ te(x,y,k=3),data=local_maxima,weights=z)# s(best_maxima$x,k=3)+s(best_maxima$y,k=5),weights=best_maxima$z)
  local_maxima$xtrue2 <- predict(gx,newdata =local_maxima)

  #*#
  #mat[,x(iter+1)] <- predict(gx,newdata =mat)
  mat[,"tempx"] <- predict(gx,newdata =mat)
  
  gy <- gam(ytrue ~ te(y,x,k=3),data=local_maxima,weights=z)# s(best_maxima$x,k=3)+s(best_maxima$y,k=5),weights=best_maxima$z)
  local_maxima$ytrue2 <- predict(gy,newdata =local_maxima,weights=z)
  #*# mat[,y(iter+1)] <- predict(gy,newdata =mat)
  mat[,"tempy"] <- predict(gy,newdata =mat)
  return(mat)
  
  
  # vis.gam(gx)
  # vis.gam(gy)
  
  # ggplot()+geom_point(aes(mydata$xtrue2,mydata$ytrue2,size=mydata$length),alpha=0.25,col="black")+
  #   geom_point(aes(x=local_maxima$xtrue2,y=local_maxima$ytrue2,size=local_maxima$z),col="red")+
  #   coord_fixed()+geom_density2d(aes(x=mydata[,"xtrue2"],y=mydata[,"ytrue2"]))+ylim(-0.1,5)+xlim(-0.1,3)

  mydata$dist_to_grid <- ((mydata$xtrue2 - round(mydata$xtrue2))^2 +(mydata$ytrue2 - round(mydata$ytrue2))^2)^0.5
  
  # g1 <- ggplot()+geom_density(aes(abs(mydata$xtrue2 - round(mydata$xtrue2)),col=factor(round(log(mydata$tcnNbrOfHets)/3))))
  # g2 <- ggplot()+geom_density(aes(abs(mydata$ytrue2 - round(mydata$ytrue2)),col=factor(round(log(mydata$tcnNbrOfHets)/3))))
  
  # grid.arrange(g1,g2)
  # ggplot(mydata)+geom_point(aes(x=log(tcnNbrOfHets),y=dist_to_grid))+
    # geom_density2d(aes(x=log(tcnNbrOfHets),y=dist_to_grid))
  
  # trick to enhance heterogeneity:
  mydata$names <- paste(round(mydata$xtrue2),round(mydata$ytrue2))
  mydata$density <- mydata$tcnNbrOfHets/mydata$length 
  
  tempdata <- mydata[which(mydata$xtrue2<4 & mydata$ytrue2 < 4),]
  # ggplot(tempdata)+geom_density(aes(abs(abs(ytrue2 - round(ytrue2))-abs(xtrue2 - round(xtrue2))),col=factor(round(log(tcnNbrOfHets)/3))))+xlim(0,1)+facet_wrap(~factor(names))
  
  
  # ggplot(tempdata)+geom_density(aes(abs(abs(ytrue2 - round(ytrue2))-abs(xtrue2 - round(xtrue2))),col=factor(round(log(length)/3))))+xlim(0,0.5)+facet_wrap(~factor(names))
  
  
  # ggplot(tempdata)+geom_density(aes(abs(abs(ytrue2 - round(ytrue2))-abs(xtrue2 - round(xtrue2))),col=factor(round(log(1/density)))))+xlim(0,0.5)+facet_wrap(~factor(names))
  
  # ggplot(mydata)+geom_point(aes(x=log(density),y=dist_to_grid))+
    # geom_density2d(aes(x=log(density),y=dist_to_grid))
  mat
}



find_opt_starting_values <- function(mat,iter,ADJUST,changes=T){
  G=10
  spread = 10
  valid = mat[1:(nrow(mat)-1),"chr"] == mat[2:nrow(mat),"chr"]
  # print(valid)
  # print(head(mat))
  best_area <- 0
  best_area_changes <- 0
  for(theta1 in pi/spread*seq(-G,G)/G){
    for(theta2 in pi/spread*seq(-G,G)/G){
      # print(theta1)
      # print(theta2)
      # Then you have two lines with gradients tan(theta1) and tan(theta2)
      T1 <- matrix(c(1,tan(theta1),tan(theta2),1),nrow=2,ncol=2)
      new <- as.data.frame(as.matrix(mat[,c(x(iter),y(iter))]) %*% solve(t(T1)))
      # print(head(new))
      # print(c(mean(mat[,x(iter)],na.rm=T),mean(new[,1],na.rm=T),mean(mat[,y(iter)],na.rm=T),mean(new[,2],na.rm=T)))
      new[,1] <- new[,1] * mean(mat[,x(iter)],na.rm=T)/mean(new[,1],na.rm=T)
      new[,2] <- new[,2] * mean(mat[,y(iter)],na.rm=T)/mean(new[,2],na.rm=T)
      # print(head(new))
      
      if(changes){
        new[,3] <- c(NA,new[1:(nrow(new)-1),1] - new[2:nrow(new),1])
        new[,4] <- c(NA,new[1:(nrow(new)-1),2] - new[2:nrow(new),2])
      } else{
        new[,3] <- new[,1]
        new[,4] <- new[,2]
      }
      # print(head(new))
      new <- new[which(valid),]
      new <- as.matrix(new)
      #area <- sum((density(new[,1],adjust=ADJUST)$y)^2) + sum((density(new[,2],adjust=ADJUST)$y)^2)
      as <- density(new[2:nrow(new),3],adjust=ADJUST)
      bs <- density(new[2:nrow(new),4],adjust=ADJUST)
      
      area <- sum((as$y/length(as$x))^2) + sum((bs$y/length(bs$x))^2)
      
      if(area > best_area){
        # print(c(area,theta1,theta2))
        best_theta1 <- theta1
        best_theta2 <- theta2
        best_area <- area
      }
    }
  }
  optimal <- list("best_theta1" = best_theta1,"best_theta2" = best_theta2)#,
                  # "best_theta1_c" = best_theta1_c,"best_theta2_c" = best_theta2_c)
  return(optimal)
}

# a funciton to set the number of xticks when using ggplot
BASE_FUNCTIONS = "/Users/lmcintosh/Documents/GRIDITH_AUTOMATION"
#BASE_FUNCTIONS = "."
# load a file containing the centromeres:
load(paste(BASE_FUNCTIONS,"/centromeres.RData",sep=""))
colnames(centromeres) <- c("chr","x")

Mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux)); 
  ux[tab == max(tab)] 
}

DC <- function(X,Y,L,maxiter=400,burnin=200,M=2){
  L <- ceil(L^0.5)
  # L <- 1
  index = rep(c(1:length(X)),L)
  
  X <- rep(X,L)
  Y <- rep(Y,L)

  Xc <- round(X)
  Yc <- round(Y)
  # Xc and Yc tell you the closest point and Xnc, Ync likewise the next closest
  Xnc <- Xc
  Ync <- Yc
  Xnc <- Xc + sign(X-Xc)
  Ync <- Yc + sign(Y-Yc)
  
  # great now calculate the tumour cellular fraction, say that you are going up if you are outside the 00 to 11 square
  Xup <- apply(cbind(Xc,Xnc),1,max)
  Yup <- apply(cbind(Yc,Ync),1,max)
  Xdown <- apply(cbind(Xc,Xnc),1,min)
  Ydown <- apply(cbind(Yc,Ync),1,min)
  
  Xparent <- Xdown
  Yparent <- Ydown
  Xparent[which(X<0)] <- Xup[which(X<0)]
  Yparent[which(Y<0)] <- Yup[which(Y<0)]
  
  xBAF <- abs(X-Xparent)
  yBAF <- abs(Y-Yparent)
  BAF <- c(xBAF,yBAF)
  
  # want to do a rotation on the baf so that we do not cut a cluster in half
  # want to find the smallest part of dens
  dens <- density(BAF)
  rotation_fac <- dens$x[which(dens$y == min(dens$y))]
  # ok now rotate by this guy
  BAF <- (BAF-rotation_fac)%%1

  rec <- burnin;
  ngrid <- 100
  res <- dpmixsim(BAF,maxiter=maxiter, rec=rec,nclinit=1,M=M,upalpha=0)#a0=1,b0=2)
  z <-  postdpmixciz(BAF, res=res, rec=rec, ngrid=ngrid, plot=T)
  cx   <- postdataseg(BAF, z, ngrid=ngrid)
  
  BAFx=cx[1:(length(cx)/2)]
  BAFy=cx[(length(cx)/2+1):length(cx)]
  
  main_cluster = Mode(cx)
  BAF = BAFx
  BAF[which(BAFy != main_cluster)] <- BAFy[which(BAFy != main_cluster)]
  BAF[which(BAFx != main_cluster)] <- BAFx[which(BAFx != main_cluster)]
  # BAF[which(BAFy != main_cluster & BAFx != main_cluster)] <- NA
  # BAF[which(BAFy == BAFx )] <- BAFx[which(BAFy == BAFx )]
  
  t <- data.frame("old"=unique(BAF),"new"=1:length(unique(BAF)))
  CLUST <- sapply(BAF,function(x) t$new[which(x == t$old)])
  return(CLUST[which(!duplicated(index))])
}

# find all the files
# dir = "/Users/lmcintosh/chris_CBSfits_example_copy/"
# full_files <- sapply(list.files(path=dir, pattern=".RData$"), function(x) paste(dir,x[[1]],sep=''))
# mat <- readCBS(full_files[[1]])
# str(mat)
# sum(mat$output$tcnNbrOfSNPs^0.5)
# set the first file to the first
# Define server logic for slider examples
shinyServer(function(input,output,session) {
  # Reactive expression to compose a data frame containing all of the values
  v <- reactiveValues(file = "",folder = "")
  
  update_iter <- function(){
    req(v$mat,cancelOutput = TRUE)
    if(input$automatically_update){
      v$iter <- v$iter + 1
      v$iter_max = max(v$iter,v$iter_max)    
      v$mat[,x(v$iter)] <- v$mat$tempx
      v$mat[,y(v$iter)] <- v$mat$tempy
    }
  }
  
  observe({
    new_file = req(input$file1)
    if(v$file != new_file$name){ # if this new file is different to the old file update the state.
      v$file = new_file$name
      v$data = readCBS(new_file$datapath)
      v$mat = v$data$output
      v$iter <- get_max_iter(v$mat)
      v$iter_max <- v$iter_min <- 1
      v$manual_update=F
      v$iter_prev <- v$iter
      values = c(firstPercentile(v$mat,x(v$iter)),ninetyNinthPercentile(v$mat,x(v$iter)))
      print(values)
      updateSliderInput(session, "xlimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
      print(values)
      values = c(firstPercentile(v$mat,y(v$iter)),ninetyNinthPercentile(v$mat,y(v$iter)))
      updateSliderInput(session, "ylimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
    }
  })
  
  observe({
    folder = req(input$folder)
    if(v$folder != folder){
      v$folder = folder
      v$file_num = 1
      v$files = sapply(list.files(path=folder, pattern=".RData$"), function(x) paste(dir,x[[1]],sep=''))
      v$mat = readCBS(v$files[[1]])
      v$iter <- v$iter_prev <- 1
      v$iter_max <- v$iter_min <- 1
      values = c(firstPercentile(v$mat,x(v$iter)),ninetyNinthPercentile(v$mat,x(v$iter)))
      updateSliderInput(session, "xlimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
      values = c(firstPercentile(v$mat,y(v$iter)),ninetyNinthPercentile(v$mat,y(v$iter)))
      updateSliderInput(session, "ylimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
    }
  })
  
  observe({
    req(v$file, cancelOutput = TRUE)
    if(input$draw_rect & !is.null(input$butterflybrush)){
      isolate({
        xmin = input$butterflybrush$xmin
        xmax = input$butterflybrush$xmax
        ymin = input$butterflybrush$ymin
        ymax = input$butterflybrush$ymax
        
        v$purity1 = ((xmax-xmin)*(ymax-ymin)/(ymin*xmax))^0.5
        
        print(input$butterflybrush)
        # Q: which point is 11? A: either (xmin,ymax) or (xmax,ymin)
        if(abs(xmin/ymax-1) < abs(xmax/ymin-1)){
          xmintemp = xmin
          xmaxtemp = xmax
          xmin = ymin
          xmax = ymax
          
          ymin = xmintemp
          ymax = xmaxtemp
        }
        # this ensures that we are looking at a rectangle in the upper diagonal
        # hence we have the mapping
        # (xmin,ymin) to (0,1)
        # (xmax,ymin) to (1,1)
        # (xmin,ymax) to (0,2)
        # (xmax,ymax) to (1,2)
        # the form of the transformation is
        # (xnew,ynew) = A%*% ((xold,yold) + B -(1,1)) +(1,1)
        # hence B is (1-xmax,1-ymin)
        # and A is the amount we want to strech each coordinate by so A is c(1/(xmax-xmin),0,0,1/(ymax-ymin))
        #*# v$mat[,x(v$iter+1)] = (v$mat[,x(v$iter)] -xmax) /(xmax-xmin) +1
        #*# v$mat[,y(v$iter+1)] = (v$mat[,y(v$iter)] -ymin) /(ymax-ymin) +1
        v$mat[,"tempx"] = (v$mat[,x(v$iter)] -xmax) /(xmax-xmin) +1
        v$mat[,"tempy"] = (v$mat[,y(v$iter)] -ymin) /(ymax-ymin) +1
        update_iter()
        values = c(firstPercentile(v$mat,x(v$iter)),ninetyNinthPercentile(v$mat,x(v$iter)))
        updateSliderInput(session, "xlimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
        values= c(firstPercentile(v$mat,x(v$iter)),ninetyNinthPercentile(v$mat,y(v$iter)))
        updateSliderInput(session, "ylimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
        
        updateCheckboxInput(session, "draw_rect", value = FALSE)
      })
    }
  })
  
  observe({
    updateSliderInput(session, "min_density_height", max = as.numeric(input$max_min_density_height),step=as.numeric(input$max_min_density_height)/100)
  })

  observeEvent(input$normalise_via_gam,{
    req(v$file, cancelOutput = TRUE)
    isolate({
      v$mat <- use_labels_with_gam(v$mat,v$iter,v$local_maxima,v$z)
      update_iter()
    })
  })
  
  observeEvent(input$cluster,{
    req(v$file, cancelOutput = TRUE)
    isolate({
      v$mat$CLUSTER <- as.factor(DC(v$mat[,x(v$iter)],v$mat[,y(v$iter)],v$mat[,"tcnNbrOfSNPs"],M=input$alpha))
      updateCheckboxInput(session, "plot_clusters", value = TRUE)
    })
  })
  
  observe({
    if(input$equal_scales){
      updateSliderInput(session, "ylimits", value = input$xlimits,min = input$xlimits[1], max = roundedMax(v$mat,v$iter))
    }
  })
  
  observe({
    if(input$update_scale_max){
      values = c(firstPercentile(v$mat,x(v$iter)),ninetyNinthPercentile(v$mat,x(v$iter)))
      updateSliderInput(session, "xlimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
      if(input$equal_scales){
        values = c(firstPercentile(v$mat,x(v$iter)),ninetyNinthPercentile(v$mat,x(v$iter)))
        updateSliderInput(session, "ylimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
      }else{
        values = c(firstPercentile(v$mat,y(v$iter)),ninetyNinthPercentile(v$mat,y(v$iter)))
        updateSliderInput(session, "ylimits", value = values,min = values[1], max = roundedMax(v$mat,v$iter))
        
      }
    }
  })
  
  # observeEvent(input$next_file,{
  #   v$file_num = v$file_num + 1
  #   # have some sort of error message to show that this is the last file in the folder
  #   v$mat = readCBS(v$files[[v$file_num]])
  #   v$iter <- 1
  #   updateSliderInput(session, "xmax", value = ninetyNinthPercentile(v$mat),
  #                     min = 0, max = roundedMax(v$mat))
  # })
  # observeEvent(input$previous_file,{
  #   v$file_num = v$file_num - 1
  #   # have some sort of error message to show that this is the last file in the folder
  #   v$mat = readCBS(v$files[[v$file_num]])
  #   #work out how to find the largest iter this file is up to. 
  #   v$iter <- get_max_iter(v$mat) #should make a diaglogue box appear here. 
  #   updateSliderInput(session, "xmax", value = ninetyNinthPercentile(v$mat),
  #                     min = 0, max = roundedMax(v$mat))
  # })
  # 
  # observeEvent(input$save_file,{
  #   v$file_num = v$file_num - 1
  #   # have some sort of error message to show that this is the last file in the folder
  #   v$mat = readCBS(v$files[[v$file_num]])
  #   #work out how to find the largest iter this file is up to. 
  #   v$iter <- 1
  #   updateSliderInput(session, "xmax", value = ninetyNinthPercentile(v$mat),
  #                     min = 0, max = roundedMax(v$mat))
  # })
  # 
  observeEvent(input$undo,{
    if(v$iter > v$iter_min){
      v$iter_prev = v$iter
      v$iter <- v$iter - 1
      v$mat$xtemp <- v$mat[,x(v$iter)] 
      v$mat$ytemp <- v$mat[,y(v$iter)] 
    }
  })
  
  observeEvent(input$redo,{
    if(v$iter < v$iter_max){
      update_iter()
      
      # v$iter_prev = v$iter
      # v$iter <- v$iter + 1
      # v$iter_max = max(v$iter,v$iter_max)
    }
  })
  # 
  #  
  # observeEvent(input$AB_normalise,{
  #   v$iter <- v$iter + 1
  # })
  # 
  observeEvent(input$shear_on_density,{
    req(v$file, cancelOutput = TRUE)
    opt <- find_opt_starting_values(v$mat,v$iter,ADJUST=input$shear_adjust,changes=F)
    updateSliderInput(session, "angle1", value = opt$best_theta1+input$angle1)
    updateSliderInput(session, "angle2", value = opt$best_theta2+input$angle2)
    T1 <- matrix(c(1,tan(opt$best_theta1),tan(opt$best_theta2),1),nrow=2,ncol=2)
    #*# v$mat[,c(x(v$iter+1),y(v$iter+1))] <- as.data.frame(as.matrix(v$mat[,c(x(v$iter),y(v$iter))]) %*% solve(t(T1)))
    v$mat[,c("tempx","tempy")] <- as.data.frame(as.matrix(v$mat[,c(x(v$iter),y(v$iter))]) %*% solve(t(T1)))
    
    update_iter()
    # v$iter <- v$iter+1
    # v$iter_max = max(v$iter,v$iter_max)
  })

  observeEvent(input$Do_IMBA,{
    #Find index of points to move and theta
    #den = density(atan(mat$a1/mat$a2), bw = 0.05)
    den = density(atan(v$mat[,x(v$iter)]/v$mat[,y(v$iter)]), bw = 0.05)
    #plot(den, main = 'Angle theta')
    der = (den$y[-1] - den$y[-length(den$y)])/(den$x[2]-den$x[1])
    #plot(den, main = 'Angle theta')
    pos = which(der >= 0)
    neg = which(der < 0)
    cut = neg[neg<pos[length(pos)]]
    cut = cut[length(cut)]
    #abline(v = den$x[cut], col = 'red')
    theta = den$x[pos[length(pos)]]
    #abline(v = theta, lty = 'dotted')
    #moveix = which(atan(mat$a1/mat$a2) > den$x[cut])
    moveix = which(atan(v$mat[,x(v$iter)]/v$mat[,y(v$iter)]) > den$x[cut])
    
    g = c(1, tan(theta))
    ng = sqrt(sum(g * g))
    g = g/ng
    #p.prime = mat$a2 * g[1] + mat$a1 * g[2]
    p.prime = v$mat[,y(v$iter)] * g[1] + v$mat[,x(v$iter)] * g[2]
    x.prime = p.prime * cos(theta)
    y.prime = p.prime * sin(theta)
    #lims = c(0, max(mat$a1, mat$a2))
    #plot(mat$a2, mat$a1, xlim = lims, ylim = lims,
    #main = 'Red points move to blue', xlab = 'Major array CN',
    #ylab = 'Minor array CN', pch = 16, col = col.transp('grey'))
    #abline(0,1)
    #lines(c(0,100*g[1]), c(0, 100*g[2]), col = 'red')
    #points(mat$a2[moveix], mat$a1[moveix], pch = 16,
    #col = col.transp('red'))
    d = (x.prime - y.prime)/2
    #points((mat$a2-d)[moveix], (mat$a1+d)[moveix],
    #col = col.transp('blue'), pch = 16)
    #if(!is.null(printout)) dev.print(png, file = printout, width= width, height = height)
    

    # fill in the new data structure for undo
    req(v$file, cancelOutput = TRUE)
    #Rearrange a1 and a2 according to size
    v$mat[,c(x(v$iter+1),y(v$iter+1))] <- v$mat[,c(x(v$iter),y(v$iter))]
    #mat$a1[moveix] = pmin((mat$a2-d)[moveix], (mat$a1+d)[moveix])
    v$mat[moveix,x(v$iter+1)] = pmin((v$mat[,y(v$iter)]-d)[moveix], (v$mat[,x(v$iter)]+d)[moveix])
    #mat$a2[moveix] = pmax((mat$a2-d)[moveix], (mat$a1+d)[moveix])
    v$mat[moveix,y(v$iter+1)] = pmax((v$mat[,y(v$iter)]-d)[moveix], (v$mat[,x(v$iter)]+d)[moveix])
    #***NO TRACKING HERE FOR NOW ****
      #mat$fixedIMBA = TRUE  
      #mat$allelic_balance = FALSE
      #mat$allelic_balance[moveix] = TRUE
    #par(mfrow = c(1,1))
    #obj
    update_iter()
    # 
    # v$iter <- v$iter+1
    # v$iter_max = max(v$iter,v$iter_max)
    
  })
  
  # observeEvent(input$slide_on_density,{
  #   req(v$file, cancelOutput = TRUE)
  #   
  #   dforig <- df <- v$mat[,c(x(v$iter),y(v$iter))]
  #   colnames(df) <- c("xnoisy","ynoisy")
  #   colnames(dforig) <- c("xnoisy","ynoisy")
  #   
  #   for( i in 1:10){
  #     lm = lm(data=df,ynoisy~xnoisy)
  #     print(coef(lm))
  #     
  #     df <- df[which(df$ynoisy < df$xnoisy*coef(lm)[2] + coef(lm)[1] +input$slider_bonus),]
  #     if(nrow(df) < num_points*0.05) break
  #   }
  #   
  #   # df is only points on the line
  #   twoddensity <- kde2d(x=df$xnoisy,y=df$ynoisy)
  #   
  #   # find the local maxima
  #   z <- as.matrix(twoddensity$z)
  #   local_maxima <- as.data.frame(which(z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[1:(nrow(z)-2),2:(ncol(z)-1)] &
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[3:(nrow(z)),2:(ncol(z)-1)] &
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[2:(nrow(z)-1),1:(ncol(z)-2)] &
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[2:(nrow(z)-1),3:(ncol(z))]&
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[3:(nrow(z)),3:(ncol(z))]&
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[1:(nrow(z)-2),1:(ncol(z)-2)]&
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[1:(nrow(z)-2),3:(ncol(z))]&
  #                                         z[2:(nrow(z)-1),2:(ncol(z)-1)] > z[3:(nrow(z)),1:(ncol(z)-2)],
  #                                       arr.ind = T))
  #   # array indicies are off by 1,1 because which reports the bottom corner of the square that we are comparing.
  #   local_maxima[,1] <- local_maxima[,1] +1
  #   local_maxima[,2] <- local_maxima[,2] +1
  #   local_maxima$z <-  z[as.matrix(local_maxima)]
  #   
  #   oneone <- local_maxima[which(local_maxima$z == max(local_maxima$z)),]
  #   oneonex <- twoddensity$x[oneone$row] 
  #   oneoney <- twoddensity$y[oneone$col] 
  #   
  #   df$xcorr <- df$xnoisy /oneonex
  #   df$ycorr <- df$ynoisy /oneoney
  #   
  #   lm = lm(data=df,ycorr~xcorr)
  # 
  #   dforig$xcorr <- dforig$xnoisy/oneonex
  #   dforig$ycorr <- dforig$ynoisy/oneoney
  #   dforig$xcorr2 <- (dforig$xcorr-1)*coef(lm)[[2]]+1
  #   dforig$ycorr2 <- dforig$ycorr
  #   
  #   ##
  #   v$mat[,c(x(v$iter+1),y(v$iter+1))] <- dforig[,c("xcorr2","ycorr2")]
  #   v$iter <- v$iter+1
  #   v$iter_max = max(v$iter,v$iter_max)
  # })
  
  observeEvent(input$shear_on_changes,{
    req(v$file, cancelOutput = TRUE)
    opt <- find_opt_starting_values(v$mat,v$iter,input$shear_adjust)
    updateSliderInput(session, "angle1", value = opt$best_theta1+input$angle1)
    updateSliderInput(session, "angle2", value = opt$best_theta2+input$angle2)
    T1 <- matrix(c(1,tan(opt$best_theta1),tan(opt$best_theta2),1),nrow=2,ncol=2)

    #*# v$mat[,c(x(v$iter+1),y(v$iter+1))] <- as.data.frame(as.matrix(v$mat[,c(x(v$iter),y(v$iter))]) %*% solve(t(T1)))
    v$mat[,c("tempx","tempy")] <- as.data.frame(as.matrix(v$mat[,c(x(v$iter),y(v$iter))]) %*% solve(t(T1)))
    update_iter()
    # v$iter <- v$iter+1
    # v$iter_max = max(v$iter,v$iter_max)
  })
  
  observeEvent(input$resize_by_change_ratio,{
    
    #need to know what one one is.
    
    horizontal <- abs(c(NA,v$mat[1:(nrow(v$mat)-1),x(v$iter)] - v$mat[2:nrow(v$mat),x(v$iter)]))
    vertical <- abs(c(NA,v$mat[1:(nrow(v$mat)-1),y(v$iter)] - v$mat[2:nrow(v$mat),y(v$iter)]))
    choice <- vertical > horizontal
    strech_ratio <- median(horizontal[!choice],na.rm=T)/median(vertical[choice],na.rm=T)
    
    #*#
    v$mat[,"tempx"] <- v$mat[,x(v$iter)]*strech_ratio
    v$mat[,"tempy"] <- v$mat[,y(v$iter)]
    update_iter()
    
    # v$iter <- v$iter+1
    # v$iter_max = max(v$iter,v$iter_max)
  })
  
  observeEvent(input$resize_by_ab_ratio,{
    ab_points <- v$mat[which(v$mat$abCall),c(x(v$iter),y(v$iter))]
    
    # lm  = lm(ab_points[,x(v$iter)] ~ ab_points[,y(v$iter)])
    # coeff(lm)
    ab_ratio <- mean(abs(ab_points[,2]/ab_points[,1]),na.rm=T)
    
    #*# v$mat[,x(v$iter+1)] <- v$mat[,x(v$iter)]*ab_ratio
    #*# v$mat[,y(v$iter+1)] <- v$mat[,y(v$iter)]
    v$mat[,"tempx"] <- v$mat[,x(v$iter)]*ab_ratio
    v$mat[,"tempy"] <- v$mat[,y(v$iter)]
    update_iter()
    
    # v$iter <- v$iter+1
    # v$iter_max = max(v$iter,v$iter_max)
    
  })
  
  observeEvent(input$phase,{
    mat[,c(x(v$iter+1),y(v$iter+1))] <- as.data.frame(as.matrix(mat[,c(x(v$iter),y(v$iter))]) %*% solve(t(T1)))
    valid = mat[1:(nrow(mat)-1),"chr"] == mat[2:nrow(mat),"chr"]
    horizontal <- c(NA,v$mat[1:(nrow(v$mat)-1),x(v$iter)] - v$mat[2:nrow(v$mat),x(v$iter)])
    vertical <- c(NA,v$mat[1:(nrow(v$mat)-1),y(v$iter)] - new[2:nrow(v$mat),y(v$iter)])
    
    # THIS IS NOT FINISHED!!!!
    # v$iter_max = max(v$iter,v$iter_max)
    
  })

## Functions for gridITH 
  dist2 = function(x, y){
    if (is.null(dim(x))) x = matrix(x, ncol = 2)
    if (is.null(dim(y))) y = matrix(y, ncol = 2)
    sqrt((x[,1] - y[,1])**2 + (x[,2] - y[,2])**2)
  }
  mahalanobis.w = function(X, nu, w){
    require(MASS)
    c = cov.trob(X, center = c(0,0), nu = nu)$cov
    D = numeric(nrow(X))
    for (j in 1:nrow(X)){
      D[j] = mahalanobis(X[j,], center = c(0,0), cov = c/w[j])
    } 
    D
  }
  gridITH = function(mat){
    v$mat$types = 'A'
    v$mat$types[v$mat$D>input$Mahalanobis_cut] = 'B'
    v$mat$types[v$mat[,y(v$iter)]>input$Mahalanobis_cut_right] = 'D'
    v$mat$types[v$mat[,x(v$iter)]<input$Mahalanobis_cut_bottom]= 'C'
    v$mat$types[v$mat[,x(v$iter)]>input$Mahalanobis_cut_top] = 'D'
    #list(ith = ith, mat = mat)
  }
  
##end functions for gridITH 

Calculate_Mahalanobis<- function(){
    #Overlay lattice points into X
    maxCN = 24
    w = 1-exp(-v$mat$length/500000)
    if (length(w)==1) w = rep(w, nrow(v$mat[,x(v$iter)])) 
    E2 = 0:maxCN
    E1 = matrix(E2, ncol = maxCN + 1, nrow = maxCN + 1)
    E2 = t(E1)
    E1 = c(E1)
    E2 = c(E2)
    
    ix = which(v$mat[,y(v$iter)]<= input$Mahalanobis_cut_top & v$mat[,y(v$iter)]>=input$Mahalanobis_cut_bottom & v$mat[,x(v$iter)]<=input$Mahalanobis_cut_right)
    x = cbind(v$mat[,x(v$iter)], v$mat[,y(v$iter)])
    alldist = apply(x, 1, dist2, y = cbind(E1, E2))
    closest = apply(alldist, 2, which.min)
    e = cbind(E1, E2)[closest,]
    di = dist2(x, e)
    xcoord = x[,2]-e[,2]
    ycoord = x[,1]-e[,1]
    X = cbind(ycoord, xcoord)
    
    #Multivariate t robust covariance qq=plots
    nu = 2
    D = mahalanobis.w(X, nu, w)  
    v$mat$X <- X
    v$mat$D <- NA
    v$mat$D[ix] <- D[ix]
    #mahanalobis_cut <- cut
    gridITH()
} 

  observeEvent(input$Mahalanobis,{
    Calculate_Mahalanobis()
    })
  # observeEvent(input$Mahalanobis_cut_top,{
   #   Calculate_Mahalanobis()
   # })
  # observeEvent(input$Mahalanobis_cut_bottom,{
  #   Calculate_Mahalanobis()
  # })
  # observeEvent(input$Mahalanobis_cut_right,{
  #   Calculate_Mahalanobis()
  # })
  # observeEvent(input$Mahalanobis_cut,{
  #   Calculate_Mahalanobis(v)
  # })

  observe({
    input$length
    isolate({
      invalidated = which(req(v$mat, cancelOutput = TRUE)$length >= quantile(v$mat$length,input$length[[1]]/100)   & v$mat$length <= quantile(v$mat$length,input$length[[2]]/100) )
      v$mat[-invalidated,"invalid_length"] <- TRUE
      v$mat[invalidated,"invalid_length"] <- FALSE
    })
  })
  
  observe({
    input$number
    isolate({
      invalidated = which(req(v$mat, cancelOutput = TRUE)$tcnNbrOfSNPs >= quantile(v$mat$tcnNbrOfSNPs,input$number[[1]]/100)   & v$mat$tcnNbrOfSNPs <= quantile(v$mat$tcnNbrOfSNPs,input$number[[2]]/100) )
      v$mat[-invalidated,"invalid_number"] <- TRUE
      v$mat[invalidated,"invalid_number"] <- FALSE
    })
  })
  
  observe({
    input$strechx
    isolate({
      v$mat[,x(v$iter+1)] <- req(v$mat[,x(v$iter)],cancelOutput = T)*input$strechx
      v$mat[,y(v$iter+1)] <- v$mat[,y(v$iter)]
      v$iter <- v$iter+1
    })
  })
  
  observe({
    input$strechy
    isolate({
      v$mat[,y(v$iter+1)] <- req(v$mat[,y(v$iter)],cancelOutput = T)*input$strechy
      v$mat[,x(v$iter+1)] <- v$mat[,x(v$iter)]
      v$iter <- v$iter+1
    })
  })

  observe({
    output$purity_out <- renderText(paste("purity_1 =",req(v$purity1,cancelOutput = T),sep=" " ))

  })
  
  # #   v$mat
  # # })
  # 
  # # observe({
  # #   files =  list.files(path=input$input_folder, pattern=".RData$")
  # #   updateSelectInput(session, "file",
  # #                         choices =  sapply(files, function(x) paste(dir,x[[1]],sep='')),
  # #                         selected = NULL ) #head(sapply(list.files(path=input$input_folder, pattern=".RData$"), function(x) paste(dir,x[[1]],sep='')), 1))
  # # })
  # 
  # observe({
  #   req(userFile, cancelOutput = TRUE)
  #   req(v$file,cancelOutput = TRUE)
  #   req(file.exists(v$file),cancelOutput = TRUE)
  #   eps = 0.01
  #   max_copy = round(max(v$mat[,pairs()]),2)
  #   updateSliderInput(session, "xmax", value = max_copy+eps,
  #                       min = 0, max = max_copy+eps)
  # })
  # 
  # # observe({
  # #
  # #
  # #
  
  

  observe({
    alpha  = input$adjsame #(10-1)/33 # proportion of adjacent segs that are the same in truth
    beta = input$oneadjsame# (64-1)/66 # proportion of adjacent segs where only one of the two estimates change
    gamma = input$allelicbal #(7-1)/33 # proportion in allelic balance
    #update_iter()
    isolate({
      req(v$mat,cancelOutput = T)
      x<-v$mat[,x(v$iter)]
      y<-v$mat[,y(v$iter)]
      data = data.frame('x'=x,'y'=y)
      data$weight <- 1
      data$xnext = c(x[2:length(x)],x[1])
      data$ynext = c(y[2:length(y)],y[1])
      data$chr <- v$mat$chr #c(rep(1,15),rep(2,18))
      # now we want to find new locations for each thing, so lets move everything around
      # things that are close we want to move together
      # things that are far we want to stay the same distance
      # things that are adjacen we want be vertical or horizontal

      data$distx <- abs(data$x - data$xnext)
      data$disty <- abs(data$y - data$ynext)
      data$dist <-(data$distx^2 + data$disty^2)^0.5
      data$horizontal <- data$distx < data$disty # if true horizontal dist is an error
      data$smallerdist <- data$distx
      data$smallerdist[which(!data$horizontal)] <- data$disty[which(!data$horizontal)]
      data$imbalance <- abs(data$x - data$y)

      # delta = 0.1 # proportion same estimate
      # epsilon = 0.1 # proportion in allelic balance

      (threshalpha <- quantile(data$dist,alpha)) # both adj same
      (threshbeta <- quantile(data$smallerdist,beta)) # one adj same
      (threshgamma <- quantile(data$imbalance,gamma)) # both same
      # (threshdelta <- quantile(data$?,delta)) # both globally same
      # (threshepsilon <- quantile(data$?,epsilon)) # one globally same
      # worry about the deltas and epsilons later on

      # it might actually make sense to get rid of alpha - since beta does this for us
      # for delta and epsilon the best thing to do is let the overall clustering algorithm resolve that

      # we need some way to find cliques:
      A <- matrix(0,nrow=2*nrow(data),ncol=2*nrow(data))
      for(i in 1:(nrow(data)-1)){
        if(data[i,"chr"] == data[i+1,"chr"]){
          if(data[i,"dist"] < threshalpha){
            A[2*i-1,2*(i+1)-1] = 1
            A[2*i,2*(i+1)] = 1
          }
          if(data[i,"distx"] < threshbeta){
            A[2*i-1,2*(i+1)-1] = 1
          }
          if(data[i,"disty"] < threshbeta){
            A[2*i,2*(i+1)] = 1
          }
        }
        if(data[i,"imbalance"] < threshgamma){
          A[2*i-1,2*i] = 1
        }
      }
      A <- A+t(A)
      sum(A)/2
      # install.packages("igraph")
      # library(igraph)
      similaritygraph <- graph_from_adjacency_matrix(A, mode = "undirected")
      (components <- components(similaritygraph)$membership)
      data[,c("componentx","componenty")] <- matrix(components,nrow=nrow(data),ncol=2,byrow=T)
      # x=1
      means <- sapply(unique(components),function(x)
        (sum(data[which(x == data$componentx),"x"]*data[which(x == data$componentx),"weight"])+
           sum(data[which(x == data$componenty),"y"]*data[which(x == data$componenty),"weight"]))/
          (sum(data[which(x == data$componentx),"weight"])+sum(data[which(x == data$componenty),"weight"])))
      data$x2 <- means[data$componentx]
      data$y2 <- means[data$componenty]

      v$mat[,"tempx"] <- data$x2
      v$mat[,"tempy"] <- data$y2
    })
  })
  
  observe({
    req(v$mat,cancelOutput = TRUE)
    v$iter
    if(input$KDE){
      # shouldn't need this bit here, move it later
      density <- get_density(v$mat,v$iter,input$bandwidth,input$num_grid_points,input$min_density_height)
      v$local_maxima <- density$local_maxima
      v$z <- density$z
    }
  })
  
  observeEvent(input$accept,{
    input$accept
    req(v$mat,cancelOutput = TRUE)
    v$mat[,x(v$iter)] <- v$mat$tempx
    v$mat[,y(v$iter)] <- v$mat$tempy
  })
  
  
  
  output$butterflyplot <- renderPlot({
    mat <- req(v$mat,cancelOutput = TRUE)
    var1 = x(v$iter)
    var2 = y(v$iter)
    var1b = "tempx"
    var2b = "tempy"
    
    
    shadeselected = min(0.2*input$shade_points,1)
    mat$invalid <- mat$invalid_length | mat$invalid_number
    imp <- which(!mat$invalid) # important
    if(length(imp) > 1){
      rimp <- imp[2:length(imp)] # left important
      limp <- imp[1:(length(imp)-1)] # right important
    } else{
      rimp <- limp <- c()
    }
    
    mat[limp,right(var1)] <- mat[rimp,var1]
    mat[limp,right(var2)] <- mat[rimp,var2]
    
    mat[limp,right(var1b)] <- mat[rimp,var1b]
    mat[limp,right(var2b)] <- mat[rimp,var2b]
    
    mat[rimp,right_minus_left(var1)] <- mat[rimp,var1] - mat[limp,var1]
    mat[rimp,right_minus_left(var2)] <- mat[rimp,var2] - mat[limp,var2]
    
    mat[rimp,"invalid_change"] = mat[rimp,"chr"] != mat[limp,"chr"]
    mat[which(mat$invalid_change),c(right_minus_left(var1),right_minus_left(var2),right(var1),right(var2))] <- NA
    
    v$rimp <- rimp
    v$lim <- limp
    
    mat$size_points <- mat$length^(input$size_points)
    mat$size_lines <- mat$length^(input$size_lines)
    mat$size_density <- mat$length^(input$size_density)
    
    if(input$error_bars){
      mat$xmin <- mat[,var1]+mat[,"c1_2.5%"]-mat[,"a1"]
      mat$xmax <- mat[,var1]+mat[,"c1_97.5%"]-mat[,"a1"]
      mat$ymin <- mat[,var2]+mat[,"c2_2.5%"]-mat[,"a2"]
      mat$ymax <- mat[,var2]+mat[,"c2_97.5%"]-mat[,"a2"]
    }
    
    g <- ggplot(data = mat[which(!mat$invalid),])
    
    if(!is.null(v$mat$D)){ # draw the mahanalobis selection box here
      g <- g + geom_rect(xmin = -Inf,    xmax = Inf, ymin = -Inf,    ymax = Inf, fill = "grey")
      #g <- g + geom_rect(xmin = input$Mahalanobis_cut_bottom,    xmax = input$Mahalanobis_cut_right, ymin = input$Mahalanobis_cut_bottom,    ymax = input$Mahalanobis_cut_top, fill = "white")
      g <- g + geom_rect(xmin = as.numeric(input$Mahalanobis_cut_bottom),    xmax = as.numeric(input$Mahalanobis_cut_right), ymin = as.numeric(input$Mahalanobis_cut_bottom),    ymax = as.numeric(input$Mahalanobis_cut_top), fill = "white")
    }

      #Chris change for the ultimate butterfly eye effects + move this to the bottom layer of drawing
    if(input$contours){
      if(input$contours_shade){
        g <- g + stat_density2d(aes_string(x=var1,y=var2,fill="..level..",col='"original_contours"',size="size_density"),geom="polygon",n=input$num_grid_points4contours,h=input$bandwidth4contours)+scale_fill_gradient(low="grey", high="black") # ,size=length^(1/input$size_density)
        g <- g + stat_density2d(aes_string(x=var1b,y=var2b,fill="..level..",col='"original_contours_2"',size="size_density"),geom="polygon",n=input$num_grid_points4contours,h=input$bandwidth4contours)+scale_fill_gradient(low="grey", high="black") # ,size=length^(1/input$size_density)
        
        if(input$mirror_on){
          g <- g + stat_density2d(aes_string(x=var2,y=var1,fill="..level..",col='"flipped_contours"',size="size_density"),geom="polygon")+scale_fill_gradient(low="grey", high="black") 
          g <- g + stat_density2d(aes_string(x=var2b,y=var1b,fill="..level..",col='"flipped_contours_2"',size="size_density"),geom="polygon")+scale_fill_gradient(low="grey", high="black") 
        }
      }
      else{
        g <- g + geom_density2d(aes_string(x=var1,var2,col='"original_contours"',size="size_density"))
        g <- g + geom_density2d(aes_string(x=var1b,var2b,col='"original_contours_2"',size="size_density"))
        
        if(input$mirror_on){
          g <- g + geom_density2d(aes_string(x=var2,y=var1,col='"flipped_contours"',size="size_density"))  
          g <- g + geom_density2d(aes_string(x=var2b,y=var1b,col='"flipped_contours_2"',size="size_density"))  
        }
      }
    }

    if(input$error_bars){
      g <- g+geom_errorbar(aes_string(x=var1,ymin = "ymin",ymax="ymax"),alpha=min(0.01*input$shade_points,1))+
        geom_errorbarh(aes_string(x=var1,y=var2,xmin = "xmin",xmax="xmax"),alpha=min(0.01*input$shade_points,1))
      if(input$mirror_on){
        g <- g+geom_errorbar(aes_string(x=var2,ymin = "xmin",ymax="xmax"),alpha=min(0.01*input$shade_points,1))+
          geom_errorbarh(aes_string(x=var2,y=var1,xmin = "ymin",xmax="ymax"),alpha=min(0.01*input$shade_points,1))
      }
    }
    
    if(input$lines){
      g <- g + geom_segment(aes_string(x=var1,xend=right(var1),y=var2,yend=right(var2),col='"original_changes"',size="size_lines"),alpha=min(0.01*input$shade_lines,1))
      g <- g + geom_segment(aes_string(x=var1b,xend=right(var1b),y=var2b,yend=right(var2b),col='"original_changes_2"',size="size_lines"),alpha=min(0.01*input$shade_lines,1))
      if(input$mirror_on){
        g <- g + geom_segment(aes_string(x=var2,xend=right(var2),y=var1,yend=right(var1),col='"flipped_changes"',size="size_lines"),alpha=min(0.01*input$shade_lines,1)) 
        g <- g + geom_segment(aes_string(x=var2b,xend=right(var2b),y=var1b,yend=right(var1b),col='"flipped_changes_2"',size="size_lines"),alpha=min(0.01*input$shade_lines,1)) 
      }
    }
    
    if(input$points){
      g <- g + geom_point(aes_string(x=var1,y=var2,col='"original"',size="size_points"),alpha=min(0.01*input$shade_points,1))
      g <- g + geom_point(aes_string(x=var1b,y=var2b,col='"original_2"',size="size_points"),alpha=min(0.01*input$shade_points,1))
      if(input$mirror_on){
        g <- g + geom_point(aes_string(x=var2,y=var1,col='"flipped"',size="size_points"),alpha=min(0.01*input$shade_points,1))
        g <- g + geom_point(aes_string(x=var2b,y=var1b,col='"flipped_2"',size="size_points"),alpha=min(0.01*input$shade_points,1))
      }
    }
    
    g <- g + xlim(input$xlimits[[1]],input$xlimits[[2]])+ylim(input$ylimits[[1]],input$ylimits[[2]])+theme_bw()+ # +coord_cartesian()
      theme(legend.position="none") #,panel.spacing = unit(, "lines")
    
    if(input$selected_points & nrow(brushedtrack()) >0){
      g <- g + geom_point(data = brushedtrack(),aes_string(x=var1,y=var2,col='"selected"',size="size_points"),alpha=min(0.01*input$shade_selected_points,1))
      if(input$mirror_on){
        g <- g+geom_point(data=brushedtrack(),aes_string(x=var2,y=var1,col='"selected_flipped"',size="size_points"),alpha=min(0.01*input$shade_selected_points,1)) # ,size=length^(1/input$size_points)
      }
    }
    
    if(input$selected_points & nrow(brushedbutterfly()) >0){
      g <- g + geom_point(data = brushedbutterfly(),aes_string(x=var1,y=var2,col='"selected"',size="size_points"),alpha=min(0.01*input$shade_selected_points,1))
      if(input$mirror_on){
        g <- g+geom_point(data=brushedbutterfly(),aes_string(x=var2,y=var1,col='"selected_flipped"',size="size_points"),alpha=min(0.01*input$shade_selected_points,1)) # ,size=length^(1/input$size_points)
      }
    }
    
    if(input$plot_clusters){
      # need to check if we have already clustered...
      if(!("CLUSTER" %in% colnames(mat))){
        updateCheckboxInput(session, "plot_clusters", value = FALSE)
      }
      updateCheckboxInput(session, "subclone_barchart", value = TRUE)
      g <- g + geom_point(aes_string(x=var1,y=var2,col="CLUSTER",size="length"))#+scale_fill_gradient2()
      if(input$mirror_on){
        g <- g + geom_point(aes_string(x=var2,y=var1,col="CLUSTER",size="length"))#+scale_fill_gradient2()
      }
       if(input$save_plot){
        #ggsave(filename="/Users/lmcintosh/myPlot.pdf", plot=g,width=2000,height=2000,limitsize=FALSE,scale=5)
        ggsave(filename="./plots/myPlot.pdf", plot=g,width=2000,height=2000,limitsize=FALSE,scale=5)
      }
   }
    
    #chris 8/5/17 move this to last from before last so that labels are always printed on top
    if(input$KDE & !is.null(v$local_maxima)){
      print(v$Local_maxima)
      g <- g +
        geom_point(data=v$local_maxima,aes(x=x,y=y,size=z))+
        geom_text(data=v$local_maxima,aes(x=x,y=y,label=name),size=7)
    }

        return(g)
  })
  
  
  
  output$trackplot <- renderPlot({
    CN_track <- get_CN_track(req(v$mat,cancelOutput = TRUE),input$xlimits,var1 = x(v$iter),var2 = y(v$iter),centromeres,num_rows = input$num_rows)
    if(input$selected_points & nrow(brushedbutterfly()) > 0){
      CN_track <- CN_track +  geom_rect(data = brushedbutterfly(),aes(xmin = Start, ymin = 0, xmax = End, ymax = Inf),size=2,alpha = 0.2)+
        geom_rect(data = brushedbutterfly(),aes(xmin = Start, ymin = 0, xmax = End, ymax = Inf),size=1,alpha = input$user_shade*0.01)
    }
    if(input$selected_points & nrow(brushedtrack()) > 0){
      CN_track <- CN_track +  geom_rect(data = brushedtrack(),aes(xmin = Start, ymin = 0, xmax = End, ymax = Inf),size=2,alpha = 0.2)+
        geom_rect(data = brushedtrack(),aes(xmin = Start, ymin = 0, xmax = End, ymax = Inf),size=1,alpha = input$user_shade*0.01)
    }
    return(CN_track)
  })
  output$trackplotui <- renderUI({
    if(input$track_plot){
      plotOutput("trackplot",hover = hoverOpts(id = "trackhover",clip = TRUE, nullOutside = TRUE,delayType = "debounce",delay=2000),
                 brush = brushOpts(id = "trackbrush", resetOnNew = FALSE, direction = "x",delayType = "debounce",delay=2000,opacity=input$user_shade*0.01),
                 height=300*input$num_rows,width=800)      
    }

  })
  
  output$trackplot_histogram <- renderPlot({
    summary_histogram(v$mat,input$xlimits,var1 = x(v$iter),var2 = y(v$iter))
  })
  
  output$trackplot_histogramui <- renderUI({
    if(input$bar_plots){
      plotOutput("trackplot_histogram",height=600)
    }
  })
  
  output$subclone_barchart <- renderPlot({
    # proportion of cells vs proportion of genome
    # and as a normal curve
    req(v$mat,cancelOutput = T)
    X <- v$mat[,x(v$iter)]
    Y <- v$mat[,y(v$iter)]
    
    Xc <- round(X)
    Yc <- round(Y)
    # Xc and Yc tell you the closest point and Xnc, Ync likewise the next closest
    Xnc <- Xc
    Ync <- Yc
    Xnc <- Xc + sign(X-Xc)
    Ync <- Yc + sign(Y-Yc)
    
    # great now calculate the cellular fraction, say that you are going up if you are outside the 00 to 11 square
    Xup <- apply(cbind(Xc,Xnc),1,max)
    Yup <- apply(cbind(Yc,Ync),1,max)
    Xdown <- apply(cbind(Xc,Xnc),1,min)
    Ydown <- apply(cbind(Yc,Ync),1,min)
    
    Xparent <- Xdown
    Yparent <- Ydown
    Xparent[which(X<0)] <- Xup[which(X<0)]
    Yparent[which(Y<0)] <- Yup[which(Y<0)]
    
    xBAF <- abs(X-Xparent)
    yBAF <- abs(Y-Yparent)
    
    BAF = xBAF
    BAF[which(xBAF > yBAF)] <- yBAF[which(xBAF > yBAF)]
    # above copied from the DC function
    
    v$mat$BAF <- BAF
    # v$mat$CLUSTER # gives the cluster
    # v$mat$tcnNbrOfSNPs
    # v$mat$length    
    v$mat$prop_genome_SNPS <- v$mat$tcnNbrOfSNPs / sum(v$mat$tcnNbrOfSNPs)
    v$mat$prop_genome_length <- v$mat$length / sum(v$mat$length)
    
    # fragmentwise plot:
    g1 <- ggplot() + geom_point(aes(x=v$mat$BAF, y=v$mat$prop_genome_SNPS,col=v$mat$CLUSTER,))+xlim(0,1)+theme_bw()
    #g1+geom_label_repel(aes(label = v$mat$prop_genome_SNPS))
                          
    g2 <- ggplot() + geom_point(aes(x=v$mat$BAF, y=v$mat$prop_genome_length,col=v$mat$CLUSTER))+xlim(0,1)+theme_bw()
    
    # global plot
    v$mat$BAF_mirrored <- sapply(v$mat$BAF,function(x) min(x,1-x))
    mat_BAF <- aggregate(x = list(BAF = v$mat$BAF_mirrored),by=list(CLUSTER = factor(v$mat$CLUSTER,levels = sort(unique(v$mat$CLUSTER)))),FUN="mean")
    print(mat_BAF)
    mat_prop <- aggregate(v$mat[,c("prop_genome_SNPS","prop_genome_length")],by=list(CLUSTER = factor(v$mat$CLUSTER,levels = sort(unique(v$mat$CLUSTER)))),FUN="sum")
    print(mat_prop)
    g3 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_SNPS,col=mat_BAF$CLUSTER))+xlim(0,0.5)+
      geom_text(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_SNPS,label = sapply(mat_prop$prop_genome_length,function(x){sprintf("%.2f",x)})),nudge_y = 0.03)+
      #guides(fill = guide_legend(title = NULL))+
      #scale_fill_discrete(name="ITH_length")+
      theme_bw()
    g4 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,col=mat_BAF$CLUSTER))+xlim(0,0.5)+
      geom_text(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,label = sapply(mat_prop$prop_genome_SNPS,function(x){sprintf("%.2f",x)})),nudge_y = 0.03)+
      #guides(fill = guide_legend(title = "ITH_length"))+
      #scale_fill_discrete(name="ITH_length")+
      theme_bw()
    #g4 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length))+xlim(0,0.5)
    #g4+scale_fill_discrete(name="Quality of the Cut",breaks=c("Fair", "Good", "Very Good", "Premium", "Ideal"),labels=c("Fair Stones", "Good Stones", "Very Good Stones", "Premium Stones", "Ideal Stones"))
    #g4 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,col= mat_BAF$CLUSTER, labs = sapply(mat_prop$prop_genome_length,function(x){sprintf("%.2f",x)})))+xlim(0,0.5)+theme_bw()
    #g4 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,col=factor( mat_BAF$CLUSTER, labels = sapply(mat_prop$prop_genome_length,function(x){sprintf("%.2f",x)}))))+xlim(0,0.5)+theme_bw()
    #g4 + scale_fill_discrete(name="ITH",breaks=mat_BAF$CLUSTER,labels=sprintf("%.2f",mat_prop$prop_genome_length))+theme_bw()
    #g4 + scale_fill_discrete(name="ITH",breaks=mat_BAF$CLUSTER,labels=sprintf("%.2f",mat_prop$prop_genome_length))
    
    #g4 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,col= factor( mat_BAF$CLUSTER, labels = sprintf("%.2f",mat_prop$prop_genome_length))))+xlim(0,0.5)+theme_bw()
    #g4 <- ggplot() + geom_point(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,color = factor( mat_BAF$CLUSTER, labels = apply(mat_prop$prop_genome_length,function(x){sprintf("%2d",x)})))+xlim(0,0.5)+theme_bw()
    #g4+geom_text(aes(x=mat_BAF$BAF, y=mat_prop$prop_genome_length,label =mat_prop$prop_genome_length))
    grid.arrange(g1,g2,g3,g4)
  })
  
  output$subclone_barchartui <- renderUI({
    if(input$subclone_barchart){
      plotOutput("subclone_barchart",height=800,width=800)
    }
  })
  
  
  
  output$butterflyplotui <- renderUI({
    if(input$butterfly_plot){
      plotOutput("butterflyplot",hover = hoverOpts(id = "butterflyhover",clip = TRUE, nullOutside = TRUE,delayType = "debounce",delay=2000),
                 brush = brushOpts(id = "butterflybrush", resetOnNew = FALSE, direction = "xy",delayType = "debounce",delay=2000,opacity=input$user_shade*0.01),height=800,width=800)      
    }

  })
  
  
  addPopover(session, "butterflyplotui", "Butterfly plot", content = butterfly_plot_description, trigger = 'hover',placement = "left")
  addPopover(session, "trackplot_histogramui", "CN track plot UI", content = CN_track_description, trigger = 'hover',placement = "left")
  addPopover(session, "trackplotui", "CN track plot", content = CN_track_description, trigger = 'hover',placement = "left")
  addPopover(session, "subclone_barchartui", "INFO", content = bar_plot_description, trigger = 'hover',placement = "left")
  addPopover(session, "densityplotui", "INFO", content = density_plots_description, trigger = 'hover',placement = "left")
  addPopover(session, "densitychangeplotui", "INFO", content = density_change_plots_description, trigger = 'hover',placement = "left")
  addPopover(session, "Mahanalobis_control_plotui", "INFO", content = Mahanalobis_control_plotui_description, trigger = 'hover',placement = "left")
  addPopover(session, "GridIT_dialPlotui", "INFO", content = GridIT_dialPlotui_description, trigger = 'hover',placement = "left")
  addPopover(session, "ITH_res_text", "INFO", content = subclone_barchart_description, trigger = 'hover',placement = "left")
  # suposedly it is much better to use tippify or popify with these ui elements.
  
  
  output$densityplots <- renderPlot({
    densityx <- get_density_plot(mat = req(v$mat,cancelOutput = TRUE),allele = "x",miniter=1,maxiter=v$iter,density_adjust = input$shear_adjust,limits = input$xlimits)
    densityy <- get_density_plot(mat = req(v$mat,cancelOutput = TRUE),allele = "y",miniter=1,maxiter=v$iter,density_adjust = input$shear_adjust,limits = input$ylimits)
    
    return(grid.arrange(densityx,densityy,ncol=2))
  })
  
  output$densityplotui <- renderUI({
    if(input$density_plots){
      plotOutput("densityplots",hover = hoverOpts(id = "densityhover",clip = TRUE, nullOutside = TRUE,delayType = "debounce",delay=2000),
               brush = brushOpts(id = "densitybrush", resetOnNew = FALSE, direction = "x",delayType = "debounce",delay=2000,opacity=input$user_shade*0.01),height=350)
    }
  })
  
  output$change_densityplots <- renderPlot({
    densityx <- get_change_density_plot(mat = req(v$mat,cancelOutput = TRUE),allele = "x",miniter=1,maxiter=v$iter,density_adjust = input$shear_adjust,limits = input$xlimits)
    densityy <- get_change_density_plot(mat = req(v$mat,cancelOutput = TRUE),allele = "y",miniter=1,maxiter=v$iter,density_adjust = input$shear_adjust,limits = input$ylimits)
    
    return(grid.arrange(densityx,densityy,ncol=2))
  })
  
  output$densitychangeplotui <- renderUI({
    if(input$density_change_plots){
      plotOutput("change_densityplots",hover = hoverOpts(id = "densityhover",clip = TRUE, nullOutside = TRUE,delayType = "debounce",delay=2000),
                 brush = brushOpts(id = "densitybrush", resetOnNew = FALSE, direction = "x",delayType = "debounce",delay=2000,opacity=input$user_shade*0.01),height=350)      
    }

  })

  
  output$Mahanalobis_control_plotui <- renderUI({
    if(input$Mahanalobis_control_plot) {
      if (is.null(v$mat$D) == FALSE){
        plotOutput("Mahanalobis_control_plot", height=600, width=600)
        }
    }
  })

  output$GridIT_dialPlotui <- renderUI({
      if (is.null(v$mat) == FALSE){
        plotOutput("GridIT_dialPlot")#, height=1000, width=1000)
      }
  })

    colors.by = function(x, pal = brewer.pal(9,'YlOrRd'), span = range(na.omit(x)),
                       alpha = .5){
    
    cols = colorRampPalette(pal)(100)
    colbreaks = seq(min(span), max(span), length = length(cols))
    delta = colbreaks[2] - colbreaks[1]
    colind = pmax((x-colbreaks[1]) %/% delta + 1, 1)
    colind[colind>length(cols)] = length(cols)
    mycols = col.transp(cols[colind], alpha = alpha)
    list(mycols=mycols, cols = cols)
  }
  col.transp <- function(col1, alpha=0.5){
    out = NULL
    for (j in 1:length(col1)){
      tmp=as.list(t(col2rgb(col1[j])))
      tmp$alpha=alpha*255
      tmp$maxColorValue=255
      out = c(out, do.call('rgb',tmp))
    }
    out
  }
  plot.cn = function(x, y, main = NULL, xlim = NULL, xlab = NULL, ylab = NULL,
                     ylim = NULL, maxCN = 24, mycols = NULL, grid = T){
    E1 = E2 = 0:maxCN
    
    #Plot frame
    if (is.null(xlim)) xlim = c(0, max(c(x, E2), na.rm = T))
    if (is.null(ylim)) ylim = c(0, max(y, na.rm = T)*1.2)
    type = ifelse(is.null(mycols), 'p', 'n')
    plot(x, y, pch = 16, col = col.transp('blue', .3), type = type,
         xlim = xlim, ylim = ylim, xlab = '', ylab = '')
    if (is.vector(mycols)){
      points(x, y, pch = 16, col = mycols, cex = 1)    
    }
    if (is.list(mycols)) {
      points(x, y, pch = 16, col = mycols$mycols, cex = 1)
      
      #Colour scale
      tmp = c(xlim[2]/2, xlim[2])
      tmp = seq(tmp[1], tmp[2], length=length(mycols$cols))
      points(tmp, rep(ylim[1], length(tmp)), col = mycols$cols, pch = 16)
      
    }
    #Gridlines
    if (grid){
      EEs = unique(E1)
      for (j in 1:length(EEs)){
        lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
        lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
      }
    }
    abline(0,1)
    title(main = main, outer = T, line = -1)
    title(xlab = xlab, ylab = ylab, line = 2)
  }
  
  output$Mahanalobis_control_plot <- renderPlot({

  #Colour scale for figures
  pal = brewer.pal(11, 'Spectral')
  pal = rev(pal[c(1:4, 9:11)])
  
  par(mfrow = c(2,2), mar = c(4, 4, 3, .1))
  #Multivariate t robust covariance qq=plots
  cols = (colors.by(v$mat$D, pal, span = c(0,10))$mycols)[order(v$mat$D)]
  qqplot(qexp(ppoints(length(v$mat$D)), rate = .5), v$mat$D, col = cols,
         main = 'Exponential QQ-plot',ylim = c(0,20),
         ylab = 'Mahalanobis squared distance',
         xlab = 'Exponential quantiles')
  abline(0, 1, col = 'gray')
  abline(h = input$Mahalanobis_cut, lty = 'dotted')
  #Overlay lattice points
  col1 = colors.by(v$mat$D, pal, span = c(0,10))
  plot(v$mat$X[,2], v$mat$X[,1], xlim = c(-.75, .75), col = col1$mycols, pch = 16,
       ylim = c(-.75, .75), main = '', cex.main = 1,
       xlab = 'Distance to closest major integer CN',
       ylab = 'Distance to closest minor integer CN')
  mtext('Colouring by squared Mahalanobis distance to lattice point', cex = .8)
  abline(v=0, h=0)
  rect(xleft = -.5, ybottom = -.5, xright = .5, 
       ytop = .5, lty = 'dotted')
  #Colour scale
  tmp = c(0.1, .75)
  tmp = seq(tmp[1], tmp[2], length=length(col1$cols))
  points(tmp, rep(-.75, length(tmp)), col = col1$cols, pch = 16)
  
  #Subclonal CN segments given cutoff with 1-pexp(cut, .5) sign.level
  cols = rep('blue', nrow(v$mat[,cbind(x(v$iter),y(v$iter))]))
  cols[v$mat$D > input$Mahalanobis_cut] = 'red'
  cols = col.transp(cols)
  plot(v$mat$X[,2], v$mat$X[,1], xlim = c(-.75, .75), col = cols, pch = 16,
       ylim = c(-.75, .75), main = paste('Subclonal CNs with cutoff =', 
                                         input$Mahalanobis_cut), cex.main = 1,
       xlab = 'Distance to closest major integer CN',
       ylab = 'Distance to closest minor integer CN')
  abline(v=0, h=0)
  rect(xleft = -.5, ybottom = -.5, xright = .5, 
       ytop = .5, lty = 'dotted')
  legend('topleft', bty = 'n', pch = 16, col = c('blue','red'), 
         legend = c('Segment with integer CN','Segment with subclonal CN'))
  
  #Choose how extreme CNs can be judged to be subclonal: CN subclonality defined here
  xlim = c(0, quantile(v$mat[,y(v$iter)], probs = c(.95))*c(1.2))
  ylim = c(0, quantile(v$mat[,x(v$iter)], probs = c(.95))*c(1.2))
  plot.cn(x = v$mat[,y(v$iter)], y = v$mat[,x(v$iter)], mycols = cols, xlim = xlim, ylim = ylim,
          xlab = 'Major CN, unscaled', ylab = 'Minor CN, unscaled')
  })
  
  
  output$ITH_res_text <- renderUI({
    if(input$Mahanalobis_control_plot) {
      if (!is.null(v$mat$types)){
        str1 = paste('A:', round(sum(v$mat$W[v$mat$type == 'A'])*100), '% of genome')
        str2 = paste('B:', round(sum(v$mat$W[v$mat$type == 'B'])*100), '% of genome')
        str3 = paste('C:', round(sum(v$mat$W[v$mat$type == 'C'])*100), '% of genome')
        str4 = paste('D:', round(sum(v$mat$W[v$mat$type == 'D'])*100), '% of genome') 
        str_ITH = paste('ITH:', round (sum(v$mat$W[v$mat$types == 'B'])/(sum(v$mat$W[v$mat$types == 'A']) + 
                                                            sum(v$mat$W[v$mat$types == 'B'])) * 100) , '% of genome') 
        
        #legend('bottomright', legend = c('A','B','C','D'),pch = 16, col = c('blue','magenta','green','red'), bty = 'n')
        #HTML(paste(str1, str2,str3,str4,str_ith, sep = '<br/>'))
        HTML(paste(str1, str2,str3,str4,str_ITH, sep = '<br/>'))
      }
    }
  })
  
  # # id	Input value name. For example, if the value is "plot_brush", then the coordinates will be available as input$plot_brush. Multiple imageOutput/plotOutput calls may share the same id value; brushing one image or plot will cause any other brushes with the same id to disappear.
  # # fill	Fill color of the brush.
  # # stroke	Outline color of the brush.
  # # opacity	Opacity of the brush
  # # delay	How long to delay (in milliseconds) when debouncing or throttling, before sending the brush data to the server.
  # # delayType	The type of algorithm for limiting the number of brush events. Use "throttle" to limit the number of brush events to one every delay milliseconds. Use "debounce" to suspend events while the cursor is moving, and wait until the cursor has been at rest for delay milliseconds before sending an event.
  # # clip	Should the brush area be clipped to the plotting area? If FALSE, then the user will be able to brush outside the plotting area, as long as it is still inside the image.
  # # direction	The direction for brushing. If "xy", the brush can be drawn and moved in both x and y directions. If "x", or "y", the brush wil work horizontally or vertically.
  # # resetOnNew	When a new image is sent to the browser (via renderImage), should the brush be reset? The default, FALSE, is useful if you want to update the plot while keeping the brush. Using TRUE is useful if you want to clear the brush whenever the plot is updated.
  
  output$GridIT_dialPlot <- renderPlot({
    
    # x    <- faithful[, 2]
    # bins <- seq(min(x), max(x), length.out = 5)
    # # draw the histogram with the specified number of bins
    # hist(x, breaks = bins, col = 'darkgray', border = 'white')

    obj = list(fileName = 'blank')
    obj$mat = v$mat
    # plot(obj$mat$a2, obj$mat$a1, pch = 16, 
    #      xlab = 'Major intensity ratio a2', 
    #      ylab = 'Minor intensity ratio a1',
    #      main = 'Minor and major CN intensity')
    # abline(0,1, col = 'red')
    
    obj$mystart = start.alphamax.f(mat = obj$mat, colbychrom = TRUE, xlim = c(0,1.5), ylim = c(0, 1.5),
                                    dx.eq.dy = TRUE, force.diag = TRUE)
    # plot.transformed(obj$mat, obj$mystart, xlim = c(0,5), ylim = c(0,5), color.by = obj$mat$W,
    #                  printout = NULL)
                    #printout = file.path(getwd(),"fixIMBA", paste(obj$fileName,'_2_startValues.png',   sep = '')))
  
  })      
  

  output$contents <- renderDataTable({
    if(input$data_table) v$mat
  })
  
  
  # # In server.R:
  # output$downloadData <- downloadHandler(
  #   filename = function(){
  #     if(is.null(input$outputfolder)){
  #       paste(dirname(req(v$file)),basename(v$file),"-processed",Sys.Date(),sep="")
  #     } else{
  #       paste(dirname(req(input$outputfolder)),basename(req(v$file)),"-processed",Sys.Date(),sep="")
  #     }
  #   },
  # 
  #   content = function(con) {
  #     write.csv(v$mat,con)
  #   }
  # )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(dirname(req(v$file)),basename(v$file),"-processed",Sys.Date(),".RData",sep="")
    },
    content = function(con) {
      v$data$output <- v$mat
      saveObject(v$data,con)
    }
  )
  
  # output$GridIT_plotui<-renderPlot({
  #   require(ggplot2)
  #   ggplot(data = mat[which(!mat$invalid),])
  # })
  
  # # output$downloadImage <- downloadHandler(
  # #   filename = function() {
  # #     paste(basename(userFile),"_transformed_plot","-",Sys.Date(),".png",sep="")
  # #   },
  # #
  # #   content = function(con) {
  # #     g1 <- ggplot(v$mat)+
  # #       geom_segment(aes(x=x(v$iter),xend=x(v$iter)_next,y=y(v$iter),yend=y(v$iter)_next,size=length,col="original changes"),alpha=0.05)+
  # #       geom_segment(aes(x=y(v$iter),xend=y(v$iter)_next,y=x(v$iter),yend=x(v$iter)_next,size=length,col="flipped changes"),alpha=0.05)+
  # #       geom_point(aes(x=x(v$iter),y=y(v$iter),size=length,col="original"),alpha=0.05)+
  # #       geom_point(aes(x=y(v$iter),y=x(v$iter),size=length,col="flipped"),alpha=0.05)+
  # #       xlim(0,input$xlimits)+ylim(0,input$xlimits)+coord_fixed()
  # #     ggsave(g1,con)
  # #
  # #   }
  # # )
  # #
  # # output$downloadInfo <- downloadHandler(
  # #   filename = function() {
  # #     paste(basename(userFile),"_final_parameters","-",Sys.Date(),"",sep="")
  # #   },
  # #
  # #   content = function(con) {
  # #     save(v$mat,input$angle1, input$angle2, input$ratioxy, input$strech,file=con)
  # #   }
  # # )
  # #
  
  brushedtrack <- reactive({
    BP <- brushedPoints(req(v$mat,cancelOutput = T), input$trackbrush)
    
    BP$size_points <- BP$length^(input$size_points)
    BP$size_lines <- BP$length^(input$size_lines)
    BP$size_density <- BP$length^(input$size_density)
    
    BP
  })
  
  brushedbutterfly <- reactive({
    BP <- unique(rbind(
      brushedPoints(req(v$mat,cancelOutput = T), input$butterflybrush, xvar = x(v$iter), yvar = y(v$iter)),
      brushedPoints(v$mat, input$butterflybrush , xvar = y(v$iter), yvar = x(v$iter))
    ))
    
    BP$size_points <- BP$length^(input$size_points)
    BP$size_lines <- BP$length^(input$size_lines)
    BP$size_density <- BP$length^(input$size_density)
    
    BP
  })
  
  # nearCNV <- reactive({
  #   if(!is.null(input$plot_hover1)){
  #     nearPoints(MAT(), input$plot_hover1, threshold = 10, maxpoints = 1,addDist = TRUE, xvar = "location",yvar="height")    
  #   }
  #   # else{ 
  #     # MAT()
  #   # }
  # })  
  # nearGRID <- reactive({
  #   if(!is.null(input$plot_hover)){
  #     nearPoints(MAT(), input$plot_hover, threshold = 10, maxpoints = 1,addDist = TRUE , xvar = "x(v$iter)", yvar = "y(v$iter)")
  #   } 
  #   # else{ 
  #     # MAT()
  #   # }
  # })
})


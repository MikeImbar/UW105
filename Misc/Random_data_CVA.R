Proc_two_d <- two.d.array(Proc$rotated)

Prb_data <- length(Proc_two_d[gp=="Prob",])
EHomo_data <- length(Proc_two_d[gp=="Early Homo",])
Aus_data <- length(Proc_two_d[gp=="Aus",])

EHomo_data/435

sample_size <- sum(Prb_data/435 + EHomo_data/435 + Aus_data/435)
(Prb_data/435)/sample_size*100

data1 <- matrix(rnorm(Prb_data, mean = 1, sd  = 1), ncol = 435, nrow = Prb_data/435)
data2 <- matrix(rnorm(EHomo_data, mean = 1, sd  = 1.1), ncol = 435, nrow = EHomo_data/435)
data3 <-  matrix(rnorm(Aus_data, mean = 1, sd  = 1.2), ncol = 435, nrow = Aus_data/435)

data1 <- matrix(rnorm(20000, mean = 1, sd  = 1), ncol = 500, nrow = 40)
data2 <- matrix(rnorm(20000, mean = 1, sd  = 1), ncol = 500, nrow = 40)
data3 <-  matrix(rnorm(20000, mean = 1, sd  = 1), ncol = 500, nrow =40)

data1 <- matrix(rnorm(20000, mean = 1, sd  = 1), ncol = 500, nrow = 40)
data2 <- matrix(rnorm(20000, mean = 1, sd  = 1), ncol = 500, nrow = 40)
data3 <-  matrix(rnorm(20000, mean = 1, sd  = 1), ncol = 500, nrow =40)

data_gp <- factor(c(rep("A", Prb_data/435),
             rep("B", EHomo_data/435),
             rep("C", Aus_data/435)))

data_gp <- factor(c(rep("A", 40), 
                    rep("B", 40),
                    rep("C", 40)))

data_rand <- data.frame(rbind(data1,
                              data2,
                              data3),
                        group = data_gp)


library(dplyr)

data_rand %>%
  group_by(group) %>%
  summarise(
    mean_all = mean(unlist(across(where(is.numeric))), na.rm = TRUE),
    sd_all   = sd(unlist(across(where(is.numeric))),  na.rm = TRUE),
    min_all  = min(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_all  = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    median_all = median(unlist(across(where(is.numeric))), na.rm = TRUE)
  )



rand_data_PCA <- prcomp(data_rand[, 1:ncol(data_rand)-1])
rand_data_PCA_df <- data.frame(rand_data_PCA$x, group = data_rand$group)

library(dplyr)

data_rand %>%
  group_by(group) %>%
  summarise(
    mean_all = mean(unlist(across(where(is.numeric))), na.rm = TRUE),
    sd_all   = sd(unlist(across(where(is.numeric))),  na.rm = TRUE),
    min_all  = min(unlist(across(where(is.numeric))), na.rm = TRUE),
    max_all  = max(unlist(across(where(is.numeric))), na.rm = TRUE),
    median_all = median(unlist(across(where(is.numeric))), na.rm = TRUE)
  )

# summary statistics
tapply(data_rand, data_rand$group, summary)


yrange<-sum(abs(range(rand_data_PCA$x[,2])));ymin<-min(rand_data_PCA$x[,2]);ymax<-max(rand_data_PCA$x[,2])
Xrange<-sum(abs(range(rand_data_PCA$x[,1])));Xmin<-min(rand_data_PCA$x[,1]);Xmax<-max(rand_data_PCA$x[,1])

plot(cbind(rand_data_PCA$x[,1],rand_data_PCA$x[,2]),type="n",asp=1,cex=1,ylim=c(ymin-0.25*yrange,ymax+0.25*yrange),xlab=paste("CV 1 (",round(CVA_res$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA_res$Var[2,2]),"%)",sep=""))
#plot(cbind(CVA$CVscores[,1],CVA$CVscores[,2]),type="n",asp=1,cex=1,xlim=c(Xmin-0.55*Xrange,Xmax+0.35*Xrange),xlab=paste("CV 1 (",round(CVA$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA$Var[2,2]),"%)",sep=""))
for(gruppe in 1:length(levels(data_gp))){
  sub<-data_gp==levels(data_gp)[gruppe]
  tr <- NULL
  try(tr<-tri.mesh(x=rand_data_PCA$x[sub,1],y=rand_data_PCA$x[sub,2],duplicate = "error"))
  if(!is.null(tr)){
    polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[gruppe], alpha.f = 0.5)),border=farbe[gruppe])
  } else if(sum(sub)==2){
    lines(rand_data_PCA$x[sub,1],rand_data_PCA$x[sub,2],col=farbe[gruppe],lwd = 2)
  }}
points(rand_data_PCA$x[,c(1,2)],col=farbe[gp_class],pch=19)
# text(rand_data_PCA$x[,c(1,2)],label=classifier$Name[classified],col=farbe[classifier$Class2[classified]],pos=c(1,2),cex=0.6,offset=0.5)
# points(cbind(projected[,1],projected[,2]),pch=19)
text(cbind(projected[,1],projected[,2]), label = rownames(projected), cex = 0.6, pos = c(1,2))                                   

par3d(FOV = 0, windowRect = c(30, 30, 1770, 1000))

clear3d("shapes");plot.range<-1.02*max(abs(cbind(rand_data_PCA$x[,1],rand_data_PCA$x[,2],rand_data_PCA$x[,3])))
plot3d(cbind(rand_data_PCA$x[,1],rand_data_PCA$x[,2],rand_data_PCA$x[,3]),type="p",lwd=12,box=FALSE,size=10,col=farbe[gp],aspect =T,axes=T,
       xlim = c(min(rand_data_PCA$x[,1])+min(rand_data_PCA$x[,1])*0.10,max(rand_data_PCA$x[,1])+max(rand_data_PCA$x[,1])*0.10),
       ylim = c(min(rand_data_PCA$x[,2])+min(rand_data_PCA$x[,2])*0.10,max(rand_data_PCA$x[,2])+max(rand_data_PCA$x[,2])*0.10),
       zlim = c(min(rand_data_PCA$x[,3])+min(rand_data_PCA$x[,3])*0.10,max(rand_data_PCA$x[,3])+max(rand_data_PCA$x[,3])*0.10),
       xlab = paste("PC1 (", round(Proc$Variance[1,2],1), "%)", sep = ""),
       ylab = paste("PC2 (", round(Proc$Variance[2,2],1), "%)", sep = ""),
       zlab = paste("PC3 (", round(Proc$Variance[4,2],1), "%)", sep = ""))
# texts3d(cbind(rand_data_PCA$x[,1],rand_data_PCA$x[,2],rand_data_PCA$x[,3]),texts = classifier$Name[classified],font=2,adj=1,cex=0.7)
legend3d("topleft", levels(data_gp), pch = 16, col = farbe, cex = 3)

for(a in 1:length(levels(data_gp))){
  sub<-data_gp==levels(data_gp)[a]
  PCsub<-cbind(rand_data_PCA$x[sub,1],rand_data_PCA$x[sub,2],rand_data_PCA$x[sub,3])
  if(length(PCsub[,1])>3){
    hull<-convhulln(PCsub)
    for(b in 1:length(hull[,1])){
      sub1<-PCsub[hull[b,],]
      triangles3d(sub1[,1],sub1[,2],sub1[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
    } 
  } else if(length(PCsub[,1])==3){
    triangles3d(PCsub[,1],PCsub[,2],PCsub[,3],col=farbe[a],alpha=0.3,lit=F,fog=F)
  } else if(length(PCsub[,1])==2){
    lines3d(PCsub,col=farbe[a])
  }
}

# points3d(cbind(projected[,1],projected[,2],projected[,3]), size =10)
# texts3d(cbind(projected[,1],projected[,2],projected[,3]),texts = rownames(projected),font=2,adj=1,cex=0.7)
legend3d("topleft", levels(gp), pch = 16, col = farbe, cex = 3)


#### CVA ####
summary(rand_data_PCA)

par(mfrow = c(3, 3))

# start PC = 80%, end PC = 95%
start_PC <- 100
end_PC <- 120

for(i in start_PC:end_PC){
  print(i)
  CVA_res<-CVA(dataarray=rand_data_PCA$x[,1:i],groups=data_gp,weighting = TRUE, cv = TRUE, prior=rep(1/nlevels(data_gp),nlevels(data_gp)))

  #CVA 1 v 2 (labels)
  yrange<-sum(abs(range(CVA_res$CVscores[,2])));ymin<-min(CVA_res$CVscores[,2]);ymax<-max(CVA_res$CVscores[,2])
  Xrange<-sum(abs(range(CVA_res$CVscores[,1])));Xmin<-min(CVA_res$CVscores[,1]);Xmax<-max(CVA_res$CVscores[,1])
  
  #pdf(paste("../Figures/All PCAs/Mar23 - Nat_comm_R1_Sang/CVA/",position,"_",number,"PCs_CVA_labels.pdf",sep=""),height=7,width=width)
  plot(cbind(CVA_res$CVscores[,1],CVA_res$CVscores[,2]),type="n",asp=1,cex=1,ylim=c(ymin-0.25*yrange,ymax+0.25*yrange),xlab=paste("CV 1 (",round(CVA_res$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA_res$Var[2,2]),"%)",sep=""))
  title(main = paste0(i ," PCs"))
  #plot(cbind(CVA$CVscores[,1],CVA$CVscores[,2]),type="n",asp=1,cex=1,xlim=c(Xmin-0.55*Xrange,Xmax+0.35*Xrange),xlab=paste("CV 1 (",round(CVA$Var[1,2]),"%)",sep=""),ylab=paste("CV 2 (",round(CVA$Var[2,2]),"%)",sep=""))
  for(gruppe in 1:length(levels(gp_class))){
    sub<-gp_class==levels(gp_class)[gruppe]
    tr <- NULL
    try(tr<-tri.mesh(x=CVA_res$CVscores[sub,1],y=CVA_res$CVscores[sub,2],duplicate = "error"))
    if(!is.null(tr)){
      polygon(convex.hull(tr)$x,convex.hull(tr)$y,col=(adjustcolor(farbe[gruppe], alpha.f = 0.5)),border=farbe[gruppe])
    } else if(sum(sub)==2){
      lines(CVA_res$CVscores[sub,1],CVA_res$CVscores[sub,2],col=farbe[gruppe],lwd = 2)
    }}
  points(CVA_res$CVscores[,c(1,2)],col=farbe[gp_class],pch=19)
  # text(CVA_res$CVscores[,c(1,2)],label=classifier$Name[classified],col=farbe[classifier$Class2[classified]],pos=c(1,2),cex=0.6,offset=0.5)
  # points(cbind(projected[,1],projected[,2]),pch=19)
  # text(cbind(projected[,1],projected[,2]), label = rownames(projected), cex = 0.6, pos = c(1,2))
  
  
}


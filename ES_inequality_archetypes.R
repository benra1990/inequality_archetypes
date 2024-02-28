##--------------
##
##Script name: Ecosystem services and (in)equality archetypes
##
##Purpose of script: To calculate archetypes using SOM clustering techniques for ES supply (eight ES) and distributional inequalities.  
##
##Author: 
##Affiliation:
##
##Copyright
##
##--------------
##
##Notes:Input data is made available under: https://github.com/benra1990/inequality_archetypes.git
##
##Scrips is provided to develop the practical processes described in the paper 
##
##Description of the input data required to replicate these analyses can be found at [] doi TBC
##
##


# Load needed libraries

library(rgdal)
library(raster)
library(rgeos)
library(readxl)
library(foreign)
library(kohonen)
library(som)
library(missSOM)
library(clusterSim)
library(reshape2)
library(openxlsx)
library(dplyr)
library(data.table)
library(scales)
library(ggplot2)
library(gridExtra)


####LOAD DATA####

db_inequities<-read.xlsx("Z:/Paper_inequalities/Datos/db_inequities.xlsx")#load database

glimpse(db_inequities)#glimpse list of variables
data.table(colnames(db_inequities))#list of variables with line numbers
db_inequities<-db_inequities[,c(1:18)]#select ES supply variables and inequality variables (gini coefficient)

####ES supply data####

# create placeholder variable and select variables to be used
db_es <- db_inequities[,c(11:18)] #select just ES variables variables for developing SOM
glimpse(db_es)#view of selected database

#check number of Na's in each column and in each row
colSums(is.na(db_es))#nr of Na within rows
rowSums(is.na(db_es))#nr of Na within columns
db_es["tot_supl_nntimber.y"][is.na(db_es["tot_supl_nntimber.y"])]<-0 #replace NAs by 0 in nntimber column

### 1 PREPARE DATA FOR SOMS: Z-TRANSFORM + MATRIX #####

db_es_scaled<-scale(db_es, center=TRUE, scale=TRUE)#z-scale data and assign new name

db_es_matrix <- as.matrix(db_es_scaled)#dataframe to matrix

### 2. SOMS #####
### 2.1 SET UP PARAMETERS #####
som.dims.es <- as.data.frame(matrix(ncol=2,nrow=16))
colnames(som.dims.es) <- c("som_cols","som_rows")
som.dims.es[,1] <- c(2,3,2,3,4,3,5,4,5,4,6,5,7,6,5,6)
som.dims.es[,2] <- c(1,1,2,2,2,3,2,3,3,4,3,4,3,4,5,5)

set.seed(2022)

col_sel_es <- c(1:length(db_es))

som_test_es <- list()

for(i in 1:nrow(som.dims.es)){
  temp.name.es <- paste0("SOM_cols",som.dims.es[i,1],"_rows",som.dims.es[i,2])
  temp.som.es <- supersom(db_es_matrix[,col_sel_es],
                          grid=somgrid(som.dims.es[i,1],som.dims.es[i,2],topo="hexagonal"),
                          rlen = 1000,
                          alpha = c(0.05,0.01),
                          keep.data=T,
                          maxNA.fraction=0.90)
  #between <- sum(rdist(as.matrix(temp.som$grid$pts)))
  #within <- sum(temp.som$distances)
  #som.homogeneity <- within
  #som.variance <- between/(between+within)
  som.dist.mean.es <- mean(temp.som.es$distances)
  som.dist.sd.es <- sd(temp.som.es$distances)
  temp.DB.es <- index.DB(x=as.data.frame(temp.som.es$data),
                         cl=temp.som.es$unit.classif,
                         centrotypes="centroids")
  som_test_es[[i]] <- list(temp.name.es,temp.som.es,som.dist.mean.es,som.dist.sd.es,temp.DB.es)
}

length(som_test_es)
length(som_test_es[[14]])

### 2.1 PERFORMANCE ####

som.perf.df.es <- as.data.frame(matrix(ncol=6,nrow=length(som_test_es)))
colnames(som.perf.df.es) <- c("name","cluster_combo","clus","DB","mean","sd")
som.perf.df.es[,2] <- c("2x1","3x1","2x2","3x2","4x2","3x3","5x2","4x3","5x3","4x4","6x3","5x4","7x3","6x4","5x5", "6X5")
som.perf.df.es[,3] <- c(2,3,4,6,8,9,10,12,15,16,18,20,21,24,25,30)

for(i in 1:length(som_test_es)){
  som.perf.df.es[i,1] <- som_test_es[[i]][[1]]   # name
  som.perf.df.es[i,4] <- som_test_es[[i]][[5]]$DB # DB
  som.perf.df.es[i,5] <- som_test_es[[i]][[3]]    # mean
  som.perf.df.es[i,6] <- som_test_es[[i]][[4]]    # sd
}

# plot

par(mar = c(5, 4, 4, 5))
plot(som.perf.df.es[c(1:16),3],som.perf.df.es[c(1:16),5],pch=2,col=4,
     ylim=range(pretty(c(min(som.perf.df.es[c(1:16),5]),max(som.perf.df.es[c(1:16),5])))),
     ylab="Mean distance to cluster centroid",xlab="Archetype size ES supply")
lines(som.perf.df.es[c(1:16),3],som.perf.df.es[c(1:16),5],col=4)
par(new=T)
plot(som.perf.df.es[c(1:16),3],som.perf.df.es[c(1:16),4],type="p",axes=F,bty="n",xlab="",ylab="",pch=3,col=2,
     ylim=range(pretty(c(-3,3))))
lines(som.perf.df.es[c(1:16),3],som.perf.df.es[c(1:16),4],col=2)
axis(side=4)
mtext(4, text = "DB Index", line =3)
legend(10,5,legend=c("Mean dist","DB index"),pch=c(2,3),col=c(4,2),box.col="white",cex=1)


### 3. INTERPRETATION ####

som_select_es <- 4 # make it flexible! Set the optimal parameterisation for SOM clusters. Important: the selected option ( in this case 4) does not represent the number of clusters but the fouth option( see"som.dims").

final.som.es <- som_test_es[[som_select_es]] # select final clusters

### get vector codebook (= cluster centroid values for each indicator)

code_vec_es <- as.data.frame(final.som.es[[2]]$codes)

### get min/max values for plot xlim and ylim
glob_max <- ceiling(max(code_vec_es))
glob_min <- floor(min(code_vec_es))

### prepare data for plotting
code_vec_es <- t(code_vec_es)
colnames(code_vec_es) <- rep(paste0("Cluster",1:ncol(code_vec_es)))

### get x- and y-dim
ydim <- som.dims.es[som_select_es,1]
xdim <- som.dims.es[som_select_es,2]

### setup colour scheme
pal_setup <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pal_sample <- seq(1,length(pal_setup),length(pal_setup)/16)
pal_sample <- round(pal_sample,0)


####dataframe with municipality and cluster number####
df.exp.es <- cbind(db_inequities[,2],
                   as.data.frame(final.som.es[[2]][[2]])) # classification units (clusters)
colnames(df.exp.es)<- c("comuna","cluster")


#plot Counts of selected SOM (Fig. S2.XX?)
plot(final.som.es[[2]], type="counts")


#Description of the clusters with barplot figure (Figs. 3-4)

barplot_db_es<-as.data.frame(cbind(db_es, comuna=db_inequities$comuna))#add municipality name to database

barplot_db_es<-as.data.frame(cbind(barplot_db_es, cluster=df.exp.es$cluster))#add cluster column
str(barplot_db_es)

barplot_db_es<-barplot_db_es%>%#summarize clusters through columns
  group_by(cluster)%>%
  dplyr::summarise_at(vars(tot_water_sup.y:tot_supl_recreation.y), mean, na.rm=TRUE) 

#export total values (Table S2.xx?)

write.xlsx( barplot_db_es, "Z:/Paper_inequalities/Datos/total_values_cluster_es.xlsx", rownames=TRUE )

#scale data for plotting

barplot_db_es [2:9]<-scale(barplot_db_es[2:9], center=TRUE, scale=TRUE)#z-score data
barplot_db1_es<-as.data.frame(barplot_db_es)

#rearrage data for plotting and annotations

barplot_db1_es<-barplot_db1_es[,-1]
barplot_db1_es<-barplot_db1_es[,c(1,7,6,2,3,4,5,8)]
barplot_db1_es<-as.data.frame(t(barplot_db1_es))
colnames(barplot_db1_es)[1:6]<-c("cluster1","cluster2", "cluster3", "cluster4", "cluster5", "cluster6")
barplot_db1_es$ES<-c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation")

#absolute values for annotation

glimpse(barplot_db1_es)
barplot_es_annotation<-barplot_db1_es[,-7]
rownames(barplot_es_annotation)<-c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation")
barplot_es_annotation<-data.frame(t(barplot_es_annotation))
barplot_es_annotation<-abs(barplot_es_annotation)

annotation_es<-barplot_es_annotation%>%rowwise%>%
  mutate(provisioning=sum(c(water.supply, native.timber,non.native.timber)))

annotation_es<-annotation_es%>%rowwise%>%
  mutate(regulating=sum(c(water.regulation, carbon.sequestration, carbon.storage, erosion.prevention)))

annotation_es<-annotation_es%>%rowwise%>%
  mutate(cultural=c(recreation))

annotation_final_es<-annotation_es%>% 
  mutate(total = sum(c_across(provisioning:cultural))) %>% 
  ungroup() %>% 
  mutate(across(provisioning:cultural, ~ . *100/total))

colnames(annotation_final_es)[1:6]<-c("cluster1","cluster2", "cluster3", "cluster4", "cluster5", "cluster6")


##barplot figures of each archetype with GGPLOT

ESA1<-ggplot(barplot_db1_es, aes(x=factor(ES, levels=ES), y=cluster1, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("ESA-1") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")

ESA1<-ESA1 +annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "39.1")+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "29.8")+
  
  annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "31.1")
ESA1


ESA2<-ggplot(barplot_db1_es, aes(x=factor(ES, levels=ES), y=cluster2, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("ESA-2") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")

ESA2<-ESA2 +annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "54.2")+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "33.4")+
  
  annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "12.4")
ESA2



ESA3<-ggplot(barplot_db1_es, aes(x=factor(ES, levels=ES), y=cluster3, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("ESA-3") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")

ESA3<-ESA3 +annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "43.4")+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "46.9")+
  
  annotate("rect", xmin=6.8,xmax=7.8,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7.3,y=1.85, label = "9.7")
ESA3

ESA4<-ggplot(barplot_db1_es, aes(x=factor(ES, levels=ES), y=cluster4, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("ESA-4") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")

ESA4<-ESA4 +annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "38")+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "48.2")+
  
  annotate("rect", xmin=5.4,xmax=6.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=5.9,y=1.85, label = "13.8")
ESA4


ESA5<-ggplot(barplot_db1_es, aes(x=factor(ES, levels=ES), y=cluster5, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("ESA-5") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")
ESA5<-ESA5 +annotate("rect", xmin=3.2,xmax=4.2,ymin=1.9,ymax=2.2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=2.05, label = "74.2")+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "17")+
  
  annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "8.8")
ESA5


ESA6<-ggplot(barplot_db1_es, aes(x=factor(ES, levels=ES), y=cluster6, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("ESA-6") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
ESA6<-ESA6 +annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "53.2")+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "34.8")+
  
  annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "12")
ESA6


#put all graphs together with gridExtra

grid.arrange(ESA1,ESA2,ESA3,ESA4,ESA5,ESA6, nrow=3, ncol=2, top="")

#spatialize archetypes

all_municipalities<-readOGR(dsn = "Z:/Paper_inequalities/Datos/espaciales", layer="todas_comunas_sort")#import reference shapefile

all_municipalities@data$NOMBRE<-tolower(all_municipalities@data$NOMBRE)
all_municipalities@data$ID<-1:nrow(all_municipalities@data)#uptade row ID for latter matching

#modifiy ID column (municipality name) so that there is a match between spatial and non-spatial data

all_municipalities@data[84,4]<-"maullin"
all_municipalities@data[167,4]<-"trehuaco"
all_municipalities@data[57,4]<-"hualane"
all_municipalities@data[170,4]<-"vichuquen"
all_municipalities@data[56,4]<-"hualaihue"

#merge
spdf_es<-sp::merge(all_municipalities,df.exp.es, by.x="NOMBRE",by.y="comuna", sort=FALSE)
data.table(colnames(spdf_es@data))##check if columns (variables were merged)
spdf_es@data<-spdf_es@data[,-c(2:10)]##eliminate columns from all_municipalities shape because they are no necessary an remain with the ones of "db_inequities" which will be used in the analysis

id_merge_es<-as.data.frame(rbind(spdf_es@data$NOMBRE,df.exp.es$comuna))#corroborate order of municipalities

writeOGR(spdf_es,dsn = "Z:/Paper_inequalities/Datos/espaciales", layer="clusters_ES_SOM", driver = "ESRI Shapefile", overwrite_layer = TRUE, verbose = TRUE)##Fig. 2 Panel a



####Inequality data####

glimpse(db_inequities)#view of database
data.table(colnames(db_inequities))#column numbers
# create placeholder variable and select variables to be used
db_ineq <- db_inequities [,c(3:10)] #select just Gini variables for developing SOM
 
### 1 PREPARE DATA FOR SOMS: Z-TRANSFORM + MATRIX #####

db_ineq_scaled<-scale(db_ineq, center=TRUE, scale=TRUE)#z-scale data and assign new name

db_ineq_matrix <- as.matrix(db_ineq_scaled)#dataframe to matrix

### 2. SOMS #####
### 2.1 SET UP PARAMETERS #####
som.dims <- as.data.frame(matrix(ncol=2,nrow=16))
colnames(som.dims) <- c("som_cols","som_rows")
som.dims[,1] <- c(2,3,2,3,4,3,5,4,5,4,6,5,7,6,5,6)
som.dims[,2] <- c(1,1,2,2,2,3,2,3,3,4,3,4,3,4,5,5)

set.seed(2022)

col_sel <- c(1:length(db_ineq))

som_test <- list()

for(i in 1:nrow(som.dims)){
  temp.name <- paste0("SOM_cols",som.dims[i,1],"_rows",som.dims[i,2])
  temp.som <- supersom(db_ineq_matrix[,col_sel],
                  grid=somgrid(som.dims[i,1],som.dims[i,2],topo="hexagonal"),
                  rlen = 1000,
                  alpha = c(0.05,0.01),
                  keep.data=T,
                  maxNA.fraction=0.90)
  #between <- sum(rdist(as.matrix(temp.som$grid$pts)))
  #within <- sum(temp.som$distances)
  #som.homogeneity <- within
  #som.variance <- between/(between+within)
  som.dist.mean <- mean(temp.som$distances)
  som.dist.sd <- sd(temp.som$distances)
  temp.DB <- index.DB(x=as.data.frame(temp.som$data),
                      cl=temp.som$unit.classif,
                      centrotypes="centroids")
  som_test[[i]] <- list(temp.name,temp.som,som.dist.mean,som.dist.sd,temp.DB)
}

length(som_test)
length(som_test[[14]])

### 2.1 PERFORMANCE ####

som.perf.df <- as.data.frame(matrix(ncol=6,nrow=length(som_test)))
colnames(som.perf.df) <- c("name","cluster_combo","clus","DB","mean","sd")
som.perf.df[,2] <- c("2x1","3x1","2x2","3x2","4x2","3x3","5x2","4x3","5x3","4x4","6x3","5x4","7x3","6x4","5x5", "6X5")
som.perf.df[,3] <- c(2,3,4,6,8,9,10,12,15,16,18,20,21,24,25,30)

for(i in 1:length(som_test)){
  som.perf.df[i,1] <- som_test[[i]][[1]]   # name
  som.perf.df[i,4] <- som_test[[i]][[5]]$DB # DB
  som.perf.df[i,5] <- som_test[[i]][[3]]    # mean
  som.perf.df[i,6] <- som_test[[i]][[4]]    # sd
}

#visualize and select the number of archetypes

par(mar = c(5, 4, 4, 5))
plot(som.perf.df[c(1:16),3],som.perf.df[c(1:16),5],pch=2,col=4,
     ylim=range(pretty(c(min(som.perf.df[c(1:16),5]),max(som.perf.df[c(1:16),5])))),
     ylab="Mean distance to cluster centroid",xlab="Archetype size for (in)equality")
lines(som.perf.df[c(1:16),3],som.perf.df[c(1:16),5],col=4)
par(new=T)
plot(som.perf.df[c(1:16),3],som.perf.df[c(1:16),4],type="p",axes=F,bty="n",xlab="",ylab="",pch=3,col=2,
     ylim=range(pretty(c(1,3))))
lines(som.perf.df[c(1:16),3],som.perf.df[c(1:16),4],col=2)
axis(side=4)
mtext(4, text = "DB Index", line =3)
legend(10,5,legend=c("Mean dist","DB index"),pch=c(2,3),col=c(4,2),box.col="white",cex=1)



### 3. INTERPReTATION ####

som_select <- 7# make it flexible! Set the optimal parameterisation for SOM clusters. Important: the selected option ( in this case 4) does not represent the number of clusters but the fouth option( see"som.dims").

final.som <- som_test[[som_select]] # select final clusters

code_vec <- as.data.frame(final.som[[2]]$codes)

### get min/max values for plot xlim and ylim
glob_max <- ceiling(max(code_vec))
glob_min <- floor(min(code_vec))

### prepare data for plotting
code_vec <- t(code_vec)
colnames(code_vec) <- rep(paste0("Cluster",1:ncol(code_vec)))

### get x- and y-dim
ydim <- som.dims[som_select,1]
xdim <- som.dims[som_select,2]

### setup colour scheme
pal_setup <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pal_sample <- seq(1,length(pal_setup),length(pal_setup)/16)
pal_sample <- round(pal_sample,0)


###dataframe with municipality and cluster number ####
df.exp <- cbind(db_inequities[,2],
                as.data.frame(final.som[[2]][[2]])) # classification units (clusters)
colnames(df.exp)<- c("comuna","cluster")


#plot Counts of selected SOM (Fig. S2.XX?)

plot(final.som[[2]], type="counts")


#Description of the clusters with barplot figure (Figs. 3-4)

barplot_db<-as.data.frame(cbind(db_ineq, comuna=db_inequities$comuna))#add municipality name to database

barplot_db<-as.data.frame(cbind(barplot_db, cluster=df.exp$cluster))#add cluster column


barplot_db<-barplot_db%>%
  group_by(cluster)%>%
  dplyr::summarise_at(vars(gini_new_ws:gini_new_recre), mean, na.rm=TRUE)

#export total values (Table S2.xx?)

write.xlsx( barplot_db, "Z:/Paper_inequalities/Datos/total_values_cluster_gini.xlsx", rownames=TRUE, overwrite = TRUE )

#scale data for plotting

barplot_db [2:9]<-scale(barplot_db[2:9], center=TRUE, scale=TRUE)#scale with z-score data
barplot_db1<-as.data.frame(barplot_db)

#rearrage data for plotting and annotations


barplot_db1<-barplot_db1[,-1]
barplot_db1<-barplot_db1[,c(1,7,6,2,3,4,5,8)]
barplot_db1<-as.data.frame(t(barplot_db1))
colnames(barplot_db1)[1:10]<-c("cluster1","cluster2", "cluster3", "cluster4", "cluster5", "cluster6","cluster7","cluster8","cluster9","cluster10")
barplot_db1$ES<-c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation")


#absolute values for annotation

glimpse(barplot_db1)
barplot_annotation<-barplot_db1[,-11]
rownames(barplot_annotation)<-c("water supply","native timber","non-native timber","water regulation","carbon sequestration","carbon storage","erosion prevention","recreation")
barplot_annotation<-data.frame(t(barplot_annotation))
barplot_annotation<-abs(barplot_annotation)

annotation<-barplot_annotation%>%rowwise%>%
  mutate(provisioning=sum(c(water.supply, native.timber,non.native.timber)))

annotation<-annotation%>%rowwise%>%
  mutate(regulating=sum(c(water.regulation, carbon.sequestration, carbon.storage, erosion.prevention)))

annotation<-annotation%>%rowwise%>%
  mutate(cultural=c(recreation))

annotation_final<-annotation%>% 
  mutate(total = sum(c_across(provisioning:cultural))) %>% 
  ungroup() %>% 
  mutate(across(provisioning:cultural, ~ . *100/total))

annotation_final<-t(annotation_final)
colnames(annotation_final)[1:10]<-c("cluster1","cluster2", "cluster3", "cluster4", "cluster5", "cluster6","cluster7","cluster8","cluster9","cluster10")

##barplot figures of each archetype with GGPLOT

GiniA1<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster1, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-1") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")
GiniA1<-GiniA1 + annotate("rect", xmin=5.4,xmax=6.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=5.9,y=1.85, label = "47.5", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "41", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "11.5", size=3)
GiniA1



GiniA2<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster2, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-2") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
GiniA2<-GiniA2 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "46.8", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "41.3", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "11.9", size=3)
GiniA2



GiniA3<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster3, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-3") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
GiniA3<-GiniA3 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "26.2", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "61", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "12.8", size=3)
GiniA3




GiniA4<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster4, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-4") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")
GiniA4<-GiniA4 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "25.5", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "59.7", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "14.8", size=3)
GiniA4




GiniA5<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster5, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-5") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
GiniA5<-GiniA5 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.5,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.75, label = "56.3", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.5,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.75, label = "43.2", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.5,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.75, label = "0.5", size=3)
GiniA5



GiniA6<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster6, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,3,0.5), labels=seq(-2,3,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-6") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
GiniA6<-GiniA6 + annotate("rect", xmin=6.3,xmax=7.3,ymin=2.5,ymax=2.8, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=6.8,y=2.65, label = "27.5", size=3)+
  
  annotate("rect", xmin=3.7,xmax=4.7,ymin=2.5,ymax=2.8, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=4.2,y=2.65, label = "53.4", size=3)+
  
  annotate("rect", xmin=0.5,xmax=1.5,ymin=2.5,ymax=2.8, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=1,y=2.65, label = "19.1", size=3)
GiniA6



GiniA7<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster7, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-7") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")
GiniA7<-GiniA7 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "16", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "70.8", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "13.2", size=3)
GiniA7



GiniA8<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster8, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-8") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
GiniA8<-GiniA8 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "31.9", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "49.9", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "18.2", size=3)
GiniA8



GiniA9<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster9, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-9") +
  theme_minimal()+
  theme(legend.position="none", axis.text.y = element_blank()) +
  labs(x="", y ="z-score")
GiniA9<-GiniA9 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "33.9", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "54.5", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "11.6", size=3)
GiniA9



GiniA10<-ggplot(barplot_db1, aes(x=factor(ES, levels=ES), y=cluster10, color=ES))+ geom_bar(stat="identity",position="dodge", fill="white", lwd=1.2) +
  coord_flip()+
  scale_y_continuous(breaks=seq(-2,2,0.5), labels=seq(-2,2,0.5))+
  scale_x_discrete(limits=rev)+
  scale_color_manual(breaks=c("Gini water supply", "Gini native timber", "Gini non-native timber","Gini water regulation", "Gini carbon sequestration", "Gini carbon storage", "Gini erosion prevention", "Gini recreation"), values=c("red","red", "red","purple","purple","purple","purple","green"))+
  ggtitle("InqA-10") +
  theme_minimal()+
  theme(legend.position="none",axis.text.y = element_text(size=11)) +
  labs(x="", y ="z-score")
GiniA10<-GiniA10 + annotate("rect", xmin=6.5,xmax=7.5,ymin=1.7,ymax=2, alpha=.5, fill="white", color="red",lwd=1)+
  annotate("text", x=7,y=1.85, label = "33.9", size=3)+
  
  annotate("rect", xmin=3.2,xmax=4.2,ymin=1.7,ymax=2, alpha=.5, fill="white", color="purple",lwd=1)+
  annotate("text", x=3.7,y=1.85, label = "54.5", size=3)+
  
  annotate("rect", xmin=0.4,xmax=1.4,ymin=1.7,ymax=2, alpha=.5, fill="white", color="green",lwd=1)+
  annotate("text", x=0.9,y=1.85, label = "11.6", size=3)
GiniA10


#put all graphs together with gridExtra

grid.arrange(GiniA1,GiniA2,GiniA3,GiniA4,GiniA5,GiniA6,GiniA7,GiniA8,GiniA9,GiniA10, nrow=4, ncol=3, top="")


#spatialize archetypes

#merge
spdf<-sp::merge(all_municipalities,df.exp, by.x="NOMBRE",by.y="comuna", sort=FALSE)
data.table(colnames(spdf@data))##check if columns (variables were merged)
spdf@data<-spdf@data[,-c(2:10)]##eliminate columns from all_municipalities shape because they are no necessary an remain with the ones of db_inequities which will be used in the analysis
id_merge<-as.data.frame(rbind(spdf@data$NOMBRE,df.exp$comuna))#corroborate order

writeOGR(spdf,dsn = "Z:/Paper_inequalities/Datos/espaciales", layer="clusters_gini_SOM", driver = "ESRI Shapefile", overwrite_layer = TRUE, verbose = TRUE)#(Fig. 4. Panel a)




#Overlap analyis and generation of Fig. 5

#libraries
library(lattice)
library(extrafont)

#import data
table_ES_Gini <- read.xlsx("Z:/Paper_inequalities/Datos/clusters_ES_gini.xlsx")
table_Gini_ES<- read.xlsx("Z:/Paper_inequalities/Datos/clusters_Gini_ES.xlsx")

ES <-table_ES_Gini[,-3]#delete column with area
Gini <-table_Gini_ES[,-3]#delete column with area

#font_import("") #only execute the first time
loadfonts(device = "win", quiet = TRUE)

p1 <- ggplot(ES, aes(colour=PERCENTAGE,x=cluster_1, y=cluster, size = PERCENTAGE)) + 
  geom_point() +
  scale_x_continuous(breaks=seq(1,12,1),labels=c("InqA-1", "InqA-2","InqA-3","InqA-4","InqA-5","InqA-6","InqA-7","InqA-8","InqA-9","InqA-10","InqA-11","InqA-12"))+
  scale_y_continuous(breaks=seq(1,12,1), labels=c("ESA-1","ESA-2", "ESA-3", "ESA-4","ESA-5","ESA-6","ESA-7","ESA-8","ESA-9","ESA-10","ESA-11","ESA-12"))+
  scale_color_continuous(breaks = seq(from=0,to=90,by=10),low="lightblue", high="darkblue") +
  theme(panel.background = element_rect(fill = "white")) +
  labs(size=expression(atop("InqA contribution","per ESA [%]")), colour=expression(atop("InqA contribution","per ESA [%]"))) +
  xlab("") +
  ylab("") +
  scale_size_continuous(range=c(1,10),breaks = seq(from=0,to=90,by=10)) + 
  theme(legend.background = element_rect(fill="white"),
        plot.margin=margin (0.5,0.5,0.5,0.5, "cm"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_line(colour="grey80"),
        legend.key = element_rect(colour = "white",fill="white"),
        text = element_text(family="Times New Roman", size=10),
        axis.text = element_text(colour = "grey10"),
        axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(colour = guide_legend(), size = guide_legend()) 



p2 <- ggplot(Gini, aes(colour=PERCENTAGE,x=cluster, y=cluster_1, size = PERCENTAGE)) +
  geom_point() +
  scale_x_continuous(breaks=seq(1,12,1),labels=c("InqA-1", "InqA-2","InqA-3","InqA-4","InqA-5","InqA-6","InqA-7","InqA-8","InqA-9","InqA-10","InqA-11","InqA-12"))+
  scale_y_continuous(breaks=seq(1,12,1), labels=c("ESA-1","ESA-2", "ESA-3", "ESA-4","ESA-5","ESA-6","ESA-7","ESA-8","ESA-9","ESA-10","ESA-11","ESA-12"))+
  scale_colour_continuous(breaks = seq(from=0,to=90,by=10),low="lightblue", high="darkblue")+
  theme(panel.background = element_rect(fill = "white")) +
  labs(size=expression(atop("ESA contribution","per InqA [%]")), colour=expression(atop("ESA contribution","per InqA [%]"))) +
  xlab("") +
  ylab("") +
  scale_size_continuous(range=c(1,10), breaks = seq(from=0,to=90,by=10)) +
  theme(legend.background = element_rect(fill="white"),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_line(colour="grey80"),
        legend.key = element_rect(colour = "white",fill="white"),
        plot.margin=margin (0.5,0.5,0.5,0.5, "cm"),
        text = element_text(family="Times New Roman", size=10),
        axis.text = element_text(colour = "grey10"),
        axis.text.x = element_text(angle=90, vjust=0.5)) +
  guides(colour = guide_legend(),size = guide_legend()) #Adaptar seg?n como queramos que se vea el plot.

grid.arrange(p1,p2,nrow=2)#final bubble plot figure for the overlapping analysis


#End of the script#




       












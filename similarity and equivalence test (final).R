library(dismo)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(usdm)
library(foreign)
library(spocc)
library(corrplot)
library(usdm)
library(XML)
library(dplyr)
library(raster)
library(dismo)  
library(SDMTools)   # varios analisis relacionados con modelos de nicho; prueba de Warren et al. (2008)
library(ecospat)    # comparacion de nichos en espacio ambiental sensu Broenimann et al. (2012)


##Read .DBF file that contain the occurrence records for both Piaya lineages.
piaya_mex <- read.csv("C:/project_sigs/Piaya/piaya_mexicana2.csv", header = T, sep = ",")
piaya_terh <- read.csv("C:/project_sigs/Piaya/piaya_therma2.csv", header = T, sep = ",")

# The layers have to be cut by the relevant geographic area for the comparison of niches
setwd("C:/project_sigs/Piaya/capas_clima/ascii/present/") 
pca_path <- list.files(".",pattern = "*.asc$",full.names = T)###crea el stack 
varclim <- stack(pca_path)

# Generate a table that for each pixel (coordinate x, y) with the values of all environmental variables
climpunto <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)

# Extract environmental data for each layer in varclim of each coordinate in climpoint
clim <- extract(varclim, climpunto)

# Format the clim table to be a normal table with two first columns x and y
clim <- data.frame(coordinates(climpunto),clim)
clim <- subset(clim, !is.na(bio_02) & !is.na(bio_07) & !is.na(bio_08) & !is.na(bio_12) & !is.na(bio_13) 
               & !is.na(bio_14) & !is.na(bio_15) & !is.na(bio_17) & !is.na(bio_19))


# Delete  points that are less than 1 km away
# If 1 grade = 111 km , so 1 km = 0.041666669
occ.sp1 <- piaya_mex[2:3]
occ.sp2 <- piaya_terh[2:3]

#Make a dataframe with the geographic points and enviromental values
occ_sp1<-na.exclude(ecospat.sample.envar(dfsp=occ.sp1,colspxy=1:2, 
                                         colspkept=1:2,dfvar=clim,
                                         colvarxy=1:2,colvar="all",resolution=0.041666669))

occ_sp2<-na.exclude(ecospat.sample.envar(dfsp=occ.sp2,colspxy=1:2,
                                         colspkept=1:2,dfvar=clim, 
                                         colvarxy=1:2,colvar="all",resolution=0.041666669))

#Iterations 
iterations<-1000
# Enviromental resolution for the density presences 
R=500

########################################################################
#                                Enviromental PCA                      #
########################################################################
#PCA data for study area and species presence for two species
data<-rbind(clim[,3:11],occ_sp1[,3:11],occ_sp2[,3:11])

data <- subset(data, !is.na(bio_02) & !is.na(bio_07) & !is.na(bio_08) & !is.na(bio_12) & !is.na(bio_13) 
               & !is.na(bio_14) & !is.na(bio_15) & !is.na(bio_17) & !is.na(bio_19))

# 0 to ocurrences and 1 for study area
w<-c(rep(1,nrow(clim)),rep(0,nrow(occ_sp1)),rep(0,nrow(occ_sp2)))

# PCA with all the data for study area, the presence data will be calculate but not to calibrate PCA.
pca.cal <-dudi.pca(data, row.w = w, center = T, scale = T, scannf = F, nf = 2)

# Data clim and data species:
row.clim <- 1:nrow(clim)
row.sp1 <-  (1+nrow(clim)):(nrow(clim) + nrow(occ_sp1))
row.sp2 <-  (1 + nrow(clim) + nrow(occ_sp1)) : (nrow(clim) + nrow(occ_sp1) + nrow(occ_sp2))

summary(pca.cal)

# PCA coordinates for study area and species
scores.clim <- pca.cal$li[row.clim,] #all rowns
scores.sp1  <- pca.cal$li[row.sp1,]   #sp1
scores.sp2  <- pca.cal$li[row.sp2,]   #sp2


#Variable contribution for each variable to each PCA component
ecospat.plot.contrib(contrib=pca.cal$co, eigen=pca.cal$eig)

###Density ocurrence####
#Presence points  in the enviromental space
z1<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp1, R=300) #Spp 1
z2<-ecospat.grid.clim.dyn (scores.clim, scores.clim, scores.sp2, R=300) #Spp 2

#Overlap statistics (D -  Schoener, I -  Warren)
D.overlap <- ecospat.niche.overlap (z1=z1, z2=z2, cor=T)$D
D.overlap

#Graphical density ocurrence in enviromental space
windows() # en mac o windows() en PC
par(mfrow=c(1,2))
ecospat.plot.niche (z1, title="Piaya cayana mexicana", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z2, title="Piaya cayana thermophila", name.axis1="PC1", name.axis2="PC2", cor=F)


###############################################################
#               Equivalency niche test                        #
###############################################################
a.dyn<-ecospat.niche.equivalency.test(z1=z1 , z2=z2, rep=1000)
#Niche equivalence 
windows() # en mac o windows() en PC
par(mfrow=c(1,2))
ecospat.plot.overlap.test(a.dyn,"D","Equivalency")
ecospat.niche.overlap(z1, z2, cor= F)

###############################################################
#               Similarity niche test                         #
###############################################################
b.dyn_phor_phai<-ecospat.niche.similarity.test(z1=z1 , z2=z2, rep=1000, alternative = "greater",
                                               rand.type=2)
b.dyn_phor_phai2<-ecospat.niche.similarity.test(z1=z2 , z2=z1, rep=1000, alternative = "greater",
                                                rand.type=2)
#Similarity test with "D" value
windows() # en mac o windows() en PC
par(mfrow=c(1,2))
ecospat.plot.overlap.test(b.dyn_phor_phai,"D","Similarity - P.c. mexicana vs. P.c. thermophila")
ecospat.plot.overlap.test(b.dyn_phor_phai2,"D","Similarity - P.c. thermophila vs. P.c. mexicana")

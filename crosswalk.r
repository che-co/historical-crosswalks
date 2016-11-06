#!/bin/R

library(rgdal)
library(rgeos)
library(maptools)
library(RColorBrewer)
colors <- colorRampPalette(c("red","grey"))(4)

#### Function to add w/o NA ####

`%*na%` <- function(a, b){ifelse(is.na(a), b, ifelse(is.na(b), a, a*b))}

#### Paths for Shape Files ####

dsn1860 <- "slavery/shapefiles/historical_counties/nhgis0001_shape/tl2000"
dsn2015 <- "slavery/shapefiles/current_counties"
layer1860 <- "US_county_1860"
layer2015 <- "cb_2015_us_county_20m"
proj <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
tol <- 1e-7

#### Clean Historical Data ####

nhgis1860.slv <- read.csv("./slavery/data/nhgis0003_csv/nhgis0003_ds14_1860_county.csv")
nhgis1860.slv <- nhgis1860.slv[, c("GISJOIN", "AH3001", "AH3002", "AH3003")]
names(nhgis1860.slv)[2:4] <- c("n_wht", "n_free_blk", "n_slv")

nhgis1860.pop <- read.csv("./slavery/data/nhgis0004_csv/nhgis0004_ds14_1860_county.csv")
names(nhgis1860.pop)[length(names(nhgis1860.pop))] <- "t_pop"
nhgis1860.pop <- data.frame(nhgis1860.pop, nhgis1860.slv[match(nhgis1860.pop$GISJOIN, nhgis1860.slv$GISJOIN), names(nhgis1860.slv)[2:4]])

nhgis1860.pop$n_blk <- nhgis1860.pop$n_free_blk + nhgis1860.pop$n_slv
nhgis1860.pop$n_free <- rowSums(cbind(nhgis1860.pop$t_pop - nhgis1860.pop$n_slv))
nhgis1860.pop$frac_free <- nhgis1860.pop$n_free / nhgis1860.pop$t_pop
nhgis1860.pop$frac_free_blk <- nhgis1860.pop$n_free_blk / nhgis1860.pop$n_blk
nhgis1860.pop[is.nan(nhgis1860.pop$frac_free_blk), "frac_free_blk"] <- NA 
nhgis1860.pop$frac_slv <- nhgis1860.pop$n_slv / nhgis1860.pop$t_pop  

#### Calculate Spatial Adjustments ####

us.cty1860 <- readOGR(dsn=dsn1860, layer=layer1860)
us.cty2015 <- readOGR(dsn=dsn2015, layer=layer2015)
us.cty1860 <- spTransform(us.cty1860, CRS(proj))
us.cty2015 <- spTransform(us.cty2015, CRS(proj))
us.cty2015.ak <- us.cty2015[us.cty2015$STATEFP=="02", ]
us.cty2015.ak <- elide(us.cty2015.ak, rotate=-50)
us.cty2015.ak <- elide(us.cty2015.ak, scale=max(apply(bbox(us.cty2015.ak), 1, diff))/2.3)
us.cty2015.ak <- elide(us.cty2015.ak, shift=c(-2100000, -2500000))
proj4string(us.cty2015.ak) <- proj4string(us.cty2015)
us.cty2015.hi <- us.cty2015[us.cty2015$STATEFP=="15", ]
us.cty2015.hi <- elide(us.cty2015.hi, rotate=-35)
us.cty2015.hi <- elide(us.cty2015.hi, shift=c(5400000, -1400000))
proj4string(us.cty2015.hi) <- proj4string(us.cty2015)
us.cty2015 <- us.cty2015[!us.cty2015$STATEFP %in% c("02","15","72"),]

us.cty1860@data <- data.frame(us.cty1860, nhgis1860.pop[match(us.cty1860@data[, "GISJOIN"], nhgis1860.pop[, "GISJOIN"]), names(nhgis1860.pop)[8:length(names(nhgis1860.pop))]])

us.cty.instersection <- gIntersection(us.cty2015, us.cty1860, byid=T)
us.cty.area.inter <- gArea(us.cty.instersection, byid=T)

tmp <- as.character(scan(text=names(us.cty.area.inter)))
tmp <- matrix(tmp, nrow=length(us.cty.area.inter), ncol=2, byrow=T, dimnames=list(c(NULL), c("cty2015", "cty1860")))
tmp <- as.data.frame(tmp, stringAsFactors=F)
overlap <- cbind(tmp, us.cty.area.inter)
rownames(overlap) <- NULL
us.cty1860.area <- gArea(us.cty1860, byid=T)
us.cty1860.area <- data.frame(cty1860=names(us.cty1860.area), area1860=us.cty1860.area, stringsAsFactors=F)
overlap <- data.frame(overlap, area1860=us.cty1860.area[match(overlap$cty1860, us.cty1860.area$cty1860), "area1860"])
names(overlap)[3] <- "area_inter"
overlap$frac_1860 <- overlap$area_inter / overlap$area1860
tmp <- aggregate(frac_1860 ~ cty1860, data=overlap, sum)
names(tmp)[2] <- "sum_frac_1860"
overlap <- data.frame(overlap, sum_frac_1860=tmp[match(overlap$cty1860, tmp$cty1860), "sum_frac_1860"])
overlap$scaled_frac_1860 <- overlap$frac_1860 * (1/overlap$sum_frac_1860)

vars <- c("t_pop", "n_wht", "n_free_blk", "n_slv", "n_blk", "n_free")

us.cty2015@data[, vars] <- 0 

for(i in levels(overlap$cty2015)) us.cty2015@data[i, vars] <- colSums(overlap$scaled_frac_1860[which(overlap$cty2015==i)] * us.cty1860@data[overlap$cty1860[which(overlap$cty2015==i)], vars])

vars_wona <- paste0(vars, "_wona")

us.cty2015@data[, vars_wona] <- 0 

for(i in levels(overlap$cty2015)) us.cty2015@data[i, vars_wona] <- colSums(overlap$scaled_frac_1860[which(overlap$cty2015==i)] * us.cty1860@data[overlap$cty1860[which(overlap$cty2015==i)], vars], na.rm=T)


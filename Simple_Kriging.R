library(sp)
library(gstat)
library(automap)
library(maptools)
library(rgdal)
library(raster)
library(sf)
library(mgcv)

setwd("E:/Masters_Program/Semester 3/Geostatistics/Practicals/Shapefiles/Datasets")
data=read.csv("Ozone_samples.csv", header = T)
hillshad=raster("ca_hillshade.tif")
dist_city=raster("EucDist_cities.tif")


# convert simple data frame into spatial point data frame and specify correct EPSG
proj4string <- "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +
ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
statPoints <- SpatialPointsDataFrame(coords      = data[,c("X","Y")], 
                                     data        = data,
                                     proj4string = CRS(proj4string))


# Compute variogram parameters automatically that will be required later.
ok_param = autofitVariogram(OZONE~1, 
                            statPoints, 
                            model = c("Sph", "Exp", "Gau"),
                            kappa = c(0.05, seq(0.2, 2, 0.1,0.5,0.001), 5, 10,15,20,25))
summary(ok_param)
plot(ok_param)


##Create experimental variogram object based on sample points
variog=variogram(OZONE~1, statPoints) 


## create variogram model by passing appropriate variogram parameters - Equivalent to model training
est <- vgm(0.0005844913,"Sph",238410.1,0.0000000000)

## fit Variogram model in to the variogram object
fitted <- fit.variogram(variog, est)
plot(variog, model=fitted)


Grid_locs=read.csv("Grid_loc.csv", header = T)

# convert simple data frame into spatial point data frame and specify correct EPSG
proj4string <- "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +
y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
pred_grid=SpatialPointsDataFrame(coords      = Grid_locs[,c("X","Y")], 
                                    data        = Grid_locs,
                                    proj4string = CRS(proj4string))

head(pred_grid)


# interpolate using simple kriging method on unsampled grid points
b=mean(data$OZONE)
s.k = krige(OZONE~1, statPoints, pred_grid, fitted,beta=b)
summary(s.k)



# write the regression results as csv table
write.csv(s.k, file = "Interpolated_ozone.csv")
interpolated_residuals = read.csv("Interpolated_ozone.csv",header = T)
print(interpolated_residuals)


ext=extent(hillshad)

## create raster object using same raster parameters as control raster
raster_object = raster(ncol=hillshad@ncols, nrow=hillshad@nrows, 
                       xmn=ext@xmin, xmx=ext@xmax, ymn=ext@ymin, ymx=ext@ymax)


# retrieve predicted and variance data
Predicted_values_df=interpolated_residuals[,c("X","Y","var1.pred")]
Variance_values_df=interpolated_residuals[,c("X","Y","var1.var")]
# will need to rename colnames for raster
colnames(Predicted_values_df) <- c('x', 'y', 'vals')
colnames(Variance_values_df) <- c('x', 'y', 'vals')

# use rasterize to create predicted raster
Predicted_rasterize <- rasterize(x=Predicted_values_df[, 1:2], # lon-lat data
                                y=raster_object, # raster object
                                field=Predicted_values_df[, 3], # vals to fill raster with
                                fun=mean) # aggregate function

# use rasterize to create Variance raster
Variance_rasterize <- rasterize(x=Variance_values_df[, 1:2], # lon-lat data
                                y=raster_object, # raster object
                                field=Variance_values_df[, 3], # vals to fill raster with
                                fun=mean) # aggregate function

# write the raster data as tiff
Predicted_raster<-writeRaster(Predicted_rasterize,'Predicted_Ozone.tiff',overwrite=TRUE)
Variance_raster<-writeRaster(Variance_rasterize,'Ozone_Variance.tiff',overwrite=TRUE)


spplot(Predicted_raster, main = "Simple kriging predictions")
spplot(Variance_raster,  main = "Simple kriging variance")

# validation using measured and predicted ozone conc
# convert simple data frame into spatial point data frame and specify correct EPSG

## Split sampled data into training and testing data
set.seed(101)
n=nrow(data)
train_index=sample(1:n,size = round(0.6*n),replace=FALSE)
train=data[train_index,]
test=data[-train_index,]


proj4string <- "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +
ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
validation <- SpatialPointsDataFrame(coords      = test[,c("X","Y")], 
                                     data        = test,
                                     proj4string = CRS(proj4string))

head(validation)

valid <- krige(OZONE~1, statPoints, validation, fitted,beta=b)
write.csv(valid,'predicted_actual_ozone_conc_points_2.csv')
summary(valid)
diff <- validation$OZONE - valid$var1.pred 
summary(diff)

sum(diff)/length(diff) # mean error (bias)
sqrt(sum(diff^2)/length(diff)) # RMSE (precision)
median(validation$OZONE) # median error



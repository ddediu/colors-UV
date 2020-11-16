# Preprocess various types of covariates for the main analyses
# both for the languages and for their families
# Dan Dediu, 2020

##
## Load the language data ####
##

d_colors <- read.table("./input_files/data_colors.csv", header=TRUE, sep=",", quote='"', stringsAsFactors=FALSE); # comma-separated double-quoted CVS file
# Alternative view of the globe: all longitudes > 180 are "flipped" to negative values:
d_colors$longitude_180        <- ifelse(d_colors$longitude <= 180,        d_colors$longitude,        d_colors$longitude - 360);
d_colors$longitude_180_family <- ifelse(d_colors$longitude_family <= 180, d_colors$longitude_family, d_colors$longitude_family - 360);

# Extract the geographic coordinates:
# For entries that have precisely the same geographic coordinates as other entries (if any), add 0.01 degrees (on the order of hundreds of meters/the scale of towns or villages):
same_coords <- duplicated(d_colors[, c("latitude", "longitude_180")]); 
d_colors$latitude[ same_coords ] <- d_colors$latitude[ same_coords ] + 0.01; d_colors$longitude_180[ same_coords ] <- d_colors$longitude_180[ same_coords ] + 0.01; 

# Assemble together the location of the languages and of their families:
d_coords <- d_colors[,c("glottocode_family", "latitude_family", "longitude_family", "longitude_180_family")];
names(d_coords) <- c("glottocode", "latitude", "longitude", "longitude_180");
d_coords <- rbind(d_colors[,c("glottocode", "latitude", "longitude", "longitude_180")], d_coords);
d_coords <- unique(d_coords);


##
## Climate data from WorldClim ####
##

# General climate data from the World Clim Database...
if( !file.exists("./input_files/data_climate.tsv") )
{
  # This is based on the code in:
  # Bentz, C., Dediu, D., Verkerk, A., & Jäger, G. (2018). The evolution of language families is shaped by the environment beyond neutral drift. Nature Human Behaviour, 2(11), 816. https://doi.org/10.1038/s41562-018-0457-6
  # Download data from the World Clim Database and compute the first two PCs
  
  # Variables used are:
  wolrdclim.bio.vars <- read.table(text="
VAR   = DESCRIPTION
BIO1  = Annual Mean Temperature
BIO2  = Mean Diurnal Range (Mean of monthly (max temp - min temp))
BIO3  = Isothermality (BIO2/BIO7) (* 100)
BIO4  = Temperature Seasonality (standard deviation *100)
BIO5  = Max Temperature of Warmest Month
BIO6  = Min Temperature of Coldest Month
BIO7  = Temperature Annual Range (BIO5-BIO6)
BIO8  = Mean Temperature of Wettest Quarter
BIO9  = Mean Temperature of Driest Quarter
BIO10 = Mean Temperature of Warmest Quarter
BIO11 = Mean Temperature of Coldest Quarter
BIO12 = Annual Precipitation
BIO13 = Precipitation of Wettest Month
BIO14 = Precipitation of Driest Month
BIO15 = Precipitation Seasonality (Coefficient of Variation)
BIO16 = Precipitation of Wettest Quarter
BIO17 = Precipitation of Driest Quarter
BIO18 = Precipitation of Warmest Quarter
BIO19 = Precipitation of Coldest Quarter
                                 ", 
                                   header=TRUE, sep="=", strip.white=TRUE, stringsAsFactors=FALSE, quote="");
  
  library(rgdal);
  library(raster);
  
  # Load the various data (this is downloaded only 1st time and cached locally):
  bio_raster <- getData('worldclim', var='bio', res=5, download=TRUE, path="./input_files/");
  
  # Extract the data for the give coordinates:
  bio_data <- extract(bio_raster, as.matrix(d_coords[,c("longitude_180","latitude")])); # don't forget that latitude is y and longitude is x!
  bio_data <- cbind(d_coords[,c("glottocode", "latitude", "longitude")], as.data.frame(bio_data));
  rownames(bio_data) <- bio_data$glottocode;

  # PCA:
  library(fpc);
  library(factoextra);
  PCs <- prcomp( ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19, 
                 data=as.data.frame(bio_data), scale=TRUE, center=TRUE);
  summary(PCs); # PC1 (49.7%), PC2 (24.7%), PC3 (8.6%)
  climate_data <- data.frame("glottocode"=rownames(PCs$x), PCs$x[,1:3]);
  names(climate_data) <- c("glottocode", paste0("clim_PC",1:(ncol(climate_data)-1)));
  climate_data <- merge(climate_data, d_coords[,c("glottocode", "latitude", "longitude")], by="glottocode", all.x=TRUE, all.y=FALSE);

  # Plots for checking:
  if( FALSE )
  {
    library(ggplot2);
    library(ggrepel);
    mapWorld <- map_data("world", wrap=c(-20,340), ylim=c(-70,100));
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=climate_data, aes(x = longitude, y = latitude, color=clim_PC1)) +
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=climate_data, aes(x = longitude, y = latitude, color=clim_PC2)) +
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=climate_data, aes(x = longitude, y = latitude, color=clim_PC3)) +
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
  }
  
  # Save it:
  write.table(climate_data, file="./input_files/data_climate.tsv", quote=FALSE, sep="\t", row.names=FALSE);
  
  # Clean up the downloaded files:
  unlink("./input_files/wc5/*.bil", recursive=TRUE);
  unlink("./input_files/wc5/*.hdr", recursive=TRUE);
}


##
## Humidity data from the NOAA ####
##

if( !file.exists("./input_files/data_humidity.tsv") )
{
  library(ncdf4);
  library(raster);
  library(lubridate);
  library(dplyr);
  
  # Download the specific humidity data from the NOAA:
  # http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Diagnostic/.above_ground/.qa/datafiles.html
  if( !file.exists("./input_files/noaa-humidity.nc") )
  {
    if( download.file("http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Diagnostic/.above_ground/.qa/data.nc", 
                      destfile="./input_files/noaa-humidity.nc") != 0 )
    {
      stop("Cannot download specific retrive humidity data from the NOAA (http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Diagnostic/.above_ground/.qa/datafiles.html)!\n");
    }
  }
  
  # Info about the netCDF file:
  ncdf4::nc_open("./input_files/noaa-humidity.nc");
  
  # Load the NOAA humidity netCDF file as a brick:
  noaa_humidity <- raster::brick("./input_files/noaa-humidity.nc");
  
  # Summarize humidity at a given location:
  humidity_summary <- function(latitude, longitude)
  {
    # Extract the whole time series for the location:
    humidity_ts <- extract(noaa_humidity, matrix(c(ifelse(longitude < 0, 360+longitude, longitude), # noaa_humidity logitudes are between 0 and 360
                                                   latitude),ncol=2));
    humidity_ts <- data.frame("humidity"=t(humidity_ts)[,1],
                              "date"=as.Date("1960-01-01", format="%Y-%m-%d") + 
                                months(round(vapply(substring(colnames(humidity_ts),2), 
                                                    function(s) ifelse(substring(s,1,1)==".", -0.5-as.numeric(substring(s,2)), as.numeric(as.character(s))), 
                                                    numeric(1))))); 
    
    # Compute the yearly means, medians, sds and IQRs:
    by_year <- humidity_ts %>%
      group_by(lubridate::year(humidity_ts$date)) %>%
      summarise("annual_mean"=mean(humidity, na.rm=TRUE),
                "annual_median"=median(humidity, na.rm=TRUE),
                "annual_sd"=sd(humidity, na.rm=TRUE),
                "annual_IQR"=IQR(humidity, na.rm=TRUE));
    
    # Compute and return the overall mean, median, sd and IQR:
    return (data.frame("overall_mean"=mean(humidity_ts$humidity, na.rm=TRUE),
                       "overall_median"=median(humidity_ts$humidity, na.rm=TRUE),
                       "overall_sd"=sd(humidity_ts$humidity, na.rm=TRUE),
                       "overall_IQR"=IQR(humidity_ts$humidity, na.rm=TRUE),
                       "mean_annual_mean"=mean(by_year$annual_mean, na.rm=TRUE),
                       "mean_annual_median"=mean(by_year$annual_median, na.rm=TRUE),
                       "mean_annual_sd"=mean(by_year$annual_sd, na.rm=TRUE),
                       "mean_annual_IQR"=mean(by_year$annual_IQR, na.rm=TRUE)));
    
  }
  humidity_data <- do.call(rbind, lapply(1:nrow(d_coords), function(i) humidity_summary(d_coords$latitude[i], d_coords$longitude_180[i])));
  
  humidity_data <- cbind(d_coords[,c("glottocode", "latitude", "longitude")], humidity_data);
  names(humidity_data) <- c("glottocode", "latitude", "longitude", 
                            "humidity_overall_mean", "humidity_overall_median", "humidity_overall_sd", "humidity_overall_IQR", 
                            "humidity_mean_annual_mean", "humidity_mean_annual_median", "humidity_mean_annual_sd", "humidity_mean_annual_IQR");
  
  # Plots for checking:
  if( FALSE )
  {
    library(ggplot2);
    library(ggrepel);
    mapWorld <- map_data("world", wrap=c(-20,340), ylim=c(-70,100));
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=humidity_data, aes(x = longitude, y = latitude, color=humidity_mean_annual_median)) +
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=humidity_data, aes(x = longitude, y = latitude, color=humidity_mean_annual_IQR)) +
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
  }

  # Save the data:
  write.table(humidity_data, file="./input_files/data_humidity.tsv", quote=FALSE, sep="\t", row.names=FALSE);
}


##
## Distance to water ####
##

# Adapted from the code in Bentz, C., Dediu, D., Verkerk, A., & Jäger, G. (2018). The evolution of language families is shaped by the environment beyond neutral drift. Nature Human Behaviour, 2(11), 816. https://doi.org/10.1038/s41562-018-0457-6
if( !file.exists("./input_files/data_dist2water.tsv") )
{
  library(raster);
  library(rgdal);
  library(geosphere);
  library(parallel);
  
  # Unfortunately, the primary data from OpenStreetMap (http://openstreetmapdata.com/), the "Reduced waterbodies as raster masks" (http://openstreetmapdata.com/data/water-reduced-raster)
  # seems to no longer be available for download separately as of July 2020, but we include here, for reproducibility, the original data as downloaded in March 2018 by Dan Dediu
  # in the folder ./input_files/openstreet/ :

  extract.dist2water.for.coords <- function(coords, # a data.frame or matrix with two columns, the 1st being the longitude, and the 2nd the latitude
                                            type=c("ocean", "lakes", "river")[1],
                                            zoom=(0:6)[3], # the zoom level -- the higher the better the precision but more computationally expensive
                                            no.cores=1 # number of cores for mclapply
                                           )
  {
    cat("Processing ",type," ...\n");
    
    if( !(type %in% c("ocean", "lakes", "river")) )
    {
      warning('type must be "ocean", "lakes" or "river"!\n');
      return (NULL);
    }
    
    if( !(zoom %in% 0:6) )
    {
      warning('zoom must be 0 to 6!\n');
      return (NULL);
    }
    
    cat("Extracting the required water type at the required zoom...\n");
    unzip(paste0("./input_files/openstreet/", type, "-raster-reduced-3857.zip"), 
          files=paste0(type, "-raster-reduced-3857/", type, "_raster_z", zoom, ".tif"),
          overwrite=TRUE, exdir="./input_files");
    
    cat("Loading the OSM data...\n");
    os.data <- raster(paste0("./input_files/", type, "-raster-reduced-3857/", type, "_raster_z", zoom, ".tif")); # load the original OopenStreetMaps data
    #crs(os.data);
    
    cat("Converting it to the WGS84 projection...\n");
    os.data.proj <- projectRaster(os.data, crs=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"), method="ngb", over=TRUE); # convert it to the WGS84 projection
    
    cat("Extracting the points that are marked as being water...\n");
    os.data.rp <- rasterToPoints(os.data.proj, fun=function(x){x>0}, spatial=TRUE); # extract the points that are marked as water
    rm(os.data, os.data.proj); gc(); # free up memory
    unlink(paste0("./input_files/", type, "-raster-reduced-3857/"), recursive=TRUE, force=TRUE); # remove the unpacked files...
  
    cat("Computing the minimum distance to water for each of the given points\n");
    d <- unlist(mclapply(1:nrow(coords), 
                         function(i) min(distGeo(coords[i,], os.data.rp), na.rm=TRUE), mc.cores=no.cores)); # compute the minimum distance from every point to water
    
    return (d/1000); # return the minimum distances in Km
  }

  # Compute the distances:
  m <- as.matrix(d_coords[,c("longitude_180","latitude")]); colnames(m) <- c("longitude","latitude"); rownames(m) <- d_coords$glottocode;
  tmp <- extract.dist2water.for.coords(coords=m, type="ocean", zoom=2); # oceans are big, lower zoom should not be a problem
  d_colors_water <- cbind(d_coords[,c("glottocode", "longitude", "latitude")], "dist2ocean"=tmp);
  tmp <- extract.dist2water.for.coords(coords=m, type="lakes", zoom=4); # but for lakes and rivers, higher zoom should capture smaller ones
  d_colors_water <- cbind(d_colors_water, "dist2lakes"=tmp);
  tmp <- extract.dist2water.for.coords(coords=m, type="river", zoom=4);
  d_colors_water <- cbind(d_colors_water, "dist2rivers"=tmp);
  d_colors_water <- cbind(d_colors_water, "dist2water"=pmin(d_colors_water$dist2ocean, d_colors_water$dist2lakes, d_colors_water$dist2rivers, na.rm=TRUE));
  
  # Plots for checking:
  if( FALSE )
  {
    library(ggplot2);
    library(ggrepel);
    mapWorld <- map_data("world", wrap=c(-20,340), ylim=c(-70,100));
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_colors_water, aes(x = longitude, y = latitude, color=dist2ocean)) +
      geom_label_repel(data=d_colors_water, aes(x = longitude, y = latitude, label=round(dist2ocean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_colors_water, aes(x = longitude, y = latitude, color=dist2lakes)) +
      geom_label_repel(data=d_colors_water, aes(x = longitude, y = latitude, label=round(dist2lakes,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_colors_water, aes(x = longitude, y = latitude, color=dist2rivers)) +
      geom_label_repel(data=d_colors_water, aes(x = longitude, y = latitude, label=round(dist2rivers,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_colors_water, aes(x = longitude, y = latitude, color=dist2water)) +
      geom_label_repel(data=d_colors_water, aes(x = longitude, y = latitude, label=round(dist2water,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
  }
  
  # Remove temporary files:
  unlink("./input_files/lakes-raster-reduced-3857", recursive=TRUE);
  unlink("./input_files/river-raster-reduced-3857", recursive=TRUE);
  unlink("./input_files/ocean-raster-reduced-3857", recursive=TRUE);
  
  # Save the data:
  write.table(d_colors_water, file="./input_files/data_dist2water.tsv", quote=FALSE, sep="\t", row.names=FALSE);
}



##
## UV incidence ####
##

# Extract & parse the UVB data from NASA Goddard Space Flight Center's "New UV Irradiance and Exposure Data Product"
# TOMS (Total Ozone Mapping Spectrometer): http://toms.gsfc.nasa.gov/ery_uv/
# available in raw format at: http://toms.gsfc.nasa.gov/ery_uv/new_uv/
#
# Original file format is described in http://toms.gsfc.nasa.gov/ery_uv/new_uv/1README.UV
#
# These are different wavelengths (305, 310, 320, and 380 nm) for UV as explain on the page http://toms.gsfc.nasa.gov/ery_uv/:
# "[...] Ultraviolet radiation is at shorter wavelengths than the visible spectrum (400 to 700 nm) and is divided into three 
# components: UV-A (315 to 400 nm), UV-B (280 to 315 nm) and UV-C (less than 280 nm). The shorter wavelengths that comprise
# UV-B are the most dangerous portion of UV radiation that can reach ground level [...]"
#
# Unfortunately, as of July 2020, this data is not longer available online (seems to have been repackaged as per
# https://disc.gsfc.nasa.gov/datacollection/TOMSEPL3mery_008.html), but we provide the originally processed data for reproducibility
# in the toms_nasa_uv folder.
if( !file.exists("./input_files/data_UV_incidence.tsv") )
{
  library(dplyr);
  
  # Load the data for the year 1998:
  d_uv_1998 <- rbind(read.table(xzfile("./input_files/toms_nasa_uv/UV-for-year-1998-freq-305.csv.xz"), header=TRUE, sep="\t"),
                     read.table(xzfile("./input_files/toms_nasa_uv/UV-for-year-1998-freq-310.csv.xz"), header=TRUE, sep="\t"),
                     read.table(xzfile("./input_files/toms_nasa_uv/UV-for-year-1998-freq-325.csv.xz"), header=TRUE, sep="\t"),
                     read.table(xzfile("./input_files/toms_nasa_uv/UV-for-year-1998-freq-380.csv.xz"), header=TRUE, sep="\t"));
  
  # For a given location, obtain the summaries:
  get_uv_data_for_location <- function(longitude, latitude)
  {
    # Find the closest datapoint to the location:
    longitudes <- as.numeric(substring(names(d_uv_1998)[6:ncol(d_uv_1998)], 5)) - 0.5;
    longitude_hit <- longitudes[ which.min(abs(longitudes - longitude)) ] + 0.5;
    
    latitudes <- sort(unique(as.numeric(d_uv_1998$latitude)));
    latitude_hit <- latitudes[ which.min(abs(latitudes - latitude)) ];
    
    # Extract the data:
    d <- d_uv_1998[ d_uv_1998$latitude == latitude_hit, c("year", "month", "day", "uv", "latitude", paste0("Deg.",longitude_hit)) ];
    names(d)[ncol(d)] <- "incidence";
    if( nrow(d) == 0 ) stop(paste0("Cannot find UV data for location (",longitude,", ",latitude,")."));
    
    # And return summaries:
    return (data.frame(# for each frequency band
                       "UV_305_mean"=mean(d$incidence[ d$uv == 305 ], na.rm=TRUE),
                       "UV_305_sd"  =sd(d$incidence[   d$uv == 305 ], na.rm=TRUE),
                       "UV_310_mean"=mean(d$incidence[ d$uv == 310 ], na.rm=TRUE),
                       "UV_310_sd"  =sd(d$incidence[   d$uv == 310 ], na.rm=TRUE),
                       "UV_325_mean"=mean(d$incidence[ d$uv == 325 ], na.rm=TRUE),
                       "UV_325_sd"  =sd(d$incidence[   d$uv == 325 ], na.rm=TRUE),
                       "UV_380_mean"=mean(d$incidence[ d$uv == 380 ], na.rm=TRUE),
                       "UV_380_sd"  =sd(d$incidence[   d$uv == 380 ], na.rm=TRUE),
                       # for UV-A
                       "UV_A_mean"  =mean(d$incidence[ d$uv >= 315 & d$uv < 400 ], na.rm=TRUE),
                       "UV_A_sd"    =sd(d$incidence[   d$uv >= 315 & d$uv < 400 ], na.rm=TRUE),
                       # for UV-B
                       "UV_B_mean"  =mean(d$incidence[ d$uv >= 280 & d$uv < 315 ], na.rm=TRUE),
                       "UV_B_sd"    =sd(d$incidence[   d$uv >= 280 & d$uv < 315 ], na.rm=TRUE),
                       # for all UVs
                       "UV_mean"    =mean(d$incidence, na.rm=TRUE),
                       "UV_sd"      =sd(d$incidence, na.rm=TRUE)));
  }
  
  # Get this data for all the populations:
  d_uv <- do.call(rbind, 
                  lapply(1:nrow(d_coords), function(i) cbind(d_coords[i, c("glottocode", "longitude", "latitude")], 
                                                             get_uv_data_for_location(d_coords$longitude[i], d_coords$latitude[i]))));
  
  # Plots for checking:
  if( FALSE )
  {
    library(ggplot2);
    library(ggrepel);
    mapWorld <- map_data("world", wrap=c(-20,340), ylim=c(-70,100));
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_305_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_305_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_305_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_305_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;

    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_310_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_310_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_310_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_310_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;

    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_325_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_325_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_325_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_325_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;

    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_380_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_380_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_380_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_380_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;

    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_A_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_A_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_A_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_A_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;

    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_B_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_B_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_B_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_B_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;

    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_mean)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_mean,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_uv, aes(x = longitude, y = latitude, color=UV_sd)) +
      geom_label_repel(data=d_uv, aes(x = longitude, y = latitude, label=round(UV_sd,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
  }
  
  # Save the data:
  write.table(d_uv, file="./input_files/data_UV_incidence.tsv", quote=FALSE, sep="\t", row.names=FALSE);
}


##
## Altitude ####
##

if( !file.exists("./input_files/data_elevation.tsv") )
{
  # Use the elevatr package to get the Mapzen elevation data, currently (July 2020) still available on https://registry.opendata.aws/terrain-tiles/
  # see https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html#get_raster_elevation_data
  library(raster);
  library(elevatr);

  # Get the elevation data:
  elevation_tiles <- list();
  elevation <- rep(NA, nrow(d_coords));
  for( i in 1:nrow(d_coords) )
  {
    cat(paste0("Obtaining elevation for ",d_coords$glottocode[i]," (",i," of ",nrow(d_coords),"):\n"));
    # Try several levels of zoom, until one succeeds:
    for( z in 1:14) 
    {
      cat(paste0("... trying zoom level ",z,"...\n"));
      e <- NULL;
      try(e <- get_elev_raster(d_coords[i,c("longitude_180", "latitude"), drop=FALSE], 
                               prj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", z=z, expand=ifelse(z < 5, 10, ifelse(z < 10, 2, 1))), # adapt expansion to the zoom level
          silent=TRUE);
      if( !is.null(e) && inherits(e, "RasterLayer") )
      {
        # Succeeded:
        elevation_tiles[[i]] <- list("z"=z, "raster"=e);
        elevation[i] <- extract(e, d_coords[i,c("longitude_180", "latitude"), drop=FALSE], method="simple");
        #plot(e, main=d_coords$glottocode[i]); points(d_coords$longitude_180[i], d_coords$latitude[i]);
        break;
      }
    }
    if( is.null(e) ){ warning("Error obtaining elevation info..."); elevation_tiles[[i]] <- list(); elevation[i] <- NA; }
  }
  names(elevation_tiles) <- names(elevation) <- d_coords$glottocode;
  d_elevation <- cbind(d_coords[,c("glottocode", "longitude", "latitude")], "elevation"=round(elevation,1));
  
  # Plots for checking:
  if( FALSE )
  {
    library(ggplot2);
    library(ggrepel);
    mapWorld <- map_data("world", wrap=c(-20,340), ylim=c(-70,100));
    
    ggplot() + theme_bw() +
      theme_bw() +
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group) ,fill = "grey") +
      geom_point(data=d_elevation, aes(x = longitude, y = latitude, color=elevation)) +
      geom_label_repel(data=d_elevation, aes(x = longitude, y = latitude, label=round(elevation,0)), alpha=0.5, fill="lightyellow") + 
      theme(legend.position = c(0.75, 0.5), 
            legend.justification = c(1, 1), 
            legend.title = element_text(size = 9), 
            legend.text = element_text(size = 10)) +
      NULL;
  }  
  
  # Save the tiles for later use as an xz-compressed tar archive:
  dir.create("./input_files/elevation_Mapzen_tiles", showWarnings=FALSE);
  for(i in seq_along(elevation_tiles) )
  {
    cat(paste0("Saving tile for ",names(elevation_tiles)[i]," (",i," of ",length(elevation_tiles),"):\n"));
    writeRaster(floor(elevation_tiles[[i]]$raster), # the floating-point precision is not useful and wastes disk space
                paste0("./input_files/elevation_Mapzen_tiles/raster_for_",names(elevation_tiles)[i],"_z_",elevation_tiles[[i]]$z,".tif"), 
                format="GTiff", options=c("COMPRESS=LZW", "PREDICTOR=2"), overwrite=TRUE); # compress it as per https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
  }
  tar("./input_files/elevation_Mapzen_tiles.txz", "./input_files/elevation_Mapzen_tiles/", compression="xz", compression_level=9);
  unlink("./input_files/elevation_Mapzen_tiles/", recursive=TRUE);

  # Save the data:
  write.table(d_elevation, file="./input_files/data_elevation.tsv", quote=FALSE, sep="\t", row.names=FALSE);
}


##
## Genetic distances ####
##

if( !file.exists("./input_files/data_genetics.tsv") )
{
  # Load the genetic distances with ultrametric imputation:
  d_gen_ult <- read.table("./input_files/distmat_gen_ult.csv", header=TRUE, sep=","); rownames(d_gen_ult) <- d_gen_ult[,1]; d_gen_ult <- d_gen_ult[,-1]; d_gen_ult <- as.matrix(d_gen_ult);

  # Let's pick the first 10 MDS dimensions:
  k <- 10; x_ult <- cmdscale(d_gen_ult, k=k, eig=TRUE, add=TRUE);
  
  d_genetics <- merge(d_colors[,c("glottocode"), drop=FALSE], data.frame("glottocode"=rownames(x_ult$points), x_ult$points), by="glottocode", all.x=TRUE, all.y=FALSE);
  names(d_genetics)[ (ncol(d_genetics)-k+1):ncol(d_genetics) ] <- paste0("gen_D",1:k);

  # Save the data:
  write.table(d_genetics, file="./input_files/data_genetics.tsv", quote=FALSE, sep="\t", row.names=FALSE);
}





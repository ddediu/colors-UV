# Compute the geographic distances between populations

# Load the language data:
d_colors <- read.table("./input_files/data_colors.csv", header=TRUE, sep=",", quote='"', stringsAsFactors=FALSE); # comma-separated double-quoted CVS file
# Alternative view of the globe: all longitudes > 180 are "flipped" to negative values:
d_colors$longitude_180 <- ifelse(d_colors$longitude <= 180, d_colors$longitude, d_colors$longitude - 360);

# Extract the geographic coordinates:
# For entries that have precisely the same geographic coordinates as other entries (if any), add 0.01 degrees (on the order of hundreds of meters/the scale of towns or villages):
same_coords <- duplicated(d_colors[, c("latitude", "longitude_180")]); 
d_colors$latitude[ same_coords ] <- d_colors$latitude[ same_coords ] + 0.01; d_colors$longitude_180[ same_coords ] <- d_colors$longitude_180[ same_coords ] + 0.01; 

d_coords <- as.matrix(d_colors[,c("longitude_180", "latitude")]); rownames(d_coords) <- as.character(d_colors$glottocode);
d_all_pairs <- expand.grid("p1"=1:nrow(d_coords), "p2"=1:nrow(d_coords)); d_all_pairs <- d_all_pairs[ d_all_pairs$p1 < d_all_pairs$p2, ];


##
## Climate data ####
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
  bio_data <- extract(bio_raster, as.matrix(d_colors[,c("longitude_180","latitude")])); # don't forget that latitude is y and longitude is x!
  rownames(bio_data) <- as.character(d_colors$glottocode);

  # PCA:
  library(fpc);
  library(factoextra);
  PCs <- prcomp( ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19, 
                 data=as.data.frame(bio_data), scale=TRUE, center=TRUE);
  summary(PCs); # PC1 (49.4%), PC2 (25.5%), PC3 (8.3%)
  
  # Save it:
  climate_data <- data.frame("glottocode"=rownames(PCs$x), PCs$x[,1:3]);
  names(climate_data) <- c("glottocode", paste0("clim_PC",1:(ncol(climate_data)-1)));
  write.table(climate_data, file="./input_files/data_climate.tsv", quote=FALSE, sep="\t", row.names=FALSE);
  
  # Clean up the downloaded files:
  unlink("./input_files/wc5/*.bil", recursive=TRUE);
  unlink("./input_files/wc5/*.hdr", recursive=TRUE);
}


# Humidity data from the NOAA:
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
      stop("Cannot specific retrive humidity data from the NOAA (http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCEP-NCAR/.CDAS-1/.MONTHLY/.Diagnostic/.above_ground/.qa/datafiles.html)!\n");
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
  humidity_data <- do.call(rbind, lapply(1:nrow(d_colors), function(i) humidity_summary(d_colors$latitude[i], d_colors$longitude_180[i])));
  
  humidity_data <- cbind(d_colors[,c("glottocode", "latitude", "longitude")], humidity_data);
  names(humidity_data) <- c("glottocode", "latitude", "longitude", 
                            "humidity_overall_mean", "humidity_overall_median", "humidity_overall_sd", "humidity_overall_IQR", 
                            "humidity_mean_annual_mean", "humidity_mean_annual_median", "humidity_mean_annual_sd", "humidity_mean_annual_IQR");
  
  # Plot:
  #ggplot(humidity_data, aes(x=longitude, y=latitude, color=humidity_mean_annual_IQR)) + geom_point()

  # Save the data:
  write.table(humidity_data, file="./input_files/data_humidity.tsv", quote=FALSE, sep="\t", row.names=FALSE);
}






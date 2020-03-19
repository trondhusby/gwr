# Use rocker/geospatial as the base image
# This will inherit everyhing that was installed into base image already
# Documented at https://hub.docker.com/r/rocker/geospatial/~/dockerfile/
FROM rocker/geospatial

# Install any dependencies required for the R packages
#RUN  apt-get update \
#  && apt-get install -y --no-install-recommends \
#  libxml2-dev \
#  libcurl4-openssl-dev \
#  libssl-dev

# Install the R Packages from CRAN
RUN Rscript -e 'install.packages(c("MASS", "GWmodel", "data.table", "rgenoud"))'

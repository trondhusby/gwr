#!/bin/bash

# create docker image
sudo docker build docker_gwr -t docker_gwr/custom-geospatial

# rename tag to include user name on docker hub
sudo docker image tag docker_gwr/custom-geospatial:latest husbyt/custom-geospatial:latest

# push to docker hub
sudo docker push husbyt/custom-geospatial

# end of script


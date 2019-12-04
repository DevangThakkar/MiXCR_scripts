################################################################################
##################### Set Inital Image to work from ############################
FROM ubuntu:16.04

################################################################################
##################### Add Container Labels #####################################
LABEL "Description"="Software toolkit for T/B cell repertoire"

################################################################################
############## set environment variables #######################################

ENV version 3.0.3

################################################################################
##################### Install System Dependencies ##############################
RUN apt-get update -y && apt-get install -y \
    git \
    build-essential \
    openjdk-8-jre \
    wget \
    unzip

################################################################################
##################### Install MIXCR    #########################################

# set inital WORKDIR
WORKDIR /mixcr/

# download the source
RUN wget https://github.com/milaboratory/mixcr/releases/download/v${version}/mixcr-${version}.zip

# unzip the source
RUN unzip mixcr-${version}.zip

# download BCR script
RUN wget https://raw.githubusercontent.com/DevangThakkar/MiXCR_scripts/master/BCR.sh

# download TCR script
RUN wget https://raw.githubusercontent.com/DevangThakkar/MiXCR_scripts/master/TCR.sh

################################################################################
###################### set environment path    #################################

# add mixcr executable to path
ENV PATH="/mixcr/mixcr-${version}/:${PATH}"
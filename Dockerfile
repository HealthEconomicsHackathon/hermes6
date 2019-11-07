## Start with the tidyverse docker image
FROM rocker/tidyverse:latest


## Add latex and magick (optional)
RUN apt-get update -y && \
apt-get install -y \
texlive-latex-recommended \
texlive-fonts-extra \
texinfo \
libqpdf-dev \
libmagick++-dev \
&& apt-get clean

## Add everything to a folder
ADD . /home/rstudio/hermes6

## Install dev deps from description
RUN Rscript -e 'devtools::install_dev_deps("home/rstudio/hermes6")'

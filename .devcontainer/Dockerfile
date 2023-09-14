# Pre-built Dev Container Image for R. More info: https://github.com/rocker-org/devcontainer-images/pkgs/container/devcontainer%2Ftidyverse
# Available R version: 4, 4.1, 4.0
ARG VARIANT="4.2"
FROM ghcr.io/rocker-org/devcontainer/tidyverse:${VARIANT}

RUN install2.r --error --skipinstalled -n -1 \
        statip \
        patchwork \
        paletteer \
        here \
        doParallel \
        janitor \
        vip \
        ranger \
        palmerpenguins \
        skimr \
        nnet \
        kernlab \
        plotly \
        factoextra \
        cluster \
        tidymodels \
        markdown \
        ottr \
    && rm -rf /tmp/downloaded_packages \
    && R -q -e 'remotes::install_github("https://github.com/dcomtois/summarytools/tree/0-8-9")'

# Install Python packages
COPY requirements.txt /tmp/pip-tmp/
RUN python3 -m pip --disable-pip-version-check --no-cache-dir install -r /tmp/pip-tmp/requirements.txt

RUN apt-get update
RUN apt-get install bzip2 libbz2-dev
RUN python3 -m pip install --upgrade pip setuptools wheel
RUN python3 -m pip install rpy2

# COPY phyloModels-clorton-relax-qt5-version.zip /tmp/phyloModels-clorton-relax-qt5-version.zip
# WORKDIR "/tmp"
# RUN unzip phyloModels-clorton-relax-qt5-version.zip
# WORKDIR "/tmp/phyloModels-clorton-relax-qt5-version"
# RUN python3 -m pip install pyqt5
# RUN python3 -m pip install .

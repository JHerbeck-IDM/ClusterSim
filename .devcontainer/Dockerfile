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

COPY rpy2-RELEASE_3_5_13.tar.gz /tmp/rpy2/
WORKDIR "/tmp/rpy2"
RUN tar -xzf rpy2-RELEASE_3_5_13.tar.gz
WORKDIR "/tmp/rpy2/rpy2-RELEASE_3_5_13"
RUN python3 setup.py build install
COPY phyloModels-master.zip /tmp/phyloModels
WORKDIR "/tmp/phyloModels"
RUN unzip phyloModels-master.zip
WORKDIR "/tmp/phyloModels/phyloModels-master"
RUN python3 -m pip install .

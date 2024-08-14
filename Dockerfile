

wget -qO - https://gitlab.com/paulcarroty/vscodium-deb-rpm-repo/raw/master/pub.gpg \
    | gpg --dearmor \
    | sudo dd of=/usr/share/keyrings/vscodium-archive-keyring.gpg
# Fetch base image for R version 4.4.1
FROM rocker/rstudio:4.4.1

# Copy repository into the container image
COPY --chown=rstudio:rstudio . /opt/pkg

# Install Java
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-21-jdk

# Let R automatically find Java files
RUN R CMD javareconf

# Install devtools
RUN Rscript -e 'install.packages("devtools")'

# Install IntegratedLearner
RUN R -e 'devtools::install(pkg = "/opt/pkg", dependencies = TRUE, upgrade = "always")'

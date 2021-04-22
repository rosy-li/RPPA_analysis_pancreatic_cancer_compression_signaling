# hash:sha256:d28d4099e750d8f668ba11a37f101be780c224efc61cd87579e6c3f02d30ee81
FROM registry.codeocean.com/codeocean/r-studio:1.2.5019-r4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("Hmisc", "4.5-0")' \
    && Rscript -e 'remotes::install_version("ggrepel", "0.9.1")' \
    && Rscript -e 'remotes::install_version("gplots", "3.1.1")' \
    && Rscript -e 'remotes::install_version("hrbrthemes", "0.8.0")' \
    && Rscript -e 'remotes::install_version("plyr", "1.8.6")' \
    && Rscript -e 'remotes::install_version("rafalib", "1.0.0")'

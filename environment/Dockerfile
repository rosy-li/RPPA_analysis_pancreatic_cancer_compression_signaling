# hash:sha256:c969b56d364a01053ddb13a482a3aa07821c59c7f25e2fc73c2259a2bb6d700a
FROM registry.codeocean.com/codeocean/r-studio:1.2.5019-r4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("Hmisc", "4.5-0")' \
    && Rscript -e 'remotes::install_version("ggrepel")' \
    && Rscript -e 'remotes::install_version("gplots", "3.1.1")' \
    && Rscript -e 'remotes::install_version("hrbrthemes", "0.8.0")' \
    && Rscript -e 'remotes::install_version("rafalib", "1.0.0")'

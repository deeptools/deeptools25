FROM mcr.microsoft.com/devcontainers/anaconda:1-3

# System-wide dependencies
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
    libclang-dev \
    liblzma-dev

# Conda package-management
RUN /opt/conda/bin/conda create -n snkmk -c conda-forge -c bioconda \
    --strict-channel-priority --override-channels \
    bedtools \
    matplotlib \
    snakemake \
    ucsc-bigwiginfo

# Ad-hoc solutions
ARG USERNAME
RUN useradd -m -s /bin/bash ${USERNAME} \
    && chown -R ${USERNAME}:${USERNAME} /opt/conda/pkgs
USER ${USERNAME}

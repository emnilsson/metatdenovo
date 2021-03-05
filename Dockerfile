FROM nfcore/base:1.12.1
LABEL authors="daniel.lundin@lnu.se" \
      description="Docker image containing all software requirements for the nf-core/metatdenovo pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-metatdenovo-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-metatdenovo-1.0dev > nf-core-metatdenovo-1.0dev.yml

# Install eggnog-mapper via pip, which requires gcc
RUN apt-get update
RUN apt-get -y install gcc
RUN pip install eggnog-mapper==2.0.8.post2

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

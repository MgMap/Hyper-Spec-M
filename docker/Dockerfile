FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

# Update and install mini-conda
RUN apt-get update && apt-get install -y time wget g++ libboost-all-dev git && \
        apt-get clean

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/anaconda.sh && \
        /bin/bash ~/anaconda.sh -b -p /opt/conda && \
        rm ~/anaconda.sh && \
        ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
        echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        find /opt/conda/ -follow -type f -name '*.a' -delete && \
        find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
        /opt/conda/bin/conda clean -afy

# set path to conda
ENV PATH /opt/conda/bin:$PATH

# setup conda virtual environment
COPY ./docker/requirements.yaml /tmp/requirements.yaml
RUN conda update conda \
        && conda env create --name hyper-spec -f /tmp/requirements.yaml

RUN echo "conda activate hyper-spec" >> ~/.bashrc
ENV PATH /opt/conda/envs/hyper-spec/bin:$PATH
ENV CONDA_DEFAULT_ENV $hyper-spec

# Install Nsight Compute
RUN wget https://developer.nvidia.com/downloads/assets/tools/secure/nsight-compute/2024_3_1/nsight-compute-linux-2024.3.1.2-34702747.run && \
    chmod +x nsight-compute-linux-2024.3.1.2-34702747.run && \
    ./nsight-compute-linux-2024.3.1.2-34702747.run && \
    rm nsight-compute-linux-2024.3.1.2-34702747.run

# Add Nsight Compute to PATH
ENV PATH=$PATH:/usr/local/NVIDIA-Nsight-Compute

# Clone the code
RUN git clone https://github.com/MgMap/Hyper-Spec-M.git

# set working dir
WORKDIR /Hyper-Spec-M

# Install and build  
RUN bash ./install.sh 

# This container installs pyVCF and the code in phyloflow/vcf_transform/code
# Conda is used to manage the python env(s) inside the container:
# Use 'conda run ..' with 'vcf-transform'


FROM ubuntu:20.04

ENV HOME=/home/dockeruser
RUN mkdir /home/dockeruser && \
    touch /home/dockeruser/.bashrc && \
    chmod -R 777 /home/dockeruser 

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get -qq update && \
    apt-get -y install --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh  -O ~/miniconda.sh

RUN /bin/bash ~/miniconda.sh -b -p /opt/conda
RUN rm ~/miniconda.sh
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/dockeruser/.bashrc
RUN echo "conda activate base" >> /home/dockeruser/.bashrc

RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda create -n moss
RUN conda install moss pandas pysam scikit-allel pyyaml

#try to reduce the size of the image
RUN conda clean --all --yes

#system wide conda status
RUN conda info

#this command prefix will run what comes in the suffix
#in a bash shell with all the paths and dependencies installed
#in the 'base' conda env. Single quotes will be
#needed around the command suffix
ENV RUN_IN_CONDA="bash -i -c "


#this will print all installed packages and their versions
RUN $RUN_IN_CONDA 'conda list'

COPY ./scripts /moss_scripts

RUN ls -al /
RUN ls -al /moss_scripts/
RUN moss -h

# ENTRYPOINT [ "python", "/moss_scripts/run_moss.py" ]

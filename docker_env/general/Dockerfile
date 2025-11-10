FROM oprgroup/nextflow_conda:1.0.0

WORKDIR /app
COPY centrifuge.yaml .
COPY datasets.yaml .

RUN conda env create -f centrifuge.yaml
RUN conda env create -f datasets.yaml


# Make subsequent RUNs use bash -lc
SHELL ["/bin/bash", "-lc"]

# Initialize conda for bash
RUN conda init bash
RUN echo 'source /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc

# Optional: also set PATH so non-login invocations still see the env
ENV PATH=/opt/conda/envs/env/bin:$PATH

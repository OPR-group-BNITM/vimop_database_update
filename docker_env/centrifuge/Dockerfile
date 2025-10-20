FROM continuumio/miniconda3

WORKDIR /app
COPY env.yaml .

RUN conda env create -f env.yaml

# Make subsequent RUNs use bash -lc
SHELL ["/bin/bash", "-lc"]

# Initialize conda for bash and auto-activate the env in login shells
RUN conda init bash
RUN echo 'source /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc \
 && echo 'conda activate env' >> ~/.bashrc

# Optional: also set PATH so non-login invocations still see the env
ENV PATH=/opt/conda/envs/env/bin:$PATH

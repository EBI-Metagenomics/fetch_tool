FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Microbiome Informatics"
LABEL version="1.0.0"
LABEL description="EBI Fetch Tool."

COPY --chown=$MAMBA_USER:$MAMBA_USER conda_environment.yml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

WORKDIR /opt

COPY . .

ENV PATH="/opt/fetchtool:$PATH"
ENV PYTHONPATH="/opt/:$PYTHONPATH"

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

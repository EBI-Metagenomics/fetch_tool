FROM mambaorg/micromamba:1.5.8

LABEL maintainer="Microbiome Informatics"
LABEL version="1.0.4"
LABEL description="EBI Fetch Tool."

COPY --chown=$MAMBA_USER:$MAMBA_USER conda_environment.yml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/fetch-tool-src

WORKDIR /opt/fetch-tool-src

RUN pip install . --no-cache-dir

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

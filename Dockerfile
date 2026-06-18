FROM continuumio/miniconda3:24.1.2-0

LABEL version="1.0.0"
LABEL description="FUNGUS-SV: Triangulation-based SV discovery for haploid fungi"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl git procps \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/fungus-sv
WORKDIR /opt/fungus-sv

RUN conda env create -f workflow/envs/alignment.yaml -n sv_align && conda clean -afy
RUN conda env create -f workflow/envs/kmers.yaml -n sv_kmers && conda clean -afy
RUN conda env create -f workflow/envs/lar.yaml -n sv_lar && conda clean -afy
RUN conda env create -f workflow/envs/validation.yaml -n sv_valid && conda clean -afy
RUN conda env create -f workflow/envs/sv_calling.yaml -n sv_call && conda clean -afy

ENV PYTHONPATH=/opt/fungus-sv

RUN conda run -n sv_align minimap2 --version && \
    conda run -n sv_align samtools --version && \
    conda run -n sv_call which sniffles && \
    conda run -n sv_call which cuteSV && \
    conda run -n sv_call which svim && \
    conda run -n sv_valid python -c "import numpy,scipy,pandas,pysam,yaml; print('sv_valid OK')" && \
    conda run -n sv_lar flye --version && \
    conda run -n sv_lar which miniasm && \
    conda run -n sv_lar which racon

HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
  CMD conda run -n sv_align minimap2 --version > /dev/null && \
      conda run -n sv_call which sniffles > /dev/null && \
      conda run -n sv_valid python -c "import pysam" && \
      echo "Healthy" || exit 1

CMD ["/bin/bash"]

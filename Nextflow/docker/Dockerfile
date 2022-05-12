FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-compomics"

COPY environment.yml /
RUN conda install mamba -n base -c conda-forge
RUN mamba env create -f /environment.yml && conda clean -a;
# conda package is outdated, therefore direct install
#RUN /opt/conda/envs/nf-core-compomics-workflow/bin/R -e 'install.packages("MSqRob", repos = "http://cran.us.r-project.org")'

ENV PATH /opt/conda/envs/nf-core-compomics-workflow/bin:$PATH

RUN chmod -R a+rwx /opt/conda/envs/nf-core-compomics-workflow/share/searchgui-*/resources
RUN chmod -R a+rwx /opt/conda/envs/nf-core-compomics-workflow/share/peptide-shaker-*/resources

USER root 

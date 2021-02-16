FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for nf-core-proline-labelfree"


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a;


#SHELL ["conda", "run", "-n", "moFF", "/bin/bash", "-c"]


ENV PATH /opt/conda/envs/nf-core-moff-labelfree/bin:$PATH


RUN chmod -R a+rwx /opt/conda/envs/nf-core-moff-labelfree/share/searchgui-*/resources
RUN chmod -R a+rwx /opt/conda/envs/nf-core-moff-labelfree/share/peptide-shaker-*/resources


#ENTRYPOINT ["conda", "run", "-n", "moFF", "/bin/bash", "-c"]

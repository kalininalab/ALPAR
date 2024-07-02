FROM continuumio/miniconda3

ENV PATH /opt/conda/bin:$PATH

COPY environment.yml /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml

RUN echo "source activate alpar" > ~/.bashrc

RUN pip install panacota

ENTRYPOINT ["bash", "-c", "source activate alpar && exec bash"]

ENTRYPOINT ["alpar"]

CMD ["bash"]
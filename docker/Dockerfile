FROM bioconductor/bioconductor_docker:RELEASE_3_16

RUN mkdir /r_user_lib/ && chown -c rstudio /r_user_lib/
ENV R_LIBS_USER=/r_user_lib/

COPY rserver.conf /etc/rstudio/rserver.conf
COPY rprofile.R /home/rstudio/.Rprofile
COPY setup.R /r_user_lib/setup.R
RUN chown -c rstudio /r_user_lib/setup.R && chown -c rstudio /home/rstudio/.Rprofile

USER rstudio
RUN Rscript --vanilla /r_user_lib/setup.R

USER root
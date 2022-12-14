FROM opencpu/ubuntu-20.04:v2.2.8-1


RUN R -e 'install.packages("renv")'

# fix clang-10 span header (needed for iLincsCor)
RUN apt-get update && apt-get install -y clang-10
RUN add-apt-repository --update -y ppa:ubuntu-toolchain-r/test
RUN apt-get update -y
RUN apt-get -y --fix-broken install gcc-11 g++-11

RUN mkdir ~/.R
RUN echo "CXX20=clang++-10" > ~/.R/Makevars

# preload iLincsCorSrv (which will load signature libraries to iLincsCor class)
RUN sed -i 's/"]/","iLincsCorSrv"]/g' /etc/opencpu/server.conf

COPY . /tmp/iLincsCor

RUN cd /tmp && R -e 'renv::install("/tmp/iLincsCor", library="/usr/lib/opencpu/library")'

# FIXME
RUN cd /tmp && R -e 'renv::install("/tmp/iLincsCor/iLincsCorSrv", library="/usr/lib/opencpu/library")'


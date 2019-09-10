#=========================#
#  Operating System setup #
#=========================#
from centos:centos7

# install basic build environment
RUN yum -y install epel-release && yum -y update && \
    yum -y install gcc g++ gcc-gfortran zlib-devel libxml2-devel && \
    yum -y install libcurl-devel make netcdf-devel libpng-devel && \
    yum -y install python libconfig-devel bc udunits2-devel && \
    yum -y install gsl-devel lapack-devel hdf5-devel compat-libf2c-34 && \
    rm -rf /var/lib/apt/lists/*


#install hydroclass binaries in /opt/hydroclass
COPY hca.tar /opt
RUN mkdir /opt/hydroclass && \
    tar -xf /opt/hca.tar -C /opt/hydroclass && \
    rm -rf /opt/hca.tar
COPY Radx2Grid /opt/hydroclass
COPY RadxConvert /opt/hydroclass
COPY hc_composite /opt/hydroclass
ENV PATH=/opt/hydroclass:$PATH

#install necessary libraries for hydroclass 
COPY hc_lib.tar /opt
RUN mkdir /opt/hc_lib && \
    tar -xf /opt/hc_lib.tar -C /opt/hc_lib && \
    rm -rf /opt/hc_lib.tar
ENV LD_LIBRARY_PATH=/opt/hc_lib:$LD_LIBRARY_PATH

#install netcdf2png
COPY hydroclass_netcdf2png.tar /opt
RUN tar xf /opt/hydroclass_netcdf2png.tar -C /opt && \
    cd /opt/hydroclass_netcdf2png && make && \
    cd /opt && rm /opt/hydroclass_netcdf2png.tar
ENV PATH=/opt/hydroclass_netcdf2png:$PATH

# install libconfig in /opt/libconfig
COPY libconfig-1.5.tar.gz /opt
RUN mkdir /opt/libconfig && \
    tar xzf /opt/libconfig-1.5.tar.gz -C /opt/libconfig && \
    cd /opt/libconfig/libconfig-1.5 && \
    ./configure --prefix=/usr/local/libconfig && \
    make && make install && \
    cd /opt && \
    rm -rf /opt/libconfig && \
    rm /opt/libconfig-1.5.tar.gz
ENV LD_LIBRARY_PATH=/usr/local/libconfig/lib:$LD_LIBRARY_PATH

# install jansson in /opt/jansson
COPY jansson-2.7.tar.gz /opt
RUN mkdir /opt/jansson && \
    tar xzf /opt/jansson-2.7.tar.gz -C /opt/jansson && \
    cd /opt/jansson/jansson-2.7 && \
    ./configure --prefix=/usr/local/jansson && \
    make && make install && \
    cd /opt && \
    rm -rf /opt/jansson && \
    rm /opt/jansson-2.7.tar.gz
ENV LD_LIBRARY_PATH=/usr/local/jansson/lib:$LD_LIBRARY_PATH
			    
#install hmt
COPY hmt.tar /opt
RUN tar xf /opt/hmt.tar -C /opt && \
    cd /opt/hmt && make && \
    cd /opt && rm /opt/hmt.tar
ENV PATH=/opt/hmt:$PATH


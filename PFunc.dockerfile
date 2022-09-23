FROM ubuntu:16.04
RUN apt update \
    && apt install -y r-base python3 python3-pip python3-tk \
    && rm -rf /var/lib/apt/lists/*
RUN R -e "packageurl <- 'http://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-27.tar.gz'; install.packages(packageurl, repos=NULL, type='source')"
RUN pip3 install numpy==1.14 matplotlib==2.0.1 rpy2==2.8.5 six==1.11.0 pytz==2020.1 cycler==0.10.0 pyparsing==2.4.7
RUN mkdir -p /root/PFunc/data
COPY COPYING.txt demo_data* PFunc* README.md /root/PFunc/
CMD cd /root/PFunc/ && python3 /root/PFunc/PFunc.py

# build command: sudo docker build -f PFunc.dockerfile . -t joccalor/pfunc

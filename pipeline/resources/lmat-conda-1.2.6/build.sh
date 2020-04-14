#!/bin/bash

ARCH="$(uname 2> /dev/null)"

LinuxInstallation () {
	mkdir -vp ${PREFIX}/build;
	mkdir -vp ${PREFIX}/bin;
	tar xvf LMAT-1.2.6.tgz;
    cd LMAT-1.2.6;
    # probably sketchy but hey it means it builds...
    sed -i "s/CXXFLAGS = /CXXFLAGS = -fpermissive /" Makefile.inc;
    # otherwise the conda ldflags mess up linking by using modern g++ stds
    sed -i "s/#LDFLAGS = -static/LDFLAGS = /" Makefile.inc;
    make || return 1;

	# patch urls 
	wget https://gist.githubusercontent.com/fmaguire/f81cb90b1e7690d09c8ec22974c8e01a/raw/d1d7c98c5595d682fb28bd3b6d15e08115092c23/get_db.patch && \
    patch bin/get_db.sh get_db.patch

    sed -i 's/\/usr\/bin\/python/\/usr\/bin\/env python/' bin/*.py
	chmod +x bin/*
	cp bin/* ${PREFIX}/bin
	return 0;
}

## patch get_db.sh
#RUN # configure directory structure
#RUN mkdir -p /runtime_inputs /pipeline /data
#
#WORKDIR /data/
#ENV LMAT_DIR /runtime_inputs

case ${ARCH} in
    'Linux')
        LinuxInstallation || exit 1;
	;;
    *)
        echo -e "Unsupported machine type: ${ARCH}";
        exit 1;
        ;;
esac

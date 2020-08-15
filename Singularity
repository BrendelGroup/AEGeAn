bootstrap: docker
From: fedora:32

%help
    This container provides portable & reproducible components for AEGeAn:
    ... from the Brendel Group.
    Please see https://github.com/BrendelGroup/AEGeAn for complete documentation.
    
%post
    dnf -y update
    dnf -y install bc bzip2 findutils git lftp mlocate tcsh unzip zip wget which
    dnf -y install gcc-c++ make ruby
    dnf -y install cairo-devel pango-devel zlib-devel
    dnf -y install libnsl
    dnf -y install python3-pycurl python3-pyyaml python3-pandas
    dnf -y install python3-entrypoints python3-pytest python3-pytest-cov
    dnf -y install python3-biopython
    dnf -y install pandoc
    dnf -y install parallel

    #NOTE: "python" is refered to in AEGeAn code, but may not be defined for
    #      the Fedora 32 environment as set up above. The following will work:
    ln -sf /usr/bin/python3 /usr/bin/python

    cd /usr/local/src

    echo 'Installing the GenomeTools package:'
    git clone https://github.com/genometools/genometools.git
    cd genometools
    make
    make install
    make clean
    sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/genometools-x86_64.conf'
    ldconfig
    cd ..

    echo 'Installing the cd-hit package:'
    git clone https://github.com/weizhongli/cdhit.git
    cd cdhit
    make
    make install
    cd ..

    echo 'Installing the lastz package:'
    mkdir LASTZ
    cd LASTZ
    wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
    tar -xzf lastz-1.04.03.tar.gz
    cd lastz-distrib-1.04.03/
    sed -i -e "1i LASTZ_INSTALL=/usr/local/bin" make-include.mak
    sed -i -e "s/-Wall//;" src/Makefile
    make
    make install
    cd ../..

    git clone https://github.com/BrendelGroup/AEGeAn.git
    cd AEGeAn
    make all LocusPocus
    make install install-scripts
    sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/aegean-x86_64.conf'
    ldconfig
    cd ..

    
    echo 'Installing BLAST+ version 2.10.1 from NCBI'
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.10.1+-x64-linux.tar.gz
    cd ncbi-blast-2.10.1+/bin
    cp * /usr/local/bin/
    cd ../..
    rm ncbi-blast-2.10.1+-x64-linux.tar.gz
    cd ..

    echo 'Installing MuSeqBox version 5.5 from BrendelGroup'
    wget http://www.brendelgroup.org/bioinformatics2go/Download/MuSeqBox-3-4-2020.tar.gz
    tar -xzf MuSeqBox-3-4-2020.tar.gz
    cd MUSEQBOX5.5/src/
    make linux
    make install
    make clean
    cd ../..
    
%environment
    export LC_ALL=C


%labels
    Maintainer vpbrendel
    Version v1.1.0

#!/usr/bin/env bash

mkdir -p $HOME/local/src

wget -P $HOME/local/src https://github.com/genometools/genometools/archive/v1.6.1.tar.gz
tar -xzf $HOME/local/src/v1.6.1.tar.gz --directory $HOME/local/src/
cd $HOME/local/src/genometools-1.6.1
make -j 2 prefix=$HOME/local install
cd -

make prefix=$HOME/local install install-scripts

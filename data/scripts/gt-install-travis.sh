#!/usr/bin/env bash

wget http://genometools.org/pub/genometools-1.5.3.tar.gz
tar xzf genometools-1.5.3.tar.gz
cd genometools-1.5.3
sudo make 64bit=yes install


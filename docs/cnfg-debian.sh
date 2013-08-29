# Procedure for installing prerequisites and configuring a system running a
# Debian based distribution for GenomeTools and AEGeAn installation; tested on
# Ubuntu 11.10.
echo $PATH | grep /usr/local/bin > /dev/null
if [ $? != 0 ]; then
  export PATH=/usr/local/bin:$PATH
  echo 'export PATH=/usr/local/bin:$PATH' >> /etc/bashrc
fi
grep '/usr/local/lib' /etc/ld.so.conf /etc/ld.so.conf.d/* > /dev/null
if [ $? != 0 ]; then
  echo '/usr/local/lib' >> /etc/ld.so.conf.d/genometools-x86_64.conf
  ldconfig
fi
apt-get install -y build-essential git libcairo2-dev libncurses5-dev

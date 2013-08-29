# Procedure for installing prerequisites and configuring a system running a Red
# Hat based distribution for GenomeTools and AEGeAn installation; tested on
# CentOS 5.3
echo $PATH | grep /usr/local/bin > /dev/null
if [ $? != 0 ]; then
  export PATH=/usr/local/bin:$PATH
  echo 'export PATH=/usr/local/bin:$PATH' >> /etc/bashrc
fi
grep '/usr/local/lib' /etc/ld.so.conf /etc/ld.so.conf.d/* > /dev/null
if [ $? != 0 ]; then
  echo '/usr/local/lib' >> /etc/ld.so.conf.d/genometools-x86_64.conf
  /sbin/ldconfig
fi
yum install -y git cairo-devel ncurses-devel

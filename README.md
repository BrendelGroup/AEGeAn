## [AEGeAn: analysis and evaluation of genome annotations](http://standage.github.io/AEGeAn/)
Copyright &copy; 2010-2013, Daniel S. Standage and CONTRIBUTORS.  
See [docs/LICENSE](https://github.com/standage/AEGeAn/blob/master/docs/LICENSE) and [docs/CONTRIBUTORS](https://github.com/standage/AEGeAn/blob/master/docs/CONTRIBUTORS) for details.

The AEGeAn Toolkit started as several distinct but related efforts to build tools for managing and analyzing whole-genome gene structure annotations.
AEGeAn has brought these efforts together into a single library that includes binaries as well as several data structures and modules available via a C API.
The AEGeAn Toolkit leverages a variety of parsers, data structures, and graphics capabilities available from the GenomeTools library ([http://genometools.org](http://genometools.org)).

If you have any questions regarding AEGeAn, feel free to contact the author by email ([daniel.standage@gmail.com](mailto:daniel.standage@gmail.com)). Even better, open up a thread on [AEGeAn's issue tracker](https://github.com/standage/AEGeAn/issues) so that I can share my response with others who may have the same questions or issues.

### Documentation
Full installation instructions are provided in [docs/INSTALL](https://github.com/standage/AEGeAn/blob/master/docs/INSTALL). The example below (for the impatient) shows default installation on an Ubuntu system. 

```bash
# Install pre-requisites with your package manager
apt-get install -y build-essential git libcairo2-dev libncurses5-dev
# Make sure /usr/local/bin is in your $PATH, and that /usr/local/lib is in your
# LD path

git clone https://github.com/genometools/genometools.git
cd genometools
make 64bit=yes
sudo make 64bit=yes install

cd ..
git clone https://github.com/standage/AEGeAn.git
cd AEGeAn
make 64bit=yes
sudo make 64bit=yes install
```

Comprehensive documentation for the AEGeAn Toolkit is coming soon!
In the mean time, each program provides a descriptive usage statement, and the C API is documented in the header files located in the `inc/core` directory.

###Publications
**Daniel S. Standage** and Volker P. Brendel (2012) ParsEval: parallel comparison and analysis of gene structure annotations. BMC Bioinformatics, **13**:187, [doi:10.1186/1471-2105-13-187](http://dx/doi.org/10.1186/1471-2105-13-187).
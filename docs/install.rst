Installing AEGeAn
=================


For the impatient
-----------------

.. code-block:: bash

    # Install pre-requisites with your package manager
    sudo apt-get install -y build-essential git libcairo2-dev libpango1.0-dev
    # Make sure /usr/local/bin is in your $PATH, and that /usr/local/lib is in
    # your LD path

    git clone https://github.com/genometools/genometools.git
    cd genometools
    make 64bit=yes
    sudo make 64bit=yes install

    cd ..
    git clone https://github.com/standage/AEGeAn.git
    cd AEGeAn
    make 64bit=yes
    sudo make 64bit=yes install


Prerequisites
-------------
In principle, AEGeAn should compile and run on any POSIX-compliant UNIX system
(Linux, Mac OS X, Cygwin), although in practice, it has only been tested on
Linux and Mac systems. Native Windows support is not anticipated any time soon.
We have made an effort to minimize dependency on external software, but
compiling and installing does require a few prerequisites.

* common UNIX tools for building software: specifically, GNU make and a C
  compiler with OpenMP support
* GenomeTools library (see http://genometools.org/ for download and installation
  instructions)

While not a strict requirement, fully leveraging the graphics capabilities
provided by AEGeAn requires that the system also have an installation of the
Cairo graphics library (see the GenomeTools documentation).


Downloading
-----------
The easiest way download the latest and greatest version of AEGeAn (including
the most recent updates) is to clone from Github.

.. code-block:: bash

    git clone https://github.com/standage/AEGeAn.git

Alternatively, you can download an archive of the source from
https://github.com/standage/AEGeAn.

If you would prefer to work with a thoroughly tested stable release, do one of
the following.

* If you used git to clone AEGeAn, then run ``git tag`` to list all stable
  releases. Then before compiling run ``git checkout vx.y.z`` where ``vx.y.z``
  is the latest stable release.

* Go to `AEGeAn's release page on Github
  <https://github.com/standage/AEGeAn/releases>`_ and download an archive of the
  latest release.

**Note**: AEGeAn uses `Semantic Versioning <http://semver.org>`_ for labeling
stable releases.


Compiling and installing
------------------------

It is generally assumed that the user will follow standard Linux/UNIX procedures
when downloading, compiling, and installing the GenomeTools prerequisite. The
default install directory is ``/usr/local``, which means that we assume the
GenomeTools binaries have been installed in ``/usr/local/bin`` and the shared
libraries have been installed in ``/usr/local/lib``.

If this is the case, then the default installation procedure should work just
fine. From the directory containing the AEGeAn source code distribution, enter
the following commands.

.. code-block:: bash

    make
    make install # requires root access

If you have some reason to install GenomeTools in some other directory (e.g.
if you don't have permissions to modify ``/usr/local``), this should not be a
problem for AEGeAn. You must simply tell AEGeAn where to find GenomeTools. For
example, if you installed GenomeTools in ``/home/alice/local``, then you must
indicate this fact when compiling AEGeAn.

The following make options are available if you encounter any problems during
installation or if you want to customize your installation.

* ``64bit=yes`` if you want compile for a 64-bit architecture
* ``cairo=no`` if your system does not have Cairo graphics libraries installed
  (you will not be able to generate graphics)
* ``GT_INSTALL_DIR=?`` if you did not install GenomeTools in ``/usr/local``
* ``prefix=?`` to install AEGeAn in a directory other than
  ``/usr/local`` (e.g. if you do not have admin privileges)

For example, if I need to install GenomeTools and AEGeAn in my home directory
``/home/alice``, then I could do something like this.

.. code-block:: bash

    test -d /home/alice/local/src || mkdir -p /home/alice/local/src

    cd /home/alice/local/src
    git clone https://github.com/genometools/genometools.git
    cd genometools
    make prefix=/home/alice/local 64bit=yes
    make prefix=/home/alice/local 64bit=yes install

    cd /home/alice/local/src
    git clone git://github.com/standage/AEGeAn.git
    cd AEGeAn
    make prefix=/home/alice/local GT_INSTALL_DIR=/home/alice/local
    make prefix=/home/alice/local GT_INSTALL_DIR=/home/alice/local install

Remember that if you install the GenomeTools library in a non-standard location,
you will need to make sure that the AEGeAn binariescan find that library at
runtime. This can be done on a temporary basis using the ``LD_LIBRARY_PATH``
environmental variable or on a permanent system-wide basis using the
``ldconfig`` command.

Although the ``/usr/local`` directory is the standard install location for third
party libraries, on some distributions this directory is not pre-configured.
If you use the default install locations and still run into problems, make sure
that ``/usr/local/bin`` is in your path (using the ``export`` or ``setenv``
commands) and that ``/usr/local/lib`` is in the LD path (using the ``ldconfig``
command).

If you are not familiar with system administration, see the
:ref:`appendix <appendix-config>` below, which includes instructions for
installing prerequisites and setting up system paths on a variety of operating
systems.


.. _appendix-config:

Appendix: system setup
----------------------
Below are instructions for installing prerequisites and configuring system paths
for the most common operating systems. You'll want to 

* Debian-based systems including Ubuntu, Mint/LMDE, etc (tested on Ubuntu 11.10)

  .. code-block:: bash
  
      echo $PATH | grep /usr/local/bin > /dev/null
      if [ $? != 0 ]; then
        export PATH=/usr/local/bin:$PATH
        echo 'export PATH=/usr/local/bin:$PATH' >> /etc/bashrc
      fi
      grep '/usr/local/lib' /etc/ld.so.conf /etc/ld.so.conf.d/* > /dev/null
      if [ $? != 0 ]; then
        echo '/usr/local/lib' >> /etc/ld.so.conf.d/genometools-x86_64.conf
      fi
      ldconfig
      apt-get install -y build-essential git libcairo2-dev libpango1.0-dev

* Red Hat-based systems including CentOS, Fedora, etc (tested on CentOS 5.3)

  .. code-block:: bash
  
      echo $PATH | grep /usr/local/bin > /dev/null
      if [ $? != 0 ]; then
        export PATH=/usr/local/bin:$PATH
        echo 'export PATH=/usr/local/bin:$PATH' >> /etc/bashrc
      fi
      grep '/usr/local/lib' /etc/ld.so.conf /etc/ld.so.conf.d/* > /dev/null
      if [ $? != 0 ]; then
        echo '/usr/local/lib' >> /etc/ld.so.conf.d/genometools-x86_64.conf
      fi
      /sbin/ldconfig
      yum install -y git cairo-devel pango-devel

* Mac OS X (tested on Mac OS 10.6)

  .. code-block:: bash
  
      # Download and install Git: http://git-scm.com
      # Download and install the Fink package manager: http://www.finkproject.org/download
      # Then install the following packages using Fink
      apt-get install -y cairo-devel ncurses-devel pango1-xft2-ft219-dev
      

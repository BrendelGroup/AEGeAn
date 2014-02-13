Installing AEGeAn
=================


For the impatient
-----------------

.. code-block:: bash

    # Install pre-requisite packages with your OS's package manager
    sudo apt-get install -y build-essential git libcairo2-dev libpango1.0-dev

    # Make sure that /usr/local/lib is in your LD path
    sudo sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/aegean-x86_64.conf'
    sudo ldconfig

    # Download the AEGeAn source code
    git clone https://github.com/standage/AEGeAn.git
    cd AEGeAn

    # Download the GenomeTools source code
    git submodule init
    git submodule update

    # Compile and install GenomeTools and AEGeAn
    make
    sudo make install


Prerequisites
-------------
In principle, AEGeAn should compile and run on any POSIX-compliant UNIX system
(Linux, Mac OS X, Cygwin), although in practice, it has only been tested on
Linux and Mac systems. Native Windows support is not anticipated any time soon.

We have made an effort to minimize dependency on external software. Aside from
the GenomeTools library (which is included in the AEGeAn source code
distribution), compiling and installing AEGeAn requires only GNU make and a C
compiler.

While not a strict requirement, fully leveraging the graphics capabilities
provided by AEGeAn (through GenomeTools) requires that the system also have an
installation of the Cairo graphics library (see "`Appendix: system setup`_" for
platform-specific installation instructions).


Downloading
-----------
The easiest way download the latest and greatest version of AEGeAn (including
the most recent updates) is to clone from Github.

.. code-block:: bash

    git clone https://github.com/standage/AEGeAn.git

Alternatively, you can download an archive of the source from
https://github.com/standage/AEGeAn.

If you would prefer instead to work with a more  thoroughly tested stable
release, do one of the following.

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

The instructions under the section labeled "`For the impatient`_" are a good
starting point for most users. Chances are many users will be able to
successfully install by running those commands verbatim. The first few commands
are specific to Debian-based Linux operating systems such as Ubuntu.
If you are using a different operating system, please see "`Appendix:
system setup`_"  for platform-specific instructions. After running those
commands, complete the installation by running ``make`` and ``make install``.

If you have already installed GenomeTools system-wide and do not want to
overwrite that installation, you can replace 

..code-block:: bash

    git submodule init
    git submodule update
    make
    make install

with

..code-block:: bash

    make agn
    make agn-install

See the section labeled `Compilation flags`_ for a complete description of
configurable settings for compilation and installation.

Compilation flags
~~~~~~~~~~~~~~~~~

Compilation settings can be configured using the following flags with the
``make`` command.

* ``64bit=no``: do not compile for a 64-bit architecture
* ``cairo=no``: compile without graphics support (if your system does not have
  the Cairo graphics libraries installed)
* ``prefix=$DIR``: install GenomeTools and AEGeAn in ``$DIR`` rather than the
  default directory ``/usr/local``
* ``optimize=yes``: enable performance optimization for the AEGeAn code
* ``errorcheck=no``: allow code to compile even if there are warnings
* ``debug=no``: disable debugging support
* ``clean``: remove all compiler-generated AEGeAn files
* ``clean-all``: remove all compiler generated files for AEGeAn and GenomeTools

For example, if you want to compile the code with performance optimizations
enabled and graphics support disabled, run ``make optimize=yes cairo=no``
instead of just ``make``.

Compiling without administrative privileges
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default installation location is ``/usr/local/``, which means:

* programs are installed in ``/usr/local/bin``
* header files are installed in ``/usr/local/include``
* libraries are installed in ``/usr/local/lib``
* auxiliary data files are installed in ``/usr/local/share``

If you do not have administrative privileges on your
machine, then you will not be able to run ``make install`` without specifying an
alternative installation directory with ``prefix``. Creating an installation
directory within your home directory, as shown in the following example, is
recommended.

.. code-block:: bash

  mkdir ~/local
  make prefix=~/local
  make prefix=~/local install

This will install the programs in ``~/local/bin``, the libraries in
``~/local/lib``, etc. You will probably want to add ``~/local/bin`` to your
``PATH`` environmental variable and ``~/local/lib`` to your ``LD_LIBRARY_PATH``
environmental variable (or ``DYLD_LIBRARY_PATH`` on Mac OS X).

.. _appendix-config:

Appendix: system setup
----------------------
Below are instructions for installing prerequisites and configuring system paths
for the most common operating systems. Note that running these commands requires
administrative/sudo privileges.

* Debian-based systems including Ubuntu, Mint/LMDE, etc (tested on Ubuntu 11.10)

  .. code-block:: bash
  
      sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/aegean-x86_64.conf'
      ldconfig
      apt-get install -y build-essential git libcairo2-dev libpango1.0-dev

* Red Hat-based systems including CentOS, Fedora, etc (tested on CentOS 5.3)

  .. code-block:: bash
  
      sh -c 'echo "/usr/local/lib" > /etc/ld.so.conf.d/aegean-x86_64.conf'
      /sbin/ldconfig
      yum install -y git cairo-devel pango-devel

* Mac OS X (tested on Mac OS 10.6)

  .. code-block:: bash
  
      # Download and install Git: http://git-scm.com
      # Download and install the Fink package manager: http://www.finkproject.org/download
      # Then install the following packages using Fink
      apt-get install -y cairo-devel pango1-xft2-ft219-dev
      

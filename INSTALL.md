# AEGeAn Installation and Setup

## Installation as a singularity container [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3327)

All the AEGeAn dependencies are encapsulated in a
[Singularity](https://www.sylabs.io/docs/) container available from
[Singularity Hub](https://singularity-hub.org/).

Assuming _git_ and  _singularity_ are installed on your system, you can get the
AEGeAn code from GitHub and the container from the
[Singularity Hub](https://www.singularity-hub.org/collections/3327) as follows:

```bash
git clone https://github.com/BrendelGroup/AEGeAn.git
singularity pull --name aegean.simg shub://BrendelGroup/AEGeAn
singularity exec -e -B `pwd` aegean.simg fidibus -h
```

Here the last command (_singularity exec_) will execute the _fidibus_ script
and display its help message.

For a gentle introduction to singularity, see our group
[handbook article](https://github.com/BrendelGroup/bghandbook/blob/master/doc/06.2-Howto-Singularity-run.md).


## Optional: System-wide Installation

AEGeAn use via the singularity container is highly recommended, with no known
drawbacks.
However, if desired, you can of course install all the required third party
software and our _fidibus_ Python package individually on your computer system.
The singularity [recipe file](./Singularity) in this repository should serve as
a guide to perform such an installation.
The `aegean.simg` container was built on the 
[current Fedora 30 release](https://getfedora.org/)
and thus the instructions apply to that particular Linux version.
For different Linux distributions, you will have to install the equivalent
packages using your distribution's package manager.

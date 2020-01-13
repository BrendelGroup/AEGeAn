#!/bin/sh
#
# ... use ".uninstall.sh --user" to uninstall a user local copy of LocusPocus.

python3 setup.py install $1 --record files.txt

cat files.txt | xargs rm -rf
\rm -rf LocusPocus.egg-info files.txt LocusPocus.egg-info dist build versioneer.pyc

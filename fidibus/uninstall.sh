#!/bin/sh
#
# ... use ".uninstall.sh --user" to uninstall a user local copy of fidibus.

python setup.py install $1 --record files.txt

cat files.txt | xargs rm -rf
\rm -rf fidibus.egg-info files.txt fidibus.egg-info dist build versioneer.pyc

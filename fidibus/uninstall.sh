# Next time you need to install something with python setup.py -- which should be never but things happen. 

python setup.py install --user --record files.txt

# This will cause all the installed files to be printed to that directory.
# Then when you want to uninstall it simply run; be careful with the 'sudo'

cat files.txt | xargs rm -rf

\rm -rf fidibus.egg-info files.txt fidibus.egg-info dist build versioneer.pyc

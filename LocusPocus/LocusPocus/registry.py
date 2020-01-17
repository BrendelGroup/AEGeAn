#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2015-2016   Indiana University
#
# This file is part of AEGeAn (http://github.com/BrendelGroup/AEGeAn) and is
# licensed under the ISC license: see LICENSE.
# -----------------------------------------------------------------------------

"""Module implementing a registry for handling genome configuration files."""

from __future__ import print_function
import glob
import os
import pkg_resources
import yaml
import LocusPocus
try:
    FileNotFoundError
except NameError:  # pragma: no cover
    FileNotFoundError = IOError


class Registry(object):

    def __init__(self):
        """Initialize the registry with the default LocusPocus configs."""
        genome_configs = pkg_resources.resource_filename('LocusPocus', '../genome_configs')
        self.update(genome_configs, clear=True)

    def update(self, path, clear=False):
        """
        Update the registry from the given directory path.

        The registry will attempt to load all .yml files as genome configs and
        .txt files as batch configs. If `clear` is true, any previous entries
        in the registry will be cleared before loading new entries.
        """
        if clear:
            self.genome_configs = dict()
            self.batch_configs = dict()

        if not os.path.exists(path):
            message = 'config directory "%s" does not exist' % path
            raise FileNotFoundError(message)

        for filepath in glob.glob(path + '/*.yml'):
            configs = self.parse_genome_config(filepath)
            self.genome_configs.update(configs)

        for filepath in glob.glob(path + '/*.txt'):
            batch = self.parse_batch_config(filepath)
            filename = os.path.basename(filepath)
            batch_label = os.path.splitext(filename)[0]
            self.batch_configs[batch_label] = batch

    def check(self, genomes=None, batches=None):
        if genomes:
            for label in genomes:
                assert label in self.genome_configs, 'unknown genome ' + label
        if batches:
            for label in batches:
                assert label in self.batch_configs, 'unknown batch ' + label

    def config(self, label):
        if label not in self.genome_configs:
            return None
        return self.genome_configs[label]

    def genome(self, label, workdir='.'):
        if label not in self.genome_configs:
            return None
        config = self.genome_configs[label]
        constructor = LocusPocus.dbtype[config['source']]
        db = constructor(label, config, workdir=workdir)
        return db

    def batch(self, batch_label):
        """Retrieve a batch of genome configs from the registry."""
        if batch_label not in self.batch_configs:
            return None
        return list(self.batch_configs[batch_label])

    def parse_genome_config(self, config):
        """
        Parse reference genome config in YAML format.

        If `config` is a string it is treated as a filename, otherwise as a
        file handle or similar object.
        """
        if isinstance(config, str):
            config = open(config, 'r')
        return yaml.safe_load(config)

    def parse_batch_config(self, config):
        """
        Parse a batch of reference genome config labels from a file.

        The file should contain a single genome config label per line, and
        `config` is treated as a filename.
        """
        batch = list()
        with open(config, 'r') as infile:
            for line in infile:
                label = line.strip()
                if label != '':
                    batch.append(label)
        return batch

    def list_genomes(self):
        for label in sorted(self.genome_configs):
            yield label, self.genome_configs[label]

    def list_batches(self):
        for label in sorted(self.batch_configs):
            yield label, self.batch_configs[label]

    def list(self, outstream=None):  # pragma: no cover
        print('===== Reference genomes =====', file=outstream)
        for label, config in self.list_genomes():
            info = label + '\t' + config['species']
            if 'common' in config:
                info += ' ({})'.format(config['common'])
            info += '\t' + LocusPocus.sources[config['source']]
            print(info, file=outstream)
        print('', file=outstream)

        print('===== Reference genome batches =====', file=outstream)
        for label, batch in self.list_batches():
            print(label, ','.join(batch), sep='\t', file=outstream)


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_list():
    """Registry: listing genome and batch configs"""
    registry = Registry()

    genome_labels = [x for x in registry.list_genomes()]
    ymlfiles = [x for x in glob.glob('genome_configs/*.yml')]
    assert len(genome_labels) == len(ymlfiles)
    batch_labels = [x for x in registry.list_batches()]
    txtfiles = [x for x in glob.glob('genome_configs/*.txt')]
    assert len(batch_labels) == len(txtfiles)

    registry.check(genomes=['Otau', 'Oluc'])
    registry.check(batches=['chlorophyta'])
    try:
        registry.check(batches=['BogusFooBar'])
    except AssertionError:
        pass

    registry.update('testdata/conf', clear=True)

    genome_labels = [x for x in registry.list_genomes()]
    ymlfiles = [x for x in glob.glob('testdata/conf/*.yml')]
    assert len(genome_labels) == len(ymlfiles)
    batch_labels = [x for x in registry.list_batches()]
    txtfiles = [x for x in glob.glob('testdata/conf/*.txt')]
    assert len(batch_labels) == len(txtfiles)

    try:
        registry.update('foobar/bogus')
    except FileNotFoundError:
        pass


def test_genome():
    """Registry: loading a genome db or configuration by label"""
    registry = Registry()
    db = registry.genome('Osat')
    assert 'accession' in db.config
    assert db.config['accession'] == 'GCF_001433935.1'

    registry.update('testdata/conf')
    db = registry.genome('Osat')
    assert 'accession' in db.config
    assert db.config['accession'] == 'BogusThisIsNotARealAccession'

    dbconfig = registry.config('Bvul')
    assert 'scaffolds' in dbconfig
    assert dbconfig['scaffolds'] == 'bv_ref_1.1_chrUn.fa.gz'

    assert registry.genome('Docc') is None
    assert registry.config('gnilleps') is None


def test_batch():
    """Registry: loading a batch configuration by label"""
    registry = Registry()
    labels = registry.batch('honeybees')
    assert sorted(labels) == ['Ador', 'Aflo', 'Amel']

    registry.update('testdata/conf')
    labels = registry.batch('mythical')
    assert labels == ['Bvul']

    assert registry.batch('nonexistent') is None


def test_parse_genome_config():
    """Registry: parsing genome configurations from a file"""
    registry = Registry()
    with open('genome_configs/Pbar.yml', 'r') as filehandle:
        config = registry.parse_genome_config(filehandle)
        assert len(config) == 1
        assert 'Pbar' in config
        assert config['Pbar']['species'] == 'Pogonomyrmex barbatus'

    config = registry.parse_genome_config('genome_configs/Hlab.yml')
    assert len(config) == 1
    assert 'Hlab' in config
    assert config['Hlab']['common'] == 'blueberry bee'

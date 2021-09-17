"""Ntuples for EMTF++."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

import six
from six.moves import range, zip, map, filter


# ______________________________________________________________________________
# Configs

# Parent directory
_EOS_PREFIX = 'root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_11_1_7/'


# ______________________________________________________________________________
# Functions

def load_tree(infiles,
              load_hits=True,
              load_simhits=True,
              load_tracks=True,
              load_tttracks=False,
              load_particles=True,
              use_multiprocessing=False):
  from emtf_tree import TreeChain, TreeQueue

  my_collections = []  # used for tree.define_collection()
  my_branches = []  # used for tree.activate()
  cache_size = 10000000  # 10 MB
  if load_hits:
    my_collections.append(dict(name='hits', prefix='vh_', size='vh_size'))
    my_branches.append('vh_*')
  if load_simhits:
    my_collections.append(dict(name='simhits', prefix='vc_', size='vc_size'))
    my_branches.append('vc_*')
  if load_tracks:
    my_collections.append(dict(name='tracks', prefix='vt_', size='vt_size'))
    my_branches.append('vt_*')
  if load_tttracks:
    my_collections.append(dict(name='tttracks', prefix='vd_', size='vd_size'))
    my_branches.append('vd_*')
  if load_particles:
    my_collections.append(dict(name='particles', prefix='vp_', size='vp_size'))
    my_branches.append('vp_*')

  if use_multiprocessing:
    tree = TreeQueue('ntupler/tree', infiles, cache=True, cache_size=cache_size, branches=my_branches)
  else:
    tree = TreeChain('ntupler/tree', infiles, cache=True, cache_size=cache_size, branches=my_branches)

  for coll_kwargs in my_collections:
    tree.define_collection(**coll_kwargs)
  return tree


def load_pgun_test():
  infile = '../test/ntuple_SingleMuon_Endcap_add.20210808.root'
  return load_tree(infile)


def load_pgun_batch(k):
  dataset = SingleMuon()
  my_range = np.split(np.arange(len(dataset)), 100)[k]
  infiles = dataset[my_range]
  return load_tree(infiles)


def load_pgun_displaced_batch(k):
  dataset = SingleMuonDisplaced()
  my_range = np.split(np.arange(len(dataset)), 100)[k]
  infiles = dataset[my_range]
  return load_tree(infiles)


# ______________________________________________________________________________
# Datasets

class _BaseDataset(object):
  """Abstract base class."""

  def __len__(self):
    return 0

  def __getitem__(self, index):
    raise IndexError


class SingleMuon(_BaseDataset):
  SIZE = 2000
  PREFIX_0 = _EOS_PREFIX + 'SingleMuon_PosEnd_2GeV_Phase2HLTTDRSummer20/ParticleGuns/CRAB3/210808_184102/'
  PREFIX_1 = _EOS_PREFIX + 'SingleMuon_NegEnd_2GeV_Phase2HLTTDRSummer20/ParticleGuns/CRAB3/210808_184117/'

  def __len__(self):
    return self.SIZE

  def __getitem__(self, index):
    if hasattr(index, '__iter__'):
      index_range = index
      return [self[index] for index in index_range]
    elif isinstance(index, slice):
      index_range = range(*index.indices(len(self)))
      return self[index_range]
    return self.getitem(index)

  def getitem(self, index):
    if index % 2 == 0:
      prefix = self.PREFIX_0
    else:
      prefix = self.PREFIX_1
    i = index // 2
    filename = '{0}{1:04d}/ntuple_{2}.root'.format(prefix, (i + 1) // 1000, (i + 1))
    return filename


class SingleMuonDisplaced(_BaseDataset):
  pass

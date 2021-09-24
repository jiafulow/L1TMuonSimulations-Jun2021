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
              load_hits=False,
              load_simhits=False,
              load_tracks=False,
              load_tttracks=False,
              load_particles=False,
              use_multiprocessing=False):
  from emtf_tree import TreeChain, TreeQueue

  my_collections = []  # used for tree.define_collection()
  keep_branches = []
  ignore_branches = []  # deactivate certain branches to make it faster
  cache_size = 10000000  # 10 MB
  if load_hits:
    my_collections.append(dict(name='hits', prefix='vh_', size='vh_size'))
    keep_branches.append('vh_*')
    ignore_branches.append('vh_glob_phi')
    ignore_branches.append('vh_glob_theta')
    ignore_branches.append('vh_glob_perp')
    ignore_branches.append('vh_glob_z')
    ignore_branches.append('vh_glob_time')
  if load_simhits:
    my_collections.append(dict(name='simhits', prefix='vc_', size='vc_size'))
    keep_branches.append('vc_*')
    ignore_branches.append('vc_mom_p')
    ignore_branches.append('vc_mom_phi')
    ignore_branches.append('vc_mom_theta')
    ignore_branches.append('vc_tof')
  if load_tracks:
    my_collections.append(dict(name='tracks', prefix='vt_', size='vt_size'))
    keep_branches.append('vt_*')
    ignore_branches.append('vt_pt')
    ignore_branches.append('vt_phi')
    ignore_branches.append('vt_theta')
    ignore_branches.append('vt_eta')
    ignore_branches.append('vt_invpt')
    ignore_branches.append('vt_d0')
    ignore_branches.append('vt_q')
    ignore_branches.append('vt_hw_*')
  if load_tttracks:
    my_collections.append(dict(name='tttracks', prefix='vd_', size='vd_size'))
    keep_branches.append('vd_*')
  if load_particles:
    my_collections.append(dict(name='particles', prefix='vp_', size='vp_size'))
    keep_branches.append('vp_*')
    ignore_branches.append('vp_beta')
    ignore_branches.append('vp_mass')
    ignore_branches.append('vp_event')
    ignore_branches.append('vp_genp')

  if use_multiprocessing:
    tree_cls = TreeQueue
  else:
    tree_cls = TreeChain
  tree_name = 'ntupler/tree'
  tree = tree_cls(tree_name, infiles, cache=True, cache_size=cache_size,
                  branches=keep_branches, ignore_branches=ignore_branches)
  for coll_kwargs in my_collections:
    tree.define_collection(**coll_kwargs)
  return tree


def load_pgun_test():
  infile = '../test/ntuple_SingleMuon_Endcap_Phase2HLTTDRSummer20.210922.root'
  return load_tree(infile, load_hits=True, load_simhits=True, load_tracks=True,
                   load_particles=True)


def load_pgun_batch(k):
  dataset = SingleMuon()
  my_range = np.split(np.arange(len(dataset)), 100)[k]
  infiles = dataset[my_range]
  return load_tree(infiles, load_hits=True, load_simhits=True, load_tracks=True,
                   load_particles=True)


def load_pgun_displaced_batch(k):
  dataset = SingleMuonDisplaced()
  my_range = np.split(np.arange(len(dataset)), 100)[k]
  infiles = dataset[my_range]
  return load_tree(infiles, load_hits=True, load_simhits=True, load_tracks=True,
                   load_particles=True)


# ______________________________________________________________________________
# Datasets

class _BaseDataset(object):
  """Abstract base class.

  Sub-classes must override SIZE, PREFIX, and FNAME.
  """
  SIZE = 0
  PREFIX = ''
  FNAME = 'ntuple.root'

  def __len__(self):
    return self.SIZE

  def __getitem__(self, index):
    if hasattr(index, '__iter__'):
      index_range = index
      return [self[index] for index in index_range]
    elif isinstance(index, slice):
      index_range = range(*index.indices(len(self)))
      return self[index_range]
    else:
      return self.getitem(index)

  def __iter__(self):
    for i in range(len(self)):
      yield self[i]

  def __contains__(self, value):
    for v in self:
      if v is value or v == value:
        return True
    return False

  def _filename(self, prefix, fname, i):
    if not prefix or not fname:
      raise ValueError('Incorrect prefix or fname.')
    if not fname.endswith('.root'):
      raise ValueError('fname must be a .root file.')
    return '{0}{1:04d}/{2}_{3}.root'.format(prefix, (i + 1) // 1000, fname[:-5], (i + 1))

  def getitem(self, index):
    return self._filename(self.PREFIX, self.FNAME, index)


class SingleMuon(_BaseDataset):
  SIZE = 2000
  PREFIX_0 = _EOS_PREFIX + 'SingleMuon_PosEnd_2GeV_Phase2HLTTDRSummer20/ParticleGuns/CRAB3/210922_194420/'
  PREFIX_1 = _EOS_PREFIX + 'SingleMuon_NegEnd_2GeV_Phase2HLTTDRSummer20/ParticleGuns/CRAB3/210922_194440/'

  def getitem(self, index):
    if index % 2 == 0:
      prefix = self.PREFIX_0
    else:
      prefix = self.PREFIX_1
    index = index // 2
    return self._filename(prefix, self.FNAME, index)


class SingleMuonDisplaced(_BaseDataset):
  pass


class SingleNeutrinoPU200(_BaseDataset):
  SIZE = 56 - 2
  PREFIX = _EOS_PREFIX + 'SingleNeutrino_PU200_ext1_Phase2HLTTDRSummer20/MinBias_TuneCP5_14TeV-pythia8/CRAB3/210921_154251/'

  def getitem(self, index):
    indexer = np.arange(self.SIZE, dtype=np.int32)
    indexer[[32 - 1, 35 - 1]] = [56 - 2, 56 - 1]  # fill the holes left by failed jobs
    return self._filename(self.PREFIX, self.FNAME, indexer[index])


class DoubleMuonPU200(_BaseDataset):
  SIZE = 28
  PREFIX = _EOS_PREFIX + 'DoubleMuon_PU200_Phase2HLTTDRSummer20/DoubleMuon_gun_FlatPt-1To100/CRAB3/210921_145639/'


class DoublePhotonPU200(_BaseDataset):
  SIZE = 12 + 8
  PREFIX_0 = _EOS_PREFIX + 'DoublePhoton_PU140_Phase2HLTTDRSummer20/DoublePhoton_FlatPt-1To100/CRAB3/210922_190214/'
  PREFIX_1 = _EOS_PREFIX + 'DoublePhoton_PU200_ext2_Phase2HLTTDRSummer20/DoublePhoton_FlatPt-1To100/CRAB3/210922_190246/'

  def getitem(self, index):
    if index < 12:
      prefix = self.PREFIX_0
    else:
      prefix = self.PREFIX_1
    index = index % 12
    return self._filename(prefix, self.FNAME, index)

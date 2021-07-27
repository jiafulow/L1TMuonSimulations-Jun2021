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
eos_prefix = 'root://cmsxrootd-site.fnal.gov//store/group/l1upgrades/L1MuonTrigger/P2_11_1_7/'


# ______________________________________________________________________________
# Functions

# rootpy is no longer included in 11_1_7, as supposedly it doesn't work with the newer version of root/pyroot.
# But it does at least work in 11_1_7.
try:
  import rootpy.tree
except ImportError:
  import os
  import sys
  if os.environ['CMSSW_VERSION'] == 'CMSSW_11_1_7':
    sys.path.insert(1, '/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/py2-rootpy/1.0.1-bcolbf4/lib/python2.7/site-packages/')
  import rootpy.tree


def load_tree(infiles):
  from rootpy.tree import TreeChain
  from rootpy.ROOT import gROOT
  gROOT.SetBatch(True)

  def truncate_list(lst):
    if isinstance(lst, list) and len(lst) > 10:
      return lst[:5] + ['...'] + lst[-5:]
    else:
      return lst
  print('[INFO] Opening files: {0}'.format(truncate_list(infiles)))
  tree = TreeChain('ntupler/tree', infiles, cache=True)
  tree.define_collection(name='hits', prefix='vh_', size='vh_size')
  tree.define_collection(name='simhits', prefix='vc_', size='vc_size')
  tree.define_collection(name='tracks', prefix='vt_', size='vt_size')
  tree.define_collection(name='particles', prefix='vp_', size='vp_size')
  return tree


def load_pgun_test():
  infile = '../test/ntuple_SingleMuon_Endcap_add.20210724.root'
  return load_tree(infile)


def load_pgun_batch(k):
  my_range = np.split(np.arange(1000), 100)[k]
  my_eos_prefix1 = eos_prefix + 'SingleMuon_PosEnd_2GeV_Phase2HLTTDRSummer20/ParticleGuns/CRAB3/210727_191232/'
  my_eos_prefix2 = eos_prefix + 'SingleMuon_NegEnd_2GeV_Phase2HLTTDRSummer20/ParticleGuns/CRAB3/210727_191248/'
  infiles = []
  for i in my_range:
    infiles.append(my_eos_prefix1 + '{0:04d}/ntuple_{1}.root'.format((i + 1) // 1000, (i + 1)))
    infiles.append(my_eos_prefix2 + '{0:04d}/ntuple_{1}.root'.format((i + 1) // 1000, (i + 1)))
  tree = load_tree(infiles)
  return tree


def load_pgun_displaced_batch(k):
  return NotImplemented

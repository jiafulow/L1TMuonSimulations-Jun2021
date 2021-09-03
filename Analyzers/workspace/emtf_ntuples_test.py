"""Tests for functions in emtf_ntuples.py"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pytest

import numpy as np

from emtf_ntuples import *


# ______________________________________________________________________________
def test_singlemuon():
  dataset = SingleMuon()
  assert len(dataset) != 0
  assert dataset[0] != ''

def test_singleneutrino_pu200():
  dataset = SingleNeutrinoPU200()
  assert len(dataset) != 0
  assert dataset[0] != ''

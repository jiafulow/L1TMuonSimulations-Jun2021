"""Tests for functions in emtf_utils.py"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pytest

import numpy as np

from emtf_utils import *


# ______________________________________________________________________________
def test_wrap_phi_rad():
  input_range = np.deg2rad(np.arange(-180,180)) + 1e-5
  output_range = np.deg2rad(np.arange(-180,180)) + 1e-5
  assert list(wrap_phi_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(180,180+360)) + 1e-5
  assert list(wrap_phi_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(-180-360,-180)) + 1e-5
  assert list(wrap_phi_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(180+360,180+720)) + 1e-5
  assert list(wrap_phi_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(-180-720,-180-360)) + 1e-5
  assert list(wrap_phi_rad(input_range)) == pytest.approx(list(output_range))

def test_wrap_phi_deg():
  input_range = 1.0*(np.arange(-180,180)) + 1e-5
  output_range = 1.0*(np.arange(-180,180)) + 1e-5
  assert list(wrap_phi_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(180,180+360)) + 1e-5
  assert list(wrap_phi_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(-180-360,-180)) + 1e-5
  assert list(wrap_phi_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(180+360,180+720)) + 1e-5
  assert list(wrap_phi_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(-180-720,-180-360)) + 1e-5
  assert list(wrap_phi_deg(input_range)) == pytest.approx(list(output_range))

def test_wrap_theta_rad():
  input_range = np.deg2rad(np.arange(0,90)) + 1e-5
  output_range = np.deg2rad(np.arange(0,90)) + 1e-5
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(90,180) + 1) - 1e-5
  input_range = input_range[::-1]
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(-180,-90)) + 1e-5
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(-90,0) + 1) - 1e-5
  input_range = input_range[::-1]
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(180,270)) + 1e-5
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(270,360) + 1) - 1e-5
  input_range = input_range[::-1]
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

  input_range = np.deg2rad(np.arange(360,90+360)) + 1e-5
  assert list(wrap_theta_rad(input_range)) == pytest.approx(list(output_range))

def test_wrap_theta_deg():
  input_range = 1.0*(np.arange(0,90)) + 1e-5
  output_range = 1.0*(np.arange(0,90)) + 1e-5
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(90,180) + 1) - 1e-5
  input_range = input_range[::-1]
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(-180,-90)) + 1e-5
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(-90,0) + 1) - 1e-5
  input_range = input_range[::-1]
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(180,270)) + 1e-5
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(270,360) + 1) - 1e-5
  input_range = input_range[::-1]
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

  input_range = 1.0*(np.arange(360,90+360)) + 1e-5
  assert list(wrap_theta_deg(input_range)) == pytest.approx(list(output_range))

def test_get_trigger_sector():
  sectors_me21 = [6, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6]
  sectors_me22 = [6, 6, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,]

  station, ring = 1, 1
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22
  station, ring = 1, 2
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22
  station, ring = 1, 3
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22
  station, ring = 1, 4
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22
  station, ring = 2, 1
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,18+1)] == sectors_me21
  station, ring = 2, 2
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22
  station, ring = 3, 1
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,18+1)] == sectors_me21
  station, ring = 3, 2
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22
  station, ring = 4, 1
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,18+1)] == sectors_me21
  station, ring = 4, 2
  assert [get_trigger_sector(ring, station, chamber) for chamber in np.arange(1,36+1)] == sectors_me22

def test_get_trigger_subsector():
  subsectors_me1 = [2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2,]
  subsectors_me21 = [0] * 18
  subsectors_me22 = [0] * 36

  station, ring = 1, 1
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me1
  station, ring = 1, 2
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me1
  station, ring = 1, 3
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me1
  station, ring = 1, 4
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me1
  station, ring = 2, 1
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,18+1)] == subsectors_me21
  station, ring = 2, 2
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me22
  station, ring = 3, 1
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,18+1)] == subsectors_me21
  station, ring = 3, 2
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me22
  station, ring = 4, 1
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,18+1)] == subsectors_me21
  station, ring = 4, 2
  assert [get_trigger_subsector(ring, station, chamber) for chamber in np.arange(1,36+1)] == subsectors_me22

def test_get_trigger_cscid():
  cscids_me11 = [2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1,]
  cscids_me12 = [(x + 3) for x in cscids_me11]
  cscids_me13 = [(x + 6) for x in cscids_me11]
  cscids_me21 = [3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2]
  cscids_me22 = [8, 9, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7,]

  station, ring = 1, 1
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me11
  station, ring = 1, 2
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me12
  station, ring = 1, 3
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me13
  station, ring = 1, 4
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me11
  station, ring = 2, 1
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,18+1)] == cscids_me21
  station, ring = 2, 2
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me22
  station, ring = 3, 1
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,18+1)] == cscids_me21
  station, ring = 3, 2
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me22
  station, ring = 4, 1
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,18+1)] == cscids_me21
  station, ring = 4, 2
  assert [get_trigger_cscid(ring, station, chamber) for chamber in np.arange(1,36+1)] == cscids_me22

def test_get_trigger_neighid():
  neighids_me1 = [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,]
  neighids_me21 = [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]
  neighids_me22 = neighids_me1

  station, ring = 1, 1
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me1
  station, ring = 1, 2
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me1
  station, ring = 1, 3
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me1
  station, ring = 1, 4
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me1
  station, ring = 2, 1
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,18+1)] == neighids_me21
  station, ring = 2, 2
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me22
  station, ring = 3, 1
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,18+1)] == neighids_me21
  station, ring = 3, 2
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me22
  station, ring = 4, 1
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,18+1)] == neighids_me21
  station, ring = 4, 2
  assert [get_trigger_neighid(ring, station, chamber) for chamber in np.arange(1,36+1)] == neighids_me22

def test_get_trigger_endsec():
  assert [get_trigger_endsec(endcap, sector) for endcap in [1, -1] for sector in [1, 2, 3, 4, 5, 6]] == list(np.arange(12))

def test_find_median_of_three():
  assert find_median_of_three(10, 20, 30) == 20
  assert find_median_of_three(30, 20, 10) == 20
  assert find_median_of_three(20, 30, 10) == 20
  assert find_median_of_three(10, 30, 20) == 20
  assert find_median_of_three(30, 10, 20) == 20
  assert find_median_of_three(20, 10, 30) == 20
  assert find_median_of_three(20, 20, 20) == 20
  assert find_median_of_three(999999, 20, 30) == 20
  assert find_median_of_three(10, 999999, 30) == 10
  assert find_median_of_three(10, 20, 999999) == 10
  assert find_median_of_three(999999, 30, 20) == 20
  assert find_median_of_three(30, 999999, 10) == 10
  assert find_median_of_three(20, 10, 999999) == 10
  assert find_median_of_three(999999, 999999, 30) == 30
  assert find_median_of_three(10, 999999, 999999) == 10
  assert find_median_of_three(999999, 20, 999999) == 20
  assert find_median_of_three(999999, 999999, 999999) == 999999

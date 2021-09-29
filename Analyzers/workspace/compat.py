"""Various ompatibility fixes"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import collections

from emtf_utils import save_np_arrays


def load_signal(fname):
  with np.load(fname) as loaded:
    out_part = loaded['out_part']
    out_hits_values = loaded['out_hits_values']
    out_hits_row_splits = loaded['out_hits_row_splits']
    out_hits_shape = (out_hits_row_splits.shape[0] - 1,) + (None,) + out_hits_values.shape[1:]
    out_simhits_values = loaded['out_simhits_values']
    out_simhits_row_splits = loaded['out_simhits_row_splits']
    out_simhits_shape = (out_simhits_row_splits.shape[0] - 1,) + (None,) + out_simhits_values.shape[1:]
  return (out_part, (out_hits_values, out_hits_row_splits), (out_simhits_values, out_simhits_row_splits))


def recreate_signal_210112():
  old_fname = 'signal_add.20210112.npz'
  new_fname = 'signal_add.20210112-new.npz'

  # particle info (old)
  _part_fields_old = [
    'part_invpt',
    'part_eta',
    'part_phi',
    'part_vx',
    'part_vy',
    'part_vz',
    'part_d0',
    'part_sector',
    'part_zone',
  ]
  PartFieldsOld = collections.namedtuple('PartFieldsOld', _part_fields_old)
  part_fields_old = PartFieldsOld(*range(len(_part_fields_old)))

  # sector_hits info (old)
  _sector_hits_fields_old = [
    'emtf_site',
    'emtf_host',
    'emtf_chamber',
    'emtf_segment',
    'zones',
    'timezones',
    'emtf_phi',
    'emtf_bend',
    'emtf_theta',
    'emtf_theta_alt',
    'emtf_qual',
    'emtf_qual_alt',
    'emtf_time',
    'strip',
    'wire',
    'fr',
    'detlayer',
    'bx',
  ]
  SectorHitsFieldsOld = collections.namedtuple('SectorHitsFieldsOld', _sector_hits_fields_old)
  sector_hits_fields_old = SectorHitsFieldsOld(*range(len(_sector_hits_fields_old)))

  # particle info (new)
  _part_fields = [
    'part_invpt',
    'part_eta',
    'part_phi',
    'part_vx',
    'part_vy',
    'part_vz',
    'part_d0',
    'part_bx',
    'part_sector',
    'part_zone',
  ]
  PartFields = collections.namedtuple('PartFields', _part_fields)
  part_fields = PartFields(*range(len(_part_fields)))

  # sector_hits info (new)
  _sector_hits_fields = [
    'emtf_chamber',
    'emtf_segment',
    'emtf_phi',
    'emtf_bend',
    'emtf_theta1',
    'emtf_theta2',
    'emtf_qual1',
    'emtf_qual2',
    'emtf_time',
    'zones',
    'timezones',
    'cscfr',
    'gemdl',
    'bx',
    'emtf_site',
    'emtf_host',
    'valid',
  ]
  SectorHitsFields = collections.namedtuple('SectorHitsFields', _sector_hits_fields)
  sector_hits_fields = SectorHitsFields(*range(len(_sector_hits_fields)))

  # Load data
  out_part_old, out_hits_old, out_simhits_old = load_signal(old_fname)
  #print('out_part: {} out_hits: {} out_simhits: {}'.format(
  #    out_part_old.shape, (out_hits_old[0].shape, out_hits_old[1].shape),
  #    (out_simhits_old[0].shape, out_simhits_old[1].shape)))
  assert out_part_old.ndim == 2
  assert isinstance(out_hits_old, tuple) and len(out_hits_old) == 2
  assert out_hits_old[0].ndim == 2
  assert out_hits_old[1].ndim == 1
  assert isinstance(out_simhits_old, tuple) and len(out_simhits_old) == 2
  assert out_simhits_old[0].ndim == 2
  assert out_simhits_old[1].ndim == 1

  # Map out_part_old -> out_part
  _zeros = np.zeros(out_part_old.shape[0], dtype=out_part_old.dtype)
  out_part = np.column_stack((
    out_part_old[:, part_fields_old.part_invpt],
    out_part_old[:, part_fields_old.part_eta],
    out_part_old[:, part_fields_old.part_phi],
    out_part_old[:, part_fields_old.part_vx],
    out_part_old[:, part_fields_old.part_vy],
    out_part_old[:, part_fields_old.part_vz],
    out_part_old[:, part_fields_old.part_d0],
    _zeros,
    out_part_old[:, part_fields_old.part_sector],
    out_part_old[:, part_fields_old.part_zone],
  ))
  assert out_part_old.shape[1] == len(part_fields_old)
  assert out_part.shape[1] == len(part_fields)

  # Map out_hits_old -> out_hits
  _ones = np.ones(out_hits_old[0].shape[0], dtype=out_hits_old[0].dtype)
  out_hits_values = np.column_stack((
    out_hits_old[0][:, sector_hits_fields_old.emtf_chamber],
    out_hits_old[0][:, sector_hits_fields_old.emtf_segment],
    out_hits_old[0][:, sector_hits_fields_old.emtf_phi],
    out_hits_old[0][:, sector_hits_fields_old.emtf_bend],
    out_hits_old[0][:, sector_hits_fields_old.emtf_theta],
    out_hits_old[0][:, sector_hits_fields_old.emtf_theta_alt],
    out_hits_old[0][:, sector_hits_fields_old.emtf_qual],
    out_hits_old[0][:, sector_hits_fields_old.emtf_qual_alt],
    out_hits_old[0][:, sector_hits_fields_old.emtf_time],
    out_hits_old[0][:, sector_hits_fields_old.zones],
    out_hits_old[0][:, sector_hits_fields_old.timezones],
    out_hits_old[0][:, sector_hits_fields_old.fr],
    out_hits_old[0][:, sector_hits_fields_old.detlayer],
    out_hits_old[0][:, sector_hits_fields_old.bx],
    out_hits_old[0][:, sector_hits_fields_old.emtf_site],
    out_hits_old[0][:, sector_hits_fields_old.emtf_host],
    _ones,
  ))
  out_hits = (out_hits_values, out_hits_old[1])
  assert out_hits_old[0].shape[1] == len(sector_hits_fields_old)
  assert out_hits[0].shape[1] == len(sector_hits_fields)

  # Map out_simhits_old -> out_simhits
  _ones = np.ones(out_simhits_old[0].shape[0], dtype=out_simhits_old[0].dtype)
  out_simhits_values = np.column_stack((
    out_simhits_old[0][:, sector_hits_fields_old.emtf_chamber],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_segment],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_phi],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_bend],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_theta],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_theta_alt],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_qual],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_qual_alt],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_time],
    out_simhits_old[0][:, sector_hits_fields_old.zones],
    out_simhits_old[0][:, sector_hits_fields_old.timezones],
    out_simhits_old[0][:, sector_hits_fields_old.fr],
    out_simhits_old[0][:, sector_hits_fields_old.detlayer],
    out_simhits_old[0][:, sector_hits_fields_old.bx],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_site],
    out_simhits_old[0][:, sector_hits_fields_old.emtf_host],
    _ones,
  ))
  out_simhits = (out_simhits_values, out_simhits_old[1])
  assert out_simhits_old[0].shape[1] == len(sector_hits_fields_old)
  assert out_simhits[0].shape[1] == len(sector_hits_fields)

  # Output
  outdict = {
    'out_part': out_part,
    'out_hits_values': out_hits[0],
    'out_hits_row_splits': out_hits[1],
    'out_simhits_values': out_simhits[0],
    'out_simhits_row_splits': out_simhits[1],
  }
  save_np_arrays(new_fname, outdict)
  return

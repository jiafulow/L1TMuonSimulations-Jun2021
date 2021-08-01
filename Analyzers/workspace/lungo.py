"""Data exploration."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtf_algos import *
from emtf_ntuples import *

try:
  import third_party.emtf_tree as emtf_tree
except ImportError:
  raise ImportError(
    'Could not import third_party.emtf_tree. Please run get-third-party.sh first.')


class _BaseAnalysis(object):
  """Abstract base class."""

  def run(self, *args, **kwargs):
    pass


# ______________________________________________________________________________
# Analyses

class DummyAnalysis(_BaseAnalysis):
  """Dummy analysis.

  Description.
  """

  def run(self):
    # Load tree
    tree = load_pgun_test()

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      if verbosity >= 1:
        print('evt {0} has {1} particles, {2} simhits, {3} hits and {4} tracks'.format(
          ievt, len(evt.particles), len(evt.simhits), len(evt.hits), len(evt.tracks)))

      # Particles
      part = evt.particles[0]  # particle gun
      if verbosity >= 1:
        ipart = 0
        print('.. part {0} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}'.format(
          ipart, part.pt, part.eta, part.phi, part.invpt, part.d0))

      # Sim hits
      if verbosity >= 1:
        for isimhit, simhit in enumerate(evt.simhits):
          simhit_endcap = 1 if simhit.z >= 0 else -1
          simhit_sector = get_trigger_sector(simhit.ring, simhit.station, simhit.chamber)
          simhit_bx = 0
          simhit_id = (simhit.subsystem, simhit.station, simhit.ring, get_trigger_endsec(simhit_endcap, simhit_sector),
                       simhit.chamber, simhit.layer, simhit_bx)
          print('.. simhit {0} {1} {2:.3f} {3:.3f}'.format(isimhit, simhit_id, simhit.phi, simhit.theta))

      # Trigger primitives
      if verbosity >= 1:
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.subsystem, hit.station, hit.ring, get_trigger_endsec(hit.endcap, hit.sector),
                    get_trigger_cscneighid(hit.station, hit.subsector, hit.cscid, hit.neighbor), hit.gemdl, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.subsystem == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6} {7}'.format(
            ihit, hit_id, hit.emtf_phi, hit.emtf_bend, hit.emtf_theta1, hit.emtf_theta2, hit.emtf_qual1, hit_sim_tp))

      # Tracks
      if verbosity >= 1:
        for itrk, trk in enumerate(evt.tracks):
          trk_id = (get_trigger_endsec(trk.endcap, trk.sector), trk.bx)
          trk_pt = np.reciprocal(np.abs(trk.model_invpt * np.power(2., -13)))
          print('.. trk {0} {1} {2:.3f} {3} {4} {5} {6}'.format(
            itrk, trk_id, trk_pt, trk.model_phi, trk.model_eta, trk.model_invpt, trk.model_d0))

    # End loop over events
    return


# ______________________________________________________________________________
class ZoneAnalysis(_BaseAnalysis):
  """Find zone boundaries.

  Description.
  """

  def run(self):
    # Overwrite maxevents
    maxevents = -1

    out_hits_metadata = {'subsystem': 0, 'station': 1, 'ring': 2, 'zone': 3, 'emtf_theta': 4}
    out_hits = []

    # __________________________________________________________________________
    # Load tree
    tree = load_pgun_batch(jobid)

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # Particles (pt > 4 GeV)
      part = evt.particles[0]  # particle gun
      if not (part.pt > 4):
        continue

      # Find particle zone
      zone = find_particle_zone(part.eta)

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        out_hits.append([hit.subsystem, hit.station, hit.ring, zone, hit.emtf_theta1])

    # End loop over events
    out_hits = np.asarray(out_hits, dtype=np.int32)

    # __________________________________________________________________________
    # Compute results
    out_hits_subsystem = out_hits[:, out_hits_metadata['subsystem']]
    out_hits_station = out_hits[:, out_hits_metadata['station']]
    out_hits_ring = out_hits[:, out_hits_metadata['ring']]
    out_hits_zone = out_hits[:, out_hits_metadata['zone']]
    out_hits_emtf_theta = out_hits[:, out_hits_metadata['emtf_theta']]

    out_hits_emtf_host = find_emtf_host(out_hits_subsystem, out_hits_station, out_hits_ring)
    assert (out_hits_emtf_host != -99).all()

    print('Find zone boundaries')
    n_zones = len(emtf_eta_bins) - 1
    for zone in range(n_zones):
      sel = (out_hits_zone == zone)
      out_hits_emtf_host_sel = out_hits_emtf_host[sel]
      out_hits_emtf_theta_sel = out_hits_emtf_theta[sel]
      out_hits_emtf_host_uniq = np.unique(out_hits_emtf_host_sel)

      for emtf_host in out_hits_emtf_host_uniq:
        sel_1 = (out_hits_emtf_host_sel == emtf_host)
        out_hits_emtf_theta_sel_1 = out_hits_emtf_theta_sel[sel_1]
        nentries = len(out_hits_emtf_theta_sel_1)
        if nentries > 100:
          p = np.percentile(out_hits_emtf_theta_sel_1, [1, 2, 2.5, 3, 97, 97.5, 98, 99], overwrite_input=True)
          m = [out_hits_emtf_theta_sel_1.min(), out_hits_emtf_theta_sel_1.max()]
          print('{0} {1:03d} {2:5d} {3} {4} {5}'.format(zone, emtf_host, nentries, p[:4], p[4:], m))

    # Done
    return


# ______________________________________________________________________________
class ChamberAnalysis(_BaseAnalysis):
  """Check chamber num of segments.

  Description.
  """

  def run(self):
    # Overwrite maxevents
    maxevents = -1

    out_hits_metadata = {
      'subsystem': 0, 'station': 1, 'sector': 2, 'subsector': 3, 'cscid': 4, 'neighbor': 5, 'bx': 6, 'ievt': 7}
    out_hits = []

    out_simhits_metadata = out_hits_metadata
    out_simhits = []

    # __________________________________________________________________________
    # Load tree
    tree = load_pgun_batch(jobid)

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        ohit = [
          hit.subsystem, hit.station, get_trigger_endsec(hit.endcap, hit.sector),
          hit.subsector, hit.cscid, hit.neighbor, hit.bx, ievt]
        out_hits.append(ohit)

      # Sim hits
      for isimhit, simhit in enumerate(evt.simhits):
        simhit.endcap = 1 if simhit.z >= 0 else -1  #FIXME
        # Special case for ME0 as it is a 20-deg chamber in station 1
        if simhit.subsystem == kME0:
          hack_me0_hit_chamber(simhit)

        simhit.sector = get_trigger_sector(simhit.ring, simhit.station, simhit.chamber)
        simhit.subsector = get_trigger_subsector(simhit.ring, simhit.station, simhit.chamber)
        simhit.cscid = get_trigger_cscid(simhit.ring, simhit.station, simhit.chamber)
        simhit.neighid = get_trigger_neighid(simhit.ring, simhit.station, simhit.chamber)
        simhit.bx = 0
        simhit.neighbor = 0
        osimhit = [
          simhit.subsystem, simhit.station, get_trigger_endsec(simhit.endcap, simhit.sector),
          simhit.subsector, simhit.cscid, simhit.neighbor, simhit.bx, ievt]
        out_simhits.append(osimhit)

        # If neighbor, share simhit with the neighbor sector
        if simhit.neighid == 1:
          simhit.neighbor = 1
          get_next_sector = lambda sector: (sector + 1) if sector != 6 else (sector + 1 - 6)
          simhit.sector = get_next_sector(simhit.sector)
          osimhit = [
            simhit.subsystem, simhit.station, get_trigger_endsec(simhit.endcap, simhit.sector),
            simhit.subsector, simhit.cscid, simhit.neighbor, simhit.bx, ievt]
          out_simhits.append(osimhit)

    # End loop over events
    out_hits = np.asarray(out_hits, dtype=np.int32)
    out_simhits = np.asarray(out_simhits, dtype=np.int32)

    # __________________________________________________________________________
    # Compute results (1)
    out_hits_subsystem = out_hits[:, out_hits_metadata['subsystem']]
    out_hits_station = out_hits[:, out_hits_metadata['station']]
    out_hits_sector = out_hits[:, out_hits_metadata['sector']]
    out_hits_subsector = out_hits[:, out_hits_metadata['subsector']]
    out_hits_cscid = out_hits[:, out_hits_metadata['cscid']]
    out_hits_neighbor = out_hits[:, out_hits_metadata['neighbor']]
    out_hits_bx = out_hits[:, out_hits_metadata['bx']]
    out_hits_ievt = out_hits[:, out_hits_metadata['ievt']]

    out_hits_chambers = find_emtf_chamber(
      out_hits_subsystem, out_hits_station, out_hits_cscid, out_hits_subsector, out_hits_neighbor)
    assert (out_hits_chambers != -99).all()

    sel = (out_hits_bx == 0) & (out_hits_sector == 0)  # check only one bx and one sector
    out_hits_chambers_sel = out_hits_chambers[sel]
    out_hits_ievt_sel = out_hits_ievt[sel]

    n_chambers = out_hits_chambers_sel.max() + 1
    n_events = out_hits_ievt_sel.max() + 1
    hist, _, _ = np.histogram2d(
      out_hits_chambers_sel, out_hits_ievt_sel, bins=(range(n_chambers + 1), range(n_events + 1)))
    chamber_counts = hist.astype(np.int32)

    print('Check chamber num of segments')
    chamber_counts_max = np.max(chamber_counts, axis=-1)
    for chm in range(n_chambers):
      print(chm, chamber_counts_max[chm])

    # Check num of segments
    def chm_check(p):
      if not p:
        print('[WARNING] chm {} has {} segments'.format(chm, cnt))

    for chm in range(n_chambers):
      cnt = chamber_counts_max[chm]
      if chm < 54:  # CSC
        chm_check(cnt <= 2)
      elif chm in {54, 55, 56, 63, 64, 65, 99}:  # GE1/1
        chm_check(cnt <= 8)
      elif chm in {72, 73, 74, 102}:  # GE2/1
        chm_check(cnt <= 8)
      elif chm in {81, 82, 83, 90, 91, 92, 104, 106}:  # RE3,4/1
        chm_check(cnt <= 12)
      elif chm in {57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 100, 101}:  # RE1/2
        chm_check(cnt <= 6)
      elif (54 <= chm < 108):  # RE2,3,4/2
        chm_check(cnt <= 6)
      elif (108 <= chm < 115):  # ME0
        chm_check(cnt <= 25)
      else:
        raise RuntimeError('Unexpected emtf_chamber number: {}'.format(chm))

    # __________________________________________________________________________
    # Compute results (2)
    out_hits_subsystem = out_simhits[:, out_simhits_metadata['subsystem']]
    out_hits_station = out_simhits[:, out_simhits_metadata['station']]
    out_hits_sector = out_simhits[:, out_simhits_metadata['sector']]
    out_hits_subsector = out_simhits[:, out_simhits_metadata['subsector']]
    out_hits_cscid = out_simhits[:, out_simhits_metadata['cscid']]
    out_hits_neighbor = out_simhits[:, out_simhits_metadata['neighbor']]
    out_hits_bx = out_simhits[:, out_simhits_metadata['bx']]
    out_hits_ievt = out_simhits[:, out_simhits_metadata['ievt']]

    out_hits_chambers = find_emtf_chamber(
      out_hits_subsystem, out_hits_station, out_hits_cscid, out_hits_subsector, out_hits_neighbor)
    assert (out_hits_chambers != -99).all()

    sel = (out_hits_bx == 0) & (out_hits_sector == 0)  # check only one bx and one sector
    out_hits_chambers_sel = out_hits_chambers[sel]
    out_hits_ievt_sel = out_hits_ievt[sel]

    n_chambers = out_hits_chambers_sel.max() + 1
    n_events = out_hits_ievt_sel.max() + 1
    hist, _, _ = np.histogram2d(
      out_hits_chambers_sel, out_hits_ievt_sel, bins=(range(n_chambers + 1), range(n_events + 1)))
    chamber_counts = hist.astype(np.int32)

    print('Check chamber num of simhits')
    chamber_counts_max = np.max(chamber_counts, axis=-1)
    for chm in range(n_chambers):
      print(chm, chamber_counts_max[chm])

    # Done
    return


# ______________________________________________________________________________
# Main

import os
import sys
from datetime import datetime

# Era (pick one)
era = 'default'  # phase 2
#era = 'run3'

# Analysis mode (pick one)
analysis = 'dummy'
#analysis = 'zone'
#analysis = 'chamber'

# Job id (pick an integer)
jobid = 0

# Max num of events (-1 means all events)
maxevents = 100

# Verbosity (0 means quiet)
verbosity = 1

# Condor or not
# if 'CONDOR_EXEC' is defined, override the 3 arguments (era, analysis, jobid)
use_condor = ('CONDOR_EXEC' in os.environ)
if use_condor:
  nargs = 3
  if len(sys.argv) != (nargs + 1):
    raise RuntimeError('Expected num of arguments: {}'.format(nargs))
  era = sys.argv[1]
  analysis = sys.argv[2]
  jobid = int(sys.argv[3])
  maxevents = -1
  verbosity = 0

# Logger
logger = emtf_tree.get_logger()


# Decorator
def app_decorator(fn):
  def wrapper(*args, **kwargs):
    # Begin
    start_time = datetime.now()
    logger.info('Using cmssw     : {}'.format(os.environ['CMSSW_VERSION']))
    logger.info('Using condor    : {}'.format(use_condor))
    logger.info('Using era       : {}'.format(era))
    logger.info('Using analysis  : {}'.format(analysis))
    logger.info('Using jobid     : {}'.format(jobid))
    logger.info('Using maxevents : {}'.format(maxevents))
    # Run
    fn(*args, **kwargs)
    # End
    stop_time = datetime.now()
    logger.info('Elapsed time    : {}'.format(stop_time - start_time))
    return

  return wrapper


# App
@app_decorator
def app():
  if analysis == 'dummy':
    DummyAnalysis().run()
  elif analysis == 'zone':
    ZoneAnalysis().run()
  elif analysis == 'chamber':
    ChamberAnalysis().run()
  else:
    raise RuntimeError('Could not recognize analysis: {}'.format(analysis))


# Finally
if __name__ == '__main__':
  app()

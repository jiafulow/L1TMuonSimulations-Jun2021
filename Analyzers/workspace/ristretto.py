"""Data preparation."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import collections
from multiprocessing import Process, Queue, cpu_count

from emtf_algos import *
from emtf_ntuples import *

try:
  import third_party.emtf_tree as emtf_tree
except ImportError:
  raise ImportError(
      'Could not import third_party.emtf_tree. Please run get-third-party.sh first.')

try:
  import third_party.emtf_nnet as emtf_nnet
except ImportError:
  raise ImportError(
      'Could not import third_party.emtf_nnet. Please run get-third-party.sh first.')


class _BaseAnalysis(object):
  """Abstract base class."""

  def run(self, *args, **kwargs):
    pass


# ______________________________________________________________________________
# Functions

class SectorRanking(object):
  """Select the best sector for a single muon event."""

  def __init__(self):
    self.sectors = np.zeros(num_emtf_sectors, dtype=np.int32)

  def reset(self):
    self.sectors.fill(0)

  def add(self, hit):
    emtf_host = hit.emtf_host
    endsec = get_trigger_endsec(hit.endcap, hit.sector)
    assert (endsec < num_emtf_sectors)

    valid_flag = np.zeros(8, dtype=np.bool)
    valid_flag[0] = (emtf_host == 18)                   # ME0
    valid_flag[1] = (emtf_host == 0)                    # ME1/1
    valid_flag[2] = (emtf_host in (1, 2))               # ME1/2, ME1/3
    valid_flag[3] = (emtf_host in (3, 4))               # ME2/1, ME2/2
    valid_flag[4] = (emtf_host in (5, 6))               # ME3/1, ME3/2
    valid_flag[5] = (emtf_host in (7, 8))               # ME4/1, ME4/2
    valid_flag[6] = (emtf_host in (9, 10, 11, 12, 13))  # GE1/1, RE1/2, RE1/3, GE2/1, RE2/2
    valid_flag[7] = (emtf_host in (14, 15, 16, 17))     # RE3/1, RE3/2, RE4/1, RE4/2
    rank = np.packbits(valid_flag)                      # pack 8 booleans into an uint8
    self.sectors[endsec] |= rank

  def best_sector(self):
    best_sector = np.argmax(self.sectors)
    best_sector_rank = self.sectors[best_sector]
    return (best_sector, best_sector_rank)


class ZoneRanking(object):
  """Select the best zone for a single muon event."""

  def __init__(self):
    self.zones = np.zeros(num_emtf_zones, dtype=np.int32)

  def reset(self):
    self.zones.fill(0)

  def add(self, hit):
    seg_zones = hit.zones
    assert (seg_zones < (1 << num_emtf_zones))

    self.zones[0] += bool(seg_zones & (1 << (num_emtf_zones - 1 - 0)))
    self.zones[1] += bool(seg_zones & (1 << (num_emtf_zones - 1 - 1)))
    self.zones[2] += bool(seg_zones & (1 << (num_emtf_zones - 1 - 2)))

  def best_zone(self):
    best_zone = np.argmax(self.zones)
    best_zone_rank = self.zones[best_zone]
    return (best_zone, best_zone_rank)


class ChamberCouncil(object):
  def __init__(self, is_sim=False):
    self.chambers = collections.defaultdict(list)
    self.is_sim = is_sim

  def reset(self):
    self.chambers.clear()

  def add(self, hit):
    emtf_chamber = hit.emtf_chamber
    self.chambers[emtf_chamber].append(hit)

  def _get_hits_from_chambers_sim(self):
    hits_info = []
    sorted_keys = sorted(self.chambers.keys())

    for k in sorted_keys:
      tmp_hits = self.chambers[k]

      # If more than one hit, pick the median in phi.
      ind = 0
      if len(tmp_hits) > 1:
        phi_median = np.median([hit.phi for hit in tmp_hits])
        ind = np.argmin(np.abs(hit.phi - phi_median))

      hit = tmp_hits[ind]
      hit_info = get_hit_info(hit)
      hits_info.append(hit_info)

    hits_info = np.asarray(hits_info, dtype=np.int32)
    return hits_info

  def _get_hits_from_chambers(self):
    hits_info = []
    sorted_keys = sorted(self.chambers.keys())

    for k in sorted_keys:
      tmp_hits = self.chambers[k]

      # Quick truncation
      _len = 2
      if tmp_hits[0].subsystem == kCSC:
        _len = 2
      elif tmp_hits[0].subsystem == kRPC:
        if tmp_hits[0].station >= 3 and tmp_hits[0].ring == 1:  # is_irpc
          _len = 12
        else:
          _len = 6
      elif tmp_hits[0].subsystem == kGEM:
        _len = 8
      elif tmp_hits[0].subsystem == kME0:
        _len = 25

      for _, hit in zip(range(_len), tmp_hits):
        hit_info = get_hit_info(hit)
        hits_info.append(hit_info)

    hits_info = np.asarray(hits_info, dtype=np.int32)
    return hits_info

  def get_hits(self):
    if self.is_sim:
      hits_info = self._get_hits_from_chambers_sim()
    else:
      hits_info = self._get_hits_from_chambers()
    return hits_info


def get_part_info(part, best_sector=0, best_zone=0):
  # 10 variables
  part_info = (
    part.invpt, part.eta, part.phi, part.vx, part.vy, part.vz, part.d0,
    part.bx, best_sector, best_zone,
  )
  part_info = np.asarray(part_info, dtype=np.float32)
  return part_info


def get_track_info(trk):
  # 14 variables
  if trk is not None:
    trk_patt = 0  #FIXME
    trk_col = 0
    trk_zone = 0
    trk_tzone = 0
    track_info = (
      trk.model_invpt, trk.model_phi, trk.model_eta, trk.model_d0,
      trk.model_z0, trk.model_beta, trk.model_qual, trk_patt,
      trk_col, trk_zone, trk_tzone, trk.emtf_mode_v1,
      trk.emtf_mode_v2, trk.valid,
    )
  else:
    track_info = (0,) * 14
  track_info = np.asarray(track_info, dtype=np.int32)
  return track_info


def get_hit_info(hit):
  # 17 variables
  # (some of them are not available for simhits)
  hit_emtf_segment = getattr(hit, 'emtf_segment', 0)
  hit_emtf_bend = getattr(hit, 'emtf_bend', 0)
  hit_emtf_theta2 = getattr(hit, 'emtf_theta2', 0)
  hit_emtf_qual1 = getattr(hit, 'emtf_qual1', 0)
  hit_emtf_qual2 = getattr(hit, 'emtf_qual2', 0)
  hit_emtf_time = getattr(hit, 'emtf_time', 0)
  hit_cscfr = getattr(hit, 'cscfr', 0)
  hit_gemdl = getattr(hit, 'gemdl', 0)
  hit_valid = getattr(hit, 'valid', 1)
  hit_info = (
    hit.emtf_chamber, hit_emtf_segment, hit.emtf_phi, hit_emtf_bend,
    hit.emtf_theta1, hit_emtf_theta2, hit_emtf_qual1, hit_emtf_qual2,
    hit_emtf_time, hit.zones, hit.timezones, hit_cscfr, hit_gemdl,
    hit.bx, hit.emtf_site, hit.emtf_host, hit_valid,
  )
  #hit_info = np.asarray(hit_info, dtype=np.int32)  # unnecessary
  return hit_info


# ______________________________________________________________________________
# Analyses

class SignalAnalysis(_BaseAnalysis):
  """Prepare signal data used for training.

  Description.
  """

  @classmethod
  def worker_fn(cls, signal, task_queue, done_queue):
    tree = load_tree(task_queue, use_multiprocessing=True)

    if signal == 'prompt':
      unconstrained = 0
    else:
      unconstrained = 1

    out_part = []
    out_track = []
    out_hits = []
    out_simhits = []

    sectors = SectorRanking()
    zones = ZoneRanking()
    chambers = ChamberCouncil()
    chambers_sim = ChamberCouncil(is_sim=True)

    # Loop over events
    for ievt, evt in enumerate(tree):
      # Particles
      part = evt.particles[0]  # particle gun

      # First, determine the best sector and best zone (using trigger primitives)

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        hit.emtf_site = find_emtf_site(hit.subsystem, hit.station, hit.ring)  #FIXME
        hit.emtf_host = find_emtf_host(hit.subsystem, hit.station, hit.ring)
        hit.zones = find_emtf_zones(hit.emtf_host, hit.emtf_theta1, hit.emtf_theta2)
        hit.timezones = find_emtf_timezones(hit.emtf_host, hit.bx)
        assert (0 <= hit.emtf_site < num_emtf_sites)
        assert (0 <= hit.emtf_host < num_emtf_hosts)
        assert (0 <= hit.emtf_chamber < num_emtf_chambers)
        assert (hit.emtf_phi > 0)
        assert (hit.emtf_theta1 > 0)
        assert (hit.emtf_theta2 >= 0)

        sectors.add(hit)
        zones.add(hit)

      best_sector, best_sector_rank = sectors.best_sector()
      best_zone, best_zone_rank = zones.best_zone()

      # Second, fill the chambers with trigger primitives and sim hits

      # Trigger primitives
      for ihit, hit in enumerate(evt.hits):
        if get_trigger_endsec(hit.endcap, hit.sector) == best_sector:
          chambers.add(hit)

      # Sim hits
      for isimhit, simhit in enumerate(evt.simhits):
        # Special case for ME0 as it is a 20-deg chamber in station 1
        if simhit.subsystem == kME0:
          hack_me0_hit_chamber(simhit)

        # Add some variables
        simhit.sector = get_trigger_sector(simhit.ring, simhit.station, simhit.chamber)
        simhit.subsector = get_trigger_subsector(simhit.ring, simhit.station, simhit.chamber)
        simhit.cscid = get_trigger_cscid(simhit.ring, simhit.station, simhit.chamber)
        simhit.neighid = get_trigger_neighid(simhit.ring, simhit.station, simhit.chamber)
        simhit.bx = 0
        simhit.neighbor = 0

        simhit.emtf_site = find_emtf_site(simhit.subsystem, simhit.station, simhit.ring)
        simhit.emtf_host = find_emtf_host(simhit.subsystem, simhit.station, simhit.ring)
        simhit.emtf_chamber = find_emtf_chamber(
            simhit.subsystem, simhit.station, simhit.cscid, simhit.subsector, simhit.neighbor)
        simhit.emtf_phi = calc_phi_loc_int(simhit.phi, (best_sector % 6) + 1)
        simhit.emtf_theta1 = calc_theta_int(simhit.theta, 1 if best_sector < 6 else -1)
        simhit.zones = find_emtf_zones(simhit.emtf_host, simhit.emtf_theta1)
        simhit.timezones = find_emtf_timezones(simhit.emtf_host, simhit.bx)

        # If neighbor, share simhit with the neighbor sector
        if get_trigger_endsec(simhit.endcap, simhit.sector) == best_sector:
          chambers_sim.add(simhit)
        elif simhit.neighid == 1:
          simhit.neighbor = 1
          get_next_sector = lambda sector: (sector + 1) if sector != 6 else (sector + 1 - 6)
          simhit.sector = get_next_sector(simhit.sector)
          simhit.emtf_chamber = find_emtf_chamber(
              simhit.subsystem, simhit.station, simhit.cscid, simhit.subsector, simhit.neighbor)
          if get_trigger_endsec(simhit.endcap, simhit.sector) == best_sector:
            chambers_sim.add(simhit)

      # Third, get the best track

      # Tracks
      best_track = None
      for itrk, trk in enumerate(evt.tracks):
        if (trk.valid and trk.unconstrained == unconstrained and
            get_trigger_endsec(trk.endcap, trk.sector) == best_sector):
          best_track = trk
          break

      # Fourth, check for at least 2 stations (using sim hits)

      ievt_part = get_part_info(part, best_sector=best_sector, best_zone=best_zone)
      ievt_track = get_track_info(best_track)

      def require_two_stations():
        subsystems = {
            hit.subsystem for tmp_hits in six.itervalues(chambers_sim.chambers) for hit in tmp_hits}
        stations = {
            hit.station for tmp_hits in six.itervalues(chambers_sim.chambers) for hit in tmp_hits}
        min_station = min(stations) if len(stations) else 5
        max_station = max(stations) if len(stations) else 0
        both_me0_me1 = (kME0 in subsystems) and (kCSC in subsystems)
        passed = ((min_station <= 1 and max_station >= 2) or (min_station == 2 and max_station >= 3) or
                  (max_station == 1 and both_me0_me1))
        return passed

      if require_two_stations():
        ievt_hits = chambers.get_hits()
        ievt_simhits = chambers_sim.get_hits()
      else:
        ievt_hits = np.array([], dtype=np.int32)
        ievt_simhits = np.array([], dtype=np.int32)

      # Finally, add to output collections
      out_part.append(ievt_part)
      out_track.append(ievt_track)
      out_hits.append(ievt_hits)
      out_simhits.append(ievt_simhits)

      # Reset before the next event
      sectors.reset()
      zones.reset()
      chambers.reset()
      chambers_sim.reset()

      # Debug
      if verbosity >= 2 and ievt < 100:
        # Event
        print('evt {0} has {1} particles, {2} simhits, {3} hits and {4} tracks'.format(
            ievt, len(evt.particles), len(evt.simhits), len(evt.hits), len(evt.tracks)))
        # Particles
        ipart = 0
        print('.. part {0} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}'.format(
            ipart, part.pt, part.eta, part.phi, part.invpt, part.d0))
        # Sim hits
        for isimhit, simhit in enumerate(evt.simhits):
          simhit_id = (simhit.subsystem, simhit.station, simhit.ring, get_trigger_endsec(simhit.endcap, simhit.sector),
                       simhit.chamber, simhit.layer, simhit.bx)
          print('.. simhit {0} {1} {2:.3f} {3:.3f}'.format(isimhit, simhit_id, simhit.phi, simhit.theta))
        # Trigger primitives
        for ihit, hit in enumerate(evt.hits):
          hit_id = (hit.subsystem, hit.station, hit.ring, get_trigger_endsec(hit.endcap, hit.sector),
                    get_trigger_cscneighid(hit.station, hit.subsector, hit.cscid, hit.neighbor), hit.gemdl, hit.bx)
          hit_sim_tp = hit.sim_tp1
          if (hit.subsystem == kCSC) and (hit_sim_tp != hit.sim_tp2):
            hit_sim_tp = -1
          print('.. hit {0} {1} {2} {3} {4} {5} {6} {7}'.format(
              ihit, hit_id, hit.emtf_phi, hit.emtf_bend, hit.emtf_theta1, hit.emtf_theta2, hit.emtf_qual1, hit_sim_tp))
        # Tracks
        for itrk, trk in enumerate(evt.tracks):
          trk_id = (get_trigger_endsec(trk.endcap, trk.sector), trk.bx)
          trk_pt = np.reciprocal(np.abs(trk.model_invpt * np.power(2., -13)))
          print('.. trk {0} {1} {2:.3f} {3} {4} {5} {6} {7}'.format(
              itrk, trk_id, trk_pt, trk.model_phi, trk.model_eta, trk.model_invpt, trk.model_d0, trk.unconstrained))
        # Output
        print('best_sector: {0} rank: {1} best_zone: {2} rank: {3}'.format(
            best_sector, best_sector_rank, best_zone, best_zone_rank))
        with np.printoptions(linewidth=140, threshold=1000):
          print('hits:')
          print(ievt_hits)
          print('simhits:')
          print(ievt_simhits)

    # End loop over events

    # __________________________________________________________________________
    # Output
    out_part = np.asarray(out_part)
    out_track = np.asarray(out_track)
    out_hits = emtf_nnet.ragged.create_ragged_array(out_hits)
    out_simhits = emtf_nnet.ragged.create_ragged_array(out_simhits)
    outdict = {
      'out_part': out_part,
      'out_track': out_track,
      'out_hits_values': out_hits.values,
      'out_hits_row_splits': out_hits.row_splits,
      'out_simhits_values': out_simhits.values,
      'out_simhits_row_splits': out_simhits.row_splits,
    }
    done_queue.put(outdict)
    return

  @classmethod
  def union_fn(cls, num_workers, done_queue):
    outdict = [done_queue.get() for _ in range(num_workers)]
    outdict = stack_np_arrays(outdict)
    if use_condor:
      outfile = '{0}_{1}.npz'.format(analysis, jobid)
    else:
      outfile = '{0}.npz'.format(analysis)
    save_np_arrays(outfile, outdict)
    return

  @classmethod
  def multiprocessing_fn(cls, signal='prompt', num_workers=8, num_files=None):
    if num_workers > cpu_count():
      num_workers = cpu_count()

    if signal == 'prompt':
      dataset = SingleMuon()
    elif signal == 'displaced':
      dataset = SingleMuonDisplaced()
    else:
      raise ValueError('Unexpected signal: {0}'.format(signal))

    if num_files is None:
      infiles = dataset[:]
    else:
      infiles = dataset[:num_files]
    task_queue = Queue()
    done_queue = Queue()
    sentinels = []

    for _ in range(num_workers):
      Process(target=cls.worker_fn, args=(signal, task_queue, done_queue)).start()
      sentinels.append(None)
    for infile in infiles + sentinels:
      task_queue.put(infile)
    result = cls.union_fn(num_workers, done_queue)
    return result

  def run(self, signal='prompt'):
    self.multiprocessing_fn(signal=signal)
    return


# ______________________________________________________________________________
class BkgndAnalysis(_BaseAnalysis):
  """Prepare background data used for training.

  Description.
  """

  def run(self, bkgnd='minbias'):
    pass


# ______________________________________________________________________________
# Main

import os
import sys
from datetime import datetime

# Era (pick one)
era = 'default'  # phase 2
#era = 'run3'

# Analysis mode (pick one)
analysis = 'signal'
#analysis = 'signal_dxy'
#analysis = 'bkgnd'

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

# Slurm or not
use_slurm = ('SLURM_JOB_ID' in os.environ)

# Logger
logger = emtf_tree.get_logger()


# Decorator
def app_decorator(fn):
  def wrapper(*args, **kwargs):
    # Begin
    start_time = datetime.now()
    logger.info('Using cmssw     : {}'.format(os.environ['CMSSW_VERSION']))
    logger.info('Using condor    : {}'.format(use_condor))
    logger.info('Using slurm     : {}'.format(use_slurm))
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
  if analysis == 'signal':
    SignalAnalysis().run(signal='prompt')
  elif analysis == 'signal_dxy':
    SignalAnalysis().run(signal='displaced')
  elif analysis == 'bkgnd':
    BkgndAnalysis().run(bkgnd='minbias')
  else:
    raise RuntimeError('Could not recognize analysis: {}'.format(analysis))


# Finally
if __name__ == '__main__':
  app()

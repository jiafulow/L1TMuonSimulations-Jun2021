"""Data analysis."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from emtf_algos import *
from emtf_ntuples import *
from emtf_tree import get_logger


class _BaseAnalysis(object):
  """Abstract base class."""

  def run(self, *args, **kwargs):
    pass


# ______________________________________________________________________________
# Functions

def find_sw_pt(trk):
  eps_trk_pt = np.power(2., -4)  # 0.0625
  return eps_trk_pt * trk.emtf_pt


def find_sw_phi(trk):
  return calc_phi_glob_deg_from_loc(calc_phi_loc_deg_from_int(trk.model_phi), trk.sector)


def find_sw_eta(trk):
  return calc_eta_from_theta_deg(calc_theta_deg_from_int(trk.model_eta)) * trk.endcap


def adjust_model_qual(trk):
  # Zone 2: emtf_theta >= 54 with ME1/2 or RE1/2 hit
  is_in_zone2 = (trk.model_eta >= 54) and ((trk.hitmode & (1 << 1)) or (trk.hitmode & (1 << 5)))
  # Promote Zone 2 tracks with emtf_mode_v2 >= 12 and model_qual >= 24 to model_qual >= 48
  if is_in_zone2 and (trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 24):
    model_qual = max(trk.model_qual, 48)
  else:
    model_qual = trk.model_qual
  return model_qual


def fill_highest_qual(tracks, select_fn, lst):
  highest = np.maximum.reduce(
      [trk.model_qual for trk in tracks if select_fn(trk)],
      initial=0)
  lst.append(highest)


def fill_highest_pt(tracks, select_fn, lst):
  highest = np.maximum.reduce(
      [trk.sw_pt for trk in tracks if select_fn(trk)],
      initial=-999999.)
  lst.append(highest)


def fill_efficiency_pt(part, tracks, select_part_fn, select_trk_fn, lst):
  if select_part_fn(part):
    trigger = any(select_trk_fn(trk) for trk in tracks)
    if trigger:
      lst.append(part.pt)


def fill_efficiency_eta(part, tracks, select_part_fn, select_trk_fn, lst):
  if select_part_fn(part):
    trigger = any(select_trk_fn(trk) for trk in tracks)
    if trigger:
      lst.append(part.eta)


# ______________________________________________________________________________
# Analyses

class RatesAnalysis(_BaseAnalysis):
  """Make performance plots using full-sim Min-Bias dataset.

  Description.
  """

  def run(self, pileup=200):
    # Load tree
    dataset = SingleNeutrinoPU200()
    infiles = dataset[:]
    tree = load_tree(infiles, load_tracks=True)

    rate_vs_qual = []
    rate_vs_pt = []
    rate_vs_pt_1 = []
    rate_vs_pt_2 = []
    rate_vs_pt_3 = []

    rate_vs_qual_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 20.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 0)))
    rate_vs_pt_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 0.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    rate_vs_pt_sel_1 = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 1.6) and (trk.sw_pt >= 0.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    rate_vs_pt_sel_2 = lambda trk: (
        (trk.valid != 0) and (1.6 <= abs(trk.sw_eta) <= 2.0) and (trk.sw_pt >= 0.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    rate_vs_pt_sel_3 = lambda trk: (
        (trk.valid != 0) and (2.0 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 0.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # Find sw_pt, sw_phi, sw_eta
      tracks = evt.tracks
      for trk in tracks:
        trk.sw_pt = find_sw_pt(trk)
        trk.sw_phi = find_sw_phi(trk)
        trk.sw_eta = find_sw_eta(trk)
        trk.model_qual = adjust_model_qual(trk)

      # Fill
      fill_highest_qual(tracks, rate_vs_qual_sel, rate_vs_qual)
      fill_highest_pt(tracks, rate_vs_pt_sel, rate_vs_pt)
      fill_highest_pt(tracks, rate_vs_pt_sel_1, rate_vs_pt_1)
      fill_highest_pt(tracks, rate_vs_pt_sel_2, rate_vs_pt_2)
      fill_highest_pt(tracks, rate_vs_pt_sel_3, rate_vs_pt_3)

    # End loop over events

    # __________________________________________________________________________
    # Output
    rate_vs_qual = np.asarray(rate_vs_qual, dtype=np.int32)
    rate_vs_pt = np.asarray(rate_vs_pt, dtype=np.float32)
    rate_vs_pt_1 = np.asarray(rate_vs_pt_1, dtype=np.float32)
    rate_vs_pt_2 = np.asarray(rate_vs_pt_2, dtype=np.float32)
    rate_vs_pt_3 = np.asarray(rate_vs_pt_3, dtype=np.float32)

    outdict = {
      'rate_vs_qual': rate_vs_qual,
      'rate_vs_pt': rate_vs_pt,
      'rate_vs_pt_1': rate_vs_pt_1,
      'rate_vs_pt_2': rate_vs_pt_2,
      'rate_vs_pt_3': rate_vs_pt_3,
    }
    if use_condor:
      outfile = '{}_{}.npz'.format(analysis, jobid)
    else:
      outfile = '{}.npz'.format(analysis)
    save_np_arrays(outfile, outdict)
    logger.info('Wrote to {}'.format(outfile))
    return


# ______________________________________________________________________________
class EffieAnalysis(_BaseAnalysis):
  """Make performance plots using full-sim mu+ mu- pair dataset.

  Description.
  """

  def run(self, pileup=200):
    # Load tree
    dataset = DoubleMuonPU200()
    infiles = dataset[:]
    tree = load_tree(infiles, load_tracks=True, load_particles=True)

    eff_vs_genpt_all = []
    eff_vs_genpt_l1pt5 = []
    eff_vs_genpt_l1pt10 = []
    eff_vs_genpt_l1pt20 = []
    eff_vs_genpt_l1pt30 = []
    eff_vs_genpt_l1pt40 = []
    eff_vs_genpt_l1pt50 = []

    eff_vs_geneta_lowpt_all = []
    eff_vs_geneta_lowpt_l1pt5 = []

    eff_vs_geneta_highpt_all = []
    eff_vs_geneta_highpt_l1pt20 = []

    eff_vs_genpt_part_sel = lambda part: (
        (part.bx == 0) and (1.24 <= abs(part.eta) <= 2.4))
    eff_vs_geneta_lowpt_part_sel = lambda part: (
        (part.bx == 0) and (5 <= part.pt < 20.))
    eff_vs_geneta_highpt_part_sel = lambda part: (
        (part.bx == 0) and (part.pt >= 20.))

    eff_vs_genpt_l1pt5_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 5.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_genpt_l1pt10_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 10.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_genpt_l1pt20_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 20.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_genpt_l1pt30_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 30.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_genpt_l1pt40_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 40.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_genpt_l1pt50_sel = lambda trk: (
        (trk.valid != 0) and (1.2 <= abs(trk.sw_eta) <= 2.4) and (trk.sw_pt >= 50.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_geneta_lowpt_l1pt5_sel = lambda trk: (
        (trk.valid != 0) and (0. <= abs(trk.sw_eta) <= 3.) and (trk.sw_pt >= 5.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))
    eff_vs_geneta_highpt_l1pt20_sel = lambda trk: (
        (trk.valid != 0) and (0. <= abs(trk.sw_eta) <= 3.) and (trk.sw_pt >= 20.) and
        ((trk.emtf_mode_v2 >= 12) and (trk.model_qual >= 48)))

    true_fn = lambda x: True

    # Loop over events
    for ievt, evt in enumerate(tree):
      if maxevents != -1 and ievt == maxevents:
        break

      # Expected particles in the DoubleMuon dataset
      if len(evt.particles) < 2:
        continue

      particles = evt.particles[:2]
      assert ((abs(particles[0].pdgid) == 13) and (particles[0].status == 1))
      assert ((abs(particles[1].pdgid) == 13) and (particles[1].status == 1))

      # Loop over particles
      for part in particles:
        endcap = 1 if (part.eta >= 0.) else -1

        # Find sw_pt, sw_phi, sw_eta
        tracks = list(filter(lambda trk: (trk.endcap * endcap) > 0, evt.tracks))
        for trk in tracks:
          trk.sw_pt = find_sw_pt(trk)
          trk.sw_phi = find_sw_phi(trk)
          trk.sw_eta = find_sw_eta(trk)
          trk.model_qual = adjust_model_qual(trk)

        # Fill
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, true_fn,
                           eff_vs_genpt_all)
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, eff_vs_genpt_l1pt5_sel,
                           eff_vs_genpt_l1pt5)
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, eff_vs_genpt_l1pt10_sel,
                           eff_vs_genpt_l1pt10)
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, eff_vs_genpt_l1pt20_sel,
                           eff_vs_genpt_l1pt20)
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, eff_vs_genpt_l1pt30_sel,
                           eff_vs_genpt_l1pt30)
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, eff_vs_genpt_l1pt40_sel,
                           eff_vs_genpt_l1pt40)
        fill_efficiency_pt(part, tracks, eff_vs_genpt_part_sel, eff_vs_genpt_l1pt50_sel,
                           eff_vs_genpt_l1pt50)

        fill_efficiency_eta(part, tracks, eff_vs_geneta_lowpt_part_sel, true_fn,
                            eff_vs_geneta_lowpt_all)
        fill_efficiency_eta(part, tracks, eff_vs_geneta_lowpt_part_sel, eff_vs_geneta_lowpt_l1pt5_sel,
                            eff_vs_geneta_lowpt_l1pt5)
        fill_efficiency_eta(part, tracks, eff_vs_geneta_highpt_part_sel, true_fn,
                            eff_vs_geneta_highpt_all)
        fill_efficiency_eta(part, tracks, eff_vs_geneta_highpt_part_sel, eff_vs_geneta_highpt_l1pt20_sel,
                            eff_vs_geneta_highpt_l1pt20)

    # End loop over events

    # __________________________________________________________________________
    # Output
    eff_vs_genpt_all = np.asarray(eff_vs_genpt_all, dtype=np.float32)
    eff_vs_genpt_l1pt5 = np.asarray(eff_vs_genpt_l1pt5, dtype=np.float32)
    eff_vs_genpt_l1pt10 = np.asarray(eff_vs_genpt_l1pt10, dtype=np.float32)
    eff_vs_genpt_l1pt20 = np.asarray(eff_vs_genpt_l1pt20, dtype=np.float32)
    eff_vs_genpt_l1pt30 = np.asarray(eff_vs_genpt_l1pt30, dtype=np.float32)
    eff_vs_genpt_l1pt40 = np.asarray(eff_vs_genpt_l1pt40, dtype=np.float32)
    eff_vs_genpt_l1pt50 = np.asarray(eff_vs_genpt_l1pt50, dtype=np.float32)
    eff_vs_geneta_lowpt_all = np.asarray(eff_vs_geneta_lowpt_all, dtype=np.float32)
    eff_vs_geneta_lowpt_l1pt5 = np.asarray(eff_vs_geneta_lowpt_l1pt5, dtype=np.float32)
    eff_vs_geneta_highpt_all = np.asarray(eff_vs_geneta_highpt_all, dtype=np.float32)
    eff_vs_geneta_highpt_l1pt20 = np.asarray(eff_vs_geneta_highpt_l1pt20, dtype=np.float32)

    outdict = {
      'eff_vs_genpt_all': eff_vs_genpt_all,
      'eff_vs_genpt_l1pt5': eff_vs_genpt_l1pt5,
      'eff_vs_genpt_l1pt10': eff_vs_genpt_l1pt10,
      'eff_vs_genpt_l1pt20': eff_vs_genpt_l1pt20,
      'eff_vs_genpt_l1pt30': eff_vs_genpt_l1pt30,
      'eff_vs_genpt_l1pt40': eff_vs_genpt_l1pt40,
      'eff_vs_genpt_l1pt50': eff_vs_genpt_l1pt50,
      'eff_vs_geneta_lowpt_all': eff_vs_geneta_lowpt_all,
      'eff_vs_geneta_lowpt_l1pt5': eff_vs_geneta_lowpt_l1pt5,
      'eff_vs_geneta_highpt_all': eff_vs_geneta_highpt_all,
      'eff_vs_geneta_highpt_l1pt20': eff_vs_geneta_highpt_l1pt20,
    }
    if use_condor:
      outfile = '{}_{}.npz'.format(analysis, jobid)
    else:
      outfile = '{}.npz'.format(analysis)
    save_np_arrays(outfile, outdict)
    logger.info('Wrote to {}'.format(outfile))
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
analysis = 'rates'
#analysis = 'effie'

# Job id (pick an integer)
jobid = 0

# Max num of events (-1 means all events)
maxevents = 200000

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
logger = get_logger()


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
  if analysis == 'rates':
    RatesAnalysis().run()
  elif analysis == 'effie':
    EffieAnalysis().run()
  else:
    raise RuntimeError('Could not recognize analysis: {}'.format(analysis))


# Finally
if __name__ == '__main__':
  app()

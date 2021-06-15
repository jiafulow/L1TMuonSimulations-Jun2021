#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

"""Print cmsDriver options as YAML string

Usage:

./cmsDriverOptionsToYaml.py Configuration/GenProduction/python/TSG-Phase2HLTTDRWinter20GS-00253-fragment.py \
  --python_filename TSG-Phase2HLTTDRWinter20GS-00253_1_cfg.py \
  --eventcontent RAWSIM \
  --customise Configuration/DataProcessing/Utils.addMonitoring \
  --datatier GEN-SIM \
  --fileout file:TSG-Phase2HLTTDRWinter20GS-00253.root \
  --conditions 110X_mcRun4_realistic_v3 \
  --beamspot HLLHC14TeV \
  --step GEN,SIM \
  --geometry Extended2026D49 \
  --era Phase2C9 \
  --no_exec \
  --mc \
  -n 10000

(Do not use '--option=value'. Always use '--option value'.)
"""


# Reference:
#     https://github.com/cms-sw/cmssw/blob/master/Configuration/Applications/scripts/cmsDriver.py
#     https://github.com/cms-sw/cmssw/blob/master/Configuration/Applications/python/cmsDriverOptions.py

import sys
import yaml
from Configuration.Applications.Options import parser, threeValued


def OptionsFromCommandLine():
  options = OptionsFromItems2(sys.argv[1:])
  ## memorize the command line arguments
  # options.arguments = reduce(lambda x, y: x+' '+y, sys.argv[1:])
  return options


def OptionsFromItems(items):
  # Handle three-valued options
  for (index, item) in enumerate(items):
    for (opt, value) in threeValued:
      if (str(item) in opt) and (index == len(items) - 1 or items[index + 1].startswith('-')):
        items.insert(index + 1, value)

  (options, args) = parser.parse_args(items)
  return (options, args)


def OptionsFromItems2(items):
  (options, args) = OptionsFromItems(items)

  # Insert '--args'
  if args:
    items.append('--args')
    items.append(args[0])

  # Replace '-s', '-n', '-o'
  oneDashed = [('-s', '--step'), ('-n', '--number'), ('-o', '--number_out')]
  for (index, item) in enumerate(items):
    for (opt, opt2) in oneDashed:
      if str(item) in opt:
        items[index] = opt2

  # Gather keys and values
  keys = []
  values = []
  prev_arg_is_key = False
  curr_arg_is_key = False

  for (index, item) in enumerate(items):
    if item.startswith('--'):
      curr_arg_is_key = True

    if curr_arg_is_key:
      if prev_arg_is_key:
        values.append(None)
      keys.append(item[2:])  # strip '--' from item
    else:
      if prev_arg_is_key:
        values.append(item)

    prev_arg_is_key = curr_arg_is_key
    curr_arg_is_key = False

    if index == len(items) - 1:
      if prev_arg_is_key:
        values.append(None)

  options = dict(zip(keys, values))
  return options


def run():
  options = OptionsFromCommandLine()
  print(yaml.dump(options))


if __name__ == '__main__':
  run()

#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

"""Print YAML file as cmsDriver options

Usage:

./yamlToCmsDriverOptions.py TSG-Phase2HLTTDRWinter20GS-00253.yaml | xargs cmsDriver.py
"""


# Reference:
#     https://github.com/cms-sw/cmssw/blob/master/Configuration/Applications/scripts/cmsDriver.py
#     https://github.com/cms-sw/cmssw/blob/master/Configuration/Applications/python/cmsDriverOptions.py

import sys
import yaml
from Configuration.Applications.Options import parser, threeValued


def get_cmdline_arguments(d):
  d_new = {}
  for k, v in d.iteritems():
    d_new['--' + k] = v  # add '--' to key

  # Retrieve positional arguments
  args1 = d_new.pop('--args', None)
  if args1:
    args1 = [args1]
  else:
    args1 = []

  # Retrieve keyword arguments
  def formatter_fn(item):
    (k, v) = item

    # Replace by '-s', '-n', '-o'
    oneDashed = [('-s', '--step'), ('-n', '--number'), ('-o', '--number_out')]
    for (opt, opt2) in oneDashed:
      if k == opt2:
        k = opt
        break

    # If argument is a string, enclose with double quotes
    if v:
      if k in ['--customise_commands', '--inputCommands', '--outputCommands']:
        v = '"%s"' % v

    # If argument is a comma-separated list, remove spaces
    if v:
      if k in ['--step', '--customise']:
        v = v.replace(' ', '')

    # Done
    return (k, v)

  def generator_fn(items):
    for item in items:
      (k, v) = formatter_fn(item)
      yield k
      if v:
        yield v

  items = d_new.items()
  args2 = list(generator_fn(items))

  # Combine positional and keyword arguments
  args = args1 + args2
  _ = parser.parse_args(items)  # ensure that it can be parsed
  arguments = ' '.join(args)  # concatenated as one line
  return arguments


def run():
  args = sys.argv[1:]
  if len(args) != 1:
    raise RuntimeError('Expected one argument. Got: %s' % args)

  with open(args[0]) as f:
    d = yaml.safe_load(f)

  arguments = get_cmdline_arguments(d)
  print(arguments)


if __name__ == '__main__':
  run()

"""Create ragged array."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import collections

RaggedTensorNamedTuple = collections.namedtuple(
    'RaggedTensorNamedTuple', ['values', 'row_splits'])


def create_ragged_array(pylist, dtype=None, row_splits_dtype='int32'):
  """Constructs a RaggedTensorValue from a nested Python list."""

  # Ragged rank for returned value
  ragged_rank = 1

  # Build the splits for each ragged rank, and concatenate the inner values
  # into a single list.
  nested_splits = []
  values = pylist
  for dim in range(ragged_rank):
    nested_splits.append([0])
    concatenated_values = []
    for row in values:
      nested_splits[dim].append(nested_splits[dim][-1] + len(row))
      concatenated_values.extend(row)
    values = concatenated_values

  values = np.asarray(values, dtype=dtype)
  for row_splits in reversed(nested_splits):
    row_splits = np.asarray(row_splits, dtype=row_splits_dtype)
    values = RaggedTensorNamedTuple(values, row_splits)
  return values

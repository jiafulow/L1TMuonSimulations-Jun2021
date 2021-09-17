try:
  import third_party.emtf_tree
except ImportError:
  raise ImportError(
      'Could not import third_party.emtf_tree. Please run get-third-party.sh first.')

from third_party.emtf_tree import (ROOT,
                                   ROOTError,
                                   ROOT_VERSION,
                                   Tree,
                                   TreeChain,
                                   TreeQueue,
                                   get_logger,
                                   keepalive,
                                   root_open)

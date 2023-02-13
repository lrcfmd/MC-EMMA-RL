"""Extended Module Materials Assembly (EMMA)"""

from distutils.version import LooseVersion

import ase
import numpy

from mcemma_rl.mc_emma import *
import numpy as np

__all__ = ['mc_emma','all','restart_gulp','gulp','make_random_structure']

__version__ = '4.00'

if LooseVersion(np.__version__) < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers

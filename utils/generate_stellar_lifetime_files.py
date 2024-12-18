#!/usr/bin/env python

"""Create the input Meraxes stellar lifetimes hdf5 file.

Data is taken from table 14 of:

Portinari, L., Chiosi, C. & Bressan, A.
Galactic chemical enrichment with new metallicity dependent stellar yields.
Astronomy and Astrophysics 334, 505–539 (1998).
"""

import numpy as np
import h5py as h5
from cStringIO import StringIO

__author__ = "Simon Mutch"
__date__   = "2014-07-06"

raw_data = """ M    Z=0.0004  Z=0.004   Z=0.008   Z=0.02    Z=0.05
 0.6  4.28E+10  5.35E+10  6.47E+10  7.92E+10  7.18E+10
 0.7  2.37E+10  2.95E+10  3.54E+10  4.45E+10  4.00E+10
 0.8  1.41E+10  1.73E+10  2.09E+10  2.61E+10  2.33E+10
 0.9  8.97E+09  1.09E+10  1.30E+10  1.59E+10  1.42E+10
 1.0  6.03E+09  7.13E+09  8.46E+09  1.03E+10  8.88E+09
 1.1  4.23E+09  4.93E+09  5.72E+09  6.89E+09  5.95E+09
 1.2  3.08E+09  3.52E+09  4.12E+09  4.73E+09  4.39E+09
 1.3  2.34E+09  2.64E+09  2.92E+09  3.59E+09  3.37E+09
 1.4  1.92E+09  2.39E+09  2.36E+09  2.87E+09  3.10E+09
 1.5  1.66E+09  1.95E+09  2.18E+09  2.64E+09  2.51E+09
 1.6  1.39E+09  1.63E+09  1.82E+09  2.18E+09  2.06E+09
 1.7  1.18E+09  1.28E+09  1.58E+09  1.84E+09  1.76E+09
 1.8  1.11E+09  1.25E+09  1.41E+09  1.59E+09  1.51E+09
 1.9  9.66E+08  1.23E+09  1.25E+09  1.38E+09  1.34E+09
 2.0  8.33E+08  1.08E+09  1.23E+09  1.21E+09  1.24E+09
 2.5  4.64E+08  5.98E+08  6.86E+08  7.64E+08  6.58E+08
  3   3.03E+08  3.67E+08  4.12E+08  4.56E+08  3.81E+08
  4   1.61E+08  1.82E+08  1.93E+08  2.03E+08  1.64E+08
  5   1.01E+08  1.11E+08  1.15E+08  1.15E+08  8.91E+07
  6   7.15E+07  7.62E+07  7.71E+07  7.45E+07  5.67E+07
  7   5.33E+07  5.61E+07  5.59E+07  5.31E+07  3.97E+07
  9   3.42E+07  3.51E+07  3.44E+07  3.17E+07  2.33E+07
 12   2.13E+07  2.14E+07  2.10E+07  1.89E+07  1.39E+07
 15   1.54E+07  1.52E+07  1.49E+07  1.33E+07  9.95E+06
 20   1.06E+07  1.05E+07  1.01E+07  9.15E+06  6.99E+06
 30   6.90E+06  6.85E+06  6.65E+06  6.13E+06  5.15E+06
 40   5.45E+06  5.44E+06  5.30E+06  5.12E+06  4.34E+06
 60   4.20E+06  4.19E+06  4.15E+06  4.12E+06  3.62E+06
 100  3.32E+06  3.38E+06  3.44E+06  3.39E+06  3.11E+06
 120  3.11E+06  3.23E+06  3.32E+06  3.23E+06  3.11E+06
"""


if __name__ == '__main__':
    
    # Read in the raw data
    fd = StringIO(raw_data)
    cols = fd.readline().strip().split()
    data = np.loadtxt(fd, dtype=zip(cols, ['f',]*len(cols)))
    fd.close()

    # Create the ouput file and write the data
    with h5.File("../input/stellar_evo/stellar_lifetimes.hdf5", "w") as fd:
        fd.create_dataset("lifetimes", data=data)

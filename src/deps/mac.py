# Mac library dependencies
import os

deps = {}

deps['exec'] = {
    'mpicc' : 'mpicc',
    'cc' : 'cc',
    'git' : 'git',
}

deps['gsl'] = {
    'inclp' : '/usr/local/Cellar/gsl/1.16/include',
    'libp'  : '/usr/local/Cellar/gsl/1.16/lib',
    'lib'   : ['gsl', 'gslcblas'],
}

deps['hdf5'] = {
    'inclp' : '/usr/local/Cellar/hdf5/1.8.16_1/include',
    'libp'  : '/usr/local/Cellar/hdf5/1.8.16_1/lib',
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

deps['gbpCode'] = {
    'inclp' : os.path.expanduser('~/3rd_party/gbpCode/myInclude'),
    'libp'  : os.path.expanduser('~/3rd_party/gbpCode/myLib'),
    'lib'   : 'gbpLib',
}

deps['fftw'] = {
    'inclp' : os.path.expanduser('~/3rd_party/fftw/fftw-3.3.4/include'),
    'libp'  : os.path.expanduser('~/3rd_party/fftw/fftw-3.3.4/lib'),
    'lib'   : ['fftw3f_mpi', 'fftw3f'],
}
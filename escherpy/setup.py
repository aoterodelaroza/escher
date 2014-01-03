
from distutils.core import setup
import sys

sys.path.append('escherpy')
import escherpy


setup(name='escherpy',
      version='1.0',
      author='Daniel Menendez',
      author_email='danielmail7@gmail.com',
      url='http://azufre.quimica.uniovi.es',
      #download_url='https://sourceforge.net/projects/py-googlemaps/files/',
      description='Aid tool for the computational chemistry workflow.',
      long_description=escherpy.__init__.__doc__,
      packages=['escherpy'],
      provides=['escherpy'],
      keywords='quantum chemistry physics vtk',
      license='LICENSE.txt',
      requires=['VTK',             # ==5.10.1',
                'numpy',           # ==1.6.2',
                'matplotlib',      # ==1.1.1rc2',
                'pexpect',         # ==2.4',
                'memory_profiler', # ==0.30',     # developers
                'wxPython',        # ==2.8.12.1', # developers
                'RunSnakeRun',     # ==2.0.4',    # developers
                'Cython',          # ==0.19.2',   # future
               ],
      classifiers=['Development Status :: 1 - Planning',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python :: 2.7',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Scientific/Engineering :: Visualization',
                   'Topic :: Multimedia :: Graphics :: 3D Modeling',
                  ],
      libraries = [('foo', dict(sources=['escherpy/foo.f90']))],
#      ext_modules = [Extension('escherpy.foo',
#                                sources=['escherpy/foo.pyx'],
#                                libraries= = ['foo'])]
# or better
#      ext_modules = cythonize([Extension('escherpy.foo',
#                                ['escherpy/foo.pyx',
#                                'escherpy/foo.f90',])])
     )

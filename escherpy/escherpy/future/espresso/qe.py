
from __future__ import print_function, division
import os
import sys
import glob
import shlex as sh
import shutil
import exceptions
from textwrap import dedent
from subprocess import Popen, PIPE, STDOUT
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from matplotlib import rcParams
import numpy as np
from solidpy.source import scf, bands, bandspp
from solidpy.source import plotband, bandplot, dos, dospp
from solidpy.source import dosbplt, ph, dynmat, q2r, matdyn
from solidpy import logcolors
from IPython.core.debugger import Tracer; debug_here = Tracer()

log = logcolors.qelog()
#root.setLevel(log.DEBUG)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

class QEProcessException(exceptions.Exception):
    """Raised when QEspresso does not end correctly"""
    pass

class QEData(dict):
    """

    This class is empty at the begining and it is gonna be
    fulfilled with data.
    For obtaining data dict:
        Data().__dict__
    For geting attributes:
        Data().__getattribute__

    Attributes:
        espresso -- True if it is a Quantum Espresso calculation.
        sucess -- True if JOB DONE.
        niter -- number of iterations.
        energies -- energy of each step in the iteration.
        ...
    """

    def __init__(self, **kwargs):
        """ Put calculation specifications here. """

        dict.__init__(self,kwargs)
        self.__dict__.update(kwargs)

        self.iscf = scf
        self.tmpscf = scf
        """ pw.x scf input """
        self.tmpbands = bands
        """ pw.x bands input """
        self.tmpbandspp = bandspp
        self.tmpplotband = plotband
        self.tmpbandplot = bandplot
        self.tmpdos = dos
        self.tmpdospp = dospp
        self.tmpdosbplt = dosbplt
        self.tmpph = ph
        self.tmpdynmat = dynmat
        self.tmpq2r = q2r
        self.tmpmatdyn = matdyn

        self.amin = 10.0
        self.amax = 10.7
        self.astep = 0.1
        self.z = 1.
        # &control
        self.calculation = 'scf'
        self.restart_mode = 'from_scratch'
        self.title = 'si'
        self.prefix='si'
        self.pseudo_dir = '/home/daniel/src/espresso-4.3.1/pseudo/'
        self.outdir='/home/daniel/tmp/'
        self.tprnfor = 'true'
        self.tstress = 'true'
        # &system
        self.ibrav = 2
        self.celldm1 = 10.2625
        self.nat=  2
        self.ntyp= 1
        self.ecutwfc = 60.
        self.smearing = 'gaussian'
        self.nbnd = 8
        self.occupations = 'tetrahedra'
        self.ecutrho = 720.0
        # &electrons
        self.mixing_beta = 0.7
        self.conv_thr =  '1.0d-8'
        self.diagonalization = 'cg'
        self.mixing_mode = 'plain'
        # atomic species
        self.atomic_species = dedent("""
         Ge  72.63  Ge.pbe-paw_kj.UPF
         """)
        # ATOMIC_POSITIONS
        self.atomic_positions = dedent("""\
         crystal
         Si 0.125 0.125 0.125
         Si -0.125 -0.125 -0.125
         """)
        # K_POINTS
        self.k_points = dedent("""\
         automatic
         3 3 3 1 1 1
        """)
        self.k_pointspath = dedent("""\
         tpiba_b
         5
         0.5 0.5 0.5 40
         0.0 0.0 0.0 40
         0.5 0.5 0. 40
         0.75 0.5 0.25 40
         0.0 0.0 0.0 1
        """)
        """
         0.0 0.0 0.0 40
         0.5 0.5 0. 40
         0.75 0.5 0.25 40
         0.5 0.5 0.5 40
         0.0 0.0 0.0 1
         0.5 0.5 0.5 40
         0.0 0.0 0.0 40
         0.0 0.0 1.0 40
         0.0 1.0 1.0 40
         0.0 0.0 0.0 1
         tpiba_b
         5
         0.5 0.5 0.5 100
         0.0 0.0 0.0 100
         0.5 0.5 0. 100
         0.75 0.5 0.25 100
         0.0 0.0 0.0 1
         para diamante
         tpiba_b
         5
         0.5 0.5 0.5 40
         0.0 0.0 0.0 40
         0.0 0.0 1.0 40
         0.0 1.0 1.0 40
         0.0 0.0 0.0 1
         mia
         tpiba_b
         5
         0.5 0.5 0.5 40
         0.0 0.0 0.0 40
         0.5 0.5 0. 40
         0.75 0.5 0.25 40
         0.0 0.0 0.0 1
         """
        self.cell_parameters = dedent("""\
         cubic
         0.125 0.125 0.125
         0.125 0.625 0.625
         0.625 0.125 0.625
         0.625 0.625 0.125
         0.875 0.875 0.875
         0.875 0.375 0.375
         0.375 0.875 0.375
         0.375 0.375 0.875
        """)
        # &inputpp
        self.Emin = -9.0
        self.Emax = 26.0
        self.DeltaE = 0.1
        self.lsym = 'true'
        # &inputph
        self.tr2_ph = '1d-14'
        self.ldisp = 'true'
        self.nq1 = 4
        self.nq2 = 4
        self.nq3 = 4
        self.amass1 = 28.0855
        # &input
        self.asr = 'simple'
        self.q1 = 0.0
        self.q2 = 0.0
        self.q3 = 0.0
        self.zasr = 'simple'
        self.dos = 'true'
        self.nk1 = 50
        self.nk2 = 50
        self.nk3 = 50
        # plot
        self.kticks = [0.0, 0.86, 1.57, 1.92, 2.9]
        self.klabels = [r'L', r'\textit{\Gamma}', r'X', r'W', r'\textit{\Gamma}']

        shutil.rmtree(self.outdir)
        os.mkdir(self.outdir)
        self.subs()
        self.check()

    def abrir(self, logfile):
        """
        Opens an output file and returns lines without whitespaces.
        """
        with open(logfile, 'r') as file:
            text = file.read()

        lines = [line for line in text.splitlines() if line]
        lines = map(lambda s: s.strip(), lines)
        return lines

    def scftake(self):
        """
        Extracts data from {prefix}.scf.out
        """
        log.info('Extracting from {}.scf.out...'.format(self.prefix))
        lines = self.abrir('{}.scf.out'.format(self.prefix))
        ks = {}
        for n,i in enumerate(lines):
            if 'Quantum ESPRESSO suite' in i:
                self.espresso = True
                log.info('Quantum ESPRESSO calculation.')
            if i.startswith('JOB DONE'):
                self.sucesscf = True
                log.info('JOB DONE calculation.')
            if i.startswith('convergence has been achieved'):
                self.niter = int(sh.split(i)[-2])
                log.info('Converged.')

            if i.startswith('!'):
                self.energies = []
                self.energies.append(float(sh.split(i)[-2]))
                self.energies_units = sh.split(i)[-1]
                log.debug('Energy.')
            if 'bravais-lattice' in i:
                self.ibrav = int(sh.split(i)[-1])
            if i.startswith('PWSCF'):
                times = [sh.split(i)[x] for x in [2, 4]]
                self.times = map(lambda s : float(s.strip('s')), times)
            if i.startswith('unit-cell volume'):
                self.volume = float(sh.split(i)[-2])
                self.volume_units = sh.split(i)[-1]
            if i.startswith('number of atoms/cell'):
                self.natom = int(i.rsplit(' ', 1)[1])
            if i.startswith('number of atomic types'):
                self.ntype = int(i.rsplit(' ', 1)[1])
            if i.startswith('number of electrons'):
                self.nelec = float(i.rsplit(' ', 1)[1])

            # Cell parameters
            if i.startswith('lattice parameter'):
                self.alat = float(sh.split(i)[-2])
                self.alat_units = sh.split(i)[-1]
            if i.startswith('celldm(1)'):
                alats = [float(sh.split(i)[x]) for x in [1,3,5]]
            if i.startswith('celldm(4)'):
                alatsd = [float(sh.split(i)[x]) for x in [1,3,5]]

            # Metrics matrix
            if i.startswith('a(1)'):
                v1 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('a(2)'):
                v2 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('a(3'):
                v3 = [float(sh.split(i)[x]) for x in [3,4,-2]]

            # Metrics inverse matrix
            if i.startswith('b(1)'):
                w1 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('b(2)'):
                w2 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('b(3'):
                w3 = [float(sh.split(i)[x]) for x in [3,4,-2]]

            # K points
            if i.startswith('number of k points'):
                self.nkpoints = int(sh.split(i)[-1])
                n = n + 1
                for k in range(self.nkpoints):
                    n = n + 1
                    ks['{}'.format(k)] = \
                        [sh.split(lines[n])[y] for y in [-6, -5, -4, -1]]
                    ks['{}'.format(k)] = \
                            map(lambda s : float(s.strip('),')),
                                ks['{}'.format(k)])
                    self.ks = ks
                log.debug('K points.')

        self.celldm = alats + alatsd
        log.debug('Cell parameters.')
        self.R = np.array([v1,v2,v3])
        log.debug('Metrics matrix.')
        self.rR = np.array([v1,v2,v3])
        log.debug('Metrics inverse matrix.')

    def dostake(self):
        """
        Extracts data form {prefix}.dos.out
        """
        dosout = '{}.dos.out'.format(self.prefix)
        log.info('Extracting from: {}...'.format(dosout))
        if not os.path.isfile(dosout):
            return 'No efermi taken still. Run dosbands(), then plotbandx()'

        lines = self.abrir(dosout)
        for i in lines:
            if i.startswith('the Fermi energy is '):
                self.efermi = float(sh.split(i)[-2])
                self.efermi_units = sh.split(i)[-1]
                log.debug('Efermi.')

    def qe_funciona(self, input, task, name):
        """
        Pipe for any Quantum Espresso task.
        """

        proc = Popen(
            'pw.x < si.scf.in',
            shell=True,
            stdin=PIPE,
            stdout=PIPE  )


        while True:
            next_line = proc.stdout.readline()
            if next_line == '' and proc.poll() != None:
                break
            sys.stdout.write(next_line)
            sys.stdout.flush()

    def qe(self, input, task, name):
        """
        Pipe for any Quantum Espresso task.
        """

        # TODO iname oname
        log.debug('Saving input: {}.in'.format(name))
        print(input, file=open('{}.in'.format(name), 'w'))
        proc = Popen(
            '{} < {}.in'.format(task, name),
            shell=True,
            stdin=PIPE,
            stdout=PIPE)

        if os.path.isfile('{}.out'.format(name)):
            log.debug('Removing previus: {}.out'.format(name))
            os.remove('{}.out'.format(name))

        log.debug('Running {} < {}.in...'.format(task, name))
        while True:
            next_line = proc.stdout.readline()
            print(next_line, file=open('{}.out'.format(name), 'a'))
            # TODO scftake, ... handle calc types
            if next_line == '' and proc.poll() != None:
                break
            #sys.stdout.write(next_line)
            #sys.stdout.flush()

        proc.communicate()[0]
        exitCode = proc.returncode
        if (exitCode == 0):
            log.info('Sucessful calculation.')
            return
        else:
            raise QEProcessException(task, exitCode, '{}.out'.format(name))

    def __qe(self, input, task, name):
        """
        Pipe for any Quantum Espresso task.
        """

        #print('Running {} < {}.in...'.format(task, name))
        log.debug('Running {} < {}.in...'.format(task, name))
        print(input, file=open('{}.in'.format(name), 'w'))
        #proc = Popen(task + ' ' + input,
        #proc = Popen(task,
        proc = Popen(task + ' {}.in'.format(name),
                shell = True,
                stdin = PIPE,
                stdout=PIPE)
                #stderr=STDOUT)

        #proc.stdin.write(input)
        #proc.stdin.write('')
        while True:
            log.debug('entering while')
            nextline = proc.stdout.readline()
            log.debug('nextline assigned')
            if nextline == '' and proc.poll() != None:
                break
            sys.stdout.write(nextline)
            sys.stdout.flush()

        output=proc.communicate()[0]
        exitCode = proc.returncode
        #print(output, file=open('{}.out'.format(name), 'w'))
        #print('....................................done.')
        if (exitCode == 0):
            return output
        else:
            raise ProcessException(command, exitCode, output)
    def _qe(self, input, task, name):
        """
        Pipe for any Quantum Espresso task.
        """

        print('Running {} < {}.in...'.format(task, name))
        print(input, file=open('{}.in'.format(name), 'w'))
        pscf = Popen(task,
                stdin = PIPE,
                stdout=PIPE,
                stderr=STDOUT)

        output=pscf.communicate(input)[0]
        print(output, file=open('{}.out'.format(name), 'w'))
        #print('....................................done.')

    def subs(self):
        """
        Proper substitutions in templates using Data.
        """
        data = self.__dict__
        self.iscf = self.tmpscf.format(q=data)
        self.ibands = self.tmpbands.format(q=data)
        self.bandspp = self.tmpbandspp.format(q=data)

        self.idos = self.tmpdos.format(q=data)
        self.dospp = self.tmpdospp.format(q=data)

        self.ph = self.tmpph.format(q=data)
        self.dynmat = self.tmpdynmat.format(q=data)
        self.q2r = self.tmpq2r.format(q=data)
        self.matdyn = self.tmpmatdyn.format(q=data)

    def check(self):
        """
        Initial check.
        """
        log.debug(os.getcwd() + ': starting')
        if not os.path.isdir(self.pseudo_dir):
            log.debug('El directorio de pseudo no existe.')
        if not os.path.isdir(self.outdir):
            log.debug('El directorio de outdir no existe.')

        log.debug('Progress:')

    def scf(self):
        """
        Runs pw.x SCF
        """
        self.qe(self.iscf, 'pw.x',
                name='{}.scf'.format(self.prefix))


    def bands(self):
        """
        Runs pw.x BANDS
        """
        self.scf()
        log.debug('Cleaning any previous bands coordinates xmgr...')
        map(lambda s : os.remove(s), glob.glob('*bands.xmgr*'))
        self.qe(self.ibands, 'pw.x',
                name='{}.bands'.format(self.prefix))
        self.qe(self.bandspp, 'bands.x',
                name='{}.bands.pp'.format(self.prefix))

    def dosbands(self):
        """
        Runs dos.x DOS_BANDS
        """
        self.scf()
        self.qe(self.idos, 'pw.x',
                name='{}.dos'.format(self.prefix))
        self.qe(self.dospp, 'dos.x',
                name='{}.dos.pp'.format(self.prefix))

    def plotbandx(self):
        """
        Runs plotband.x and band_plot.x
        """

        self.dostake()
        data = self.__dict__
        self.plotband = self.tmpplotband.format(q=data)
        self.bandplot = self.tmpbandplot.format(q=data)
        #self.dosbplt = self.dosbplt.format(q=data)

        self.qe(self.plotband, 'plotband.x',
                name='{}.plotband'.format(self.prefix))
        self.qe(self.bandplot, 'band_plot.x',
                name='{}.bandplot'.format(self.prefix))
        #qe(dosbplt, 'gnuplot', name='si.dos.plt')
        #pbs.gnuplot('si_bands.gnu')

    def plt2bands(self):
        """
        Plots DOS and bands with matplotlib.
        """

        self.dostake()
        bands = glob.glob('*bands.xmgr*')
        fig = plt.figure()
        ax = fig.add_subplot(1,2,1)

        for i in bands:
            bandsdata = mlab.csv2rec(i, delimiter=' ')
            ax.plot(zip(*bandsdata)[0],
                     map(lambda s: (s -self.efermi), zip(*bandsdata)[1]),
                     '.', color='red')

        ax.set_ylabel(r'\textit{Energy (eV)}', fontsize=28)
        plt.xticks(self.kticks, self.klabels, fontsize=28)
        ax.set_xlim([0.0,2.9])
        ax.set_ylim([-15, 7.])
        plt.yticks(fontsize=18)
        ax.grid(True)

        ax2 = fig.add_subplot(1, 2, 2)
        b = np.genfromtxt('{}.dos'.format(self.prefix), unpack = True)
        b0 = map(lambda s: (s -self.efermi), b[0])
        ax2.plot(b[1], b0, color = 'red')
        ax2.fill_between(b[1], 200, b0, color='red')
        ax2.set_ylim([-15, 7.])
        ax2.set_xlim([0.0,2.5])
        ax2.set_xlabel(r'\textit{DOS}', fontsize = 28)
        ax2.get_yaxis().set_visible(False)
        plt.xticks(fontsize=18)
        plt.savefig('plot.pdf')
        plt.show()

    def various(self):
        """ For various calculations """

        print('#Calculation on Si fcc',
              '# Volume bohr^3',
              '# Energy ry',
              '# z 1', sep = '\n', file=open('evdata.txt', 'w'))

        for i in np.arange(self.amin, self.amax, self.astep):
            log.info('a = {}'.format(i))
            self.celldm1 = i
            old  = self.prefix
            self.prefix = '{}.{}'.format(self.prefix, i)
            self.subs()
            self.qe(self.iscf, 'pw.x',
                    name='{}.scf'.format(self.prefix))
            self.prefix  = old

        es = []
        for i in np.arange(self.amin, self.amax, self.astep):
            log.info('a = {}'.format(i))
            self.celldm1 = i
            old  = self.prefix
            self.prefix = '{}.{}'.format(self.prefix, i)
            self.scftake()
            es.append(self.energies[-1])
            print(self.volume/self.z, self.energies[-1],
                    sep=' ', end= '\n', file=open('evdata.txt', 'a'))
            self.prefix  = old

        with open('evdata.txt', 'r') as text:
            text = text.read()
            proc1 = Popen(['asturfit', 'evdata.txt'], stdout=PIPE, stderr=STDOUT)
            result = proc1.communicate()[0]
            print(result)


        plt.grid(True)
        plt.plot(list(np.arange(self.amin, self.amax, self.astep)), es)
        plt.savefig('ev.png')

if __name__ == '__main__':
    a = QEData()

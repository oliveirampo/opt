"""Configuration object for all main options of the pipeline.

Classes:
    Conf
    PlotConf
"""

import numpy as np
import sys

import ChargeDistribution
from matrix import Matrix
import combiningRule
import IO


class Conf:
    """Configuration object for all main options of the pipeline.

    Attributes:
    -----------
        _it: (int) Iteration number.
        _nJobs: (it) Number of simulation jobs.
        _scl_dns: (float) Scaling factor for density.
        _scl_hvp: (float) Scaling factor for vaporization enthalpy.
        _opt_nit: (it) Number of steps for the search of optimal parameters.
        _rng_scl: (float) Range of parameter variation.
        _eq_stp: (int) Number of equilibration steps during the simulation.
        _eq_frq: (int) Frequency that output is printed during equilibration step.
        _prd_stp: (int) Number of production steps during the simulation.
        _prd_frq: (int) Frequency that output is printed during production steps.
        _wall_time: (str, int) Simulation wall time.
        _cr: (combiningRule) Combining rule.
        _matrix: (Matrix) Matrix with usage of C12(II) parameter.
        _charge_distribution_method: (ChargeDistributionMethod) Method for charge distribution.
        _charge_group_type: (Charge_group_type) Type of charge group.
        _kappa: (float) Parameter for charge distribution calculation.
        _lam: (float) Parameter for charge distribution calculation.
        _scl_sig_NEI: (float) Scaling factor for 1-4 sigma.
        _scl_eps_NEI: (float) Scaling factor for 1-4 epsilon.
        _ignoreIAC: (list) List of atom types indexes which are ignored during charge distribution.
        _plotConf: (PlotConf) Configuration object for plots.

        _inpDir: (str) Directory of input files.
        _molDataFile: (str) Input file with list of molecules ("mol.dat").
        _atomListFile: (str) Input file with list of atoms ("listAtom.dat").
        _bondListFile: (str) Input file with list of bonds ("listBond.dat").
        _angListFile: (str) Input file with list of angles ("listAng.dat").
        _matrixFile: (str) Input file with information about the usage of C12(II) parameters ("matrix.dat").
        _symSigFile: (str) Input file with list of symmetrical atom types in terms of sigma ("symmetry_sig.dat").
        _symEpsFile: (str) Input file with list of symmetrical atom types in terms of epsilon ("symmetry_eps.dat").
        _prmNeiFile: (str) Input file with 1-4 interaction parameters ("prm_NEI.dat").
        _samTemplateFile: (str) Input file with template to run SAMOS ("model.sam").
        _prmFile: (str) Input file with parameters at given iteration ("prm_[0-9]*.dat").
        _ifpFile: (str) IFP input file with parameters from 2016H66 force field ("2016H66_upd.ifp").
        _vdWFile: (str) File with van der Waals radii.
        _prmCrFile: (str) File with alpha parameter for linear combination of combining rules.

        _optDir: (str) Output directory.
        _outPrmFile: (str) Output file with optimized parameters at given iteration ("prm_[0-9]*.dat").
        _optOutFile: (str) Output file with results from optimization.

        self._inpFiles: (list) List of most used input files.

    Methods:
    --------
        ignoreAtom(iac): Checks if given atom type should be ignored in the charge distribution scheme.
        __str__(): Returns a string of this object.
    """

    def __init__(self, fileName, it):
        """Constructs all the necessary attributes for the configuration object from input file
        and checks if input files exist.

        :param fileName: (str) Name of configuration input file.
        :param it: (int) Iteration number.
        """

        self._it = int(it)
        self._nJobs = 4
        self._scl_dns = 20.0
        self._scl_hvp = 1.0
        self._opt_nit = 5000
        self._rng_scl = 5.0
        self._eq_stp = 150000
        self._eq_frq = 2500
        self._prd_stp = 150000
        self._prd_frq = 2500
        self._wall_time = 12
        self._cr = combiningRule.GeometricCR()
        self._matrix = Matrix()
        self._charge_distribution_method = ChargeDistribution.NONE()
        self._charge_group_type = ChargeDistribution.Atomic()
        self._kappa = 1.0
        self._lam = 1.0
        self._scl_sig_NEI = 1.0
        self._scl_eps_NEI = 1.0
        self._ignoreIAC = []
        self._plotConf = PlotConf

        # input files
        self._inpDir = '00_inp'
        self._molDataFile = '{}/{}'.format(self._inpDir, "mol.dat")
        self._atomListFile = '{}/{}'.format(self._inpDir, "listAtom.dat")
        self._bondListFile = '{}/{}'.format(self._inpDir, "listBond.dat")
        self._angListFile = '{}/{}'.format(self._inpDir, "listAng.dat")
        self._matrixFile = '{}/{}'.format(self._inpDir, "matrix.dat")
        self._symSigFile = '{}/{}'.format(self._inpDir, "symmetry_sig.dat")
        self._symEpsFile = '{}/{}'.format(self._inpDir, "symmetry_eps.dat")
        self._prmNeiFile = '{}/{}'.format(self._inpDir, "prm_NEI.dat")
        self._samTemplateFile = '{}/{}'.format(self._inpDir, "model.sam")
        self._prmFile = '{}/{}_{}.dat'.format(self._inpDir, "prm", self._it)
        self._ifpFile = '{}'.format("prm/2016H66_upd.ifp")
        self._vdWFile = '{}/{}'.format(self._inpDir, "vdw.dat")
        self._prmCrFile = '{}/{}_cr_{}.dat'.format(self._inpDir, "prm", self._it)

        self._optDir = 'opt_{}'.format(self._it + 1)
        self._outPrmFile = '{}/{}_{}.dat'.format(self._inpDir, "prm", self._it + 1)
        self._optOutFile = '{}/{}_{}.out'.format(self._optDir, "opt", self._it + 1)

        self._inpFiles = [self._molDataFile, self._atomListFile, self._bondListFile, self._angListFile,
                          self._matrixFile, self._ifpFile, self._prmNeiFile, self._prmFile]

        IO.readConf(self, fileName)
        # TODO
        # IO.checkInpDir(self)
        IO.checkInpFiles(self)

    @property
    def it(self):
        return self._it

    @property
    def nJobs(self):
        return self._nJobs

    @property
    def scl_dns(self):
        return self._scl_dns

    @property
    def scl_hvp(self):
        return self._scl_hvp

    @property
    def opt_nit(self):
        return self._opt_nit

    @property
    def rng_scl(self):
        return self._rng_scl

    @property
    def eq_stp(self):
        return self._eq_stp

    @property
    def eq_frq(self):
        return self._eq_frq

    @property
    def prd_stp(self):
        return self._prd_stp

    @property
    def prd_frq(self):
        return self._prd_frq

    @property
    def wall_time(self):
        return self._wall_time

    @property
    def cr(self):
        return self._cr

    @property
    def matrix(self):
        return self._matrix

    @property
    def charge_distribution_method(self):
        return self._charge_distribution_method

    @property
    def charge_group_type(self):
        return self._charge_group_type

    @property
    def kappa(self):
        return self._kappa

    @property
    def lam(self):
        return self._lam

    @property
    def scl_sig_NEI(self):
        return self._scl_sig_NEI

    @property
    def scl_eps_NEI(self):
        return self._scl_eps_NEI

    @property
    def ignoreIAC(self):
        return self._ignoreIAC

    @property
    def plotConf(self):
        return self._plotConf

    @property
    def inpFiles(self):
        return self._inpFiles

    @property
    def inpDir(self):
        return self._inpDir

    @property
    def molDataFile(self):
        return self._molDataFile

    @property
    def atomListFile(self):
        return self._atomListFile

    @property
    def bondListFile(self):
        return self._bondListFile

    @property
    def angListFile(self):
        return self._angListFile

    @property
    def matrixFile(self):
        return self._matrixFile

    @property
    def symSigFile(self):
        return self._symSigFile

    @property
    def symEpsFile(self):
        return self._symEpsFile

    @property
    def prmNeiFile(self):
        return self._prmNeiFile

    @property
    def prmFile(self):
        return self._prmFile

    @property
    def ifpFile(self):
        return self._ifpFile

    @property
    def vdWFile(self):
        return self._vdWFile

    @property
    def prmCrFile(self):
        return self._prmCrFile

    @property
    def samTemplateFile(self):
        return self._samTemplateFile

    @property
    def optDir(self):
        return self._optDir

    @property
    def outPrmFile(self):
        return self._outPrmFile

    @property
    def optOutFile(self):
        return self._optOutFile

    @it.setter
    def it(self, n):
        self._it = n

    @nJobs.setter
    def nJobs(self, n):
        try:
            n = int(n)
            self._nJobs = n

        except ValueError as err:
            print('\nWrong number of jobs')
            print(err)
            sys.exit('\n')

    @scl_dns.setter
    def scl_dns(self, n):
        self._scl_dns = float(n)

    @scl_hvp.setter
    def scl_hvp(self, n):
        self._scl_hvp = float(n)

    @opt_nit.setter
    def opt_nit(self, n):
        self._opt_nit = int(n)

    @rng_scl.setter
    def rng_scl(self, n):
        n = float(n)
        self._rng_scl = float(n)

    @eq_stp.setter
    def eq_stp(self, n):
        self._eq_stp = n

    @eq_frq.setter
    def eq_frq(self, n):
        self._eq_frq = n

    @prd_stp.setter
    def prd_stp(self, n):
        self._prd_stp = n

    @prd_frq.setter
    def prd_frq(self, n):
        self._prd_frq = n

    @wall_time.setter
    def wall_time(self, n):
        self._wall_time = n

    @cr.setter
    def cr(self, n):
        cR = combiningRule.CR.get_object(n)
        self._cr = cR

    @matrix.setter
    def matrix(self, n):
        self._matrix = n

    @charge_distribution_method.setter
    def charge_distribution_method(self, n):
        method = ChargeDistribution.ChargeDistributionMethod.get_object(n)

        if isinstance(method, ChargeDistribution.BondChargeDistributionMethod):
            if not isinstance(self._charge_group_type, ChargeDistribution.Atomic):
                print(
                    'This charge distribution method ({}) is not compatible with ({}) charge group.\n'
                    'Try ({}) charge group, for instance.'.format(type(method), type(self._charge_group_type),
                                                                  type(ChargeDistribution.Atomic())))

                sys.exit(123)

        self._charge_distribution_method = method

    @charge_group_type.setter
    def charge_group_type(self, n):
        method = ChargeDistribution.Charge_group_type.get_object(n)
        self._charge_group_type = method

    @kappa.setter
    def kappa(self, n):
        self._kappa = float(n)

    @lam.setter
    def lam(self, n):
        self._lam = float(n)

    @scl_sig_NEI.setter
    def scl_sig_NEI(self, n):
        n = float(n)
        self._scl_sig_NEI = n

    @scl_eps_NEI.setter
    def scl_eps_NEI(self, n):
        n = float(n)
        self._scl_eps_NEI = n

    @ignoreIAC.setter
    def ignoreIAC(self, n):
        if isinstance(n, str):
            self._ignoreIAC = [int(n)]
        elif isinstance(n, list):
            self._ignoreIAC = [int(i) for i in n]
        else:
            if np.isnan(n):
                self._ignoreIAC = []
            else:
                sys.exit('ERROR: wrong type for ignoreIAC')

    @plotConf.setter
    def plotConf(self, plt):
        self._plotConf = plt

    @inpFiles.setter
    def inpFiles(self, n):
        self._inpFiles = n

    def ignoreAtom(self, iac):
        if iac in self._ignoreIAC:
            return True
        return False

    def __str__(self):
        s = '\t{}:00 {} {}'.format(self._wall_time, self._prd_frq, self._prd_stp)
        return s


class PlotConf:
    """Configuration object for plots.

    Attributes:
    -----------
        _map_iac_name: (Dict) Dictionary which maps iac code to atom name.
        _map_cod_family: (Dict) Dictionary which maps molecule cod to family code.
        _map_cod_color: (Dict) Dictionary which maps molecule code to colors.
        _map_cod_marker: (Dict) Dictionary which maps molecule cod to marker.
        _settings: (Dict) Dictionary with plot settings.
        _plotDir: (str) Directory of plots.
    """

    def __init__(self):
        """Constructs all the necessary attributes for the plot configuration object."""

        self._map_iac_name = {}
        self._map_cod_family = {}
        self._map_cod_color = {}
        self._map_cod_marker = {}
        self._settings = {}
        self._plotDir = ''

    @property
    def map_iac_name(self):
        return self._map_iac_name

    @property
    def map_cod_family(self):
        return self._map_cod_family

    @property
    def map_cod_color(self):
        return self._map_cod_color

    @property
    def map_cod_marker(self):
        return self._map_cod_marker

    @property
    def settings(self):
        return self._settings

    @property
    def plotDir(self):
        return self._plotDir

    @map_iac_name.setter
    def map_iac_name(self, d):
        self._map_iac_name = d

    @map_cod_family.setter
    def map_cod_family(self, d):
        self._map_cod_family = d

    @map_cod_color.setter
    def map_cod_color(self, d):
        self._map_cod_color = d

    @map_cod_marker.setter
    def map_cod_marker(self, d):
        self._map_cod_marker = d

    @settings.setter
    def settings(self, d):
        self._settings = d

    @plotDir.setter
    def plotDir(self, plotDir):
        self._plotDir = plotDir

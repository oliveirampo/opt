"""Implementation of the possible options that the pipeline can execute.

Classes:
--------
	Action
	Gen(Action)
	Ana(Action)
	Plot(Action)
	Opt(Action)
	Sub(Action)
	GROMOSCoTo(Action)
"""

from abc import ABC, abstractmethod
import shutil
import os

from scr.base import IO
from scr.base import iac
from scr.base import ana
from scr.base import plot
from scr.base import check
from scr.base import optimize
from scr.base import gromos_mtb
from scr.base import myExceptions
from scr.base import createNewIFP
from scr.base import prepareFiles
from scr.base import writeOutFiles
from scr.base import molecules_utils


class Action(ABC):
	"""Base class for the implementation of the possible options that the pipeline can execute.
	The options are:
	GEN, ANA, OPT, SUB, PLOT.
	See their respective implementations for more details.

	Attributes
	----------
	:parameter it: (int) Iteration number.

	Methods
	-------
	get_choices(): Returns list of possible choices of actions.
	get_object(runType, it): Returns object of class that implements base class Action.
	read_inp_files(conf): Reads input files.
	run(conf, molecules, atomTypes): Executes option.
	"""

	def __init__(self, it):
		"""Constructs all the necessary attributes for the given action.

		:param it: (int) Iteration number.
		"""
		self.it = it

	@staticmethod
	def get_choices():
		"""Returns list of possible choices of actions."""

		return ['PREP', 'GEN', 'ANA', 'OPT', 'SUB', 'PLOT', 'MTB', 'IFP', 'CHECK']

	@staticmethod
	def get_object(runType, it):
		"""Returns object of class that implements base class Action.

		:param runType: (str) Code for class that implements base class Action.
		:parameter it: (str) Iteration number.
		:return: One of the classes that implements Action.
		"""

		classes = {"PREP": Prep, "GEN": Gen, "ANA": Ana, "PLOT": Plot, "SUB": Sub, "OPT": Opt,
				   "MTB": MTB, "IFP": IFP, "CHECK": Check}

		if runType not in classes:
			raise myExceptions.ClassNotImplemented(runType, 'Action')

		act = classes[runType](it)

		return act

	@staticmethod
	def read_inp_files(conf):
		"""Reads input files specified in the configuration file

		:param conf: (configuration.Conf) Configuration object read from configuration input file.
		:return:
			molecules: (collections.OrderedDict) Ordered dictionary of molecules.
			atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.

		1 - Reads molecule file.
		2 - Reads prm file.
		3 - Reads prm_nei file.
		4 - Reads matrix file.
		"""

		# read molecule file
		molecules = IO.readMolData(conf)

		# read prm file
		atomTypes = IO.readPrm(conf.prmFile)

		# read prm_nei file
		IO.readPrmNei(atomTypes, conf.prmNeiFile)

		# read matrix file
		conf.matrix = IO.readMatrix(conf.matrixFile)

		return molecules, atomTypes

	@staticmethod
	def read_extra_inp_files(conf, molecules, atomTypes):
		"""Reads additional input files specified in the configuration file

		:param conf: (configuration.Conf) Configuration object read from configuration input file.
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		:return:
			molecules: (collections.OrderedDict) Ordered dictionary of molecules.
			atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
			crPrms: (dict) Parameters to do linear combination of combining rules.
		"""

		IO.readListAtom(conf, molecules, conf.atomListFile)
		bnd, ang = IO.readRefBondAndAngle(conf.ifpFile)

		vdw = IO.readVdWRadii(conf.vdWFile)
		iac.add_vdw_radii(atomTypes, vdw)

		crPrms = IO.readPrmCrFile(conf.prmCrFile)

		IO.readListBond(bnd, molecules, conf.bondListFile)
		IO.readListAngle(ang, molecules, conf.angListFile)
		IO.readSymmetry(atomTypes, "sig", conf.symSigFile)
		IO.readSymmetry(atomTypes, "eps", conf.symEpsFile)
		molecules_utils.createEffectivePrms(atomTypes, molecules, conf.charge_group_type)

		return molecules, atomTypes, crPrms

	@abstractmethod
	def run(self, conf, molecules, atomTypes):
		"""Executes the implemented action.

		:param conf: (configuration.Conf) Configuration object.
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		"""
		pass


class Prep(Action):
	"""Setups up default input files."""

	def run(self, conf, molecules, atomTypes):
		"""Setups up default input files.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		"""

		prepareFiles.run(molecules)


class Gen(Action):
	"""Generates input files to run simulations."""

	def run(self, conf, molecules, atomTypes):
		"""Generates input files to run simulations.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.

		The following output files are created for each molecule:
			mol file: mol.mol
			param file: para.mod
			sam file: mol_[0-3].sam
			coordinate files (copied from cfg/ directory to the same directory)
			submission files: mol_[0-3]_sub.sh

		The number of files for each type is determined by the variable nJobs in the conf.
		"""

		molecules, atomTypes, crPrms = Action.read_extra_inp_files(conf, molecules, atomTypes)

		kap = conf.kappa
		lam = conf.lam

		molecules_utils.computeChargeDistribution(conf.charge_distribution_method, molecules, atomTypes, kap, lam)
		molecules_utils.computeCR(conf.cr, crPrms, conf.scl_sig_NEI, conf.scl_eps_NEI, molecules, atomTypes, conf.matrix)

		writeOutFiles.writeParamMod(self.it, molecules)

		samTemplate = IO.readSamTemplateFile(conf.samTemplateFile)
		writeOutFiles.writeSamFile(conf, molecules, samTemplate)
		writeOutFiles.copyCOTO(molecules)


class Ana(Action):
	"""Performs analysis of the simulation results."""

	def run(self, conf, molecules, atomTypes):
		"""Perform analysis of the simulation results.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.

		See ana.runAna method for more information.
		"""

		molecules, atomTypes, crPrms = Action.read_extra_inp_files(conf, molecules, atomTypes)

		anaDir = 'ana_{}'.format(conf.it)
		if os.path.exists(anaDir):
			shutil.rmtree(anaDir)
		os.makedirs(anaDir)

		ana.runAna(conf, molecules, anaDir)


class Plot(Action):
	"""Plots results from optimization and from simulations."""
	def run(self, conf, molecules, atomTypes):
		plotDir = 'plot_{}'.format(conf.it)
		conf.plotConf.plotDir = plotDir

		if not os.path.exists(plotDir):
			os.makedirs(plotDir)

		plot.run(conf)


class Sub(Action):
	"""Submits jobs to run the simulations."""
	def run(self, conf, molecules, atomTypes):
		molecules, atomTypes, crPrms = Action.read_extra_inp_files(conf, molecules, atomTypes)
		writeOutFiles.writeSubScript(conf, molecules)


class Opt(Action):
	"""Optimize force-field parameters."""

	def __init__(self, it):
		"""Constructs all the necessary attributes for the given action.

		:param it: (int) Iteration number.
		"""

		Action.__init__(self, it - 1)

	def run(self, conf, molecules, atomTypes):
		"""Optimize force-field parameters.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.

		See optimize.runOptimization method for more information.
		"""

		molecules, atomTypes, crPrms = Action.read_extra_inp_files(conf, molecules, atomTypes)
		optimize.runOptimization(conf, crPrms, molecules, atomTypes)


class MTB(Action):
	"""Updates GROMOS mtb files with the final force-field parameters."""

	def __init__(self, it):
		"""Constructs all the necessary attributes for the given action.

		param it: (int) Iteration number.
		"""

		Action.__init__(self, it)

	def run(self, conf, molecules, atomTypes):
		"""Reads parameters, calculates effective parameters to each molecule and updates mtb file.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		:return:
		"""

		molecules, atomTypes, crPrms = Action.read_extra_inp_files(conf, molecules, atomTypes)

		mtbDir = 'mtb'
		mtbGROMOSDir = 'mtb_GROMOS'

		kap = conf.kappa
		lam = conf.lam

		molecules_utils.computeChargeDistribution(conf.charge_distribution_method, molecules, atomTypes, kap, lam)
		molecules_utils.computeCR(conf.cr, crPrms, conf.scl_sig_NEI, conf.scl_eps_NEI, molecules, atomTypes, conf.matrix)

		# create new IFP file with N new dummy atom types.
		newIacDict = createNewIFP.create(conf, atomTypes, crPrms)

		gromos_mtb.updateMTB(newIacDict, atomTypes, molecules, mtbDir, mtbGROMOSDir)


class IFP(Action):
	"""Extends GROMOS ifp file with dummy atoms types."""

	def __init__(self, it):
		"""Constructs all the necessary attributes for the given action.

		param it: (int) Iteration number.
		"""

		Action.__init__(self, it)

	def run(self, conf, molecules, atomTypes):
		"""Reads parameters, calculates effective parameters to each molecule and updates mtb file.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		:return:
		"""

		# molecules, atomTypes, crPrms = Action.read_extra_inp_files(conf, molecules, atomTypes)
		crPrms = IO.readPrmCrFile(conf.prmCrFile)

		# create new IFP file with N new dummy atom types.
		_ = createNewIFP.create(conf, atomTypes, crPrms)


class Check(Action):
	"""Perform a few checks on the simulation results"""

	def run(self, conf, molecules, atomTypes):
		"""Check parameter ratio.

		:param conf: (configuration.Conf) Configuration object
		:param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
		:param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
		"""

		molecules, atomTypes, _ = Action.read_extra_inp_files(conf, molecules, atomTypes)

		check.prm_ratio(molecules, atomTypes)

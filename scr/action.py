from scr import IO, writeOutFiles, molecules_utils


class Action:
	def __init__(self, it):
		self.it = it
		#self.read_inp_files(it)
		#self.check_inp_files(it)

	def read_inp_files(self, conf):
		# read molecule file
		molecules = IO.readMolData(conf.molDataFile)
		# read prm file
		atomTypes = IO.readPrm(conf.prmFile)
		# read prm_nei file
		IO.readPrmNei(atomTypes, conf.prmNeiFile)
		# read matrix file
		matrix = IO.readMatrix(conf.matrixFile)
		IO.readListAtom(conf, molecules, conf.atomListFile)
		bnd, ang = IO.readRefBondAndAngle(conf.ifpFile)
		IO.readListBond(bnd, molecules, conf.bondListFile)
		IO.readListAngle(ang, molecules, conf.angListFile)
		IO.readSymmetry(atomTypes, "sig", conf.symSigFile)
		IO.readSymmetry(atomTypes, "eps", conf.symEpsFile)

		molecules_utils.createEffectivePrms(atomTypes, molecules, conf.eem, matrix)

		return molecules, atomTypes

	def run(self, conf):
		print("TODO: writeSam")


class Gen(Action):
	def __init__(self, it):
		Action.__init__(self, it)

	def run(self, conf, molecules, atomTypes):
		writeOutFiles.writeSam(self.it, molecules)

		molecules_utils.computeEEM(conf.eem, molecules, atomTypes)
		# TODO - read matrix file
		molecules_utils.computeCR(conf.cr, molecules, atomTypes, conf.matrix)

		writeOutFiles.writeParamMod(self.it, molecules)


class Ana(Action):
	def __init__(self, it):
		Action.__init__(self, it)


class Opt(Action):
	def __init__(self, it):
		Action.__init__(self, it - 1)
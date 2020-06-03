import molecules_utils
import writeOutFiles
import IO


class Action:
	def __init__(self, it):
		self.it = it
		#self.read_inp_files(it)
		#self.check_inp_files(it)

	def read_inp_files(self, conf):
		molecules = IO.readMolData(conf.molDataFile)
		atomTypes = IO.readPrm(conf.prmFile)
		IO.readPrmNei(atomTypes, conf.prmNeiFile)
		matrix = IO.readMatrix(conf.matrixFile)
		IO.readListAtom(conf, molecules, conf.atomListFile)
		bnd, ang = IO.readRefBondAndAngle(conf.ifpFile)
		IO.readListBond(bnd, molecules, conf.bondListFile)
		IO.readListAngle(ang, molecules, conf.angListFile)
		IO.readSymmetry(atomTypes, "sig", conf.symSigFile)
		IO.readSymmetry(atomTypes, "eps", conf.symEpsFile)

		molecules_utils.setCG(molecules, conf.eem)
		molecules_utils.createEffectivePrms(molecules, conf.eem)

		return molecules, atomTypes

	def run(self, conf):
		print("TODO: writeSam")


class Gen(Action):
	def __init__(self, it):
		Action.__init__(self, it)

	def run(self, conf, molecules, atomTypes):
		writeOutFiles.writeSam(self.it, molecules)

		molecules_utils.computeEEM(conf.eem, molecules, atomTypes)

		writeOutFiles.writeParamMod(self.it, atomTypes, molecules)




class Ana(Action):
	def __init__(self, it):
		Action.__init__(self, it)


class Opt(Action):
	def __init__(self, it):
		Action.__init__(self, it - 1)
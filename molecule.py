import sys


from effectivePrm import EffectivePrm


class Molecule:
	def __init__(self, cod, frm, run, pre_sim, tem_sim, dns_wei, dns_ref, hvp_wei, hvp_ref, mlp_ref, blp_ref, eps_ref):
		self._cod     = cod
		self._frm     = frm
		self._run     = float(run)
		self._pre_sim = float(pre_sim)
		self._tem_sim = float(tem_sim)
		self._dns_wei = float(dns_wei)
		self._dns_ref = float(dns_ref)
		self._hvp_wei = float(hvp_wei)
		self._hvp_ref = float(hvp_ref)
		self._mlp_ref = float(mlp_ref)
		self._blp_ref = float(blp_ref)
		self._eps_ref = float(eps_ref)

		self._atoms = {}
		self._LJPairs = []
		self._CGs = []

	@property
	def bondListFile(self):
		return self._bondListFile

	@property
	def cod(self):
		return self._cod
	@property
	def frm(self):
		return self._frm
	@property
	def run(self):
		return self._run
	@property
	def pre_sim(self):
		return self._pre_sim
	@property
	def tem_sim(self):
		return self._tem_sim
	@property
	def dns_wei(self):
		return self._dns_wei 
	@property
	def dns_ref(self):
		return self._dns_ref
	@property
	def hvp_wei(self):
		return self._hvp_wei
	@property
	def hvp_ref(self):
		return self._hvp_ref
	@property
	def mlp_ref(self):
		return self._mlp_ref
	@property
	def blp_ref(self):
		return self._blp_ref
	@property
	def eps_ref(self):
		return self._eps_ref
	@property
	def atoms(self):
		return self._atoms
	@property
	def LJPairs(self):
		return self._LJPairs

	@property
	def CGs(self):
		return self._CGs

	@CGs.setter
	def CGs(self, n):
		self._CGs = n

	def addAtom(self, atom, conf):
		try:
			idx = atom.idx
		except KeyError:
			sys.exit('\tERROR: Atom has no idx')

		if conf.ignoreAtom(atom.iac):
			atom.ignore = True

		self._atoms[idx] = atom

	def getAtom(self, idx):
		idx = int(idx)
		return self._atoms[idx]

	def checkAtoms(self):
		if len(self._atoms) == 0:
			sys.exit('No atoms were found for {}'.format(self._cod))

	def createEffectivePrms(self, eem):
		self.createLJPairs()

	def createLJPairs(self):
		iacList = []
		for idx in self._atoms:
			iac = self._atoms[idx].iac
			if not iac in iacList:
				iacList.append(iac)

		iacList = sorted(iacList)
		pairs = []
		N = len(iacList)
		for i in range(N):
			for j in range(i, N):
				iac1 = iacList[i]
				iac2 = iacList[j]
				pairs.append([iac1, iac2])
		print(pairs)
		print('\n\tHERE')
		sys.exit(666)



	def writePrmMod(self, f):
		# write charges
		for idx in self._atoms:
			atom = self._atoms[idx]
			# print(atom.curChg)
			f.write('{0:4} {1:>7} {2:>3} {3:>3} {4:>5} {5:3} {6:>15.4f} {7:13.4f}\n'
			.format(idx, 'CHG_ATM', 1, idx, 'MOLEC', atom.nam, 0.0, atom.curChg))

		# write LJ
		print(self._LJPairs)


import traceback
import sys

from configuration import Conf
import myExceptions
import inpParser


def main():
	"""
	Pipeline to run optimization of force-field parameters.
	The options are:
		GEN: generates param.mod files to run simulations.
		ANA: performs analysis of the simulation results.
		OPT: optimizes force-field parameters.
		SUB: submits jobs to run on euler (https://scicomp.ethz.ch/wiki/Main_Page).
		PLOT: plots main simulation results.
		MTB: Creates GROMOS coordinate and topology (mtb) files.
	"""

	try:
		action = inpParser.parse()
		conf = Conf('00_inp/Conf.dat', action.it)

		molecules, atomTypes = action.read_inp_files(conf)
		action.run(conf, molecules, atomTypes)

	except NotImplementedError as err:
		print(err)

	except myExceptions.NoSuchFile as err:
		print(err)
		sys.exit(1)

	except myExceptions.MissingKeyWordError as err:
		print(err)
		sys.exit(1)

	except KeyError as err:
		# for more information about the lines use traceback
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=4)
		print('\n\tKeyError: {}'.format(err))
		sys.exit(1)


if __name__ == "__main__":
	main()

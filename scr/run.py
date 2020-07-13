import traceback
import sys

from scr.configuration import Conf
from scr import myExceptions, inpParser


def main():
	try:
		action = inpParser.parse()
		conf = Conf('00_inp/Conf.dat', action.it)

		molecules, atomTypes = action.read_inp_files(conf)
		action.run(conf, molecules, atomTypes)

	except (FileNotFoundError) as err:
		print(err)
		sys.exit(1)

	except (myExceptions.MissingKeyWordError) as err:
		# for more information about the lines use traceback
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=4)
		print(err)
		sys.exit(1)

	except (KeyError) as err:
		print('\n\tKeyError: {}'.format(err))
		sys.exit(1)



if __name__ == "__main__":
	main()
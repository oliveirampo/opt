"""Parse arguments given to the command line.

Methods:
	parse()
"""

import argparse

from scr.base import action


def parse():
	"""Parse arguments given to the command line.

	:return:
		act: (Action) One of the possible options that the pipeline can execute.
		See Action for options.
	"""

	choices = action.Action.get_choices()

	parser = argparse.ArgumentParser()
	parser.add_argument("--type", dest="runType", required=True, type=str, choices=choices, help="Type of task to be done")
	parser.add_argument("--it", dest="it", required=True, type=int, help="Iteration")

	args = parser.parse_args()
	runType = args.runType
	it = args.it

	act = action.Action.get_object(runType, it)

	return act

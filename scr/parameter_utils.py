"""Module for handling parameter information.

Methods:
    getType(iac, atomTypes)
"""


import sys


def getType(iac, atomTypes):
    """Returns type of parameter given atom type index.

    :param iac: (int) Atom type index.
    :param atomTypes: (collections.OrderedDict) Ordered dictionary of atom types.
    :return: (str) Type of parameter.
    """

    if iac not in atomTypes:
        print('No such IAC = {}'.format(iac))
        sys.exit(123)

    return atomTypes[iac].typ


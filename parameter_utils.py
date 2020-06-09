import sys

def getType(iac, atomTypes):
    if not iac in atomTypes:
        print('No such IAC = {}'.format(iac))
        sys.exit(123)

    return atomTypes[iac].typ


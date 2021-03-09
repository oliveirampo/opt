import shutil
import sys
import os


def run(molecules):
    """Setups up default input files.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    """

    copyFiles(molecules)

    cmd = 'bash ../scr/make_lists.sh'
    os.system(cmd)


def copyFiles(molecules):
    """Copy top and cfg files.

    :param molecules: (collections.OrderedDict) Ordered dictionary of molecules.
    """

    print('Copying top and cfg files...')

    for cod in molecules:
        short_cod = cod[:5]

        srcTop = 'top/{}.top'.format(short_cod)
        srcCfgLiq = 'cfg/00_{}_liq_ini.cfg'.format(short_cod)

        trgTop = 'top/{}.top'.format(cod)
        trgCfgLiq = 'cfg/00_{}_liq_ini.cfg'.format(cod)
        trgCfgGas = 'cfg/00_{}_gas_ini.cfg'.format(cod)

        if not os.path.exists(srcTop):
            sys.exit('\n\tNo such file: {}'.format(srcTop))
        if not os.path.exists(srcCfgLiq):
            sys.exit('\n\tNo such file: {}'.format(srcCfgLiq))

        shutil.copyfile(srcTop, trgTop)
        shutil.copyfile(srcCfgLiq, trgCfgLiq)
        shutil.copyfile(srcCfgLiq, trgCfgGas)

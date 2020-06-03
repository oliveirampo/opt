from shutil import copyfile
import sys
import os

def writeSam(it, molecules):
    samDir = 'sam_' + str(it)
    if not os.path.exists(samDir):
        os.makedirs(samDir)

    cfgDir = 'cfg/'
    if not os.path.exists(cfgDir):
        os.makedirs(cfgDir)

    for cod in molecules:
        top    = 'top/' + cod + '.top'
        cfgGas = 'cfg/00_' + cod + '_gas_ini.cfg'
        cfgLiq = 'cfg/00_' + cod + '_liq_ini.cfg'

        # copy cnf files
        copyfile(cfgGas, cfgDir + cod + '_gas_ini.cfg')
        copyfile(cfgLiq, cfgDir + cod + '_liq_ini.cfg')

        samFile = open(samDir + '/' + cod + '.mol', 'w')

        mol = molecules[cod]
        liq_pre = mol.pre_sim * 0.06022141790

        samFile.write('mol_nam = {}\n'.format(cod))
        samFile.write('prm_mod = param_{}\n'.format(cod))
        samFile.write('liq_eps = {}\n'.format(mol.eps_ref))
        samFile.write('liq_pre = {0:.6e}\n'.format(liq_pre))
        samFile.write('liq_tem = {}\n'.format(mol.tem_sim))
        samFile.write('gas_tem = {}\n'.format(mol.tem_sim))
        samFile.write('bar_pre = {0:.6e}\n'.format(liq_pre))
        samFile.write('the_tem = {}\n'.format(mol.tem_sim))

        samFile.close()

def writeParamMod(it, atomTypes, molecules):
    prmDir = 'prm_{}'.format(it)
    if not os.path.exists(prmDir):
        os.makedirs(prmDir)

    outDir = 'out_{}'.format(it)
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    for cod in molecules:
        mol = molecules[cod]

        prmModFile = '{}/param_{}.mod'.format(prmDir, cod)
        with open(prmModFile, 'w') as f:
            mol.writePrmMod(f)

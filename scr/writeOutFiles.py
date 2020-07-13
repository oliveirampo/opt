from shutil import copyfile
import sys
import os


import myExceptions


def copyCOTO(molecules):
    cfgDir = 'cfg/'
    if not os.path.exists(cfgDir):
        os.makedirs(cfgDir)

    for cod in molecules:
        cfgGas = 'cfg/00_' + cod + '_gas_ini.cfg'
        cfgLiq = 'cfg/00_' + cod + '_liq_ini.cfg'

        # copy cnf files
        copyfile(cfgGas, cfgDir + cod + '_gas_ini.cfg')
        copyfile(cfgLiq, cfgDir + cod + '_liq_ini.cfg')


def writeMolFile(it, molecules):
    samDir = 'sam_' + str(it)
    if not os.path.exists(samDir):
        os.makedirs(samDir)

    for cod in molecules:
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


def writeSamFile(conf, molecules, samTemplate):
    it = conf.it
    nJobs = conf.nJobs

    samDir = 'sam_' + str(it)
    if not os.path.exists(samDir):
        os.makedirs(samDir)

    for cod in molecules:

        for i in range(nJobs):
            outFileName = 'sam_{}/{}_{}.sam'.format(it, cod, i)

            # production
            num_stp = conf.prd_stp
            prt_res_frq = conf.prd_frq
            prt_prp_frq = conf.prd_frq

            # equilibration
            if i == 0:
                num_stp = conf.eq_stp
                prt_res_frq = conf.eq_frq
                prt_prp_frq = conf.eq_frq

            # write sam file
            with open(outFileName, 'w') as out:
                j = 0
                lines = samTemplate[:]
                while j < len(lines):
                    row = lines[j]

                    if len(row) == 0:
                        j += 1

                    elif row[0] == 'num_stp':
                        lines[j][2] = num_stp

                    elif row[0] == 'prt_res_frq':
                        lines[j][2] = prt_res_frq

                    elif row[0] == 'prt_prp_frq':
                        lines[j][2] = prt_prp_frq

                    s = ' '.join(lines[j])
                    out.write('{}\n'.format(s))

                    j += 1


def writeParamMod(it, molecules):
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


def writeSubScript(conf, molecules):
    it = conf.it
    nJobs = conf.nJobs
    wall_time = conf.wall_time

    datDir = 'dat_' + str(it)
    if not os.path.exists(datDir):
        os.makedirs(datDir)

    subDir = 'sub_' + str(it)
    if not os.path.exists(subDir):
        os.makedirs(subDir)

    # create scr files with bash script
    for cod in molecules:
        cmd = '../scr/prepare.sh {} {} {} {}' \
            .format(cod, it, nJobs, wall_time)
        os.system(cmd)

    # check if all files were created
    for cod in molecules:
        submitFile = 'sub_{}/{}_{}_sub.sh'.format(it, cod, 0)
        if not os.path.exists(submitFile):
            raise myExceptions.NoSuchFile(submitFile)

    # submit scripts
    for cod in molecules:
        submitFile = 'sub_{}/{}_{}_sub.sh'.format(it, cod, 0)
        cmd = '../scr/bsub.sh {} {} {} {}'.format(cod, it, 0, wall_time)
        #os.system(cmd)




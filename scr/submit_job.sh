#!/bin/sh

SIMULDIR=$PWD
echo "Simulation directory ${SIMULDIR}"

# first we set some variables
NAME=`whoami`

# samos dir
INPDIR=INP_DIR_MAPE
PRMDIR=PRMDIR_MAPE

CMP=CMP_MAPE
IT=IT_MAPE
JOB=JOB_MAPE

PRM_MOD_DIR=${INPDIR}/prm_${IT}
TOPDIR=${INPDIR}/top
CFGDIR=${INPDIR}/cfg

#input
TOPO=${CMP}.top
PRMMOD=param_${CMP}.mod
#MOLFILE=MOL_FILE_MAPE
SAMFILE=SAM_FILE_MAPE

#output
OUTPUT=${CMP}_${JOB}_lag.out
OUT=o_${CMP}_${JOB}
SEN=s_${CMP}_${JOB}

# create temporary directory
WORKDIR=${TMPDIR}/${CMP}

mkdir -p ${WORKDIR}
cd       ${WORKDIR}
echo "Working directory ${WORKDIR}"

mkdir -p prm/
mkdir -p top/
mkdir -p cfg/

echo ${PRMDIR}

cp ${PRMDIR}/2016H66_upd.ifp	prm/
cp ${PRM_MOD_DIR}/${PRMMOD}	prm/
cp ${TOPDIR}/${TOPO}		top/
cp ${CFGDIR}/${CMP}_gas_ini.cfg	cfg/
cp ${CFGDIR}/${CMP}_liq_ini.cfg	cfg/

#SAM=/cluster/work/igc/marinap/opt_setup_19_06/bin/sam_19_02_trc_1000
#SAM=/cluster/work/igc/marinap/opt_setup_19_06/bin/sam_19_10_trc_1000
SAM=/cluster/work/igc/marinap/opt_setup_19_06/bin/sam_20_01_thermostat_test

${SAM} ${SAMFILE} > ${OUTPUT}

# copy the file back
#cat ${OUTPUT}
#cp ${OUTPUT} ${SIMULDIR}/dat_${IT}/${CMP}_${JOB}.log
tail ${OUTPUT} > ${SIMULDIR}/out_${IT}/${CMP}_${JOB}.log
TMPFILE=${OUTPUT}.tmp

if [ ! -f  ${OUTPUT} ]; then
        echo "ERROR: no output"
	cd ${SIMULDIR};
	rm -rf ${WORKDIR}
	exit
fi

egrep '<ave{D_liq}>  |<ave{dH_vap}>  |<sen_dave{D_liq}/dP> [0-9]|<sen_dave{dH_vap}/dP> [0-9]' $OUTPUT > ${TMPFILE}

cat ${TMPFILE} |\
egrep '<ave{D_liq}>  |<ave{dH_vap}>  ' |\
sed 's/:/ /g' | awk '{printf "%-8s %8.3f %9.2f\n", $5, $3*0.002, $NF}' |\
sed 's/<ave{D_liq}>/D/g;s/<ave{dH_vap}>/H/g' > ${SIMULDIR}/dat_${IT}/${OUT}

cat ${TMPFILE} |\
egrep '<sen_dave{D_liq}/dP> [0-9]|<sen_dave{dH_vap}/dP> [0-9]' |\
sed 's/:/ /g' | awk '{printf "%-8s %8.3f ", $5, $3*0.002; for(i=7;i<=NF;i++){printf"%s ", $i};printf"\n"}' |\
sed 's/<sen_dave{D_liq}[/]dP> */D /g;s/<sen_dave{dH_vap}[/]dP> */H /g' > ${SIMULDIR}/dat_${IT}/${SEN}

#grep tem $OUTPUT > ${SIMULDIR}/dat_${IT}/tem_${CMP}_${JOB}
#grep pre $OUTPUT > ${SIMULDIR}/dat_${IT}/pre_${CMP}_${JOB}

# No output
err=`grep ERROR ${OUTPUT}`
if [ ! -z "$err" ]; then
        echo $err
	cd ${SIMULDIR};
	rm -rf ${WORKDIR}
	exit
fi

FILE=cfg/${CMP}_liq_fin.cfg
if [ ! -f "$FILE"  ]; then
        echo "ERROR: no liquid cnf"
	cd ${SIMULDIR};
	rm -rf ${WORKDIR}
	exit
fi
FILE=cfg/${CMP}_gas_fin.cfg
if [ ! -f "$FILE"  ]; then
        echo "ERROR: no gas cnf"
	cd ${SIMULDIR};
	rm -rf ${WORKDIR}
	exit
fi

# copy the files back
mv cfg/${CMP}_gas_fin.cfg ${SIMULDIR}/cfg/${CMP}_gas_ini.cfg
mv cfg/${CMP}_liq_fin.cfg ${SIMULDIR}/cfg/${CMP}_liq_ini.cfg

#mv ${CMP}.trc ${SIMULDIR}/trc/${CMP}_${JOB}.cfg

# clean up after successful run
cd ${SIMULDIR};
rm -rf ${WORKDIR}

# perform last command (usually submit next job)






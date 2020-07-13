#!/bin/bash


cmp=$1
it=$2
nJobs=$3
wallTime=$4

localInpDir=${PWD}
localInpDirSed=${localInpDir//\//\\\/}

subDir="sub_${it}"
mkdir -p $subDir

file=submit_job.sh
localInpDir=${PWD}
localInpDirSed=${localInpDir//\//\\\/}
prmDir=${PWD}/prm
prmDirSed=${prmDir//\//\\\/}

molFile=${localInpDir}/sam_${it}/${cmp}.mol
molFileSed=${molFile//\//\\\/}

samDir=${PWD}/sam_${it}


for i in `seq 0 $nJobs`; do
	f=${subDir}/${cmp}_${i}_sub.sh
	cp ../scr/$file $f

	sed -i "s/INP_DIR_MAPE/${localInpDirSed}/g"	$f
	sed -i "s/PRMDIR_MAPE/${prmDirSed}/g"		$f
	sed -i "s/CMP_MAPE/${cmp}/g"			$f
	sed -i "s/IT_MAPE/${it}/g"			$f
	sed -i "s/JOB_MAPE/${i}/g"			$f

	sed -i "s/MOL_FILE_MAPE/${molFileSed}/g"	$f

	if [ $i -ne $nJobs ]; then
		j=$((i+1))
		echo "cd $localInpDir" >> $f
		echo "bsub -W ${wallTime} -n 1 -J ${cmp}_${it}_${j} -e out_${it}/${cmp}_${j}.err -o out_${it}/${cmp}_${j}.out < sub_${it}/${cmp}_${j}_sub.sh" >> $f
	fi
done


#Generate the lists of bonds, angle and atoms including H atoms

TOPDIR=top
INPDIR=00_inp

if [ ! -d "$TOPDIR" ]; then
        echo "$TOPDIR does not exist"
fi

if [ ! -d "$INPDIR" ]; then
        echo "$INPDIR does not exist"
fi

cd $TOPDIR || exit

function mount {
  file=$1
  < "$file" column -t > tempXXX.dat
  mv tempXXX.dat "$file"

}

#---------------------------
#List Bonds
#---------------------------

echo "Creating listBond.dat"

for f in [A-Z][0-9][0-9][0-9][[:alnum:]][a-z].top; do cod=(${f//./ }); awk '/NBON:/{flag=1;next}/END/{flag=0}flag' "${cod[0]}".top | awk -v cod="${cod[0]}" '{if (NF==3) print cod, $1, $2, $3}' ; done > listBond.dat

for f in [A-Z][0-9][0-9][0-9][[:alnum:]][a-z].top; do cod=(${f//./ }); awk '/NBONH:/{flag=1;next}/END/{flag=0}flag' "${cod[0]}".top | awk -v cod="${cod[0]}" '{if (NF==3) print cod, $1, $2, $3}' ; done >> listBond.dat

sort -V listBond.dat -o listBond.dat
mount listBond.dat

#---------------------------
#List Angles
#---------------------------

echo "Creating listAng.dat"

for f in [A-Z][0-9][0-9][0-9][[:alnum:]][a-z].top; do cod=(${f//./ }); awk '/NTHE:/{flag=1;next}/END/{flag=0}flag' "${cod[0]}".top | awk -v cod="${cod[0]}" '{if (NF==4) print cod, $1, $2, $3, $4}' ; done > listAng.dat

for f in [A-Z][0-9][0-9][0-9][[:alnum:]][a-z].top; do cod=(${f//./ }); awk '/NTHEH:/{flag=1;next}/END/{flag=0}flag' "${cod[0]}".top | awk -v cod="${cod[0]}" '{if (NF==4) print cod, $1, $2, $3, $4}' ; done >> listAng.dat

sort -V listAng.dat -o listAng.dat
mount listAng.dat

#--------------------------
#List Atoms
#--------------------------

echo "Creating listAtom.dat"

for f in [A-Z][0-9][0-9][0-9][[:alnum:]][a-z].top; do cod=(${f//./ }); awk '/  INE14/{flag=1;next}/END/{flag=0}flag' "${cod[0]}".top | awk -v cod="${cod[0]}" '{if (NF>6 && $7 == 0 || $7 == 1) print cod, $1, $3, $4, $6}' ; done > listAtom.dat

sort -V listAtom.dat -o listAtom.dat
mount listAtom.dat



mv list*.dat ../$INPDIR
cd ../

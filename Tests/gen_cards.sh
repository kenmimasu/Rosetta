#!/bin/bash
EXEC=../Rosetta/bin/rosetta

for BASIS in HiggsBasis SILHBasis WarsawBasis HISZ TemplateBasis;
    do NAME=`echo ${BASIS/Basis/}  | tr "[:upper:]" "[:lower:]"`;
    for FLAV in general diagonal universal;
        do $EXEC defaultcard -w -o ../Rosetta/Cards/${BASIS}_${FLAV}.dat --flavor $FLAV $NAME;
        $EXEC  defaultcard -w -o Cards/${BASIS}_${FLAV}.dat --flavor $FLAV $NAME;
        $EXEC  defaultcard -w --value 0.001 -o Cards/${BASIS}_${FLAV}_1e-3.dat --flavor $FLAV $NAME;
        $EXEC  defaultcard -w --value "random" -o Cards/${BASIS}_${FLAV}_rand.dat --flavor $FLAV $NAME;
    done;
done

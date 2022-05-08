#!/bin/bash
# File : test.sh
# Author : Adam Švenk
# Date : 2022-05-05
# Version : 1.0

SOURCE='pro.cpp'
OUTPUT='pro'

if [ $# != 1 ]
    then
        echo 'Error: Wrong parameters'
        exit -1
    else
        INPUT=$1
        PROCESSES=${#INPUT}
        PROCESSES=$PROCESSES
        PROCESSES=$(expr $PROCESSES \* 2 - 2)
fi

#Compile pro.cpp
mpic++ --prefix /usr/local/share/OpenMPI -o $OUTPUT $SOURCE

#Launch
mpirun --prefix /usr/local/share/OpenMPI -oversubscribe -np $PROCESSES $OUTPUT $INPUT

rm $OUTPUT
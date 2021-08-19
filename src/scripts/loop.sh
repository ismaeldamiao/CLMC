#!/bin/bash

if ! [ -x ./clmc ]; then 
   echo "Nao foi possivel encontrar o executavel.\nNao te esquecestes de compilar o programa?"
   exit 1
fi

#array $cpus, o valor do número de cpus vai estar no índice 1 do array.
cpus=( $( lscpu | grep '^CPU(s):' ) )

JOBMAX=${cpus[1]}

for (( i = ${1}; i <= ${2}; ++i )); do
   for (( ; ; )); do
      JOBS=( $( jobs -p ) )
      if [ ${#JOBS[@]} -ge ${JOBMAX} ]; then
         sleep 60
      else
         time ./clmc ${i} &
         break
      fi
   done
done
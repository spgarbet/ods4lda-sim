#!/bin/bash

DIR="output"

echo "Batch 1"
for i in $(seq 1 2000)
do
    if [ ! -f "${DIR}/run-$i-1.RData"  ]; then
      echo -ne "$i,"
    fi
done

echo "Batch 2"
for i in $(seq 1 2000)
do
    if [ ! -f "${DIR}/run-$i-2.RData"  ]; then
      echo -ne "$i,"
    fi
done 

echo "Batch 3"
for i in $(seq 1 2000)
do
    if [ ! -f "${DIR}/run-$i-3.RData"  ]; then
      echo -ne "$i,"
    fi
done 

echo "Batch 4"
for i in $(seq 1 2000)
do
    if [ ! -f "${DIR}/run-$i-4.RData"  ]; then
      echo -ne "$i,"
    fi
done

echo ""

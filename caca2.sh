#!/bin/bash
echo "100000" > /tmp/como2.txt
for i in $(seq 100000)
do
	aer=$(( (i-1)*2+1 ))
	echo "$aer" >> /tmp/como2.txt
done
echo "100000" >> /tmp/como2.txt
for i in $(seq 100000)
do
	aer=$(( i*2 ))
	echo "$aer" >> /tmp/como2.txt
done

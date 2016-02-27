#!/bin/bash
echo "200000" > /tmp/como.txt
for i in $(seq 200000)
do
	echo "$i" >> /tmp/como.txt
done
echo "200000" >> /tmp/como.txt
for i in $(seq 200000)
do
	echo "$i" >> /tmp/como.txt
done

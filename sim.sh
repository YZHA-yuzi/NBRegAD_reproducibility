#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------

for((i=1;i<=100;++i))
do
	for ((j=1;j<=5;++j))
	do 
		for ((k=1;k<=2;++k))
		do
			Rscript Simulations.R "$i" "$j" "$k"
		done
	done
done
#!/bin/bash
#
#$ -l hostname=micro4.local 
#$ -N SM1b-2v-2
#$ -l h_rt=180:00:00                   #estimate max run time
#$ -q all.q
#$ -m ea
#$ -M lucasrdf@fisica.ufrn.br
#$ -cwd
#$ -o /home/rodrigues/Mathematica2/.o
#$ -e /home/rodrigues/Mathematica2/.e
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/sge/lib/lx-amd64
export CLASSPATH=/opt/sge/lib/drmaa.jar:/opt/sge/lib/juti.jar:/opt/sge/lib/jgdi.jar
#
math -script 705-SM1b.m

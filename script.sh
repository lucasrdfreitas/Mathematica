#!/bin/bash
#
#$ -l hostname=micro4.local 
#$ -N MF_K
#$ -q all.q
#$ -m ea
#$ -M lucasrdf@fisica.ufrn.br
#$ -cwd
#$ -o /home/rodrigues/Mathematica/.o
#$ -e /home/rodrigues/Mathematica/.e
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/sge/lib/lx-amd64
export CLASSPATH=/opt/sge/lib/drmaa.jar:/opt/sge/lib/juti.jar:/opt/sge/lib/jgdi.jar
#
math -script 705-MF_V0_CLUSTER.m

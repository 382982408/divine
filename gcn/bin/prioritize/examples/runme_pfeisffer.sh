#!/bin/bash -l

echo "../divine.py -q ./Pfeiffer.hpo -v ./Pfeiffer.vcf -o ./Pfeiffer -e 0 -c ../../../config/filterconf_dp10.txt --reuse"

../divine.py -q ./Pfeiffer.hpo -v ./Pfeiffer.vcf -o ./Pfeiffer -e 0 -c ../../../config/filterconf_dp10.txt --reuse

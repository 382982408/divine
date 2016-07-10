#!/bin/bash -l

echo "../divine.py -v ./Pfeiffer.vcf -o ./Pfeiffer_noHpo -e 0 -c ../../../config/filterconf_dp10.txt --reuse"

../divine.py -v ./Pfeiffer.vcf -o ./Pfeiffer_noHpo -e 0 -c ../../../config/filterconf_dp10.txt --reuse

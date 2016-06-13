#!/bin/bash -l
export SRC=$1
export DST=$2

rsync -ravP --exclude="python_libs/" --exclude="gcndb/" --exclude="gcndata/" --exclude="gcn/bin/prioritize/examples/" $SRC/divine-0.1.1 $DST

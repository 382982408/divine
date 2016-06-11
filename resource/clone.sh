#!/bin/bash -l
export SRC=/media/hong/blackedge/projects/apps
export DST=/media/hong/redfish/projects/apps

rsync -ravP $SRC/divine/gcn $DST/divine/
rsync -ravP $SRC/divine/python_libs $DST/divine/
rsync -ravP $SRC/divine/gcndb $DST/divine/
rsync -ravP $SRC/divine/gcndata/{hpo,disgenet,stringDB,snv_training} $DST/divine/gcndata/
rsync -ravzP $SRC/divine/setup.py $DST/divine/
rsync -ravzP $SRC/divine/documents $DST/divine/

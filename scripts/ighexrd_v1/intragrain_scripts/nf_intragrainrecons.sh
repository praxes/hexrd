#!/bin/sh
# -*- coding: utf-8 -*-

cd /nfs/chess/aux/user/ken38/Ti7_project/paper-near-field-run/step_1_near_field_scripts/nf_3_final_stitched_orientations/
source activate hexrd_0526_nftest

echo starting final gbg nf
#python nf_gbg_final_diffvol_1.py
#python nf_gbg_final_diffvol_2.py
python nf_gbg_final_diffvol_3.py
python nf_gbg_final_diffvol_4.py
echo Done


#!/bin/sh
# -*- coding: utf-8 -*-


cd /workingdirectory/

source activate hexrd_environment

echo starting near-field uniform orientation field

echo diff_vol_1
python near-field_uniform_diffvol_1.py

echo diff_vol_2
python near-field_uniform_diffvol_2.py

echo diff_vol_3
python near-field_uniform_diffvol_3.py

echo diff_vol_4
python near-field_uniform_diffvol_4.py

echo diff_vol_5
python near-field_uniform_diffvol_5.py

echo All done

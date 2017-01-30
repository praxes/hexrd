import os

from hexrd.imageseries import omega
from pp_dexela import PP_Dexela

CHESS_BASE = '/nfs/chess/raw/current/f2/shade-560-1/LSHR-6'
CHESS_TMPL = '%d/ff/ff2_%05d.h5'

def h5name(scan, file, base=CHESS_BASE):
    path = CHESS_TMPL % (scan, file)
    return os.path.join(base, path)

# ==================== Inputs (should not need to alter above this line)

## Room temp
#raw_scannumber = 32
#raw_filenumber = 35

# ROOM TEMP
raw_scannumber = 81
raw_filenumber = 45

## 100C 
#raw_scannumber = 82
#raw_filenumber = 46

## 300C 
#raw_scannumber = 83
#raw_filenumber = 47

flips = [('flip', 't'), ('flip', 'hv') ]

nframes = 100

ostart = 0
ostep = 0.25
fstart = 5
threshold = 150

# ==================== End Inputs (should not need to alter below this line)

input_name = h5name(raw_scannumber, raw_filenumber)
output_name = input_name.split('/')[-1].split('.')[0]

ostop = ostart + nframes*ostep
omw = omega.OmegaWedges(nframes)
omw.addwedge(ostart, ostop, nframes)

ppd = PP_Dexela(input_name, omw, flips, frame_start=fstart)
ppd.save_processed(output_name, threshold)

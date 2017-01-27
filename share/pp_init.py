from hexrd.imageseries import omega
from pp_dexla import PP_Dexela

raw_fname = '/nfs/chess/raw/current/f2/shade-560-1/LSHR-6/%d/ff/ff2_%05d.h5'
raw_scannumber = 32
raw_filenumber = 35

input_name = raw_fname %(raw_scannumber, raw_filenumber)
output_name = input_name.split('/')[-1].split('.')[0]

nframes = 1440
ostart = 0
ostep = 0.25
fstart = 5
threshold = 150

# ==================== End Inputs (should not need to alter below this line)

ostop = ostart + nframes*ostep
omw = omega.OmegaWedges(nframes)
omw.addwedge(ostart, ostop, nframes)


ppd = PP_Dexela(input_name, omw, frame_start=fstart)
ppd.save_processed(output_name, threshold)

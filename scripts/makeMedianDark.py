#!/usr/bin/env python
import sys, os

import numpy as np

from hexrd import config
from hexrd.fitgrains import get_instrument_parameters

def make_dark_frame(cfg, nframes_use=100, nbytes_header=8192, pixel_type=np.uint16):
    """
    FOR BYTE STREAM IMAGES ONLY!

    nbytes_header and pixel_bytes default for raw GE
    """
    
    nempty = cfg.image_series.images.start
    image_number = cfg.image_series.file.ids[0]

    instr_cfg = get_instrument_parameters(cfg)

    nrows = instr_cfg['detector']['pixels']['rows']
    ncols = instr_cfg['detector']['pixels']['columns']

    pixel_bytes = pixel_type().itemsize
    nbytes_frame = pixel_bytes*nrows*ncols    # uint16

    fid = open(cfg.image_series.file.stem %cfg.image_series.file.ids[0], 'rb')
    fid.seek(nbytes_header+nempty*nbytes_frame, 0)     # header plus junk frames
    frames = np.frombuffer(fid.read(nframes_use*nbytes_frame), dtype=pixel_type).reshape(nframes_use, nrows, ncols)
    fid.close()
    
    dark = np.median(frames, axis=0)
    
    # usi UNIX style path spitting...
    image_stem = cfg.image_series.file.stem.split('/')[-1].split('.')
    if len(image_stem) > 1:
        suffix = '.'+image_stem[-1]
    else:
        suffix = ''
    output_name = os.path.join(cfg.working_dir, 
                               image_stem[0] %image_number + '_medianDark' + suffix)

    fid = open(output_name, 'wb')
    fid.seek(nbytes_header)
    fid.write(dark.astype(np.uint16))
    fid.close()
    return dark.astype(np.uint16)
    
def subtract_dark_and_write(cfg, dark, filename='full_scan.ge'):
    archive = open(os.path.join(cfg.working_dir, filename), 'wb')
    archive.seek(nbytes_header)
    for f in cfg.image_series.files:
        fid = open(f, 'rb')
        fid.seek(nbytes_header, 0)     # header plus first junk frame
        for iframe in range(int(cfg.image_series.n_frames)):
            frame = np.frombuffer(fid.read(nbytes_frame), dtype=np.uint16).reshape(nrows, ncols)
            tmp = frame - dark
            tmp[tmp <= cfg.find_orientations.orientation_maps.threshold] = 0
            archive.write(tmp)
        fid.close()
    archive.close()
    return
    
if __name__ == '__main__':
    config_filename = sys.argv[1]
    cfg = config.open(config_filename)[0]
    
    if len(sys.argv) > 2:
        dark = make_dark_frame(cfg, nframes_use=int(sys.argv[2]))
    else:
        dark = make_dark_frame(cfg)
        #subtract_dark_and_write(cfg, dark)
    
    

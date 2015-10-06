#!/usr/bin/env python
import sys, os

import argparse

import yaml

import numpy as np

from hexrd import config
from hexrd.fitgrains import get_instrument_parameters

def make_dark_frame(cfg_name, nframes_use=100, nbytes_header=8192, pixel_type=np.uint16, quiet=False, output_name=None):
    """
    FOR BYTE STREAM IMAGES ONLY!

    nbytes_header and pixel_bytes default for raw GE
    """
    
    cfg = config.open(cfg_name)[0] # only first block...
    raw_cfg = [g for g in yaml.load_all(open(cfg_name, 'r'))] # take first block... maybe iterate?

    try:
        output_name = cfg.image_series.dark
    except IOError:
        if raw_cfg[0]['image_series'].has_key('dark') and output_name is None:
            output_name = raw_cfg[0]['image_series']['dark']

    if not os.path.exists(output_name):
        # create if necessary    
        open(output_name, 'w+').close()

    n_empty = cfg.image_series.images.start

    image_numbers = cfg.image_series.file.ids
    n_images = len(image_numbers)
    
    instr_cfg = get_instrument_parameters(cfg)

    nrows = instr_cfg['detector']['pixels']['rows']
    ncols = instr_cfg['detector']['pixels']['columns']

    pixel_bytes = pixel_type().itemsize
    nbytes_frame = pixel_bytes*nrows*ncols

    filename = cfg.image_series.file.stem %cfg.image_series.file.ids[0]
    file_bytes = os.stat(filename).st_size
    n_frames = getNFramesFromBytes(file_bytes, nbytes_header, nbytes_frame)
    n_frames -= n_empty
    n_frames_total = n_frames*n_images
    if nframes_use > n_frames_total:
        if not quiet:
            print "Requested %d frames, which is in image spec; defaulting to %d" %(nframes_use, n_frames_total)
        nframes_use = n_frames_total
    
    ii = 0
    jj = min(nframes_use, n_frames)

    n_frames_cum = 0
    n_frames_rem = nframes_use

    frames = np.empty((nframes_use, nrows, ncols), dtype=pixel_type)
    for i_frame in range(len(image_numbers)):
        filename = cfg.image_series.file.stem %cfg.image_series.file.ids[i_frame]

        if not quiet:
            print "Using %d frames from '%s'" %(min(n_frames, nframes_use, n_frames_rem), filename)

        nfr = jj - ii
        fid = open(filename, 'rb')
        fid.seek(nbytes_header+n_empty*nbytes_frame, 0)     # header plus junk frames
        tmp = np.frombuffer(fid.read(nfr*nbytes_frame), dtype=pixel_type).reshape(nfr, nrows, ncols)
        fid.close()
        
        # index into frames array
        frames[ii:jj] = tmp

        # increment...
        n_frames_rem -= nfr
        if n_frames_rem <= 0:
            break
        else:
            ii = jj
            jj += min(n_frames_rem, n_frames)
            if not quiet:
                print "\tframes left: %d" %n_frames_rem
        pass
    
    dark = np.median(frames, axis=0)

    if not quiet:
        print "Output file: %s" %output_name

    fid = open(cfg.image_series.dark, 'wb')
    fid.seek(nbytes_header)
    fid.write(dark.astype(pixel_type))
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


def getNFramesFromBytes(fileBytes, nbytesHeader, nbytesFrame):
    assert (fileBytes - nbytesHeader) % nbytesFrame == 0,\
        'file size not correct'
    nFrames = int((fileBytes - nbytesHeader) / nbytesFrame)
    if nFrames*nbytesFrame + nbytesHeader != fileBytes:
        raise RuntimeError, 'file size not correctly calculated'
    return nFrames

if __name__ == '__main__':
    """
    USAGE : python makeMedianDark <cfg_file> <num frames to use> <output filename>
    """
    parser = argparse.ArgumentParser(description='Make median dark from cfg file')

    parser.add_argument('cfg', metavar='cfg_filename', type=str, help='a YAML config filename')
    parser.add_argument('-n','--num-frames', help='Number of frames to use in median calculation', type=int, default=100)
    parser.add_argument('-o','--output-name', help='Output filename', type=str)

    args = vars(parser.parse_args(sys.argv[1:]))
    
    dark = make_dark_frame(args['cfg'], nframes_use=args['num_frames'], output_name=args['output_name'])
    
    

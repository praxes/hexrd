#! /usr/bin/env python
#
"""Make an imageseries from a list of image files
"""
import sys
import argparse
import logging

# Put this before fabio import and reset level if you
# want to control its import warnings.
logging.basicConfig(level=logging.DEBUG)

import numpy
import h5py
import fabio

# Error messages

ERR_NO_FILE = 'Append specified, but could not open file'
ERR_NO_DATA = 'Append specified, but dataset not found in file'
ERR_OVERWRITE = 'Failed to create new dataset. Does it already exist?'
ERR_SHAPE = 'Image shape not consistent with previous images'
ERR_NOEMPTY = 'dark-from-empty specified, but number of empty frames not given'
ERR_OMEGASPEC = 'Must specify both omega-min and omega-max if either is given'

DSetPath = lambda f, p: "%s['%s']" % (f, p)

class MakeImageSeriesError(Exception):
    """Class for MakeImageSeriesError  errors"""
    def __init__(self, message):
	self.message = message
        return

    def __str__(self):
	return self.message

    pass  # end class


def write_file(a, **kwargs):
    #
    # Get shape and dtype information from files
    #
    shp, dtp = image_info(a)
    #
    # Open file and dataset
    #
    f, ds = open_dset(a, shp, dtp)
    #
    # Image options
    #
    popts = process_img_opts(a, **kwargs)
    #
    # Now add the images
    # . empty frames only apply to multiframe images
    #
    nframes = ds.shape[0]
    nfiles = len(a.imagefiles)
    for i in range(nfiles):
        if a.max_frames and nframes >= a.max_frames:
            break
        logging.debug('processing file %d of %d' % (i, nfiles))
        popts['filenumber'] = i
        img_i = fabio.open(a.imagefiles[i])
        nfi = img_i.nframes
        for j in range(nfi):
            if a.max_frames and nframes >= a.max_frames:
                break
            logging.debug('... processing image %d of %d' % (j, img_i.nframes))
            if nfi > 1 and j < a.empty:
                logging.debug('...empty frame ... skipping')
                continue
            nframes += 1
            ds.resize(nframes, 0)
            ds[nframes - 1, :, :] = process_img(a, img_i.data, popts)
            if (j + 1) < nfi:
                img_i = img_i.next()
        pass

    add_metadata(ds, a, **kwargs)

    f.close()
    return

def open_dset(a, shp, dtp):
    """open HDF5 file and dataset"""
    #
    # If append option is true, file and target group must exist;
    # otherwise, file may exist but may not already contain the
    # target dataset.
    #
    if a.append:
        try:
            f = h5py.File(a.outfile, "r+")
        except:
            errmsg = '%s: %s' % (ERR_NO_FILE, a.outfile)
            raise MakeImageSeriesError(errmsg)

        ds = f.get(a.dset)
        if ds is None:
            errmsg = '%s: %s' % (ERR_NO_DATA, DSetPath(a.outfile, a.dset))
            raise MakeImageSeriesError(errmsg)
    else:
        f = h5py.File(a.outfile, "a")
        chsize = (1, int(numpy.floor(1e6/shp[1])), shp[1]) if shp[1] < 1.e6 else True
        try:
            ds = f.create_dataset(a.dset, (0, shp[0], shp[1]), dtp,
                                  maxshape=(None, shp[0], shp[1]), chunks=chsize,
                compression="gzip")
        except Exception as e:
            errmsg = '%s: %s\n...  exception: ' % (ERR_OVERWRITE, DSetPath(a.outfile, a.dset))
            raise MakeImageSeriesError(errmsg + str(e))

    return f, ds

def process_img_opts(a, **kwargs):
    """make dictionary to pass to process_img"""
    pdict = {}
    # dark file (possibly need to transpose [not done yet])
    if a.dark_file:
        dark = fabio.open(a.dark_file)
        pdict['dark'] = dark.data

    # dark from empty
    if a.dark_from_empty:
        if a.empty == 0:
            raise MakeImageSeriesError(ERR_NOEMPTY)
        darks = []
        for i in range(len(a.imagefiles)):
            img_i = fabio.open(a.imagefiles[i])
            drk_i = img_i.data
            for j in range(1, a.empty):
                img_i = img_i.next()
                drk_i += img_i.data
            darks += [drk_i*(1/a.empty)]
        pdict['darks'] = darks

    return pdict

def process_img(a, img, pdict):
    """process image data according to options

       * need to check on image shape in case not square
"""
    # flip: added some other option specifiers
    if a.flip in ('y','v'): # about y-axis (vertical)
        pimg = img[:, ::-1]
    elif a.flip in ('x', 'h'): # about x-axis (horizontal)
        pimg = img[::-1, :]
    elif a.flip in ('vh', 'hv', 'r180'): # 180 degree rotation
        pimg = img[::-1, ::-1]
    elif a.flip in ('t', 'T'): # transpose (possible shape change)
        pimg = img.T
    elif a.flip in ('ccw90', 'r90'): # rotate 90 (possible shape change)
        pimg = img.T[:, ::-1]
    elif a.flip in ('cw90', 'r270'): # rotate 270 (possible shape change)
        pimg = img.T[::-1, :]
    else:
        pimg = img

    # dark image(s)
    if 'dark' in pdict:
        pimg = pimg - pdict['dark']
    if 'darks' in pdict:
        fnum = pdict['filenumber']
        pimg = pimg - pdict['darks'][fnum]

    return pimg

def add_metadata(ds, a, **kwargs):
    """Add metadata. Right now, that just includes omega information."""
    omkey = 'omega'
    ominatt = 'omega_min'
    omaxatt = 'omega_max'
    hasval = lambda a, att: hasattr(a, att) and getattr(a, att) is not None
    hasone = lambda a, att1, att2: hasval(a, att1) and not hasval(a, att2)

    if (hasone(a, ominatt, omaxatt) or hasone(a, omaxatt, ominatt)):
        raise MakeImageSeriesError(ERR_OMEGASPEC)

    if not (hasval(a, ominatt) and hasval(a, omaxatt)):
        return

    omin = getattr(a, ominatt)
    omax = getattr(a, omaxatt)

    if a.append:
        om = ds.dims[0][omkey]
        n0 = om.shape[1]
        n1 = ds.shape[0]
    else:
        om = ds.parent.file.create_dataset(a.dset + '_omega', (0, 2),
                                           numpy.dtype(float),
                                           maxshape=(None, 2))
        n0 = 0
        n1 = ds.shape[0]
        ds.dims.create_scale(om, omkey)
        ds.dims[0].attach_scale(om)

    dn = n1 - n0
    ominmax = numpy.linspace(omin, omax, num=(dn + 1))

    om.resize(n1, 0)
    om[n0:n1, :] = numpy.array([ominmax[0:dn], ominmax[1:(dn+1)]]).T

    return

def image_info(a):
    """Return shape and dtype of first image

    * See process_img for options that transpose shape
"""
    img_0 = fabio.open(a.imagefiles[0])
    imgshape = img_0.data.shape
    if a.flip in ('t', 'T') + ('ccw90', 'r90') + ('cw90', 'r270'):
        imgshape = imgshape[::-1]

    return imgshape, img_0.data.dtype

def describe_imgs(a):
    print 'image files are: ', a.imagefiles
    im0 = fabio.open(a.imagefiles[0])
    print 'Total number of files: %d' % len(a.imagefiles)
    print 'First file: %s' % a.imagefiles[0]
    print '... fabio class: %s' % im0.__class__
    print '... number of frames: %d' % im0.nframes
    print '... image dimensions: %d X %d' % (im0.dim1, im0.dim2)
    print '... image data type: %s' % im0.data.dtype

    pass

def set_options():
    """Set options for command line"""
    parser = argparse.ArgumentParser(description="imageseries builder")

    parser.add_argument("-i", "--info", help="describe the input files and quit",
                        action="store_true")

    # file options
    parser.add_argument("-o", "--outfile", help="name of HDF5 output file",
                        default="imageseries.h5")
    parser.add_argument("-a", "--append",
                        help="append to the dataset instead of making a new one",
                        action="store_true")

    help_d = "path to HDF5 data set"
    parser.add_argument("-d", "--dset", help=help_d, default="/imageseries")

    parser.add_argument("imagefiles", nargs="+", help="image files")

    # image processing options
    parser.add_argument("--flip",
                        help="reorient the image according to specification",
                        metavar="FLIPARG", action="store", default=None)

    parser.add_argument("--empty", "--blank",
                        help="number of blank frames in beginning of file",
                        metavar="N", type=int, action="store", default=0)

    parser.add_argument("--dark-file", help="name of file containing dark image")

    parser.add_argument("--dark-from-empty",
                        help="use empty frames to build dark image",
                        action="store_true")
    # metadata
    parser.add_argument("--omega-min",
                        help="minimum omega for this series of images",
                        type=float, action="store")
    parser.add_argument("--omega-max",
                        help="minimum omega for this series of images",
                        type=float, action="store")
    parser.add_argument("--max-frames",
                        help="maximum number of frames in file (for testing)",
                        metavar="N", type=int, action="store", default=0)

    return parser

def execute(args, **kwargs):
    """Main execution

    * kwargs added to allow passing further options when not called from command line
    """
    p = set_options()
    a = p.parse_args(args)
    logging.info(str(a))

    if a.info:
        describe_imgs(a)
        return

    write_file(a, **kwargs)

    return

if __name__ == '__main__':
    #
    #  run
    #
    execute(sys.argv[1:])


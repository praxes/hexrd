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

class ImageFiles(object):
    """List of image files in sequence"""
    def __init__(self, a):
        """a is a namespace from parser"""
        self.files = a.imagefiles
        self.nfiles = len(self.files)
        self.nempty = a.empty
        self.maxframes = a.max_frames

        self._info()


        self.ntowrite = numpy.min((self.maxframes, self.nframes))\
          if self.maxframes > 0 else self.nframes

        self.outfile = a.outfile
        self.dgrppath = a.dset
        self.dsetpath = '/'.join((a.dset, 'images'))

    @staticmethod
    def _checkvalue(v, vtest, msg):
        """helper: ensure value set conistently"""
        if v is None:
            val = vtest
        else:
            if vtest != v:
                raise MakeImageSeriesError(msg)
            else:
                val = v

        return val

    def _info(self):
        """basic info: dtype, shape, nframes, and verify consistency"""
        cn = None
        shp = None
        dtp = None

        nf = 0
        for imgf in self.files:
            img = fabio.open(imgf)
            dat = img.data
            shp = self._checkvalue(shp, dat.shape, "inconistent image shapes")
            dtp = self._checkvalue(dtp, dat.dtype, "inconistent image dtypes")
            cn = self._checkvalue(cn, img.classname, "inconistent image types")
            if img.nframes >= self.nempty:
                nf += img.nframes - self.nempty
            else:
                raise MakeImageSeriesError("more empty frames than images")

        self.nframes = nf
        self.shape = shp
        self.dtype = dtp
        self.imagetype = cn

    def describe(self):
        print('Number of Files: %d' % self.nfiles)
        print('... image type: %s' % self.imagetype)
        print('... image dimensions: %d X %d' % self.shape)
        print('... image data type: %s' % self.dtype)
        print('... empty frames per file: %d' % self.nempty)
        print('... number of nonempty frames: %d' % self.nframes)
        maxf = self.maxframes if self.maxframes > 0 else 'unlimited'
        print('... max frames requested: %s' % maxf)
        print('... will write: %d' % self.ntowrite)

    def opendset(self):
        """open the HDF5 data set"""
        # note: compression implies chunked storage
        msg = 'writing to file/path: %s:%s' % (self.outfile, self.dgrppath)
        logging.info(msg)
        f = h5py.File(self.outfile, "a")
        try:
            shp = (self.ntowrite,) + self.shape
            ds = f.create_dataset(self.dsetpath, shp, self.dtype,
                                  compression="gzip")
        except Exception as e:
            errmsg = '%s: %s\n...  exception: ' % \
              (ERR_OVERWRITE, DSetPath(self.outfile, self.dsetpath))
            raise MakeImageSeriesError(errmsg + str(e))

        return f, ds

    def write(self):
        """write to HDF5 file"""
        f, ds = self.opendset()
        #
        # Now add the images
        #
        nframes = 0 # number completed
        for i in range(self.nfiles):
            if nframes >= self.ntowrite: break

            logging.debug('processing file %d of %d' % (i+1, self.nfiles))
            img_i = fabio.open(self.files[i])
            nfi = img_i.nframes
            for j in range(nfi):
                msg = '... file %d/image %d' % (i, j)
                logging.debug(msg)
                if j < self.nempty:
                    logging.debug('... empty frame ... skipping')
                else:
                    ds[nframes, :, :] = img_i.data
                    nframes += 1
                    logging.debug('... wrote image %s of %s' %\
                                  (nframes, self.ntowrite))
                    if nframes >= self.ntowrite:
                        logging.debug('wrote last frame: stopping')
                        break
                if j < nfi - 1:
                    # on last frame in file, fabio will look for next file
                    img_i = img_i.next()

        f.close()

def set_options():
    """Set options for command line"""
    parser = argparse.ArgumentParser(description="imageseries builder")

    parser.add_argument("-i", "--info", help="describe the input files and quit",
                        action="store_true")

    # file options
    parser.add_argument("-o", "--outfile", help="name of HDF5 output file",
                        default="imageseries.h5")
    help_d = "path to HDF5 data set"
    parser.add_argument("-d", "--dset", help=help_d, default="/imageseries")

    # image options
    parser.add_argument("imagefiles", nargs="+", help="image files")

    parser.add_argument("--empty", "--blank",
                        help="number of blank frames in beginning of file",
                        metavar="N", type=int, action="store", default=0)
    parser.add_argument("--max-frames",
                        help="maximum number of frames in file (for testing)",
                        metavar="N", type=int, action="store", default=0)

    return parser

def execute(args, **kwargs):
    """Main execution

    * kwargs added to allow passing further options when not called from
      command line
    """
    p = set_options()
    a = p.parse_args(args)
    logging.info(str(a))

    ifiles = ImageFiles(a)

    if a.info:
        ifiles.describe()
    else:
        ifiles.write()

    # write_file(a, **kwargs)

if __name__ == '__main__':
    #
    #  run
    #
    execute(sys.argv[1:])

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 18:43:52 2018

@author: s1iduser
"""
from skimage import io

# %%
img_stem = 'MapZBeam_5mmTo9mm_%s__%06d.tif'

im_or = ['m90deg', '0deg', '90deg']
im_idx = 31

img_list = [io.imread(img_stem % (i, im_idx)) for i in im_or]

# %%
for i, img in enumerate(img_list):
    fig, ax = plt.subplots()

    ax.imshow(img, cmap=plt.cm.inferno, vmin=np.percentile(img, 50))
    fig.suptitle("%s" % im_or[i])
    ax.axis('normal')

# %%
x_range = [0, 2047]  # eyeballed
y_range = [1990, 2040] # eyeballed

beam_img_list = [img[np.ix_(range(*y_range), range(*x_range))] for img in img_list]

sinogram = np.vstack([np.sum(bimg, axis=0) for bimg in beam_img_list])


# %%
pix_size = 0.00148

left = pix_size*np.r_[361, 2015]
right = pix_size*np.r_[1669, 2029]

diff = right - left
incl = np.degrees(np.arctan2(diff[1], diff[0]))
# Copyright (c) 2016 by Mike Jarvis and the other collaborators on GitHub at
# https://github.com/rmjarvis/Piff  All rights reserved.
#
# Piff is free software: Redistribution and use in source and binary forms
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.

from __future__ import print_function
import galsim
import numpy
import piff
import os
import subprocess
import yaml
import fitsio

from test_helper import get_script_name

def test_Gaussian():
    """This is about the simplest possible model I could think of.  It just uses the
    HSM adaptive moments routine to measure the moments, and then it models the
    PSF as a Gaussian.
    """

    # Here is the true PSF
    sigma = 1.3
    g1 = 0.23
    g2 = -0.17
    psf = galsim.Gaussian(sigma=sigma).shear(g1=g1, g2=g2)

    # Draw the PSF onto an image.  Let's go ahead and give it a non-trivial WCS.
    wcs = galsim.JacobianWCS(0.26, 0.05, -0.08, -0.29)
    image = galsim.Image(64,64, wcs=wcs)
    # This is only going to come out right if we (unphysically) don't convolve by the pixel.
    psf.drawImage(image, method='no_pixel')

    # Make a StarData instance for this image
    stardata = piff.StarData(image, image.trueCenter())
    star = piff.Star(stardata, None)

    # Fit the model from the image
    model = piff.Gaussian()
    fit = model.fit(star).fit

    print('True sigma = ',sigma,', model sigma = ',fit.params[0])
    print('True g1 = ',g1,', model g1 = ',fit.params[1])
    print('True g2 = ',g2,', model g2 = ',fit.params[2])

    # This test is pretty accurate, since we didn't add any noise and didn't convolve by
    # the pixel, so the image is very accurately a sheared Gaussian.
    true_params = [ sigma, g1, g2 ]
    numpy.testing.assert_almost_equal(fit.params[0], sigma, decimal=7)
    numpy.testing.assert_almost_equal(fit.params[1], g1, decimal=7)
    numpy.testing.assert_almost_equal(fit.params[2], g2, decimal=7)
    numpy.testing.assert_almost_equal(fit.params, true_params, decimal=7)

    # Now test running it via the config parser
    config = {
        'model' : {
            'type' : 'Gaussian'
        }
    }
    logger = piff.config.setup_logger()
    model = piff.process_model(config, logger)
    fit = model.fit(star).fit

    # Same tests.
    numpy.testing.assert_almost_equal(fit.params[0], sigma, decimal=7)
    numpy.testing.assert_almost_equal(fit.params[1], g1, decimal=7)
    numpy.testing.assert_almost_equal(fit.params[2], g2, decimal=7)
    numpy.testing.assert_almost_equal(fit.params, true_params, decimal=7)


def test_Mean():
    """For the interpolation, the simplest possible model is just a mean value, which barely
    even qualifies as doing any kind of interpolating.  But it tests the basic glue software.
    """
    import numpy
    # Make a list of paramter vectors to "interpolate"
    numpy.random.seed(123)
    nstars = 100
    vectors = [ numpy.random.random(10) for i in range(nstars) ]
    mean = numpy.mean(vectors, axis=0)
    print('mean = ',mean)

    # Make some dummy StarData objects to use.  The only thing we really need is the properties,
    # although for the Mean interpolator, even this is ignored.
    target_data = [
            piff.StarData.makeTarget(x=numpy.random.random()*2048, y=numpy.random.random()*2048)
            for i in range(nstars) ]
    fit = [ piff.StarFit(v) for v in vectors ]
    stars = [ piff.Star(d, f) for d,f in zip(target_data,fit) ]

    # Use the piff.Mean interpolator
    interp = piff.Mean()
    interp.solve(stars)

    print('True mean = ',mean)
    print('Interp mean = ',interp.mean)

    # This should be exactly equal, since we did the same calculation.  But use almost_equal
    # anyway, just in case we decide to do something slightly different, but equivalent.
    numpy.testing.assert_almost_equal(mean, interp.mean)

    # Now test running it via the config parser
    config = {
        'interp' : {
            'type' : 'Mean'
        }
    }
    logger = piff.config.setup_logger()
    interp = piff.process_interp(config, logger)
    interp.solve(stars)
    numpy.testing.assert_almost_equal(mean, interp.mean)


def test_single_image():
    """Test the simple case of one image and one catalog.
    """
    # Make the image
    image = galsim.Image(2048, 2048, scale=0.26)

    # Where to put the stars.  Include some flagged and not used locations.
    x_list = [ 123.12, 345.98, 567.25, 1094.94, 924.15, 1532.74, 1743.11, 888.39, 1033.29, 1409.31 ]
    y_list = [ 345.43, 567.45, 1094.32, 924.29, 1532.92, 1743.83, 888.83, 1033.19, 1409.20, 123.11 ]
    flag_list = [ 0, 0, 12, 0, 0, 1, 0, 0, 0, 0 ]
    use_list = [ 1, 1, 1, 1, 1, 0, 1, 1, 0, 1 ]

    # Draw a Gaussian PSF at each location on the image.
    sigma = 1.3
    g1 = 0.23
    g2 = -0.17
    psf = galsim.Gaussian(sigma=sigma).shear(g1=g1, g2=g2)
    for x,y,flag,use in zip(x_list, y_list, flag_list, use_list):
        bounds = galsim.BoundsI(int(x-31), int(x+32), int(y-31), int(y+32))
        offset = galsim.PositionD( x-int(x)-0.5 , y-int(y)-0.5 )
        psf.drawImage(image=image[bounds], method='no_pixel', offset=offset)
        # corrupt the ones that are marked as flagged
        if flag:
            print('corrupting star at ',x,y)
            ar = image[bounds].array
            im_max = numpy.max(ar) * 0.2
            ar[ar > im_max] = im_max

    # Write out the image to a file
    image_file = os.path.join('data','simple_image.fits')
    image.write(image_file)

    # Write out the catalog to a file
    dtype = [ ('x','f8'), ('y','f8'), ('flag','i2'), ('use','i2') ]
    data = numpy.empty(len(x_list), dtype=dtype)
    data['x'] = x_list
    data['y'] = y_list
    data['flag'] = flag_list
    data['use'] = use_list
    cat_file = os.path.join('data','simple_cat.fits')
    fitsio.write(cat_file, data, clobber=True)

    # Use InputFiles to read these back in
    input = piff.InputFiles(image_file, cat_file)
    assert input.image_files == [ image_file ]
    assert input.cat_files == [ cat_file ]
    assert input.x_col == 'x'
    assert input.y_col == 'y'

    # Check image
    input.readImages()
    assert len(input.images) == 1
    numpy.testing.assert_equal(input.images[0].array, image.array)

    # Check catalog
    input.readStarCatalogs()
    assert len(input.cats) == 1
    numpy.testing.assert_equal(input.cats[0]['x'], x_list)
    numpy.testing.assert_equal(input.cats[0]['y'], y_list)

    # Repeat, using flag and use columns this time.
    input = piff.InputFiles(image_file, cat_file, flag_col='flag', use_col='use', stamp_size=48)
    assert input.flag_col == 'flag'
    assert input.use_col == 'use'
    input.readImages()
    input.readStarCatalogs()
    assert len(input.cats[0]) == 7

    # Make star data
    orig_stars = input.makeStarData()
    assert len(orig_stars) == 7
    assert orig_stars[0].image.array.shape == (48,48)

    # Process the star data
    model = piff.Gaussian()
    interp = piff.Mean()
    fitted_stars = [ model.fit(piff.Star(star,None)) for star in orig_stars ]
    interp.solve(fitted_stars)
    print('mean = ',interp.mean)

    # Check that the interpolation is what it should be
    target_data = piff.StarData.makeTarget(x=1024, y=123) # Any position would work here.
    target = piff.Star(target_data, None)
    true_params = [ sigma, g1, g2 ]
    test_star = interp.interpolate(target)
    numpy.testing.assert_almost_equal(test_star.fit.params, true_params, decimal=5)

    # Now test running it via the config parser
    config = {
        'input' : {
            'images' : image_file,
            'cats' : cat_file,
            'flag_col' : 'flag',
            'use_col' : 'use',
            'stamp_size' : 48
        }
    }
    logger = piff.config.setup_logger()
    orig_stars = piff.process_input(config, logger)

    # Use the PSF builder to process the stars data this time.
    psf = piff.PSF.build(orig_stars, model, interp, logger)
    test_star = psf.interp.interpolate(target)
    numpy.testing.assert_almost_equal(test_star.fit.params, true_params, decimal=5)

    # Round trip to a file
    psf_file = os.path.join('output','simple_psf.fits')
    psf.write(psf_file, logger)
    psf = piff.PSF.read(psf_file, logger)
    assert type(psf.model) is piff.Gaussian
    assert type(psf.interp) is piff.Mean
    test_star = psf.interp.interpolate(target)
    numpy.testing.assert_almost_equal(test_star.fit.params, true_params, decimal=5)

    # Do the whole thing with the config parser
    os.remove(psf_file)
    config['model'] = { 'type' : 'Gaussian' }
    config['interp'] = { 'type' : 'Mean' }
    config['output'] = { 'file_name' : psf_file }

    piff.piffify(config)
    psf = piff.PSF.read(psf_file)
    test_star = psf.interp.interpolate(target)
    numpy.testing.assert_almost_equal(test_star.fit.params, true_params, decimal=5)

    # Test using the piffify executable
    os.remove(psf_file)
    with open('simple.yaml','w') as f:
        f.write(yaml.dump(config, default_flow_style=False))
    piffify_exe = get_script_name('piffify')
    p = subprocess.Popen( [piffify_exe, 'simple.yaml'] )
    p.communicate()
    psf = piff.PSF.read(psf_file)
    test_star = psf.interp.interpolate(target)
    numpy.testing.assert_almost_equal(test_star.fit.params, true_params, decimal=5)

    # Test that we can make rho statistics
    min_sep = 1
    max_sep = 100
    bin_size = 0.1
    stats = piff.RhoStatistics(psf, orig_stars, min_sep=min_sep, max_sep=max_sep, bin_size=bin_size)

    rhos = [stats.rho1, stats.rho2, stats.rho3, stats.rho4, stats.rho5]
    import numpy as np
    for rho in rhos:
        # Test the range of separations
        radius = np.exp(rho.logr)
        # last bin can be one bigger than max_sep
        np.testing.assert_array_less(radius, np.exp(np.log(max_sep) + bin_size))
        np.testing.assert_array_less(min_sep, radius)
        np.testing.assert_array_almost_equal(np.diff(rho.logr), bin_size, decimal=5)

        # Test that the max absolute value of each rho isn't crazy
        np.testing.assert_array_less(np.abs(rho.xip), 1)

        # Check that each rho isn't precisely zero. This means the sum of abs > 0
        np.testing.assert_array_less(0, np.sum(np.abs(rho.xip)))

    # Test the plotting and writing
    rho_psf_file = os.path.join('output','simple_psf_rhostats.png')
    stats.write(rho_psf_file)

    # Test that we can make summary shape statistics, using HSM
    shapeStats = piff.ShapeStatistics(psf, orig_stars)

    # test their characteristics
    np.testing.assert_array_almost_equal(sigma, shapeStats.T, decimal=4)
    np.testing.assert_array_almost_equal(sigma, shapeStats.T_model, decimal=3)
    np.testing.assert_array_almost_equal(g1, shapeStats.g1, decimal=4)
    np.testing.assert_array_almost_equal(g1, shapeStats.g1_model, decimal=3)
    np.testing.assert_array_almost_equal(g2, shapeStats.g2, decimal=4)
    np.testing.assert_array_almost_equal(g2, shapeStats.g2_model, decimal=3)

    shape_psf_file = os.path.join('output','simple_psf_shapestats.png')
    shapeStats.write(shape_psf_file)

    # Test that we can use the config parser for both RhoStatistics and ShapeStatistics
    config['stats'] = [{'type': 'ShapeStatistics',
                        'output': {'file_name': shape_psf_file}},
                       {'type': 'RhoStatistics',
                        'output': {'file_name': rho_psf_file}}]

    os.remove(psf_file)
    os.remove(rho_psf_file)
    os.remove(shape_psf_file)
    piff.piffify(config)

    # Test using the piffify executable
    os.remove(psf_file)
    os.remove(rho_psf_file)
    os.remove(shape_psf_file)
    with open('simple.yaml','w') as f:
        f.write(yaml.dump(config, default_flow_style=False))
    piffify_exe = get_script_name('piffify')
    p = subprocess.Popen( [piffify_exe, 'simple.yaml'] )
    p.communicate()

if __name__ == '__main__':
    test_Gaussian()
    test_Mean()
    test_single_image()


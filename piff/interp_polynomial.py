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
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

"""
.. module:: interp_polynomial
"""

from __future__ import print_function
from scipy.optimize import curve_fit

from .interp import Interp

class Polynomial(Interp):
    """Polynomial interpolation of the data as a function of whatever coordinates
    """
    def __init__(self, degree):
        self.degree = degree
        self.beta   = [] # polynomial coefficients, list (for each CCD) of lists 
        self.powers = None
        
    def fitData(self, data, pos):
        """Fit for the polynomial coefficients given some data.

        :param data:        A list of lists of data vectors (numpy arrays) for each star
        :param pos:         A list of lists of positions (or whatever additional parameters for interpolation) of the stars
        """
        import numpy as np
        import itertools
        
        nlists = len(data)
        if(len(pos)!=nlists):
          raise Exception('data list and position list are not of same dimension')
          
        self.beta = []
        
        ndim=len(pos[0][0])
        generators = [_basis_vector(ndim+1,i for i in range(ndim+1)] 
        # generator matrix
        
        self.powers = map(sum, itertools.combinations_with_replacement(generators, self.degree))
        # a list of combinations of powers
          
        for i in range(nlists): # go over one CCD at a time
          thisdata=np.asarray(data[i])
          nstars=thisdata.shape[0]
          
          thispos=np.asarray(pos[i])
          if(ndim!=thispos.shape[1]):
            raise Exception('coordinates have inconsistent dimensionality')
          
          thispos = np.hstack((np.ones((nstars, 1), dtype=thispos.dtype) , thispos)) 
          # prepend ones for constant
           
          A = np.hstack(np.asarray([_as_tall((thispos**p).prod(1)) for p in powers]))
          # a matrix of star coordinates to the combinations of powers
          
          self.beta.append(np.linalg.lstsq(A, thisdata)[0])
          # solve

    def interpolate(self, image_num, pos):
        """Perform the interpolation to find the interpolated data vector at some position in
        some image.

        :param image_num:   The index of the image in the original list of data vectors.
        :param pos:         The position to which to interpolate.

        :returns: the data vector (a numpy array) interpolated to the given position.
        """
        
        if(image_num>=len(beta)):
          raise Exception('you are asking me to interpolate in an image_num I have not fitted')

        return sum([b * ((np.asarray((1,) + pos))**p).prod() # pre-pend constant term
                             for p, b in zip(self.powers, beta[image_num])])
        
        
    @staticmethod
    def _basis_vector(n, i):
        """ Return an array like [0, 0, ..., 1, ..., 0, 0]
        >>> basis_vector(3, 1)
        array([0, 1, 0])
        >>> basis_vector(5, 4)
        array([0, 0, 0, 0, 1])
        """
        x = np.zeros(n, dtype=int)
        x[i] = 1
        return x
    
    @staticmethod    
    def _as_tall(x):
        """ Turns a row vector into a column vector """
        return x.reshape(x.shape + (1,))

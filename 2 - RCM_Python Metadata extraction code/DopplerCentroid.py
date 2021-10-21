#------------------------------------------------------------------------------
# Copyright (c) Her majesty the Queen in right of Canada as represented
# by the Minister of National Defence, 2015.
#------------------------------------------------------------------------------
# MODULE
#     DopplerCentroid.py
# INPUT
#     Doppler Centroid Estimates
# OUTPUT
#     Doppler Centroid
#------------------------------------------------------------------------------

#*********************************************************************************
# Author             Version        Date               Modifications
#---------------------------------------------------------------------------------
# Adam Klein         1.0            Feb 2015           initial development 
#*********************************************************************************

import calendar
import numpy
from scipy.interpolate import interp1d


class DopplerCentroidCalculator:
    def __init__(self, dopplerCentroidEstimates):
        self.dopplerCentroidEstimates = dopplerCentroidEstimates
        self.dopplerCentroidEstimates.sort(key=lambda dc: dc.azimuthTime)

        self.__useInterpolation = len(self.dopplerCentroidEstimates) > 1
        if self.__useInterpolation:
            self.__createInterpolationFunctions()

    def __createInterpolationFunctions(self):
        azimuthTimes = [self.__datetimeTotalSeconds(dce.azimuthTime) for dce in self.dopplerCentroidEstimates]
        slantRangeReferenceTimes = [dce.slantRangeReferenceTime for dce in self.dopplerCentroidEstimates]

        numPolynomialCoefficients = len(self.dopplerCentroidEstimates[0].polynomialCoefficients)
        polynomialCoefficients = []
        for i in range(numPolynomialCoefficients):
            polynomialCoefficients.append([dce.polynomialCoefficients[i] for dce in self.dopplerCentroidEstimates])

        self.__coefficientInterpFuncs = [interp1d(azimuthTimes, pc) for pc in polynomialCoefficients]
        self.__refTimeInterpFunc = interp1d(azimuthTimes, slantRangeReferenceTimes)

    def __datetimeTotalSeconds(self, dt):
        return calendar.timegm(dt.utctimetuple())+dt.microsecond/1e6

    def getDopplerCentroid(self, azimuthTime, slantRangeTime):
        if self.__useInterpolation:
            dcEstimate = self.__getDopplerCentroidEstimate(azimuthTime)
        else:
            dcEstimate = self.dopplerCentroidEstimates[0]

        return self.__evalDopplerCentroid(slantRangeTime, dcEstimate.slantRangeReferenceTime, dcEstimate.polynomialCoefficients)

    def __evalDopplerCentroid(self, slantRangeTime, dopplerCentroidRefTime, dopplerCentroidCoeffs):
        '''evaluate the doppler centroid given the parameters
        slantRangeTime is the 2 way slant range time to the point in seconds
        dopplerCentroidRefTime is the doppler centroid reference time in seconds
        dopplerCentroidCoeffs is a list of the doppler centroid estimation coefficients
        returns the dopplerCentroid in Hz.'''

        #doppler centroid = d0 + d1(tSR - t0) +d2(tSR - t0)^2 +d3(tSR - t0)^3 + d4(tSR - t0)^4
        #where d0~d5 are the doppler centroid coefficients
        #tSR is the two way slant range time to the point in question
        #t0 is the doppler centroid reference time

        slantRangeTimeDelta = slantRangeTime - dopplerCentroidRefTime   #(tSR - t0)

        return numpy.sum( [c * slantRangeTimeDelta**i for i, c in enumerate(dopplerCentroidCoeffs)] )

    def __getDopplerCentroidEstimate(self, azimuthTime):
        try:
            return self.__interpolateDopplerCentroidEstimate(azimuthTime)
        except ValueError:
            if azimuthTime < self.dopplerCentroidEstimates[0].azimuthTime:
                return self.dopplerCentroidEstimates[0]
            elif azimuthTime > self.dopplerCentroidEstimates[-1].azimuthTime:
                return self.dopplerCentroidEstimates[-1]
            else:
                raise

    def __interpolateDopplerCentroidEstimate(self, azimuthTime):
        azimuthTimeSec = self.__datetimeTotalSeconds(azimuthTime)
        refTime = self.__refTimeInterpFunc(azimuthTimeSec)[()]
        coefficients = tuple(f(azimuthTimeSec)[()] for f in self.__coefficientInterpFuncs)

        return DopplerCentroidEstimate(azimuthTime, refTime, coefficients)


class DopplerCentroidEstimate():
    def __init__(self, azimuthTime, slantRangeReferenceTime, polynomialCoefficients):
        self.azimuthTime = azimuthTime
        self.slantRangeReferenceTime = slantRangeReferenceTime
        self.polynomialCoefficients = polynomialCoefficients
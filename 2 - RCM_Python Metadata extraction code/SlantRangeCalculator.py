#------------------------------------------------------------------------------
# Copyright (c) Her majesty the Queen in right of Canada as represented
# by the Minister of National Defence, 2015.
#------------------------------------------------------------------------------
# MODULE
#     SlantRangeCalculator.py
# INPUT
#     slantRangeToGroundRange parameters
# OUTPUT
#     Slant Range
#------------------------------------------------------------------------------

#*********************************************************************************
# Author             Version        Date               Modifications
#---------------------------------------------------------------------------------
# Adam Klein         1.0            June 2015          Refactor out of RADARSAT2Layer
#*********************************************************************************

import calendar
import scipy.interpolate


class SlantRangeCalculator:
    '''
    Calculates slant range given a list of GroundToSlantRangeEntry's.
    Performs interpolation of the ground to slant range entries.
    '''
    def __init__(self, groundToSlantRangeEntries):
        self.groundToSlantRangeEntries = groundToSlantRangeEntries
        self.groundToSlantRangeEntries.sort(key = lambda entry: entry.timestamp)

    def getSlantRange(self, groundRange, timestamp):
        entry = self.__interpolateGroundToSlantRangeEntry(timestamp)
        return self.__calculateSlantRange(groundRange, entry)

    def __interpolateGroundToSlantRangeEntry(self, timestamp):
        numEntries = len(self.groundToSlantRangeEntries)
        if numEntries == 1:
            return self.groundToSlantRangeEntries[0]
        elif timestamp < self.groundToSlantRangeEntries[0].timestamp:
            return self.groundToSlantRangeEntries[0]
        elif timestamp > self.groundToSlantRangeEntries[numEntries-1].timestamp:
            return self.groundToSlantRangeEntries[numEntries-1]

        interpolation = self.__getGroundToSlantRangeInterpolationFunctions()

        unixTimestamp = calendar.timegm(timestamp.utctimetuple())+timestamp.microsecond/1e6
        firstSampleSlantRange = interpolation.firstSampleSlantRange(unixTimestamp)
        groundRangeOrigin = interpolation.groundRangeOrigin(unixTimestamp)
        groundToSlantRangeCoefficients = tuple(c(unixTimestamp) for c in interpolation.coefficients)

        entry = GroundToSlantRangeEntry(timestamp, firstSampleSlantRange, groundRangeOrigin, groundToSlantRangeCoefficients)
        return entry

    def __getGroundToSlantRangeInterpolationFunctions(self):
        try:
            return self.__groundToSlantRangeInterpolationFunctions
        except:
            interpolation = self.__calculateGroundToSlantRangeInterpolationFunctions(self.groundToSlantRangeEntries)
            self.__groundToSlantRangeInterpolationFunctions = interpolation
            return interpolation

    def __calculateGroundToSlantRangeInterpolationFunctions(self, groundToSlantRangeEntries):
        times = [entry.timestamp for entry in groundToSlantRangeEntries]
        firstSampleSlantRanges = [entry.firstSampleSlantRange for entry in groundToSlantRangeEntries]
        groundRangeOrigins = [entry.groundRangeOrigin for entry in groundToSlantRangeEntries]
        coefficients = [entry.groundToSlantRangeCoefficients for entry in groundToSlantRangeEntries]

        #translate times into unix timestamps because the scipy interpolation function
        #cannot operate on datetimes
        times = [calendar.timegm(t.utctimetuple())+t.microsecond/1e6 for t in times]

        interpolation = SlantRangeCalculator.GroundToSlantRangeInterpolation()

        interpolation.firstSampleSlantRange = scipy.interpolate.interp1d(times, firstSampleSlantRanges)
        interpolation.groundRangeOrigin = scipy.interpolate.interp1d(times, groundRangeOrigins)
        interpolation.coefficients = tuple(scipy.interpolate.interp1d(times, c) for c in zip(*coefficients))

        return interpolation

    def __calculateSlantRange(self, groundRange, groundToSlantRangeEntry):
        groundRangeDelta = groundRange - groundToSlantRangeEntry.groundRangeOrigin
        slantRange = 0
        for i, coeff in enumerate(groundToSlantRangeEntry.groundToSlantRangeCoefficients):
            slantRange += coeff * groundRangeDelta**i

        return slantRange

    class GroundToSlantRangeInterpolation():
        pass


class GroundToSlantRangeEntry:

    def __init__(self, timestamp, firstSampleSlantRange, groundRangeOrigin, groundToSlantRangeCoefficients):

        self.timestamp = timestamp
        self.firstSampleSlantRange = firstSampleSlantRange
        self.groundRangeOrigin = groundRangeOrigin
        self.groundToSlantRangeCoefficients = groundToSlantRangeCoefficients


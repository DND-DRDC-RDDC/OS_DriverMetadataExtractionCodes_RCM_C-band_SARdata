#------------------------------------------------------------------------------
# Copyright (c) Her majesty the Queen in right of Canada as represented
# by the Minister of National Defence, 2015.
#------------------------------------------------------------------------------
# MODULE
#     OrbitCalculator.py
# INPUT
#     Orbit State Vector
# OUTPUT
#     Orbit Position and Velocity
#------------------------------------------------------------------------------

#*********************************************************************************
# Author             Version        Date               Modifications
#---------------------------------------------------------------------------------
# Adam Klein         1.0            June 2015          Refactor out of RADARSAT2Layer
#*********************************************************************************

import calendar


class OrbitCalculator:
    '''calculates the position and velocity of the platform given a list of state vectors.
    Linearly interpolates the state vectors to the given datetime.'''
    def __init__(self, orbitStateVectors):
        self.orbitStateVectors = orbitStateVectors
        self.orbitStateVectors.sort(key = lambda stateVector: stateVector.timestamp)

    def getPosition(self, timestamp):
        return self.__getOrbitStateVector('position', timestamp)

    def getVelocity(self, timestamp):
        return self.__getOrbitStateVector('velocity', timestamp)

    def __getOrbitStateVector(self, positionOrVelocity, timestamp):

        defaultStateVector = self.__checkStateVectorBounds(timestamp)
        if defaultStateVector is not None:
            return getattr(defaultStateVector,positionOrVelocity)

        index = {'position': 0, 'velocity': 1}[positionOrVelocity]
        interpolationFunctions = self.__getOrbitStateVectorFunctions()[index]

        return self.__evaluateStateVectorInterpolationFunctions(interpolationFunctions, timestamp)

    def __checkStateVectorBounds(self, timestamp):
        numStateVectors = len(self.orbitStateVectors)
        if numStateVectors == 1:
            return self.orbitStateVectors[0]
        elif timestamp < self.orbitStateVectors[0].timestamp:
            return self.orbitStateVectors[0]
        elif timestamp > self.orbitStateVectors[numStateVectors-1].timestamp:
            return self.orbitStateVectors[numStateVectors-1]

    def __getOrbitStateVectorFunctions(self):
        try:
            return self.__orbitStateVectorFunctions
        except AttributeError:
            self.__orbitStateVectorFunctions = self.__calculateStateVectorInterpolationFunctions()
            return self.__orbitStateVectorFunctions

    def __calculateStateVectorInterpolationFunctions(self):
        import scipy.interpolate

        times = [stateVector.timestamp for stateVector in self.orbitStateVectors]
        positions = [stateVector.position for stateVector in self.orbitStateVectors]
        velocities = [stateVector.velocity for stateVector in self.orbitStateVectors]
        x, y, z = list(zip(*positions))
        vx, vy, vz = list(zip(*velocities))

        #translate times into unix timestamps because the scipy interpolation function
        #cannot operate on datetimes
        times = [calendar.timegm(t.utctimetuple())+t.microsecond/1e6 for t in times]

        positionFunctions = tuple(scipy.interpolate.interp1d(times, var) for var in [x, y, z])

        velocityFunctions = tuple(scipy.interpolate.interp1d(times, var) for var in [vx, vy, vz])

        return positionFunctions, velocityFunctions

    def __evaluateStateVectorInterpolationFunctions(self, interpolationFunctions, timestamp):
        unixTimestamp = calendar.timegm(timestamp.utctimetuple())+timestamp.microsecond/1e6
        return tuple(f(unixTimestamp) for f in interpolationFunctions)


class OrbitStateVector:
    def __init__(self, timestamp, position, velocity):
        '''timestamp: a datetime
        position: a tuple (x,y,z)
        velocity: a tuple (x,y,z)'''
        self.timestamp = timestamp
        self.position = position
        self.velocity = velocity

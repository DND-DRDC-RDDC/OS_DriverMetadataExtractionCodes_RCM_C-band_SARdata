#------------------------------------------------------------------------------
# Copyright (c) Her majesty the Queen in right of Canada as represented
# by the Minister of National Defence, 2018.
#------------------------------------------------------------------------------

# ***********************************************************************************************
# Author                  Version        Date               Modifications
# -----------------------------------------------------------------------------------------------
# Shawn Gong, MDA          1.0          Dec 2018            initial development
# Shawn Gong, MDA          1.1          Feb 2019            SigmaNoughtLUT(offset, gains), SigmaNoughtNoise, are per band
# Shawn Gong, MDA          1.2          summer 2020         extract azimuth noise
# Shawn Gong, MDA          1.3          Dec 2020            remove CIS's dependency
# Shawn Gong, MDA          1.4          April 2021          add ScanSAR burstmap
# Shawn Gong, MDA          1.41         May 2021            noise azimuth_scaling
# Shawn Gong, MDA          1.5          May 2022            add MLC support
# Shawn Gong, MDA          1.51         May 2022            fix a bug in __applyLUT()
#
# http://www.asc-csa.gc.ca/eng/satellites/radarsat/radarsat-tableau.asp
# ***********************************************************************************************

import os
import math
import re
import datetime
import numpy
import scipy.constants
from xml.etree import ElementTree as ET
from SARLayer import SARLayer, UnknownProductTypeException, NotInZeroDopplerCoordsException
from RasterLayer import register
import SAR.DopplerCentroid as DC
from SAR.SlantRangeCalculator import SlantRangeCalculator, GroundToSlantRangeEntry
from SAR.OrbitCalculator import OrbitCalculator, OrbitStateVector
import gvutils
import ATD_Utils
from bisect import bisect
from operator import add

SCHEMA = "{rcmGsProductSchema}"
RCM_PREFIX = 'RCM: '
GDAL_DRIVER_SHORT_NAME = 'RCM'

# regular expression for parsing the dates in product.xml.  The decimal part of the seconds field is optional
# after re.match, the output tuple will be year, month, day, hour, minute, integerSecond, decimalSeconds
RCM_DATE_FORMAT = '(\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+)(\.\d+)?Z'

# RCM orbit constants
# RCM inclination angle is 97.74, = the angle between the equatorial plane and the
# satellite orbit in the counter-clockwise direction.
# Degrees from north in the clockwise direction = 360-7.74
ORBIT_INCLINATION  = 97.74
ORBIT_RATE = 2*math.pi*179/(12*24*60*60)    # rad/s   repeat cycle of 179 orbits in 12 days
ORBIT_DURATION = 96.4       # mins

# following resolutions are from RCM Product Description RCM-SP-52-9092, issue 1/11, 08-JAN-2018
# BEAM_PARAMS tuple is (resolution, non-slc number of looks, non-slc number of looks for Dual HH-VV)
BEAM_PARAMS = {
               'Ship Detection':            ('variable', '5 x 1'),
               'Low Noise':                 ('100 x 100', '4 x 2'),
               'Low Resolution 100m':       ('100 x 100', '8 x 1'),
               'Medium Resolution 50m':     ('50 x 50', '4 x 1'),
               'Medium Resolution 50m High PRF':   ('50 x 50', '4 x 1'),
               'Medium Resolution 30m':     ('30 x 30', '2 x 2', '2 x 1'),
               'Medium Resolution 16m':     ('16 x 16', '1 x 4', '1 x 2'),
               'High Resolution 5m':        ('5 x 5', '1 x 1'),
               'Very High Resolution 3m':   ('3 x 3', '1 x 1'),
               'Quad-Polarization':         ('9 x 9', '1 x 1'),
               'Spotlight':                 ('3 x 1', '1 x 1')
               }

# various product types
PRODUCT_TYPES = {#tuple (inZeroDopplerCoords, inSlantRange, isDetected)
                 'SLC': (True,  True,  False), # Single-Look Complex
                 'GRD': (True,  False, True),  # Ground range georeferenced Detected
                 'GRC': (True,  False, False), # Ground range georeferenced Complex
                 'GCD': (False, False, True),  # GeoCoded Detected
                 'GCC': (False, False, False), # GeoCoded Complex
                 'MLC': (True,  True,  False)  # Multi-Look Complex (Mixed: UInt16 for detected image, CInt16 for complex I&Q band; or Float32 for detected image, Float32 for complex I&Q band)
                 }

IN_ZERO_DOPPLER_COORDS = 0      # first column
IN_SLANT_RANGE = 1              # second column
IS_DETECTED = 2                 # third column


class RCMLayer(SARLayer):
    # flag to specify that a type search should not be performed when instantiating new RCMLayer
    TYPE_SEARCH = False

    def __init__(self, *args, **kwargs):
        SARLayer.__init__(self, *args, **kwargs)
        
        self.set_name(RCM_PREFIX + self.get_filename(vrtOriginalFilename=True))
        self.productPath, filename = os.path.split(self.get_filename())

        #try:    #don't fail just because metadata has problems
        self._getXmlLines()
        self._readMetaData()
        self._setDisplayInfo()
        
        self.readBurstAttributes()
        self.getSigmaNoughtNoiseParams()        # this is needed to populate self.sigmaNoughtNoiseParams
        if not self.geocoded:
            self._calcIncAngles()
        
        # single-file NITF needs to open the .ntf
#        if filename == 'product.xml' and self.imageFiles[0].lower().endswith('.ntf'):
#            gvutils.warning('RCM single-file NITF image: please close this product.xml and open '+self.imageFiles[0])
        
        #except:
        #    print 'Warning: Failed to read some metadata for %s' % self.get_filename()
    
    def get_filename(self, vrtOriginalFilename=False, removeCalibrationDirective=False):
        f = super(RCMLayer, self).get_filename(vrtOriginalFilename)
        
        if removeCalibrationDirective and f.startswith('RCM_CALIB:'):
            calib, f = f.split(':', 2)[-2:]
            if calib not in ['UNCALIB', 'SIGMA0', 'BETA0', 'GAMMA']:
                raise ValueError('Unknown calibration directive')
        return f
    
    def _getXmlLines(self):
        self.filename = self.get_filename(vrtOriginalFilename=True, removeCalibrationDirective=True)
        with open(self.filename) as fd:
            self.xmlLines = fd.read()
        return self.xmlLines

    def _parseProductXml(self):
        ''' This function is no longer used. It is kept just in case. Shawn Gong Dec 2020
        '''
        #need to remove the xml namespace information because etree will add the xmlns to every element
        #if the file is parsed directly, it would be necessary to search using the fully qualified name, ie:
        #self.__xmlTree.find('{http://www.rsi.ca/rs2/prod/xml/schemas}imageAttributes')
        #hopefully there is an etree workaround in a later version.
        #even with this hack, etree is still way faster than xml.dom.minidom
        self.xmlLines = self.xmlLines.replace('xmlns="rcmGsProductSchema"', "", 1)        # this is necessary
    
    def _getProductFile(self, filePath):
        if ("product.xml" not in filePath):
            return os.path.join(filePath, "product.xml")
        else:
            return filePath

    def _imageAcqusition(self):
        imageAcqusition = self.datetimestr.split('.')[0]
        imageAcqusition = imageAcqusition.replace('T', ' ')
        return imageAcqusition

    def _readMetaData(self):
        '''Reads metadata from product.xml and sets the corresponding variables'''
        
        if '_GCC' in self.filename or '_GCD' in self.filename:
            gvutils.warning("RCM geocoded product is not supported in IA Pro.\n It is for display purpose only!")
            self.geocoded = True
        else:
            self.geocoded = False
        
        root = ET.fromstring(self.xmlLines)
        
        if int([i.text for i in next(root.iter(SCHEMA + "sceneAttributes")).findall(".//" + SCHEMA + "numberOfEntries")][0]) > 1:
            gvutils.error("RCM ScanSAR SLC image chips need to be stitched together!")
            return
        
        self.productID = [i.text for i in root.iter(SCHEMA + "productId")][0]           # the same image will be given different productId for different user
        self.imageID = [i.text for i in root.iter(SCHEMA + "downlinkSegmentId")][0]     # images in the same swath will have the same downlinkSegmentId
        self.datetimestr = [i.text for i in root.iter(SCHEMA + "rawDataStartTime")][0]
        self.year = self.datetimestr[:4]
        self.month = self.datetimestr[5:7]
        self.datetime = self.__datetime_from_string(self.datetimestr)
        self.sensorString = [i.text for i in root.iter(SCHEMA + "satellite")][0]
        self.beamMode = [i.text for i in root.iter(SCHEMA + "beamMode")][0]                     # i.e. 'Ship Detection'
        self.beamModeMnemonic = [i.text for i in root.iter(SCHEMA + "beamModeMnemonic")][0]     # i.e. 'SCSDA'
        self.freqBandStr = "C-band"
        self.acquisitionType = [i.text for i in root.iter(SCHEMA + "acquisitionType")][0]
        self.dataType = [i.text for i in root.iter(SCHEMA + "sampleType")][0]
        if self.isSpotlight():
            self.sensorType = 'Spotlight'
        elif self.isScanSAR():
            self.sensorType = 'ScanSAR'
        else:
            self.sensorType = 'Stripmap'
        self.imageFiles = [i.text for i in next(root.iter(SCHEMA + "sceneAttributes")).findall(".//" + SCHEMA + "ipdf")]
        self.xsize = int([i.text for i in next(root.iter(SCHEMA + "sceneAttributes")).findall(".//" + SCHEMA + "samplesPerLine")][0])
        self.ysize = int([i.text for i in next(root.iter(SCHEMA + "sceneAttributes")).findall(".//" + SCHEMA + "numLines")][0])
        self.incAngleNear = float([i.text for i in next(root.iter(SCHEMA + "sceneAttributes")).findall(".//" + SCHEMA + "incAngNearRng")][0])
        self.incAngleFar = float([i.text for i in next(root.iter(SCHEMA + "sceneAttributes")).findall(".//" + SCHEMA + "incAngFarRng")][0])
        self.nSamples = self.xsize
        self.nLines = self.ysize
        self.passDirection = [i.text for i in root.iter(SCHEMA + "passDirection")][0]
        self.polarizations = [i.text for i in root.iter(SCHEMA + "polarizations")][0]
        self.pols = self.polarizations.split()
        self.productType = [i.text for i in root.iter(SCHEMA + "productType")][0]
        self.isComplex = (self.productType == 'SLC')
        self.lutApplied = [i.text for i in root.iter(SCHEMA + "lutApplied")][0]
        self.processingElevation = float([i.text for i in root.iter(SCHEMA + "geodeticTerrainHeight")][0])
        self.lookDirection = [i.text for i in root.iter(SCHEMA + "antennaPointing")][0]
        self.orbitDataSource = [i.text for i in root.iter(SCHEMA + "orbitDataSource")][0]
        try:
            self.orbitID = [i.text for i in root.iter(SCHEMA + "orbitDataFileName")][0].split('.DEF_ORBIT')[0].split('_')[-1]
        except:
            self.orbitID = ''
        self.rangeLooks = int([i.text for i in root.iter(SCHEMA + "numberOfRangeLooks")][0])
        self.azimuthLooks = int([i.text for i in root.iter(SCHEMA + "numberOfAzimuthLooks")][0])
        self.bitsPerSample = int([i.text for i in root.iter(SCHEMA + "bitsPerSample")][0])
        windowCoefficient = float([i.text for i in root.iter(SCHEMA + "windowCoefficient")][0])
        self.PreserveDataSpectrum = (windowCoefficient == 0)
        self.radarCenterFrequency = float([i.text for i in root.iter(SCHEMA + "radarCenterFrequency")][0])
        self.wavelength = scipy.constants.c / self.radarCenterFrequency
        self.ellipsoidName = [i.text for i in root.iter(SCHEMA + "ellipsoidName")][0]
        self.semiMajorAxis = float([i.text for i in root.iter(SCHEMA + "semiMajorAxis")][0])
        self.semiMinorAxis = float([i.text for i in root.iter(SCHEMA + "semiMinorAxis")][0])
        self.orbitAltitude = float([i.text for i in root.iter(SCHEMA + "satelliteHeight")][0])
        try:
            self.noiseSubtraction = [i.text for i in root.iter(SCHEMA + "noiseSubtractionPerformed")][0]
        except:
            self.noiseSubtraction = ''
        self.pixelSpacing = float([i.text for i in root.iter(SCHEMA + "sampledPixelSpacing")][0])
        self.lineSpacing = float([i.text for i in root.iter(SCHEMA + "sampledLineSpacing")][0])
        self.lineTimeOrdering = [i.text for i in root.iter(SCHEMA + "lineTimeOrdering")][0]
        self.pixelTimeOrdering = [i.text for i in root.iter(SCHEMA + "pixelTimeOrdering")][0]
        self.beams = [i.text for i in root.iter(SCHEMA + "beams")][0].split()
        self.numberOfBeams = len(self.beams)
        
        if not self.geocoded:
            #save the zero doppler first and last line times as a datetime
            zeroDopplerTimeFirstLineStr = [i.text for i in root.iter(SCHEMA + "zeroDopplerTimeFirstLine")][0]
            self.zeroDopplerTimeFirstLineStr = zeroDopplerTimeFirstLineStr
            self.zeroDopplerTimeFirstLine = self.__datetime_from_string(zeroDopplerTimeFirstLineStr)
            zeroDopplerTimeLastLineStr = [i.text for i in root.iter(SCHEMA + "zeroDopplerTimeLastLine")][0]
            self.zeroDopplerTimeLastLineStr = zeroDopplerTimeLastLineStr
            self.zeroDopplerTimeLastLine = self.__datetime_from_string(zeroDopplerTimeLastLineStr)
            self.starttime = self.zeroDopplerTimeFirstLineStr
            self.endtime = self.zeroDopplerTimeLastLineStr
            if self.passDirection.lower() == "ascending":   # lineTimeOrdering is Decreasing
                self.starttime, self.endtime = self.endtime, self.starttime

            self.__createDopplerCentroidCalculator(root)
            self.__createSlantRangeCalculator(root)
            self.__createOrbitCalculator(root)

        self.tiePoints_line = list(map(float, [i.text for i in next(root.iter(SCHEMA + "geolocationGrid")).findall(".//" + SCHEMA + "line")]))
        self.tiePoints_pixel = list(map(float, [i.text for i in next(root.iter(SCHEMA + "geolocationGrid")).findall(".//" + SCHEMA + "pixel")]))
        self.tiePoints_lat = list(map(float, [i.text for i in next(root.iter(SCHEMA + "geolocationGrid")).findall(".//" + SCHEMA + "latitude")]))
        self.tiePoints_lon = list(map(float, [i.text for i in next(root.iter(SCHEMA + "geolocationGrid")).findall(".//" + SCHEMA + "longitude")]))
        self.tiePoints_height = list(map(float, [i.text for i in next(root.iter(SCHEMA + "geolocationGrid")).findall(".//" + SCHEMA + "height")]))
        
        #check for two different polarizations on transmit.
        #If both H and V are transmitted, they are alternated, not transmitted at the same time, so the effective PRF per channel is 1/2 the system PRF
        self.tx_dualPol = any([p.startswith('H') for p in self.pols]) and any([p.startswith('V') for p in self.pols])
        self.compactPol = (self.pols == ['CH', 'CV'] or self.pols == ['CV', 'CH'])
        
        self.PRF = {}
        prfInformation_nodes = root.findall(".//" + SCHEMA + "prfInformation")
        if len(prfInformation_nodes) == 1:   # only one PRF
            self.PRF[self.beams[0]] = float([i.text for i in prfInformation_nodes[0].iter(SCHEMA + "pulseRepetitionFrequency")][0])
        else:       # more than one PRF
            for node in prfInformation_nodes:
                if node.attrib['beam'] not in self.PRF:   # take the first PRF value per beam
                    self.PRF[node.attrib['beam']] = float([i.text for i in node.iter(SCHEMA + "pulseRepetitionFrequency")][0])
        
        # this is copied from RADARSAT2Layer.py, without this RCM NITF getBandLabel() won't work
        #hack: since the vrt contains very little metadata, we always have to go back to the original dataset (not sure what the point of a VRT is anyways)
        #now we need to know which label goes with each band.  The GDAL driver takes the info from the product.xml tags fullResolutionImageData
        #the band numbers should be in the same order as the tags - extract the pol and add to metadata
        for idx, pol in enumerate(self.pols):
            if self.isVRT():    # RCM NITF
                try:            # check if the data exists
                    self._ds.GetRasterBand(idx+1).GetMetadata()['POLARIMETRIC_INTERP']
                except:         # if not, write it
                    self._ds.GetRasterBand(idx+1).SetMetadataItem('POLARIMETRIC_INTERP', pol)

        #add the pol as the band description so OpenEV will use it in the raster properties dialog
        #TODO: it might be better to edit the GDAL driver to set the band description
        try:
            for idx in range(self._ds.RasterCount):
                self._ds.GetRasterBand(idx+1).SetDescription(self._ds.GetRasterBand(idx+1).GetMetadata()['POLARIMETRIC_INTERP'])
        except:
            pass
    
    def _imageStartTime(self):
        return self.zeroDopplerTimeFirstLineStr.split('T')

    def _imageEndTime(self):
        return self.zeroDopplerTimeLastLineStr.split('T')

    def _imageRawStartTime(self):
        return self.zeroDopplerTimeFirstLineStr

    def _imageRawEndTime(self):
        return self.zeroDopplerTimeLastLineStr

    def _setDisplayInfo(self):
        '''Updates the values in the self.displayInfo dictionary.'''
        
        self.displayInfo['Sensor'] = self.sensorString
        self.displayInfo['Acquisition Date'] = self.datetime.strftime("%d-%b-%Y %H:%M:%S UTC")
        self.displayInfo['Frequency'] = self.freqBandStr
        self.displayInfo['Polarization'] = self.polarizations
        self.displayInfo['Beam Mode'] = self.beamMode + ' (' + self.beamModeMnemonic + ')'    # i.e. Ship Detection (SCSDA)
        self.displayInfo['ASC_DES'] = self.passDirection.capitalize()
        self.displayInfo['Incidence Angle'] = str(round(float(self.incAngleNear), 1)) + ' (Near), ' + \
                                              str(round(float(self.incAngleFar ), 1)) + ' (Far)'
        
        self.displayInfo['Resolution'] = BEAM_PARAMS[self.beamMode][0]
        if self.isComplex:
            self.displayInfo['Number of Looks'] = '1 x 1'
        else:
            self.displayInfo['Number of Looks'] = BEAM_PARAMS[self.beamMode][1]
            if self.polarizations == "HH VV" and self.beamMode in ['Medium Resolution 30m', 'Medium Resolution 16m']:
                self.displayInfo['Number of Looks'] = BEAM_PARAMS[self.beamMode][2]
        
        if not self.geocoded:
            self.displayInfo['SatelliteHeading'] = self.getTrackAngle()

            self.displayInfo['Antenna_Pointing'] = self.lookDirection
            if self.lookDirection.lower() == "left":
                self.displayInfo['SensorLookDirection'] = math.fmod(self.displayInfo['SatelliteHeading'] - 90.0, 360.)
            else:
                self.displayInfo['SensorLookDirection'] = math.fmod(self.displayInfo['SatelliteHeading'] + 90.0, 360.)
            
            self.displayInfo['ProcessingElevation'] = self.processingElevation
            self.displayInfo['lutApplied'] = self.lutApplied
        
    def __createDopplerCentroidCalculator(self, root):
        self.dopplerCentroidEstimates = self.__extractDopplerCentroidEstimates(root)
        self.__dopplerCentroidCalculator = DC.DopplerCentroidCalculator(self.dopplerCentroidEstimates)

    def __extractDopplerCentroidEstimates(self, root):
        dopplerCentroidEstimates = []
        for node in root.findall(".//" + SCHEMA + "dopplerCentroidEstimate"):
            timeOfEstimateStr = [i.text for i in node.iter(SCHEMA + "timeOfDopplerCentroidEstimate")][0]
            timeOfEstimate = self.__datetime_from_string(timeOfEstimateStr)
            refTime = float([i.text for i in node.iter(SCHEMA + "dopplerCentroidReferenceTime")][0])
            coefficients = [float(c) for c in [i.text for i in node.iter(SCHEMA + "dopplerCentroidCoefficients")][0].split()]

            dopplerCentroidEstimates.append(DC.DopplerCentroidEstimate(timeOfEstimate, refTime, coefficients))
        return dopplerCentroidEstimates

    def __createOrbitCalculator(self, root):
        orbitStateVectors = self.__extractOrbitStateVectors(root)
        self.__orbitCalculator = OrbitCalculator(orbitStateVectors)

        orbitXVel, orbitYVel, orbitZVel = self.__orbitCalculator.getVelocity(self.getTime(self.nSamples/2, self.nLines/2))
        self.platformVelocity = math.sqrt(orbitXVel**2+orbitYVel**2+orbitZVel**2)

    def __extractOrbitStateVectors(self, root):
        return [self.__getOrbitStateVector(node) for node in root.findall(".//" + SCHEMA + "stateVector")]

    def __getOrbitStateVector(self, node):
        timestampStr = [i.text for i in node.iter(SCHEMA + "timeStamp")][0]
        timestamp = self.__datetime_from_string(timestampStr)
        x = float([i.text for i in node.iter(SCHEMA + "xPosition")][0])
        y = float([i.text for i in node.iter(SCHEMA + "yPosition")][0])
        z = float([i.text for i in node.iter(SCHEMA + "zPosition")][0])
        vel_x = float([i.text for i in node.iter(SCHEMA + "xVelocity")][0])
        vel_y = float([i.text for i in node.iter(SCHEMA + "yVelocity")][0])
        vel_z = float([i.text for i in node.iter(SCHEMA + "zVelocity")][0])
        return OrbitStateVector(timestamp, (x, y, z), (vel_x, vel_y, vel_z))
        
    def getOrbitPosition(self, timestamp):
        return self.__orbitCalculator.getPosition(timestamp)
    
    def getOrbitVelocity(self, timestamp):
        return self.__orbitCalculator.getVelocity(timestamp)
        
    def __createSlantRangeCalculator(self, root):
        groundToSlantRangeEntries = []
        for node in root.findall(".//" + SCHEMA + "slantRangeToGroundRange"):
            entry = self.__getGroundToSlantRangeEntry(node)
            groundToSlantRangeEntries.append(entry)
            
        self.__slantRangeCalculator = SlantRangeCalculator(groundToSlantRangeEntries)
    
    def __getGroundToSlantRangeEntry(self, node):
        timestampString = [i.text for i in node.iter(SCHEMA + "zeroDopplerAzimuthTime")][0]
        timestamp = self.__datetime_from_string(timestampString)
        firstSampleSlantRange = float([i.text for i in node.iter(SCHEMA + "slantRangeTimeToFirstRangeSample")][0])
        groundRangeOrigin = float([i.text for i in node.iter(SCHEMA + "groundRangeOrigin")][0])
        groundToSlantRangeCoefficients = tuple(float(c) for c in [i.text for i in node.iter(SCHEMA + "groundToSlantRangeCoefficients")][0].split())
        return GroundToSlantRangeEntry(timestamp, firstSampleSlantRange, groundRangeOrigin, groundToSlantRangeCoefficients)
    
    def getSlantRange(self, pixel, line):
        timestamp = self.getTime(pixel, line)
        groundRange = self.getGroundRange(pixel)
        return self.__slantRangeCalculator.getSlantRange(groundRange, timestamp)
    
    def __datetime_from_string(self, datetimeStr):
        '''given the date time string in the standard RCM format, returns a datetime object.'''
        
        #parse out the datetime string
        datetimeList = list(re.match(RCM_DATE_FORMAT, datetimeStr).groups())
        
        #save the last element (which is the decimal part of the seconds) to a float
        try:
            datetimeList[-1] = float(datetimeList[-1])*1e6  # the decimal part is optional so it may be None
        except:
            datetimeList[-1] = 0    # exception probably due to no decimal seconds, set microseconds to zero
            
        #now convert all elements to int (from string)
        datetimeList = list(map(int, datetimeList))
        
        #conveniently, the output from parsing the text is the same order that is needed by datetime
        #list should now be all integers (year, month, day, hour, minute, second, microsecond)
        return datetime.datetime(*datetimeList)
    
    def _readNoiseLevel(self, calib):
        '''Read NoiseLevel.
        calib = 'Beta Nought', 'Sigma Nought', or 'Gamma'
        returns a list [pixelFirstNoiseValue, stepSize, noiseValues] where noiseValues is a list of noise levels.
        Section 7.7 of RCM-SP-53-0419 Image Product Format Definition 2-5.pdf
        '''
        
        self.azimuthNoiseParams = [[] for _ in xrange(self._ds.RasterCount)]
        
        noise = []
        for idx, pol in enumerate(self.pols):
            noise_suffix = "noiseLevels_" + pol.upper() + ".xml"
            noise_xml = os.path.join(self.productPath, "calibration", noise_suffix)  # RCM product.xml
            if not os.path.isfile(noise_xml):           # RCM NITF
                noise_xml = os.path.join(os.path.dirname(self.productPath), "metadata", "calibration", noise_suffix)
                if not os.path.isfile(noise_xml):
                    gvutils.error('{} is not a valid file!'.format(noise_xml))
                    return
            
            root = ET.parse(noise_xml).getroot()
            
            for node in root.iter(SCHEMA + "referenceNoiseLevel"):
                if [i.text for i in node.iter(SCHEMA + "sarCalibrationType")][0] == calib:
                    
                    pixelFirstNoiseValue = int([i.text for i in node.iter(SCHEMA + "pixelFirstNoiseValue")][0])
                    stepSize = float([i.text for i in node.iter(SCHEMA + "stepSize")][0])
                    noiseValues = [i.text for i in node.iter(SCHEMA + "noiseLevelValues")][0].split()
                    noiseValues = list(map(float, noiseValues))
                    break
           
            noise.append([pixelFirstNoiseValue, stepSize, noiseValues])
            
            # read azimuthNoiseLevelScaling if they exist, otherwise set azimuthNoiseLevelScaling to None
            if root.findall(SCHEMA +"azimuthNoiseLevelScaling"):
                for node in root.iter(SCHEMA + "azimuthNoiseLevelScaling"):
                    beam = [i.text for i in node.iter(SCHEMA + "beam")][0]
                    stepSize = float([i.text for i in node.iter(SCHEMA + "stepSize")][0])
                    numberOfNoiseLevelScalingValues = int([i.text for i in node.iter(SCHEMA + "numberOfNoiseLevelScalingValues")][0])
                    noiseValues = [i.text for i in node.iter(SCHEMA + "noiseLevelScalingValues")][0].split()
                    noiseValues = list(map(float, noiseValues))
                    try:
                        # only for Spotlight products
                        lineFirstNoiseScalingValue = int([i.text for i in node.iter(SCHEMA + "lineFirstNoiseScalingValue")][0])
                    except:
                        lineFirstNoiseScalingValue = None
                    
                    self.azimuthNoiseParams[idx].append([beam, stepSize, numberOfNoiseLevelScalingValues, noiseValues, lineFirstNoiseScalingValue])
            else:
                self.azimuthNoiseParams = None
        
        #print self.azimuthNoiseParams
        return noise
    
    def getBetaNoughtNoiseParams(self):
        '''returns a tuple (stepSize, values) where value is a list of noise levels.'''
        # if it has already been read in, return the saved value
        try:
            return self.betaNoughtNoiseParams
        except:
            self.betaNoughtNoiseParams = self._readNoiseLevel('Beta Nought')
            return self.betaNoughtNoiseParams
        
    def getSigmaNoughtNoiseParams(self):
        '''returns a tuple (stepSize, values) where value is a list of noise levels.'''
        # if it has already been read in, return the saved value
        try:
            return self.sigmaNoughtNoiseParams
        except:
            self.sigmaNoughtNoiseParams = self._readNoiseLevel('Sigma Nought')
            return self.sigmaNoughtNoiseParams
        
    def getGammaNoiseParams(self):
        '''returns a tuple (stepSize, values) where value is a list of noise levels.'''
        #if it has already been read in, return the saved value
        try:
            return self.gammaNoiseParams
        except:
            self.gammaNoiseParams = self._readNoiseLevel('Gamma')
            return self.gammaNoiseParams
        
    def getFullWidthSigmaNoughtNoiseLevel(self, line, band=None):  # band from 1 to 4
        if isinstance(band, int):
            FullWidthSigmaNoughtNoise = [self.getSigmaNoughtNoiseLevel(pixel, line=line)[band-1] for pixel in range(self.nSamples)]
            return FullWidthSigmaNoughtNoise
    
    def getSigmaNoughtNoiseLevel(self, pixel, line=None):
        '''returns interpolated noise level at (p, l) per band'''
        
        noise = []
        for band in range(self._ds.RasterCount):
            pixelFirstNoiseValue, stepSize, noiseValues = self.sigmaNoughtNoiseParams[band]
            noise.append(ATD_Utils.interpolate_values(pixelFirstNoiseValue, stepSize, noiseValues, pixel))
        
        if not self.azimuthNoiseParams:
            return noise

        # apply azimuth_scaling for ScanSAR and spotlight
        if line is None: line = self.nLines/2
        azimuth_scaling = []        # one value per band
        
        if self.isSpotlight():  # Spotlight has azimuthNoiseLevelScaling per band
            for band in range(self._ds.RasterCount):
                stepSize = self.azimuthNoiseParams[band][0][1]      # Spotlight can have multi-bands, but only one beam 
                noiseValues = self.azimuthNoiseParams[band][0][3]
                lineFirstNoiseScalingValue = self.azimuthNoiseParams[band][0][4]
                azimuth_scaling.append(ATD_Utils.interpolate_values(lineFirstNoiseScalingValue, stepSize, noiseValues, line))
            #print 'spot', azimuth_scaling
        
        else:   # ScanSAR or High Resolution products has azimuthNoiseLevelScaling per band per burst
            #try:
                _, topMinLine, _, bottomMaxLine, burst_height, _ = self.bursts[self.getBeam(pixel)]
                
                for band in range(self._ds.RasterCount):
                    for param in self.azimuthNoiseParams[band]:
                        if (not topMinLine <= line <= bottomMaxLine):   # out of burst range
                            azimuth_scaling.append(0.)
                            break
                        
                        if param[0] == self.getBeam(pixel):
                            stepSize = abs(param[1])           # not sure the negative stepSize!
                            numberOfNoiseLevelScalingValues = param[2]
                            noiseValues = param[3]
                            
                            # RCM-SP-53-0419 Image Product Format Definition page 161
                            # the trough of the azimuthNoiseLevelScaling is at the centre line of the bursts 
                            # The edges of the bursts are covered by adjacent ones, thus line_idx starts zero at the centre line, i.e., min azimuthNoiseLevelScaling
                            line_idx = (line - topMinLine) % burst_height - burst_height/2.
                            azimuth_scaling.append(ATD_Utils.interpolate_values(-(numberOfNoiseLevelScalingValues-1)*stepSize/2., stepSize, noiseValues, line_idx))
                            break
            #except:
            #    return noise    # skip adding azimuth_scaling when fails
        #print noise, azimuth_scaling, list(map(add, noise, azimuth_scaling))
        return list(map(add, noise, azimuth_scaling))
        
    def _readLUT(self, calib):
        '''Read a lookup table
        calib = 'Beta', 'Sigma', or 'Gamma'
        returns offset, gains where gains is a list of values for each pol in [HH, HV, VH, VV]
        offset = [offset1, offset2, offset3, offset4]
        gains = [[gains_1 ...] [gains_2 ...] [gains_2 ...] [gains_2 ...]]
        '''
        # Setup blank arrays
        #offset = numpy.empty((len(self.pols), 1), numpy.float32)
        #gains = numpy.empty((len(self.pols), self.xsize), numpy.float32)
        offset = []
        gains = []
        
        # read in gains and offsets
        for pol in self.pols:
            lut_suffix = "lut" + calib + "_" + pol + ".xml"
            lutFile = os.path.join(self.productPath, "calibration", lut_suffix)  # RCM product.xml
            if not os.path.isfile(lutFile):     # RCM NITF
                lutFile = os.path.join(os.path.dirname(self.productPath), "metadata", "calibration", lut_suffix)
                if not os.path.isfile(lutFile):
                    gvutils.error('{} is not a valid file!'.format(lutFile))
                    return
            
            root = ET.parse(lutFile).getroot()
            pixelFirstLutValue = int([i.text for i in root.iter(SCHEMA + "pixelFirstLutValue")][0])
            stepSize = [i.text for i in root.iter(SCHEMA + "stepSize")][0]
            offset.append(float([i.text for i in root.iter(SCHEMA + "offset")][0]))
            gg = list(map(float, [i.text for i in root.iter(SCHEMA + "gains")][0].split()))
            
            if stepSize == '1':
                gains.append(gg)
            elif stepSize == '-1':
                gg.reverse()
                gains.append(gg)
            else:
                #print('gg={}'.format(gg[:15]))
                gains_interpolate = ATD_Utils.interpolate_values(pixelFirstLutValue, float(stepSize), gg, numpy.arange(self.xsize), 'linear')
                gains.append(gains_interpolate)
        return offset, gains

    def getBetaNoughtLUT(self):
        '''returns offset, gains where gains contain an array for each band.'''
        # if it has already been read in, return the saved value
        try:
            return self.betaNoughtLUT
        except:
            self.betaNoughtLUT = self._readLUT('Beta')
            return self.betaNoughtLUT
    
    def getSigmaNoughtLUT(self):
        '''returns offset, gains where gains contain an array for each band.'''
        #if it has already been read in, return the saved value
        try:
            return self.sigmaNoughtLUT
        except:
            # return LUT for all bands
            self.sigmaNoughtLUT = self._readLUT('Sigma')
            return self.sigmaNoughtLUT
        
    def getGammaLUT(self):
        '''returns offset, gains where gains contain an array for each band.'''
        # if it has already been read in, return the saved value
        try:
            return self.gammaLUT
        except:
            self.gammaLUT = self._readLUT('Gamma')
            return self.gammaLUT
    
    def __applyLUT(self, extents, calib, band=None, data=None, oversample=False, SLCProductType=0):
        '''Returns an array of calibrated values for the given extents.  
        extents is a tuple: (startPixel, startLine, nPixels, nLines)
        offset and gains are from the desired lookup table.
        Uses data as raw DN if given, otherwise reads from the dataset.
        If data is None, all bands will be read unless band is specified.
        If data is not None, band parameter is ignored.'''
        
        # added by Shawn Gong in July 2020
        # if SLCProductType = 0 default, returns conventional Sigma0 = |DN|^2/LUT^2
        # if SLCProductType = 1, returns ComplexScatter = (DN/LUT)
        # if SLCProductType = 2, returns ComplexSigma = (DN^2/LUT^2)
        
        # offset = [offset1, offset2, offset3, offset4]
        # gains  = [[gains_1 ], [gains_2 ], [gains_3 ], [gains_4 ]]
        if calib == 'Sigma':
            offset, gains = self.getSigmaNoughtLUT()
        elif calib == 'Beta':
            offset, gains = self.getBetaNoughtLUT()
        elif calib == 'Gamma':
            offset, gains = self.getGammaLUT()
        
        # convert as numpy.array floats
        gains = numpy.asarray(gains).astype(numpy.float64)

        startPixel, startLine, nPixels, nLines = extents
        
        # read in data if not given
        if data is None:
            if isinstance(band, int) and band > 0:    # get just a single band, band ranges from 1 to 4
                data = self._ds.GetRasterBand(band).ReadAsArray(*extents)
            else:
                data = self._ds.ReadAsArray(*extents)

        # perform the error check for given band and data 
        # add by Shawnn Gong on May 29, 2022
        if isinstance(band, int) and band > 0 and data.ndim == 3:
            raise ValueError("a single RCM band is specified, but 'data' contains multiple bands!")
        if band is None:
            if data.ndim == 3 and data.shape[0] != self._ds.RasterCount:
                raise ValueError("'data' does not contain all bands.")
            if data.ndim == 2 and self._ds.RasterCount > 1:
                msg = " ** 'data' is given as a single band array, but this RCM is a multi-band image.\n"
                msg += " ** Therefore we can't determine which LUT[band] to apply, exit!"
                raise ValueError(msg)
        
        if oversample:
            import SARSpectrumAnalysis
            data = SARSpectrumAnalysis.oversample(self, startPixel, startPixel+nPixels, startLine, startLine+nLines, image_array=data, nonComplexOK=True)
        
        # more precision array, whether to convert to real array depends on SLCProductType parameter
        if numpy.iscomplexobj(data):
            if SLCProductType > 0:
                data = data.astype(numpy.complex128)
            else:
                data = numpy.absolute(data)
                data = data.astype(numpy.float64)
        else:
            data = data.astype(numpy.float64)
        
        if isinstance(band, int) and band > 0:    # read a single band (between 1 and 4)
            # extract the part of the lookup table that is needed
            LUT = gains[band-1][startPixel:startPixel+nPixels]
            if oversample or data.shape[-1] != len(LUT):
                # if oversample, interpolate the LUT
                LUT = numpy.interp(numpy.linspace(0, len(LUT)-1, num=data.shape[-1]), list(range(len(LUT))), LUT)
            
            if self.isComplex:
                if SLCProductType == 1:     # ComplexScatter
                    data /= LUT
                else:
                    #calibrated_values = data**2 / LUT**2
                    # The above line is correct but doing it in steps as below will allow numpy to perform in-place operations. This can save a lot of memory that would be wasted in an unnecessary copy - especially for a large image chip
                    data **= 2
                    LUT **= 2
                    data /= LUT
            else:
                #calibrated_values = (data**2 + offset) / LUT
                # The above line is correct but doing it in steps as below will allow numpy to perform in-place operations. This can save a lot of memory that would be wasted in an unnecessary copy - especially for a large image chip
                data **= 2
                data += offset[band-1]
                data /= LUT
        else:           # all bands read in
            if self._ds.RasterCount == 1:   # single-band: data.shape = (row, col)
                LUT = gains[0][startPixel:startPixel+nPixels]
                if oversample or data.shape[-1] != len(LUT):
                    # if oversample, interpolate the LUT
                    LUT = numpy.interp(numpy.linspace(0, len(LUT)-1, num=data.shape[-1]), list(range(len(LUT))), LUT)
                
                if self.isComplex:
                    if SLCProductType == 1:     # ComplexScatter
                        data /= LUT
                    else:
                        data **= 2
                        LUT **= 2
                        data /= LUT
                else:
                    data **= 2
                    data += offset[0]
                    data /= LUT
            else:                           # multi-band: data.shape = (band, row, col)
                for idx in range(self._ds.RasterCount):
                    LUT = gains[idx][startPixel:startPixel+nPixels]
                    if oversample or data.shape[-1] != len(LUT):
                        # if oversample, interpolate the LUT
                        LUT = numpy.interp(numpy.linspace(0, len(LUT)-1, num=data.shape[-1]), list(range(len(LUT))), LUT)
                    
                    if self.isComplex:
                        if SLCProductType == 1:
                            data[idx] /= LUT
                        else:
                            data[idx] **= 2
                            LUT **= 2
                            data[idx] /= LUT
                    else:
                        data[idx] **= 2
                        data[idx] += offset[idx]
                        data[idx] /= LUT

        return data
    
    def getBetaNought(self, extents, band=None, data=None):     # band from 1 to 4
        '''Returns an array of beta nought calibrated values for the given extents.  
        extents is a tuple (startPixel, startLine, nPixels, nLines)
        Uses data as raw dn if given, otherwise reads from the dataset.
        If data is None, all bands will be read unless band is specified.
        If data is not None, band parameter is ignored.'''
        
        return self.__applyLUT(extents, 'Beta', band=band, data=data)
    
    def getSigmaNought(self, extents, band=None, data=None, oversample=False):     # band from 1 to 4
        '''Returns an array of sigma nought calibrated values for the given extents.  
        extents is a tuple (startPixel, startLine, nPixels, nLines)
        Uses data as raw dn if given, otherwise reads from the dataset.
        If data is None, all bands will be read unless band is specified.
        If data is not None, band parameter is ignored.'''
        
        return self.__applyLUT(extents, 'Sigma', band=band, data=data, oversample=oversample)
    
    def getSigmaNoughtComplexScatter(self, extents, band=None, data=None, oversample=False):
        ''' This is for DRDC usage
        returns ComplexScatter = (DN/LUT) '''
        
        return self.__applyLUT(extents, 'Sigma', band=band, data=data, oversample=oversample, SLCProductType=1)
    
    def getComplexSigmaNought(self, extents, band=None, data=None, oversample=False):
        ''' This is for DRDC usage
        returns ComplexSigma = (DN^2/LUT^2) '''
        
        return self.__applyLUT(extents, 'Sigma', band=band, data=data, oversample=oversample, SLCProductType=2)
    
    def getGamma(self, extents, band=None, data=None):     # band from 1 to 4
        '''Returns an array of gamma calibrated values for the given extents.  
        extents is a tuple (startPixel, startLine, nPixels, nLines)
        Uses data as raw dn if given, otherwise reads from the dataset.
        If data is None, all bands will be read unless band is specified.
        If data is not None, band parameter is ignored.'''
        
        return self.__applyLUT(extents, 'Gamma', band=band, data=data)
        
    def old_getSigmaNoughtComplexScatter(self, extents, band=None, data=None, oversample=False):     # band from 1 to 4
        '''Returns an array of Sigma Nought Complex Scatter for SLC.'''
        
        startPixel, startLine, nPixels, nLines = extents
        
        offset, gains = self.getSigmaNoughtLUT()
        gains = numpy.array(gains)
        
        #read in data if not given
        if data is None:
            if isinstance(band, int) and band > 0:    # get just a single band, band ranges from 1 to 4
                data = self._ds.GetRasterBand(band).ReadAsArray(*extents)
            else:
                data = self._ds.ReadAsArray(*extents)
            
        if numpy.iscomplexobj(data):
            data = data.astype(numpy.complex128)
        else:
            data = data.astype(numpy.float64)
        
        if oversample:
            import SARSpectrumAnalysis
            data = SARSpectrumAnalysis.oversample(self, startPixel, startPixel+nPixels, startLine, startLine+nLines, image_array=data, nonComplexOK=True)
        
        if band:        # single band from 1 to 4
            gg = gains[band-1][startPixel:startPixel+nPixels]
            if oversample or data.shape[-1] != len(gg):
                gg = numpy.interp(numpy.linspace(0, len(gg)-1, num=data.shape[-1]), list(range(len(gg))), gg)

            data /= gg
        else:           # for all bands
            if self._ds.RasterCount == 1:   # data.shape = (row, col)
                gg = gains[0][startPixel:startPixel+nPixels]
                # if oversample, interpolate the LUT to fit
                if oversample or data.shape[-1] != len(gg):
                    gg = numpy.interp(numpy.linspace(0, len(gg)-1, num=data.shape[-1]), list(range(len(gg))), gg)
                
                data /= gg
            else:                           # data.shape = (band, row, col)
                for band_idx in range(self._ds.RasterCount):
                    gg = gains[band_idx, startPixel:startPixel+nPixels]
                    # if oversample, interpolate the LUT to fit
                    if oversample or data.shape[-1] != len(gg):
                        gg = numpy.interp(numpy.linspace(0, len(gg)-1, num=data.shape[-1]), list(range(len(gg))), gg)
                    
                    data[band_idx] /= gg
        return data
    
    def isScanSAR(self):
        return ("SC" in self.beamModeMnemonic and "RGSC" not in self.beamModeMnemonic)
    
    def isSpotlight(self):
        return self.beamModeMnemonic[:2] == "FS"
    
    def getPlatformVelocity(self):
        return self.platformVelocity
    
    def _calcIncAngles(self):
        '''Returns the incidence angles in degrees for all pixels.'''
        inc_ang_file = os.path.join(self.productPath, "calibration", "incidenceAngles.xml")
        if not os.path.isfile(inc_ang_file):        # RCM NITF case
            inc_ang_file = os.path.join(os.path.dirname(self.productPath), "metadata", "calibration", "incidenceAngles.xml")
            if not os.path.isfile(inc_ang_file):
                gvutils.error('{} is not a valid file!'.format(inc_ang_file))
                return
        
        root = ET.parse(inc_ang_file).getroot()
        pixelFirstAnglesValue = int([i.text for i in root.iter(SCHEMA + "pixelFirstAnglesValue")][0])
        stepSize = [i.text for i in root.iter(SCHEMA + "stepSize")][0]
        angles = list(map(float, [i.text for i in root.iter(SCHEMA + "angles")]))
        
        if stepSize == '1':
            self.incidenceAngles = angles
        elif stepSize == '-1':
            angles.reverse()
            self.incidenceAngles = angles
        else:
            stepSize = float(stepSize)
            self.incidenceAngles = ATD_Utils.interpolate_values(pixelFirstAnglesValue, stepSize, angles, numpy.arange(self.xsize), 'linear')
    
    def getIncidenceAngle(self, pixel, line=0):
        '''Returns an incidence angle in degrees at a pixel location.'''
        pixel = min(self.nSamples-1, max(0, int(pixel)))
        return self.incidenceAngles[pixel]

    def getBurstNumber(self, pixel, line):
        if pixel < 0 or pixel >= self.nSamples or line < 0 or line >= self.nLines:
            #print('Error getBurst for pixel {}, line {}'.format(pixel, line))
            return 0
        
        root = ET.fromstring(self.xmlLines)
        
        if self.grd_ScanSAR:
            node_1 = root.findall(".//" + SCHEMA + "grdBurstMap")[0]
            for node in node_1.iter(SCHEMA + "burstAttributes"):
                topLeftLine = int([i.text for i in node.iter(SCHEMA + "topLeftLine")][0])
                topLeftPixel = int([i.text for i in node.iter(SCHEMA + "topLeftPixel")][0])
                bottomRightLine = int([i.text for i in node.iter(SCHEMA + "bottomRightLine")][0])
                bottomRightPixel = int([i.text for i in node.iter(SCHEMA + "bottomRightPixel")][0])
                #print topLeftLine, topLeftPixel, bottomRightLine, bottomRightPixel
                #print (pixel-bottomRightPixel)*(pixel-topLeftPixel), (line-bottomRightLine)*(line-topLeftLine)
                if (pixel-bottomRightPixel)*(pixel-topLeftPixel) <= 0 and (line-bottomRightLine)*(line-topLeftLine) <= 0:
                    return int(node.attrib['burst'])
        elif self.mlc_ScanSAR:
            node_1 = root.findall(".//" + SCHEMA + "mlcBurstMap")[0]
            for node in node_1.iter(SCHEMA + "burstAttributes"):
                topLeftLine = int([i.text for i in node.iter(SCHEMA + "topLeftLine")][0])
                topLeftPixel = int([i.text for i in node.iter(SCHEMA + "topLeftPixel")][0])
                bottomRightLine = int([i.text for i in node.iter(SCHEMA + "bottomRightLine")][0])
                bottomRightPixel = int([i.text for i in node.iter(SCHEMA + "bottomRightPixel")][0])
                if (pixel-bottomRightPixel)*(pixel-topLeftPixel) <= 0 and (line-bottomRightLine)*(line-topLeftLine) <= 0:
                    return int(node.attrib['burst'])
        elif self.slc_ScanSAR:
            node_1 = root.findall(".//" + SCHEMA + "slcBurstMap")[0]
            for node in node_1.iter(SCHEMA + "burstAttributes"):
                topLeftLine = int([i.text for i in node.iter(SCHEMA + "lineOffset")][0])
                topLeftPixel = int([i.text for i in node.iter(SCHEMA + "pixelOffset")][0])
                bottomRightLine = topLeftLine + int([i.text for i in node.iter(SCHEMA + "numLines")][0]) - 1
                bottomRightPixel = topLeftPixel + int([i.text for i in node.iter(SCHEMA + "samplesPerLine")][0]) - 1
                if (pixel-bottomRightPixel)*(pixel-topLeftPixel) <= 0 and (line-bottomRightLine)*(line-topLeftLine) <= 0:
                    return int(node.attrib['burst'])
        return 0

    def getDopplerRateParams(self, pixel, line):
        root = ET.fromstring(self.xmlLines)

        # non-scanSAR has a single dopplerRate entry
        if int([i.text for i in next(root.iter(SCHEMA + "dopplerRate")).findall(".//" + SCHEMA + "numberOfEstimates")][0]) == 1:
            reftime = float([i.text for i in root.iter(SCHEMA + "dopplerRateReferenceTime")][0])
            coefficients = [float(c) for c in [i.text for i in root.iter(SCHEMA + "dopplerRateCoefficients")][0].split()]
            return reftime, coefficients
        
        # scanSAR or High Resolution dopplerRate varies per burst
        for node in root.iter(SCHEMA + "dopplerRateEstimate"):
            if self.getBurstNumber(pixel, line) == int(node.attrib['burst']):
                reftime = float([i.text for i in node.iter(SCHEMA + "dopplerRateReferenceTime")][0])
                coefficients = [float(c) for c in [i.text for i in node.iter(SCHEMA + "dopplerRateCoefficients")][0].split()]
                return reftime, coefficients
    
    def readBurstAttributes(self):
        self.slc_ScanSAR = False
        self.grd_ScanSAR = False
        self.mlc_ScanSAR = False
        
        if self.xmlLines.find('slcBurstMap') > -1:
            self.slc_ScanSAR = True
        elif self.xmlLines.find('grdBurstMap') > -1:
            self.grd_ScanSAR = True
        elif self.xmlLines.find('mlcBurstMap') > -1:
            self.mlc_ScanSAR = True
        else:
            self.bursts = None
            self.burst_pts_lon_lat = None
            #print("No BurstMap in product.xml")
            return None, None, None

        self.burst_pts_lon_lat = []
        self.bursts = {}
            # bursts[beam][0] = topLeftPixel
            # bursts[beam][1] = topLeftLine
            # bursts[beam][2] = bottomRightPixel
            # bursts[beam][3] = bottomRightLine
            # bursts[beam][4] = # of lines
            # bursts[beam][5] = # of bursts
        for beam in self.beams:
            self.bursts[beam] = [-1, 999999., -1, -1, 0, 0]
        
        root = ET.fromstring(self.xmlLines)
        if self.slc_ScanSAR:
            node_1 = root.findall(".//" + SCHEMA + "slcBurstMap")[0]     # just get the first pole in case there are more than 1
            self.numberBursts = int([i.text for i in node_1.iter(SCHEMA + "numberOfEntries")][0])
            self.numberBurstCycles = self.numberBursts/self.numberOfBeams
            
            #get the line/pixel coordinates of each burst from xml
            #just save the top and bottom points for each beam and assume that the beam seam is a straight line connecting them
            for node in node_1.iter(SCHEMA + "burstAttributes"):
                beam = node.attrib['beam']
                self.bursts[beam][5] += 1
                
                topLeftLine = int([i.text for i in node.iter(SCHEMA + "lineOffset")][0])
                topLeftPixel = int([i.text for i in node.iter(SCHEMA + "pixelOffset")][0])
                bottomRightLine = topLeftLine + int([i.text for i in node.iter(SCHEMA + "numLines")][0]) - 1
                bottomRightPixel = topLeftPixel + int([i.text for i in node.iter(SCHEMA + "samplesPerLine")][0]) - 1
                self.bursts[beam][4] += abs(topLeftLine - bottomRightLine)
                
                ll = []
                ll.extend(self.pixel_to_view(topLeftPixel, topLeftLine))
                ll.extend(self.pixel_to_view(bottomRightPixel, topLeftLine))
                ll.extend(self.pixel_to_view(bottomRightPixel, bottomRightLine))
                ll.extend(self.pixel_to_view(topLeftPixel, bottomRightLine))
                ll.extend(self.pixel_to_view(topLeftPixel, topLeftLine))
                self.burst_pts_lon_lat.append(ll)
                
                #if the current node contains data closer to the top or bottom, replace the topPixel or bottomPixel values
                if topLeftLine < self.bursts[beam][1]:
                    self.bursts[beam][1] = topLeftLine
                    self.bursts[beam][0] = topLeftPixel
                if bottomRightLine > self.bursts[beam][3]:
                    self.bursts[beam][3] = bottomRightLine
                    self.bursts[beam][2] = bottomRightPixel
        
        elif self.grd_ScanSAR or self.mlc_ScanSAR:
            if self.grd_ScanSAR:
                node_1 = root.findall(".//" + SCHEMA + "grdBurstMap")[0]     # just get the first pole in case there are more than 1
            elif self.mlc_ScanSAR:
                node_1 = root.findall(".//" + SCHEMA + "mlcBurstMap")[0]     # just get the first pole in case there are more than 1
            
            self.numberBursts = int([i.text for i in node_1.iter(SCHEMA + "numberOfEntries")][0])
            self.numberBurstCycles = self.numberBursts/self.numberOfBeams
            
            #get the line/pixel coordinates of each burst from xml
            #just save the top and bottom points for each beam and assume that the beam seam is a straight line connecting them
            for node in node_1.iter(SCHEMA + "burstAttributes"):
                beam = node.attrib['beam']
                self.bursts[beam][5] += 1
                
                topLeftLine = int([i.text for i in node.iter(SCHEMA + "topLeftLine")][0])
                topLeftPixel = int([i.text for i in node.iter(SCHEMA + "topLeftPixel")][0])
                bottomRightLine = int([i.text for i in node.iter(SCHEMA + "bottomRightLine")][0])
                bottomRightPixel = int([i.text for i in node.iter(SCHEMA + "bottomRightPixel")][0])
                self.bursts[beam][4] += abs(topLeftLine - bottomRightLine)
                
                ll = []
                ll.extend(self.pixel_to_view(topLeftPixel, topLeftLine))
                ll.extend(self.pixel_to_view(bottomRightPixel, topLeftLine))
                ll.extend(self.pixel_to_view(bottomRightPixel, bottomRightLine))
                ll.extend(self.pixel_to_view(topLeftPixel, bottomRightLine))
                ll.extend(self.pixel_to_view(topLeftPixel, topLeftLine))
                self.burst_pts_lon_lat.append(ll)
                
                # if the current node contains data closer to the top or bottom, replace the topPixel or bottomPixel values
                if topLeftLine < self.bursts[beam][1]:
                    self.bursts[beam][1] = topLeftLine
                    self.bursts[beam][0] = topLeftPixel
                if bottomRightLine > self.bursts[beam][3]:
                    self.bursts[beam][3] = bottomRightLine
                    self.bursts[beam][2] = bottomRightPixel
                
        for beam in self.beams:
            self.bursts[beam][4] /= float(self.bursts[beam][5])
        #print self.bursts

        if self.numberOfBeams == 1:
            return None, None, None
        else:
            # put the top and bottom coordinates in the correct order
            beams = []
            tops = []
            bottoms = []
            for beam in self.beams:
                beams.append(beam)
                topMinPixel, topMinLine, bottomMaxPixel, bottomMaxLine = self.bursts[beam][:4]
                tops.append([topMinPixel, topMinLine])
                bottoms.append([bottomMaxPixel, bottomMaxLine])
            
            self.burst_bin_tops = [[], []]  #[[beams], [left most pixel#]] sorted from left to right
            for beam, top in zip(beams, tops):
                self.burst_bin_tops[0].append(beam)
                self.burst_bin_tops[1].append(top[0])
            
            # for descending
            if self.burst_bin_tops[1][0] > self.burst_bin_tops[1][1]:
                self.burst_bin_tops[0].reverse()
                self.burst_bin_tops[1].reverse()
            
            # set left edge pixel# to 0
            self.burst_bin_tops[1][0] = 0
            
            return beams, tops, bottoms     # tops and bottoms are sorted from left to right (0 to nSamples)
    
    def getBeam(self, pixel, line=None):
        pixel = max(pixel, 0)
        pixel = min(self.nSamples-1, pixel)
        
        if self.numberOfBeams == 1:
            return self.beams[0]    # single beam
        else:
            return self.burst_bin_tops[0][bisect(self.burst_bin_tops[1], pixel)-1]   # for scanSAR
        
    def getPRF(self, pixel, line, effectivePerChannelPRF=True):
        '''Returns the pulse repetition frequency in Hz at the given image location.
        pixel, line location is required for a scanSAR image, otherwise it is not used.
        For a scanSAR image, the beam seams will first be detected if it has not already been done.  
        If effectivePerChannelPRF is false, the system PRF from product.xml will be returned.
        Otherwise factors like transmit polarizations and number of receive antennas will be taken into account.'''
        
        beam = self.getBeam(pixel)
        if not beam: return
        
        prf = self.PRF[beam]     # assuming the PRF list is in order of increasing incidence angle

        #account for other factors that can impact the effective PRF
        if effectivePerChannelPRF:
            if self.tx_dualPol:
                prf /= 2.    # per-channel PRF is half of system PRF when 2 polarizations are transmitting alternately 
            
            # PRF is effectively increased by the number of receive antennas
            # this seems to be correct for oversampling ultrafine products but
            # not for azimuth ambiguity marking.  Oversampling now just gets the
            # effective PRF by taking nLines/(endTime-startTime),
            # so remove the adjustment here.
            #prf *= self.nRecieveAntennas
            
        return prf
    
    def getProcessingTime(self):
        '''returns the datetime when the image was processed.'''
        processingTime = [i.text for i in node.iter(SCHEMA + "processingTime")][0]
        return self.__datetime_from_string(processingTime)
    
    def getTime(self, pixel, line):
        '''returns the datetime of the given line/pixel location'''
        
        if not self.isInZeroDopplerCoords():
            raise RuntimeError('getTime( ) does not support geocoded products.')

        # just do a linear interpolation between the first and last zero doppler azimuth times
        def total_seconds(td):
            return td.days*24*3600 + td.seconds + td.microseconds/10.0**6
            
        # get the timedelta - this may be negative if the line order has been reversed
        td = total_seconds(self.zeroDopplerTimeLastLine - self.zeroDopplerTimeFirstLine)
        
        return self.zeroDopplerTimeFirstLine + datetime.timedelta(seconds=line*td/self.nLines)
    
    def getDopplerCentroid(self, pixel, line):
        '''returns the doppler centroid in Hz at the given location.'''
        
        azimuthTime = self.getTime(pixel, line)
        slantRangeDistance = self.getSlantRange(pixel, line)
        slantRangeTime = 2. * slantRangeDistance/scipy.constants.c  # round trip (2 way) slant range time
        
        dopplerCentroid = self.__dopplerCentroidCalculator.getDopplerCentroid(azimuthTime, slantRangeTime)
        
        # adjust for doppler rate in spotlight products
        if self.isSpotlight():
            # same polynomial for doppler centroid and doppler rate
            params = self.getDopplerRateParams(pixel, line)
            doppler_rate = self.__evalDopplerCentroid(slantRangeTime, params[0], params[1])
            
            # get total seconds from the doppler centroid estimate time
            def total_seconds(td):
                return td.days*24*3600 + td.seconds + td.microseconds/10.0**6
            ref = self.dopplerCentroidEstimates[0].azimuthTime
            dt = total_seconds(azimuthTime-ref)

            # the doppler rate should be Hz/s
            # calculate the doppler centroid offset
            dc_offset = dt*doppler_rate
        else:
            dc_offset = 0

        return (dopplerCentroid - dc_offset)
    
    def __evalDopplerCentroid(self, slantRangeTime, dopplerCentroidRefTime, dopplerCentroidCoeffs):
        '''evaluate the doppler centroid given the parameters
        slantRangeTime is the 2 way slant range time to the point in seconds
        dopplerCentroidRefTime is the doppler centroid reference time in seconds
        dopplerCentroidCoeffs is a list of the doppler centroid estimation coefficients
        returns the dopplerCentroid in Hz.'''
        
        #doppler centroid = d0 + d1(tSR - t0) +d2(tSR - t0)^2 +d3(tSR - t0)^3 + d4(tSR - t0)^4
        #where d0~d5 are the doppler centroid coefficients given in product.xml
        #tSR is the two way slant range time to the point in question
        #t0 is the doppler centroid reference time given in product.xml
        
        slantRangeTimeDelta = slantRangeTime - dopplerCentroidRefTime   #(tSR - t0)
        return numpy.sum([coeff * slantRangeTimeDelta**i for i, coeff in enumerate(dopplerCentroidCoeffs)])
    
    def getGroundToSlantRangeCoeffs(self):
        entry = self.__getCenterGroundToSlantRangeEntry()
        return entry.groundToSlantRangeCoefficients
    
    def getGroundRangeOrigin(self):
        entry = self.__getCenterGroundToSlantRangeEntry()
        return entry.groundRangeOrigin
    
    def __getCenterGroundToSlantRangeEntry(self):
        return self.groundToSlantRangeEntries[len(self.groundToSlantRangeEntries)/2]
    
    def getGroundRange(self, pixel):
        if self.firstPixelIsNearRange():
            return pixel*self.pixelSpacing
        else:
            return (self.nSamples-pixel)*self.pixelSpacing
    
    def getGroundPixelSpacing(self, pt=None):
        if not self.isInSlantRange():
            return self.pixelSpacing
        else:
            if isinstance(pt, list) and len(pt) >= 2:
                inc_angle = self.getIncidenceAngle(pt[0], pt[1])
            else:
                inc_angle = (self.incAngleNear + self.incAngleFar)/2
            return self.pixelSpacing/math.sin(math.radians(inc_angle))
    
    def firstPixelIsNearRange(self):
        if not self.isInZeroDopplerCoords():
            raise NotInZeroDopplerCoordsException
        return self.pixelTimeOrdering.lower() == 'increasing'
        
    def firstLineIsFirstAcquired(self):
        if not self.isInZeroDopplerCoords():
            raise NotInZeroDopplerCoordsException
        return self.lineTimeOrdering.lower() == 'increasing'
    
    def getPlatformLookDirection(self, pixel=None, line_start=None, line_end=None):

        if ('left' in self.lookDirection.lower()):
            return math.fmod(self.getTrackAngle(pixel, line_start, line_end) - 90., 360.)
        else:
            return math.fmod(self.getTrackAngle(pixel, line_start, line_end) + 90., 360.)
    
    def getTrackAngle(self, pixel=None, line_start=None, line_end=None):
        """ calculate the satellite track angle on sphere surface using
        near-range's first_line/last_line lat/long pairs """
        
        # Note by Shawn Gong in Jul 2015: remove gdal.Transformer()
        #   because it fails when cross antemeridian
        # It is replaced by OpenEV's pixel_to_view()
        #   where cross antemeridian case has been updated in gvraster.c
        if pixel is None:
            try:
                if self.firstPixelIsNearRange():  # get the near-range azimuth line
                    pixel = 0
                else:
                    pixel = self.nSamples-1
            except:     # NotInZeroDopplerCoordsException
                return

        if line_start is None or line_start < 0 or line_start > self.nLines-1:
            line_start = 0
        if line_end is None or line_end > self.nLines-1 or line_end < 0:
            line_end = self.nLines-1
            
        firstline_lon, firstline_lat = self.pixel_to_view(pixel, line_start)
        lastline_lon, lastline_lat = self.pixel_to_view(pixel, line_end)

        # for ascending, swap the first and last lines
        if self.passDirection.lower() == "ascending":
            firstline_lon, lastline_lon = lastline_lon, firstline_lon
            firstline_lat, lastline_lat = lastline_lat, firstline_lat
            
        d = math.cos(math.radians((firstline_lat+lastline_lat)/2))
        xx = (lastline_lon-firstline_lon)*d
        yy = lastline_lat-firstline_lat

        return math.fmod((450.-math.degrees(math.atan2(yy, xx))), 360.)

    def isInZeroDopplerCoords(self):
        '''test whether the product is in zero doppler coordinates 
        (ie range is pixels, azimuth is lines).
        Raises an exception if the product type is unknown'''
        try:
            return PRODUCT_TYPES[self.productType][IN_ZERO_DOPPLER_COORDS]
        except:
            raise UnknownProductTypeException
    
    def isInSlantRange(self):
        '''test whether the product pixel spacing is in slant range.
        Raises an exception if the product type is unknown'''
        try:
            return PRODUCT_TYPES[self.productType][IN_SLANT_RANGE]
        except:
            raise UnknownProductTypeException
        
    def getProductType(self):
        return self.productType
    
    def isDetected(self):
        '''test whether the product samples are detected.
        Raises an exception if the product type is unknown'''
        try:
            return PRODUCT_TYPES[self.productType][IS_DETECTED]
        except:
            raise UnknownProductTypeException

    def getImageFiles(self):
        return self.imageFiles

    def getBandLabel(self, gdal_band_idx):
        '''Returns a string 'HH', 'HV', 'VH', 'VV' depending on the polarimetric interpretation of the given band index
        gdal_band_idx is the index of the band starting at 1 as reported by gdal.
        Use get_data(layerSourceIdx).get_band_number()to translate from layer source index to gdal band index.'''
        try:
            return self._ds.GetRasterBand(gdal_band_idx).GetMetadata()['POLARIMETRIC_INTERP']
        except KeyError:
            return SARLayer.getBandLabel(self, gdal_band_idx)

    def getOrbitRate(self):
        return ORBIT_RATE

    def getOrbitInclination(self):
        return ORBIT_INCLINATION


register(GDAL_DRIVER_SHORT_NAME, RCMLayer, openSubDatasets=False, handleVRT=True)

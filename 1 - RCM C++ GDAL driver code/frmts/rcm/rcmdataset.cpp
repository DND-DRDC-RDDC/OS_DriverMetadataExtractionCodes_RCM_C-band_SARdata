/******************************************************************************
*
* Project:  DRDC Ottawa GEOINT
* Purpose:  Radarsat Constellation Mission - XML Products (product.xml) driver
* Author:   Roberto Caron, MDA
*           on behalf of DRDC Ottawa
*
******************************************************************************
* Copyright (c) 2020, DRDC Ottawa
*
* Based on the RS2 Dataset Class
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
****************************************************************************/

#include <time.h>
#include <stdio.h>
#include <sstream>
//#include <conio.h>
#include "cpl_minixml.h"
#include "gdal_frmts.h"
#include "gdal_pam.h"
#include "ogr_spatialref.h"
#include "rcmdataset.h"
#include "gdal_io_error.h"

CPL_CVSID("$Id: rcmdataset.cpp 99999 2018-03-05 18:40:40Z rcaron $");

/* RCM has a special folder that contains all LUT, Incidence Angle and Noise Level files */
static char * CALIBRATION_FOLDER = "calibration";

/*** Function to test for valid LUT files ***/
static bool IsValidXMLFile(const char *pszPath, const char *pszLut)
{
	/* Return true for valid xml file, false otherwise */
	char *pszLutFile
		= VSIStrdup(CPLFormFilename(pszPath, pszLut, NULL));

	CPLXMLTreeCloser psLut(CPLParseXMLFile(pszLutFile));

	CPLFree(pszLutFile);

	if (psLut.get() == NULL)
	{
		const char msgError[] = "ERROR: Failed to open the LUT file %s";
		write_to_file_error(msgError, pszLutFile);

		CPLError(CE_Failure, CPLE_OpenFailed, msgError, pszLutFile);
	}

	return psLut.get() != NULL;
}


/*** Function to format calibration for unique identification for Layer Name ***/
/*
*  RCM_CALIB : { SIGMA0 | GAMMA0 | BETA0 | UNCALIB } : product.xml full path
*/
static CPLString FormatCalibration(const char *pszCalibName, const char *pszFilename)
{
	CPLString ptr;

	// Always begin by the layer calibrtion name
	ptr.append(szLayerCalibration);

	if (pszCalibName != NULL || pszFilename != NULL)
	{
		if (pszCalibName != NULL)
		{
			// A separator is needed before concat calibration name
			ptr.append(szLayerSeparator);
			// Add calibration name
			ptr.append(pszCalibName);
		}

		if (pszFilename != NULL)
		{
			// A separator is needed before concat full filename name
			ptr.append(szLayerSeparator);
			// Add full filename name
			ptr.append(pszFilename);
		}
	}
	else
	{
		// Always add a separator even though there are no name to concat
		ptr.append(szLayerSeparator);
	}

	/* return calibration format */
	return ptr;
}

/*** Function to concat 'metadata' with a folder separator with the filename 'product.xml'  ***/
/*
*  Should return either 'metadata\product.xml' or 'metadata/product.xml'
*/
static CPLString GetMetadataProduct()
{
	// Always begin by the layer calibrtion name
	CPLString ptr;
	ptr.append("metadata");
	ptr.append(szPathSeparator);
	ptr.append("product.xml");

	/* return metadata product filename */
	return ptr;
}

/*** check that the referenced dataset for each band has the
correct data type and returns whether a 2 band I+Q dataset should
be mapped onto a single complex band.
Returns BANDERROR for error, STRAIGHT for 1:1 mapping, TWOBANDCOMPLEX for 2 bands -> 1 complex band
*/
typedef enum
{
	BANDERROR,
	STRAIGHT,
	TWOBANDCOMPLEX
}BandMappingRCM;
static BandMappingRCM checkBandFileMappingRCM(GDALDataType dataType, GDALDataset* poBandFile, bool isNITF)
{

	GDALRasterBand* band1 = poBandFile->GetRasterBand(1);
	GDALDataType bandfileType = band1->GetRasterDataType();
	//if there is one band and it has the same datatype, the band file gets passed straight through
	if ( (poBandFile->GetRasterCount() == 1 || poBandFile->GetRasterCount() == 4) && dataType == bandfileType)
		return STRAIGHT;

	//if the band file has 2 bands, they should represent I+Q
	//and be a compatible data type
	if (poBandFile->GetRasterCount() == 2 && GDALDataTypeIsComplex(dataType))
	{
		GDALRasterBand* band2 = poBandFile->GetRasterBand(2);

		if (bandfileType != band2->GetRasterDataType())
			return BANDERROR;   //both bands must be same datatype

								//check compatible types - there are 4 complex types in GDAL
		if ((dataType == GDT_CInt16 && bandfileType == GDT_Int16) ||
			(dataType == GDT_CInt32 && bandfileType == GDT_Int32) ||
			(dataType == GDT_CFloat32 && bandfileType == GDT_Float32) ||
			(dataType == GDT_CFloat64 && bandfileType == GDT_Float64))
			return TWOBANDCOMPLEX;

		if ((dataType == GDT_CInt16 && bandfileType == GDT_CInt16) ||
			(dataType == GDT_CInt32 && bandfileType == GDT_CInt32) ||
			(dataType == GDT_CFloat32 && bandfileType == GDT_CFloat32) ||
			(dataType == GDT_CFloat64 && bandfileType == GDT_CFloat64))
			return TWOBANDCOMPLEX;
	}

	if (isNITF) {
		return STRAIGHT;
	}

	return BANDERROR;   //don't accept any other combinations
}

/************************************************************************/
/*                            RCMRasterBand                             */
/************************************************************************/

RCMRasterBand::RCMRasterBand(RCMDataset *poDSIn, int nBandIn,
								GDALDataType eDataTypeIn,
								const char *pszPole,
								GDALDataset *poBandFile,
								bool bTwoBandComplex, bool isOneFilePerPol, bool isNITF) :
	poBandFile(poBandFile),
	poRCMDataset(poDSIn),
	m_eCalib(eCalibration::Uncalib),
	m_nfTable(NULL),
	m_nTableSize(0),
	m_nfOffset(0),
	m_pszLUTFile(NULL)
{
	poDS = poDSIn;
	this->nBand = nBandIn;
	this->isOneFilePerPol = isOneFilePerPol;
	this->isNITF = isNITF;
	eDataType = eDataTypeIn;
	twoBandComplex = bTwoBandComplex;


	/*Check image type, whether there is one file per polarization or 
	 *one file containing all polarizations*/
	if (this->isOneFilePerPol)
        {
                poBand = poBandFile->GetRasterBand(1);
        }
        else
        {
                poBand = poBandFile->GetRasterBand(this->nBand);
        }

	poBand->GetBlockSize(&nBlockXSize, &nBlockYSize);

	if (pszPole != NULL && strlen(pszPole) != 0) {
		SetMetadataItem("POLARIMETRIC_INTERP", pszPole);
	}
}

double RCMRasterBand::GetLUT(int pixel)
{
	return NAN;
}

int RCMRasterBand::GetLUTsize()
{
	return 0;
}

const char *RCMRasterBand::GetLUTFilename()
{
	return NULL;
}

double RCMRasterBand::GetLUTOffset()
{
	return 0.0f;
}

bool RCMRasterBand::IsComplex()
{
	if (this->m_eType == GDT_CInt16 || this->m_eType == GDT_CInt32 || this->m_eType == GDT_CFloat32 || this->m_eType == GDT_CFloat64) {
		return true;
	}
	else {
		return false;
	}
}

bool RCMRasterBand::IsExistLUT()
{
    return false;
}

eCalibration RCMRasterBand::GetCalibration()
{
	return this->m_eCalib;
}

void RCMRasterBand::SetPartialLUT(int pixel_offset, int pixel_width)
{
	// nothing to do
}

double RCMRasterBand::GetNoiseLevels(int pixel)
{
	return 0.0f;
}

int RCMRasterBand::GetNoiseLevelsSize()
{
	return 0;
}

const char *RCMRasterBand::GetNoiseLevelsFilename()
{
	return NULL;
}

bool RCMRasterBand::IsExistNoiseLevels()
{
	return false;
}


/************************************************************************/
/*                            RCMRasterBand()                            */
/************************************************************************/

RCMRasterBand::~RCMRasterBand()

{
	if (poBandFile != NULL)
		GDALClose(reinterpret_cast<GDALRasterBandH>(poBandFile));
}

/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr RCMRasterBand::IReadBlock( int nBlockXOff, int nBlockYOff,
	                              void * pImage)

{
	int nRequestYSize;
	int nRequestXSize;

	/* -------------------------------------------------------------------- */
	/*      If the last strip is partial, we need to avoid                  */
	/*      over-requesting.  We also need to initialize the extra part     */
	/*      of the block to zero.                                           */
	/* -------------------------------------------------------------------- */
	if ((nBlockYOff + 1) * nBlockYSize > nRasterYSize)
	{
		nRequestYSize = nRasterYSize - nBlockYOff * nBlockYSize;
		memset(pImage, 0, (GDALGetDataTypeSize(eDataType) / 8) *
			nBlockXSize * nBlockYSize);
	}
	else
	{
		nRequestYSize = nBlockYSize;
	}

	/*-------------------------------------------------------------------- */
	/*      If the input imagery is tiled, also need to avoid over-        */
	/*      requesting in the X-direction.                                 */
	/* ------------------------------------------------------------------- */
	if ((nBlockXOff + 1) * nBlockXSize > nRasterXSize)
	{
		nRequestXSize = nRasterXSize - nBlockXOff * nBlockXSize;
		memset(pImage, 0, (GDALGetDataTypeSize(eDataType) / 8) *
			nBlockXSize * nBlockYSize);
	}
	else
	{
		nRequestXSize = nBlockXSize;
	}

	int dataTypeSize = GDALGetDataTypeSizeBytes(eDataType);
	GDALDataType bandFileType = poBandFile->GetRasterBand(1)->GetRasterDataType();
        int bandFileSize = GDALGetDataTypeSizeBytes(bandFileType);

	//case: 2 bands representing I+Q -> one complex band
	if (twoBandComplex && !this->isNITF)
	{
		//int bandFileSize = GDALGetDataTypeSizeBytes(bandFileType);
		//this data type is the complex version of the band file
		// Roberto: don't check that for the moment: CPLAssert(dataTypeSize == bandFileSize * 2);

		return
			//I and Q from each band are pixel-interleaved into this complex band
			poBandFile->RasterIO(GF_Read,
					nBlockXOff * nBlockXSize,
					nBlockYOff * nBlockYSize,
					nRequestXSize, nRequestYSize,
					pImage, nRequestXSize, nRequestYSize,
					bandFileType,
					2,NULL, dataTypeSize, dataTypeSize*nBlockXSize, bandFileSize, NULL);

	}
        else if (twoBandComplex && this->isNITF)
	{
		return
			poBand->RasterIO(GF_Read,
                                nBlockXOff * nBlockXSize,
                                nBlockYOff * nBlockYSize,
                                nRequestXSize, nRequestYSize,
                                pImage, nRequestXSize, nRequestYSize,
                                eDataType, 0,dataTypeSize*nBlockXSize,NULL);
	}
        
	if (poRCMDataset->IsComplexData())
	{
		//this data type is the complex version of the band file
		// Roberto: don't check that for the moment: CPLAssert(dataTypeSize == bandFileSize * 2);
		return
			//I and Q from each band are pixel-interleaved into this complex band
			poBandFile->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pImage, nRequestXSize, nRequestYSize,
				bandFileType,
				2, NULL, dataTypeSize, nBlockXSize * dataTypeSize, bandFileSize, NULL);
	}

	//case: band file == this band
	//NOTE: if the underlying band is opened with the NITF driver, it may combine 2 band I+Q -> complex band
	else if (poBandFile->GetRasterBand(1)->GetRasterDataType() == eDataType)
	{
		return
			poBand->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
                                nBlockYOff * nBlockYSize,
                                nRequestXSize, nRequestYSize,
                                pImage, nRequestXSize, nRequestYSize,
                                eDataType, 0,dataTypeSize*nBlockXSize,NULL);

	}
	else
	{
		CPLAssert(FALSE);
		return CE_Failure;
	}
}


/************************************************************************/
/*                            ReadLUT()                                 */
/************************************************************************/
/* Read the provided LUT in to m_ndTable                                */
/* 1. The gains list spans the range extent covered by all              */
/*    beams(if applicable).                                             */
/* 2. The mapping between the entry of gains                            */
/*    list and the range sample index is : the range sample             */ 
/*    index = gains entry index * stepSize + pixelFirstLutValue,        */ 
/*    where the gains entry index starts with ‘0’.For ScanSAR SLC,      */ 
/*    the range sample index refers to the index on the COPG            */
/************************************************************************/
void RCMCalibRasterBand::ReadLUT() {

	char bandNumber[6];
	//itoa(poDS->GetRasterCount() + 1, bandNumber, 10);
	sprintf(bandNumber, "%s", poDS->GetRasterCount()+1);

	CPLXMLNode *psLUT = CPLParseXMLFile(m_pszLUTFile);

	this->m_nfOffset = CPLAtof(CPLGetXMLValue(psLUT, "=lut.offset", "0.0"));

	this->pixelFirstLutValue = atoi(CPLGetXMLValue(psLUT, "=lut.pixelFirstLutValue", "0"));

	this->stepSize = atoi(CPLGetXMLValue(psLUT, "=lut.stepSize", "0"));

	this->numberOfValues = atoi(CPLGetXMLValue(psLUT, "=lut.numberOfValues", "0"));

	if (this->numberOfValues <= 0) {
		const char msgError[] = "ERROR: The RCM driver does not support the LUT Number Of Values  equal or lower than zero.";
		write_to_file_error(msgError, "");
		CPLError(CE_Failure, CPLE_NotSupported,"%s", msgError);
		return;
	}

	char **papszLUTList = CSLTokenizeString2(CPLGetXMLValue(psLUT, "=lut.gains", ""), " ", CSLT_HONOURSTRINGS);

	if (this->stepSize <= 0) {
		if (this->pixelFirstLutValue <= 0) {
			const char msgError[] = "ERROR: The RCM driver does not support LUT Pixel First Lut Value equal or lower than zero when theproduct is descending.";
			write_to_file_error(msgError, "");
			CPLError(CE_Failure, CPLE_NotSupported, "%s", msgError);
			return;
		}
	}

	/* Get the Pixel Per range */
	this->m_nTableSize = abs(this->stepSize) * abs(this->numberOfValues);

	if (this->m_nTableSize < this->m_poBandDataset->GetRasterXSize()) {
		const char msgError[] = "ERROR: The RCM driver does not support range of LUT gain values lower than the full image pixel range.";
		write_to_file_error(msgError, "");
		CPLError(CE_Failure, CPLE_NotSupported, "%s", msgError);
		return;
	}

	/* Allocate the right LUT size according to the product range pixel */
	this->m_nfTable = InterpolateValues(papszLUTList, this->m_nTableSize, this->stepSize, this->numberOfValues, this->pixelFirstLutValue);

	const size_t nLen = this->m_nTableSize * max_space_for_string; // 32 max + space 
	char *lut_gains = static_cast<char *>(CPLMalloc(nLen));
	memset(lut_gains, 0, nLen);

	for (int i = 0; i < this->m_nTableSize; i++) {
		char lut[max_space_for_string];
		// 6.123004711900930e+04  %e Scientific annotation
		sprintf(lut, "%e ", this->m_nfTable[i]);
		strcat(lut_gains, lut);

	}

#ifdef _TRACE_RCM
	write_to_file("RCM ReadLUT m_pszLUTFile=", m_pszLUTFile);
	write_to_file("   m_nfTable=", lut_gains);
#endif

	poDS->SetMetadataItem(CPLString("LUT_GAINS_").append(bandNumber).c_str(), lut_gains);
	// Can free this because the function SetMetadataItem takes a copy
	CPLFree(lut_gains);

	char snum[256];
	if (this->m_eCalib == eCalibration::Sigma0) {
		poDS->SetMetadataItem(CPLString("LUT_TYPE_").append(bandNumber).c_str(), "SIGMA0");
	}
	else if (this->m_eCalib == eCalibration::Beta0) {
		poDS->SetMetadataItem(CPLString("LUT_TYPE_").append(bandNumber).c_str(), "BETA0");
	}
	else if (this->m_eCalib == eCalibration::Gamma) {
		poDS->SetMetadataItem(CPLString("LUT_TYPE_").append(bandNumber).c_str(), "GAMMA");
	}
	sprintf(snum, "%d", this->m_nTableSize);
	poDS->SetMetadataItem(CPLString("LUT_SIZE_").append(bandNumber).c_str(), snum);
	sprintf(snum, "%f", this->m_nfOffset);
	poDS->SetMetadataItem(CPLString("LUT_OFFSET_").append(bandNumber).c_str(), snum);

	CPLDestroyXMLNode(psLUT);

	CSLDestroy(papszLUTList);
}

/************************************************************************/
/*                            ReadNoiseLevels()                         */
/************************************************************************/
/* Read the provided LUT in to m_nfTableNoiseLevels                     */
/* 1. The gains list spans the range extent covered by all              */
/*    beams(if applicable).                                             */
/* 2. The mapping between the entry of gains                            */
/*    list and the range sample index is : the range sample             */
/*    index = gains entry index * stepSize + pixelFirstLutValue,        */
/*    where the gains entry index starts with ‘0’.For ScanSAR SLC,      */
/*    the range sample index refers to the index on the COPG            */
/************************************************************************/
void RCMCalibRasterBand::ReadNoiseLevels() {

	this->m_nfTableNoiseLevels = NULL;

	if (this->m_pszNoiseLevelsFile == NULL) {
		return;
	}

	char bandNumber[6];
	//itoa(poDS->GetRasterCount() + 1, bandNumber, 10);
	sprintf(bandNumber, "%s", poDS->GetRasterCount()+1);

	CPLXMLNode *psNoiseLevels = CPLParseXMLFile(this->m_pszNoiseLevelsFile);

	// Load Beta Nought, Sigma Nought, Gamma noise levels
	// Loop through all nodes with spaces
	CPLXMLNode *psreferenceNoiseLevelNode =
		CPLGetXMLNode(psNoiseLevels,
			"=noiseLevels");

	CPLXMLNode *psNodeInc;
	for (psNodeInc = psreferenceNoiseLevelNode->psChild; psNodeInc != NULL;
		psNodeInc = psNodeInc->psNext)
	{
		if (EQUAL(psNodeInc->pszValue, "referenceNoiseLevel")) {
			CPLXMLNode *psCalibType =
				CPLGetXMLNode(psNodeInc, "sarCalibrationType");
			CPLXMLNode *psPixelFirstNoiseValue =
				CPLGetXMLNode(psNodeInc, "pixelFirstNoiseValue");
			CPLXMLNode *psStepSize =
				CPLGetXMLNode(psNodeInc, "stepSize");
			CPLXMLNode *psNumberOfValues =
				CPLGetXMLNode(psNodeInc, "numberOfValues");
			CPLXMLNode *psNoiseLevelValues =
				CPLGetXMLNode(psNodeInc, "noiseLevelValues");

			if (psCalibType != NULL && psPixelFirstNoiseValue != NULL &&
				psStepSize != NULL && psNumberOfValues != NULL &&
				psNoiseLevelValues != NULL) {
				const char * calibType = CPLGetXMLValue(psCalibType, "", "");
				this->pixelFirstLutValueNoiseLevels = atoi(CPLGetXMLValue(psPixelFirstNoiseValue, "", "0"));
				this->stepSizeNoiseLevels = atoi(CPLGetXMLValue(psStepSize, "", "0"));
				this->numberOfValuesNoiseLevels = atoi(CPLGetXMLValue(psNumberOfValues, "", "0"));
				const char * noiseLevelValues = CPLGetXMLValue(psNoiseLevelValues, "", "");
				char **papszNoiseLevelList = CSLTokenizeString2(noiseLevelValues, " ", CSLT_HONOURSTRINGS);
				/* Get the Pixel Per range */
				this->m_nTableNoiseLevelsSize = abs(this->stepSizeNoiseLevels) * abs(this->numberOfValuesNoiseLevels);

				if ( (EQUAL(calibType, "Beta Nought") && this->m_eCalib == Beta0) ||
					 (EQUAL(calibType, "Sigma Nought") && this->m_eCalib == Sigma0) ||
					 (EQUAL(calibType, "Gamma") && this->m_eCalib == Gamma) ) {
					/* Allocate the right Noise Levels size according to the product range pixel */
					this->m_nfTableNoiseLevels = InterpolateValues(papszNoiseLevelList, 
																	this->m_nTableNoiseLevelsSize, 
																	this->stepSizeNoiseLevels,
																	this->numberOfValuesNoiseLevels,
																	this->pixelFirstLutValueNoiseLevels);
				}

				CSLDestroy(papszNoiseLevelList);

				if (this->m_nfTableNoiseLevels != NULL) {
					break; // We are done
				}
			}
		}
	}

#ifdef _TRACE_RCM
	if (this->m_nfTableNoiseLevels != NULL) {
		const size_t nLen = this->m_nTableNoiseLevelsSize * max_space_for_string; // 12 max + space + 11 reserved
		char *noise_levels_values = static_cast<char *>(CPLMalloc(nLen));
		memset(noise_levels_values, 0, nLen);

		for (int i = 0; i < this->m_nTableNoiseLevelsSize; i++) {
			char lut[max_space_for_string];
			// 6.123004711900930e+04  %e Scientific annotation
			sprintf(lut, "%e ", this->m_nfTableNoiseLevels[i]);
			strcat(noise_levels_values, lut);
		}
		write_to_file("RCM ReadNoiseLevel m_pszLUTFile=", m_pszNoiseLevelsFile);
		write_to_file("   m_nfTableNoiseLevels=", noise_levels_values);

		CPLFree(noise_levels_values);
	}
#endif

}

/************************************************************************/
/*                        RCMCalibRasterBand()                          */
/************************************************************************/

RCMCalibRasterBand::RCMCalibRasterBand(
	RCMDataset *poDataset, const char *pszPolarization, GDALDataType eType,
	GDALDataset *poBandDataset, eCalibration eCalib,
	const char *pszLUT, const char *pszNoiseLevels, GDALDataType eOriginalType) :
	m_eCalib(eCalib),
	m_poRCMDataset(poDataset),
	m_poBandDataset(poBandDataset),
	m_eType(eType),
	m_eOriginalType(eOriginalType),
	m_nfTable(NULL),
	m_nTableSize(0),
	m_nfOffset(0),
	m_pszLUTFile(VSIStrdup(pszLUT)),
	m_pszNoiseLevelsFile(VSIStrdup(pszNoiseLevels))
{
	this->poDS = poDataset;

	if (pszPolarization != NULL && strlen(pszPolarization) != 0) {
		SetMetadataItem("POLARIMETRIC_INTERP", pszPolarization);
	}

    //this->eDataType = eType;

	if ((eType == GDT_CInt16) || (eType == GDT_CFloat32))
		this->eDataType = GDT_CFloat32;
	else
		this->eDataType = GDT_Float32;

	GDALRasterBand *poRasterBand = poBandDataset->GetRasterBand( 1 );
	poRasterBand->GetBlockSize(&nBlockXSize, &nBlockYSize);

	ReadLUT();
	ReadNoiseLevels();
}

double RCMCalibRasterBand::GetNoiseLevels(int pixel)
{
	return this->m_nfTableNoiseLevels[pixel];
}

int RCMCalibRasterBand::GetNoiseLevelsSize()
{
	return this->m_nTableNoiseLevelsSize;
}

const char *RCMCalibRasterBand::GetNoiseLevelsFilename()
{
	return this->m_pszNoiseLevelsFile;
}

bool RCMCalibRasterBand::IsExistNoiseLevels()
{
	if (this->m_nfTableNoiseLevels == NULL || this->m_pszNoiseLevelsFile == NULL || strlen(this->m_pszNoiseLevelsFile) == 0 || this->m_nTableNoiseLevelsSize == 0) {
		return false;
	}
	else {
		return true;
	}
}

bool RCMCalibRasterBand::IsComplex()
{
	if (this->m_eType == GDT_CInt16 || this->m_eType == GDT_CInt32 || this->m_eType == GDT_CFloat32 || this->m_eType == GDT_CFloat64) {
		return true;
	}
	else {
		return false;
	}
}

eCalibration RCMCalibRasterBand::GetCalibration()
{
	return this->m_eCalib;
}

double RCMCalibRasterBand::GetLUT(int pixel)
{
	return this->m_nfTable[pixel];
}

int RCMCalibRasterBand::GetLUTsize()
{
	return this->m_nTableSize;
}

const char *RCMCalibRasterBand::GetLUTFilename()
{
	return this->m_pszLUTFile;
}

double RCMCalibRasterBand::GetLUTOffset()
{
	return this->m_nfOffset;
}

bool RCMCalibRasterBand::IsExistLUT()
{
	if (this->m_nfTable == NULL || this->m_pszLUTFile == NULL || strlen(this->m_pszLUTFile) == 0 || this->m_nTableSize == 0) {
		return false;
	}
	else {
		return true;
	}
}

void RCMCalibRasterBand::SetPartialLUT(int pixel_offset, int pixel_width)
{
	/* Alway start from 0 */
	if (pixel_offset < 0) {
		pixel_offset = 0;
	}

	if (pixel_offset < this->GetLUTsize()) {
		/* Can only change if the starting pixel in the raster width range */
		if ((pixel_offset + pixel_width) > this->GetLUTsize() - 1) {
			/* Ya but the width is way too large based on the raster width range when beginning from a different offset
			Recalculate the true relative width
			*/
			pixel_width = this->GetLUTsize() - pixel_offset - 1;
		}

		if (pixel_width > 0) {
			/* Prepare a buffer */
			double *nfTableBuffer = static_cast<double *>(CPLMalloc(sizeof(double) * pixel_width));
			memset(nfTableBuffer, 0, sizeof(double) * pixel_width);

			/* Copy a range */
			int j = 0;
			for (int i = pixel_offset; i < (pixel_offset + pixel_width); i++) {
				nfTableBuffer[j++] = this->GetLUT(i);
			}

			const size_t nLen = pixel_width * max_space_for_string; // 12 max + space + 11 reserved
			char *lut_gains = static_cast<char *>(CPLMalloc(nLen));
			memset(lut_gains, 0, nLen);

			for (int i = 0; i < pixel_width; i++) {
				char lut[max_space_for_string];
				// 2.390641e+02 %e Scientific annotation
				sprintf(lut, "%e ", nfTableBuffer[i]);
				strcat(lut_gains, lut);
			}

			char bandNumber[6];
			//itoa(this->GetBand(), bandNumber, 10);
			sprintf(bandNumber, "%s", this->GetBand());

			poDS->SetMetadataItem(CPLString("LUT_GAINS_").append(bandNumber).c_str(), lut_gains);
			// Can free this because the function SetMetadataItem takes a copy
			CPLFree(lut_gains);

			/* and Set new LUT size */
			char snum[256];
			sprintf(snum, "%d", pixel_width);
			poDS->SetMetadataItem(CPLString("LUT_SIZE_").append(bandNumber).c_str(), snum);

			/* Change the internal value now */
			/* Free the old gains value table, not valid anymore */
			CPLFree(this->m_nfTable);

			/* New gains value table size */
			this->m_nTableSize = pixel_width;

			/* Realoocate and have new gain values from the buffer */
			this->m_nfTable = reinterpret_cast<double *>(CPLMalloc(sizeof(double) * this->m_nTableSize));
			memset(this->m_nfTable, 0, sizeof(double) * this->m_nTableSize);
			memcpy(this->m_nfTable, nfTableBuffer, sizeof(double) * this->m_nTableSize);

			/* Free our buffer */
			CPLFree(nfTableBuffer);
		}
	}
}

double * RCMCalibRasterBand::CloneLUT()
{
	double *values = NULL;

	if (this->m_nfTable != NULL) {
		values = reinterpret_cast<double *>(CPLMalloc(sizeof(double) * this->m_nTableSize));
		memcpy(values, this->m_nfTable, sizeof(double) * this->m_nTableSize);
	}

	return values;
}

double * RCMCalibRasterBand::CloneNoiseLevels()
{
	double *values = NULL;

	if (this->m_nfTableNoiseLevels != NULL) {
		values = (double *)malloc(sizeof(double) * this->m_nTableNoiseLevelsSize);
		memcpy(values, this->m_nfTableNoiseLevels, sizeof(double) * this->m_nTableNoiseLevelsSize);
	}

	return values;
}

/************************************************************************/
/*                       ~RCMCalibRasterBand()                          */
/************************************************************************/

RCMCalibRasterBand::~RCMCalibRasterBand() {
	CPLFree(m_nfTable);
	CPLFree(m_nfTableNoiseLevels);
	CPLFree(m_pszLUTFile);
	CPLFree(m_pszNoiseLevelsFile);
	GDALClose(m_poBandDataset);
}

/************************************************************************/
/*                        IReadBlock()                                  */
/************************************************************************/

CPLErr RCMCalibRasterBand::IReadBlock(int nBlockXOff, int nBlockYOff,
	void *pImage )
{
	CPLErr eErr;
	int nRequestYSize;
	int nRequestXSize;

	/* -------------------------------------------------------------------- */
	/*      If the last strip is partial, we need to avoid                  */
	/*      over-requesting.  We also need to initialize the extra part     */
	/*      of the block to zero.                                           */
	/* -------------------------------------------------------------------- */
	if ((nBlockYOff + 1) * nBlockYSize > nRasterYSize)
	{
		nRequestYSize = nRasterYSize - nBlockYOff * nBlockYSize;
		memset(pImage, 0, (GDALGetDataTypeSize(eDataType) / 8) *
			nBlockXSize * nBlockYSize);
	}
	else
	{
		nRequestYSize = nBlockYSize;
	}

	/*-------------------------------------------------------------------- */
	/*      If the input imagery is tiled, also need to avoid over-        */
	/*      requesting in the X-direction.                                 */
	/* ------------------------------------------------------------------- */
	if ((nBlockXOff + 1) * nBlockXSize > nRasterXSize)
	{
		nRequestXSize = nRasterXSize - nBlockXOff * nBlockXSize;
		memset(pImage, 0, (GDALGetDataTypeSize(eDataType) / 8) *
			nBlockXSize * nBlockYSize);
	}
	else
	{
		nRequestXSize = nBlockXSize;
	}

	if (this->m_eOriginalType == GDT_CInt16) {
		GInt16 *pnImageTmp;
		/* read in complex values */
		pnImageTmp = (GInt16 *)CPLMalloc(2 * nBlockXSize * nBlockYSize *
			GDALGetDataTypeSize(GDT_Int16) / 8);

		if (m_poBandDataset->GetRasterCount() == 2) {
			eErr = m_poBandDataset->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pnImageTmp, nRequestXSize, nRequestYSize,
				this->m_eOriginalType,
				2, NULL, 4, nBlockXSize * 4, 4, NULL);

			/*
			eErr = m_poBandDataset->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pnImageTmp, nRequestXSize, nRequestYSize,
				GDT_Int32,
				2, NULL, 4, nBlockXSize * 4, 2, NULL);
			*/
		}
		else {
			eErr =
				m_poBandDataset->RasterIO(GF_Read,
					nBlockXOff * nBlockXSize,
					nBlockYOff * nBlockYSize,
					nRequestXSize, nRequestYSize,
					pnImageTmp, nRequestXSize, nRequestYSize,
					this->m_eOriginalType,
					1, NULL, 4, nBlockXSize * 4, 0, NULL);

			/*
			eErr =
				m_poBandDataset->RasterIO(GF_Read,
					nBlockXOff * nBlockXSize,
					nBlockYOff * nBlockYSize,
					nRequestXSize, nRequestYSize,
					pnImageTmp, nRequestXSize, nRequestYSize,
					GDT_UInt32,
					1, NULL, 4, nBlockXSize * 4, 0, NULL);
			*/

#ifdef CPL_LSB
			/* First, undo the 32bit swap. */
			GDALSwapWords(pImage, 4, nBlockXSize * nBlockYSize, 4);

			/* Then apply 16 bit swap. */
			GDALSwapWords(pImage, 2, nBlockXSize * nBlockYSize * 2, 2);
#endif        
		}

#ifdef _TRACE_RCM
		char msgBlocks[256] = "";
		sprintf(msgBlocks, "IReadBlock: nBlockXOff=%d and nBlockYOff=%d ", nBlockXOff, nBlockYOff);
		write_to_file(msgBlocks, "");
#endif

		/* calibrate the complex values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				/* calculate pixel offset in memory*/
				int nPixOff = (2 * (i * nBlockXSize)) + (j * 2);
				int nTruePixOff = (i * nBlockXSize) + j;

				// Formula for Complex Q+J
				float real = (float)pnImageTmp[nPixOff];
				float img = (float)pnImageTmp[nPixOff + 1];
				float digitalValue = (real * real) + (img * img);
				float lutValue = m_nfTable[nBlockXOff*nBlockXSize + j];
				float calibValue = digitalValue / (lutValue * lutValue);

#ifdef _TRACE_RCM
				if (nBlockXOff*nBlockXSize + j >= 4961) {
					char msgBlocks[2048] = "";
					sprintf(msgBlocks, "IReadBlock: [%d,%d] real=%f img=%f digitalValue=%f lutValue[%d]=%f calibValue=%f", i, j, real, img, digitalValue, (nBlockXOff*nBlockXSize + j), lutValue, calibValue);
					write_to_file(msgBlocks, "");
				}
#endif

				((float *)pImage)[nTruePixOff] = calibValue;
			}
		}

		CPLFree(pnImageTmp);
	}

	//If the underlying file is NITF CFloat32
	else if (this->m_eOriginalType == GDT_CFloat32 || this->m_eOriginalType == GDT_CFloat64)
	{
		/* read in complex values */
		float *pnImageTmp;

		int dataTypeSize = GDALGetDataTypeSize(this->m_eOriginalType) / 8;
		GDALDataType bandFileType = this->m_eOriginalType;
		int bandFileSize = GDALGetDataTypeSize(bandFileType) / 8;

		/* read the original image complex values in a temporary image space */
		pnImageTmp = (float *)CPLMalloc(2 * nBlockXSize * nBlockYSize * bandFileSize);

		eErr =
			//I and Q from each band are pixel-interleaved into this complex band
			m_poBandDataset->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pnImageTmp, nRequestXSize, nRequestYSize,
				bandFileType,
				2, NULL, dataTypeSize, nBlockXSize * dataTypeSize, bandFileSize, NULL);

		/* calibrate the complex values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				/* calculate pixel offset in memory*/
				int nPixOff = (2 * (i * nBlockXSize)) + (j * 2);
				int nTruePixOff = (i * nBlockXSize) + j;

				// Formula for Complex Q+J
				float real = (float)pnImageTmp[nPixOff];
				float img = (float)pnImageTmp[nPixOff + 1];
				float digitalValue = (real * real) + (img * img);
				float lutValue = m_nfTable[nBlockXOff*nBlockXSize + j];
				float calibValue = digitalValue / (lutValue * lutValue);

				((float *)pImage)[nTruePixOff] = calibValue;
			}
		}

		CPLFree(pnImageTmp);
	}

	else if (this->m_eOriginalType == GDT_Float32) {
		/* A Float32 is actual 4 bytes */
		eErr = m_poBandDataset->RasterIO(GF_Read,
			nBlockXOff * nBlockXSize,
			nBlockYOff * nBlockYSize,
			nRequestXSize, nRequestYSize,
			pImage, nRequestXSize, nRequestYSize,
			GDT_Float32,
			1, NULL, 4, nBlockXSize * 4, 0, NULL);

		/* iterate over detected values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				int nPixOff = (i * nBlockXSize) + j;

				/* For detected products, in order to convert the digital number of a given range sample to a calibrated value,
				the digital value is first squared, then the offset(B) is added and the result is divided by the gains value(A)
				corresponding to the range sample. RCM-SP-53-0419  Issue 2/5:  January 2, 2018  Page 7-56 */
				float digitalValue = ((float *)pImage)[nPixOff];
				float A = m_nfTable[nBlockXOff*nBlockXSize + j];
				((float *)pImage)[nPixOff] = ((digitalValue * digitalValue) + this->m_nfOffset) / A;
			}
		}
	}

	else if (this->m_eOriginalType == GDT_Float64) {
		/* A Float64 is actual 8 bytes */
		eErr = m_poBandDataset->RasterIO(GF_Read,
			nBlockXOff * nBlockXSize,
			nBlockYOff * nBlockYSize,
			nRequestXSize, nRequestYSize,
			pImage, nRequestXSize, nRequestYSize,
			GDT_Float64,
			1, NULL, 8, nBlockXSize * 8, 0, NULL);

		/* iterate over detected values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				int nPixOff = (i * nBlockXSize) + j;

				/* For detected products, in order to convert the digital number of a given range sample to a calibrated value,
				the digital value is first squared, then the offset(B) is added and the result is divided by the gains value(A)
				corresponding to the range sample. RCM-SP-53-0419  Issue 2/5:  January 2, 2018  Page 7-56 */
				float digitalValue = ((float *)pImage)[nPixOff];
				float A = m_nfTable[nBlockXOff*nBlockXSize + j];
				((float *)pImage)[nPixOff] = ((digitalValue * digitalValue) + this->m_nfOffset) / A;
			}
		}
	}

	else if (this->m_eOriginalType == GDT_UInt16) {
		GUInt16 *pnImageTmp;
		/* read in detected values */
		pnImageTmp = (GUInt16 *)CPLMalloc(nBlockXSize * nBlockYSize *
			GDALGetDataTypeSize(GDT_UInt16) / 8);
		eErr = m_poBandDataset->RasterIO(GF_Read,
			nBlockXOff * nBlockXSize,
			nBlockYOff * nBlockYSize,
			nRequestXSize, nRequestYSize,
			pnImageTmp, nRequestXSize, nRequestYSize,
			GDT_UInt16,
			1, NULL, 2, nBlockXSize * 2, 0, NULL);

		/* iterate over detected values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				int nPixOff = (i * nBlockXSize) + j;

				float digitalValue = (float)pnImageTmp[nPixOff];
				float A = m_nfTable[nBlockXOff*nBlockXSize + j];
				((float *)pImage)[nPixOff] = ((digitalValue *	digitalValue) +	this->m_nfOffset) / A;
			}
		}
		CPLFree(pnImageTmp);
	} /* Ticket #2104: Support for ScanSAR products */

	else if (this->m_eOriginalType == GDT_Byte) {
		GByte *pnImageTmp;
		pnImageTmp = (GByte *)CPLMalloc(nBlockXSize * nBlockYSize *
			GDALGetDataTypeSize(GDT_Byte) / 8);
		eErr = m_poBandDataset->RasterIO(GF_Read,
			nBlockXOff * nBlockXSize,
			nBlockYOff * nBlockYSize,
			nRequestXSize, nRequestYSize,
			pnImageTmp, nRequestXSize, nRequestYSize,
			GDT_Byte,
			1, NULL, 1, nBlockXSize, 0, NULL);

		/* iterate over detected values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				int nPixOff = (i * nBlockXSize) + j;

				float digitalValue = (float)pnImageTmp[nPixOff];
				float A = m_nfTable[nBlockXOff*nBlockXSize + j];
				((float *)pImage)[nPixOff] = ((digitalValue * digitalValue) + this->m_nfOffset) / A;
			}
		}
		CPLFree(pnImageTmp);
	}
	else {
		CPLAssert(FALSE);
		return CE_Failure;
	}
	return eErr;
}

/************************************************************************/
/* ==================================================================== */
/*                              RCMDataset                              */
/* ==================================================================== */
/************************************************************************/

/************************************************************************/
/*                             RCMDataset()                             */
/************************************************************************/

RCMDataset::RCMDataset() :
	psProduct(NULL),
	nGCPCount(0),
	pasGCPList(NULL),
	pszGCPProjection(CPLStrdup("")),
	papszSubDatasets(NULL),
	pszProjection(CPLStrdup("")),
	pszLutApplied(CPLStrdup("")),
	bHaveGeoTransform(FALSE),
	bPerPolarizationScaling(FALSE),
	papszExtraFiles(NULL),
	m_nfIncidenceAngleTable(NULL),
	m_IncidenceAngleTableSize(0),
	isComplexData(FALSE),
	magnitudeBits(16),
	realBitsComplexData(32),
	imaginaryBitsComplexData(32)
{
	adfGeoTransform[0] = 0.0;
	adfGeoTransform[1] = 1.0;
	adfGeoTransform[2] = 0.0;
	adfGeoTransform[3] = 0.0;
	adfGeoTransform[4] = 0.0;
	adfGeoTransform[5] = 1.0;
}

/************************************************************************/
/*                            ~RCMDataset()                             */
/************************************************************************/

RCMDataset::~RCMDataset()

{
	FlushCache();

	CPLDestroyXMLNode(psProduct);
	CPLFree(pszProjection);
	CPLFree(pszGCPProjection);
	CPLFree(pszLutApplied);

	if (nGCPCount > 0)
	{
		GDALDeinitGCPs(nGCPCount, pasGCPList);
		CPLFree(pasGCPList);
	}

	CloseDependentDatasets();
	
	if (papszSubDatasets != NULL )
		CSLDestroy(papszSubDatasets);

	if (papszExtraFiles != NULL)
		CSLDestroy(papszExtraFiles);

	if (m_nfIncidenceAngleTable != NULL)
		CPLFree(m_nfIncidenceAngleTable);

	psProduct = NULL;
	pszProjection = NULL;
	pszGCPProjection = NULL;
}

/************************************************************************/
/*                      CloseDependentDatasets()                        */
/************************************************************************/

int RCMDataset::CloseDependentDatasets()
{
	int bHasDroppedRef = GDALPamDataset::CloseDependentDatasets();

	if (nBands != 0)
		bHasDroppedRef = TRUE;

	for (int iBand = 0; iBand < nBands; iBand++)
	{
		delete papoBands[iBand];
	}
	nBands = 0;

	return bHasDroppedRef;
}

/************************************************************************/
/*                            GetFileList()                             */
/************************************************************************/

char **RCMDataset::GetFileList()

{
	char **papszFileList = GDALPamDataset::GetFileList();

	papszFileList = CSLInsertStrings(papszFileList, -1, papszExtraFiles);

	return papszFileList;
}

/************************************************************************/
/*                             Identify()                               */
/************************************************************************/

int RCMDataset::Identify(GDALOpenInfo *poOpenInfo)
{
	/* Check for the case where we're trying to read the calibrated data: */
	CPLString calibrationFormat = FormatCalibration(NULL, NULL);

	if (STARTS_WITH_CI(poOpenInfo->pszFilename, calibrationFormat)) {
		return TRUE;
	}


	if (poOpenInfo->bIsDirectory)
	{
		/* Check for directory access when there is a product.xml file in the
		directory. */
		CPLString osMDFilename =
			CPLFormCIFilename(poOpenInfo->pszFilename, "product.xml", NULL);

		VSIStatBufL sStat;
		if (VSIStatL(osMDFilename, &sStat) == 0)
		{
			CPLXMLNode *psProduct = CPLParseXMLFile(osMDFilename);
			if (psProduct == NULL)
				return FALSE;

			CPLXMLNode *psProductAttributes = CPLGetXMLNode(psProduct, "=product");
			if (psProductAttributes == NULL)
			{
				CPLDestroyXMLNode(psProduct);
				return FALSE;
			}

			/* Check the namespace only, should be rs2 */
			const char *szNamespace = CPLGetXMLValue(psProductAttributes, "xmlns", "");

			if (strstr(szNamespace, "rcm") == NULL)
			{
				/* Invalid namespace */
				CPLDestroyXMLNode(psProduct);
				return FALSE;
			}

			CPLDestroyXMLNode(psProduct);
			return TRUE;
		}

		/* If not, check for directory extra 'metadata' access when there is a product.xml file in the
		directory. */

		CPLString osMDFilenameMetadata =
			CPLFormCIFilename(poOpenInfo->pszFilename, GetMetadataProduct(), NULL);

		VSIStatBufL sStatMetadata;
		if (VSIStatL(osMDFilenameMetadata, &sStatMetadata) == 0)
		{
			CPLXMLNode *psProduct = CPLParseXMLFile(osMDFilenameMetadata);
			if (psProduct == NULL)
				return FALSE;

			CPLXMLNode *psProductAttributes = CPLGetXMLNode(psProduct, "=product");
			if (psProductAttributes == NULL)
			{
				CPLDestroyXMLNode(psProduct);
				return FALSE;
			}

			/* Check the namespace only, should be rs2 */
			const char *szNamespace = CPLGetXMLValue(psProductAttributes, "xmlns", "");

			if (strstr(szNamespace, "rcm") == NULL)
			{
				/* Invalid namespace */
				CPLDestroyXMLNode(psProduct);
				return FALSE;
			}

			CPLDestroyXMLNode(psProduct);

			return TRUE;
		}

		return FALSE;
	}

	/* otherwise, do our normal stuff */
	if (strlen(poOpenInfo->pszFilename) < 11
		|| !EQUAL(poOpenInfo->pszFilename + strlen(poOpenInfo->pszFilename) - 11,
			"product.xml"))
		return FALSE;

	if (poOpenInfo->nHeaderBytes < 100)
		return FALSE;

	/* The RCM schema location is rcm_prod_product.xsd */
	if (strstr((const char *)poOpenInfo->pabyHeader, "/rcm") == NULL
		|| strstr((const char *)poOpenInfo->pabyHeader, "<product") == NULL)
		return FALSE;

	return TRUE;
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *RCMDataset::Open(GDALOpenInfo * poOpenInfo)

{
	/* -------------------------------------------------------------------- */
	/*      Is this a RCM Product.xml definition?                           */
	/* -------------------------------------------------------------------- */
	if (!RCMDataset::Identify(poOpenInfo)) {
		return NULL;
	}

	/* -------------------------------------------------------------------- */
	/*        Get subdataset information, if relevant                       */
	/* -------------------------------------------------------------------- */
	const char *pszFilename = poOpenInfo->pszFilename;
	eCalibration eCalib = None;

	CPLString calibrationFormat(FormatCalibration(NULL, NULL));

	if (STARTS_WITH_CI(pszFilename, calibrationFormat)) {
		// The calibration name and filename begins after the hard coded layer name
		pszFilename += strlen(szLayerCalibration) + 1;

		if (STARTS_WITH_CI(pszFilename, szBETA0))
		{
			eCalib = Beta0;
		}
		else if (STARTS_WITH_CI(pszFilename, szSIGMA0)) {
			eCalib = Sigma0;
		}
		else if (STARTS_WITH_CI(pszFilename, szGAMMA) || STARTS_WITH_CI(pszFilename, "GAMMA0")) // Cover both situation
		{
			eCalib = Gamma;
		}
		else if (STARTS_WITH_CI(pszFilename, szUNCALIB))
		{
			eCalib = Uncalib;
		}
		else
		{
			eCalib = None;
		}

		/* advance the pointer to the actual filename */
		while (*pszFilename != '\0' && *pszFilename != ':')
			pszFilename++;

		if (*pszFilename == ':')
			pszFilename++;

		//need to redo the directory check:
		//the GDALOpenInfo check would have failed because of the calibration string on the filename
		VSIStatBufL  sStat;
		if (VSIStatL(pszFilename, &sStat) == 0)
			poOpenInfo->bIsDirectory = VSI_ISDIR(sStat.st_mode);
	}

	CPLString osMDFilename;
	if (poOpenInfo->bIsDirectory)
	{
		/* Check for directory access when there is a product.xml file in the
		directory. */
		osMDFilename =
			CPLFormCIFilename(pszFilename, "product.xml", NULL);

		VSIStatBufL sStat;
		if (VSIStatL(osMDFilename, &sStat) != 0)
		{
			/* If not, check for directory extra 'metadata' access when there is a product.xml file in the
			directory. */
			osMDFilename =
				CPLFormCIFilename(pszFilename, GetMetadataProduct(), NULL);
		}
	}
	else
		osMDFilename = pszFilename;

	/* -------------------------------------------------------------------- */
	/*      Ingest the Product.xml file.                                    */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psProduct = CPLParseXMLFile(osMDFilename);
	if (psProduct == NULL)
		return NULL;

	/* -------------------------------------------------------------------- */
	/*      Confirm the requested access is supported.                      */
	/* -------------------------------------------------------------------- */
	if (poOpenInfo->eAccess == GA_Update)
	{
		const char msgError[] = "ERROR: The RCM driver does not support update access to existing dataset.";
		write_to_file_error(msgError, "");

		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_NotSupported, "%s", msgError);
		return NULL;
	}

	CPLXMLNode *psSceneAttributes = CPLGetXMLNode(psProduct, "=product.sceneAttributes");
	if (psSceneAttributes == NULL)
	{
		const char msgError[] = "ERROR: Failed to find <sceneAttributes> in document.";
		write_to_file_error(msgError, "");
		
		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed, "%s", msgError);
		return NULL;
	}

	CPLXMLNode *psImageAttributes = CPLGetXMLNode(psProduct, "=product.sceneAttributes.imageAttributes");
	if (psImageAttributes == NULL) {
		const char msgError[] = "ERROR: Failed to find <sceneAttributes.imageAttributes> in document.";
		write_to_file_error(msgError, "");

		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed, "%s", msgError);
		return NULL;
	}

	int numberOfEntries =
		atoi(CPLGetXMLValue(psSceneAttributes,
			"numberOfEntries",
			"0"));
	if (numberOfEntries != 1 ) {
		const char msgError[] = "ERROR: Only RCM with Complex Single-beam is supported.";
		write_to_file_error(msgError, "");
		
		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed, "%s", msgError);
		return NULL;
	}

	CPLXMLNode *psImageReferenceAttributes = CPLGetXMLNode(psProduct,
		"=product.imageReferenceAttributes");
	if (psImageReferenceAttributes == NULL) {
		const char msgError[] = "ERROR: Failed to find <imageReferenceAttributes> in document.";
		write_to_file_error(msgError, "");

		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed, "%s", msgError);
		return NULL;
	}

	CPLXMLNode *psImageGenerationParameters = CPLGetXMLNode(psProduct,
		"=product.imageGenerationParameters");
	if (psImageGenerationParameters == NULL) {
		const char msgError[] = "ERROR: Failed to find <imageGenerationParameters> in document.";
		write_to_file_error(msgError, "");

		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed, "%s", msgError);
		return NULL;
	}

	/* -------------------------------------------------------------------- */
	/*      Create the dataset.                                             */
	/* -------------------------------------------------------------------- */
	RCMDataset *poDS = new RCMDataset();

	poDS->psProduct = psProduct;

	/* -------------------------------------------------------------------- */
	/*      Get overall image information.                                  */
	/* -------------------------------------------------------------------- */
	poDS->nRasterXSize =
		atoi(CPLGetXMLValue(psSceneAttributes,
			"imageAttributes.samplesPerLine",
			"-1"));
	poDS->nRasterYSize =
		atoi(CPLGetXMLValue(psSceneAttributes,
			"imageAttributes.numLines",
			"-1"));
	if (poDS->nRasterXSize <= 1 || poDS->nRasterYSize <= 1) {
		const char msgError[] = "ERROR: Non-sane raster dimensions provided in product.xml. If this is "
								"a valid RCM scene, please contact your data provider for "
								"a corrected dataset.";
		write_to_file_error(msgError, "");

		delete poDS;
		CPLError(CE_Failure, CPLE_OpenFailed, "%s",msgError);
		return NULL;
	}

	/* -------------------------------------------------------------------- */
	/*      Check product type, as to determine if there are LUTs for       */
	/*      calibration purposes.                                           */
	/* -------------------------------------------------------------------- */

	const char *pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"generalProcessingInformation.productType", "UNK");
	poDS->SetMetadataItem("PRODUCT_TYPE", pszItem);
	const char *pszProductType = pszItem;

	pszItem = CPLGetXMLValue(psProduct, "=product.productId", "UNK");
	poDS->SetMetadataItem("PRODUCT_ID", pszItem);

	pszItem = CPLGetXMLValue(psProduct, "=product.securityAttributes.securityClassification", "UNK");
	poDS->SetMetadataItem("SECURITY_CLASSIFICATION", pszItem);

	pszItem = CPLGetXMLValue(psProduct, "=product.sourceAttributes.polarizationDataMode", "UNK");
	poDS->SetMetadataItem("POLARIZATION_DATA_MODE", pszItem);

	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"generalProcessingInformation.processingFacility", "UNK");
	poDS->SetMetadataItem("PROCESSING_FACILITY", pszItem);

	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"generalProcessingInformation.processingTime", "UNK");
	poDS->SetMetadataItem("PROCESSING_TIME", pszItem);

	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"sarProcessingInformation.satelliteHeight", "UNK");
	poDS->SetMetadataItem("SATELLITE_HEIGHT", pszItem);

	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"sarProcessingInformation.zeroDopplerTimeFirstLine", "UNK");
	poDS->SetMetadataItem("FIRST_LINE_TIME", pszItem);

	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"sarProcessingInformation.zeroDopplerTimeLastLine", "UNK");
	poDS->SetMetadataItem("LAST_LINE_TIME", pszItem);

	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"sarProcessingInformation.lutApplied", "");
	poDS->SetMetadataItem("LUT_APPLIED", pszItem);
	poDS->pszLutApplied = static_cast<char*>(VSIMalloc(strlen(pszItem) + 2));
	poDS->pszLutApplied[0] = '\0';
	strcpy( poDS->pszLutApplied, pszItem );

	/*--------------------------------------------------------------------- */
	/*     If true, a polarization dependent application LUT has been applied
	/*     for each polarization channel. Otherwise the same application LUT
	/*     has been applied for all polarization channels. Not applicable to
	/*     lookupTable = "Unity*" or if dataType = "Floating-Point".
	/*--------------------------------------------------------------------- */
	pszItem = CPLGetXMLValue(psImageGenerationParameters,
		"sarProcessingInformation.perPolarizationScaling", "false");
	poDS->SetMetadataItem("PER_POLARIZATION_SCALING", pszItem);
	if (EQUAL(pszItem, "true") || EQUAL(pszItem, "TRUE")) {
		poDS->bPerPolarizationScaling = TRUE;
	}

	/* the following cases can be assumed to have no LUTs, as per
	* RN-RP-51-2713, but also common sense
	* SLC represents a SLant range georeferenced Complex product
	* (i.e., equivalent to a Single-Look Complex product for RADARSAT-1 or RADARSAT-2).
	*  • GRD or GRC represent GRound range georeferenced Detected or Complex products (GRD is equivalent to an SGX, SCN or SCW product for RADARSAT1 or RADARSAT-2).
	*  • GCD or GCC represent GeoCoded Detected or Complex products (GCD is equivalent to an SSG or SPG product for RADARSAT-1 or RADARSAT-2).
	*/
	bool bCanCalib = false;
	if (!(STARTS_WITH_CI(pszProductType, "UNK") ||
		STARTS_WITH_CI(pszProductType, "GCD") || 
		STARTS_WITH_CI(pszProductType, "GCC")))  
	{
		bCanCalib = true;
	}

	/* -------------------------------------------------------------------- */
	/*      Get dataType (so we can recognise complex data), and the        */
	/*      bitsPerSample.                                                  */
	/* -------------------------------------------------------------------- */
	const char *pszSampleDataType =
		CPLGetXMLValue(psImageReferenceAttributes, "rasterAttributes.sampleType", 
			"");
	poDS->SetMetadataItem("SAMPLE_TYPE", pszSampleDataType);

	/* Either Integer (16 bits) or Floating-Point (32 bits) */
	const char *pszDataType =
		CPLGetXMLValue(psImageReferenceAttributes, "rasterAttributes.dataType", 
			"");
	poDS->SetMetadataItem("DATA_TYPE", pszDataType);

	/* 2 entries for complex data 1 entry for magnitude detected data */
	const char *pszBitsPerSample = CPLGetXMLValue(psImageReferenceAttributes,
		"rasterAttributes.bitsPerSample", "");
	const int nBitsPerSample =	atoi(pszBitsPerSample);
	poDS->SetMetadataItem("BITS_PER_SAMPLE", pszBitsPerSample);

	const char *pszSampledPixelSpacingTime =
		CPLGetXMLValue(psImageReferenceAttributes, "rasterAttributes.sampledPixelSpacingTime",
			"UNK");
	poDS->SetMetadataItem("SAMPLED_PIXEL_SPACING_TIME", pszSampledPixelSpacingTime);

	const char *pszSampledLineSpacingTime =
		CPLGetXMLValue(psImageReferenceAttributes, "rasterAttributes.sampledLineSpacingTime",
			"UNK");
	poDS->SetMetadataItem("SAMPLED_LINE_SPACING_TIME", pszSampledLineSpacingTime);

	GDALDataType eDataType;
	if (EQUAL(pszSampleDataType, "Complex"))
	{
		poDS->isComplexData = true;
		/* Usually this is the same bits for both */
		poDS->realBitsComplexData = nBitsPerSample;
		poDS->imaginaryBitsComplexData = nBitsPerSample;

		if (nBitsPerSample == 32)
		{
			eDataType = GDT_CFloat32; /* 32 bits, check read block */
		}
		else
		{
			eDataType = GDT_CInt16; /* 16 bits, check read block */
		}
	}
	else if (nBitsPerSample == 32 && EQUAL(pszSampleDataType, "Magnitude Detected"))
	{
		/* Actually, I don't need to test the 'dataType'=' Floating-Point', we know that a 32 bits per sample */
		eDataType = GDT_Float32;
		poDS->isComplexData = false;
		poDS->magnitudeBits = 32;
	}
	else if (nBitsPerSample == 16 && EQUAL(pszSampleDataType, "Magnitude Detected"))
	{
		/* Actually, I don't need to test the 'dataType'=' Integer', we know that a 16 bits per sample */
		eDataType = GDT_UInt16;
		poDS->isComplexData = false;
		poDS->magnitudeBits = 16;
	}
	else
	{
		char msgError[256] = "";
		sprintf(msgError, "ERROR: dataType=%s and bitsPerSample=%d are not a supported configuration.", pszDataType, nBitsPerSample);
		write_to_file_error(msgError, "");

		delete poDS;
		CPLError(CE_Failure, CPLE_AppDefined, "%s", msgError);
		return NULL;
	}

	/* Indicates whether pixel number (i.e., range) increases or decreases with range time.  
	For GCD and GCC products, this applies to intermediate ground range image data prior to geocoding. */
    const char *pszPixelTimeOrdering = CPLGetXMLValue(psImageReferenceAttributes,
			"rasterAttributes.pixelTimeOrdering", "UNK");
	poDS->SetMetadataItem("PIXEL_TIME_ORDERING", pszPixelTimeOrdering);

	/* Indicates whether line numbers (i.e., azimuth) increase or decrease with azimuth time.  
	For GCD and GCC products, this applies to intermediate ground range image data prior to geocoding. */
	const char *pszLineTimeOrdering = CPLGetXMLValue(psImageReferenceAttributes,
		"rasterAttributes.lineTimeOrdering", "UNK");
	poDS->SetMetadataItem("LINE_TIME_ORDERING", pszLineTimeOrdering);

	/* while we're at it, extract the pixel spacing information */
	const char *pszPixelSpacing = CPLGetXMLValue(psImageReferenceAttributes,
		"rasterAttributes.sampledPixelSpacing", "UNK");
	poDS->SetMetadataItem("PIXEL_SPACING", pszPixelSpacing);

	const char *pszLineSpacing = CPLGetXMLValue(psImageReferenceAttributes,
		"rasterAttributes.sampledLineSpacing", "UNK");
	poDS->SetMetadataItem("LINE_SPACING", pszLineSpacing);

	/* -------------------------------------------------------------------- */
	/*      Open each of the data files as a complex band.                  */
	/* -------------------------------------------------------------------- */
	char *pszBeta0LUT = NULL;
	char *pszGammaLUT = NULL;
	char *pszSigma0LUT = NULL;
	// Same file for all calibrations except the calibration is plit inside the XML
	char *pszNoiseLevelsValues = NULL;

	char *pszPath = CPLStrdup(CPLGetPath(osMDFilename));
	const int nFLen = static_cast<int>(osMDFilename.size());

	/* Get a list of all polarizations */
	CPLXMLNode *psSourceAttrs = CPLGetXMLNode(psProduct,
		"=product.sourceAttributes");
	if (psSourceAttrs == NULL)	{
		const char msgError[] = "ERROR: RCM source attributes is missing. Please contact your data provider for a corrected dataset.";
		write_to_file_error(msgError, "");
		
		CPLFree(pszPath);
		delete poDS;
		CPLError(CE_Failure, CPLE_OpenFailed, msgError);
		return NULL;
	}

	CPLXMLNode *psRadarParameters = CPLGetXMLNode(psProduct,
		"=product.sourceAttributes.radarParameters");
	if (psRadarParameters == NULL)	{
		const char msgError[] = "ERROR: RCM radar parameters is missing. Please contact your data provider for a corrected dataset.";
		write_to_file_error(msgError, "");

		CPLFree(pszPath);
		delete poDS;
		CPLError(CE_Failure, CPLE_OpenFailed, msgError);
		return NULL;
	}

	const char *pszPolarizations = CPLGetXMLValue(psRadarParameters,
		"polarizations", "");
	if (pszPolarizations == NULL || strlen(pszPolarizations) == 0)	{
		const char msgError[] = "ERROR: RCM polarizations list is missing. Please contact your data provider for a corrected dataset.";
		write_to_file_error(msgError, "");

		CPLFree(pszPath);
		delete poDS;
		CPLError(CE_Failure, CPLE_OpenFailed, msgError);
		return NULL;
	}
	poDS->SetMetadataItem("POLARIZATIONS", pszPolarizations);

	const char *psAcquisitionType = CPLGetXMLValue(psRadarParameters,
		"acquisitionType", "UNK");
	poDS->SetMetadataItem("ACQUISITION_TYPE", psAcquisitionType);

	const char *psBeams = CPLGetXMLValue(psRadarParameters,
		"beams", "UNK");
	poDS->SetMetadataItem("BEAMS", psBeams);

	char** papszPolarizationsGrids = CSLTokenizeString2(pszPolarizations, " ", 0);
	CPLStringList imageBandList;
	CPLStringList imageBandFileList;
	const int nPolarizationsGridCount = CSLCount(papszPolarizationsGrids);

	/* File names for full resolution IPDFs. For GeoTIFF format, one entry per pole;
	For NITF 2.1 format, only one entry. */
	bool bIsNITF = false;
	const char *pszNITF_Filename;
	int imageBandFileCount = 0;
	int imageBandCount = 0;

	/* Split the polarization string*/
	auto iss = std::istringstream((CPLString(pszPolarizations)).c_str());
        auto pol = std::string{};
	/* Count number of polarizations*/
        while(iss >> pol)
             imageBandCount++;

	CPLXMLNode *psNode = psImageAttributes->psChild;
	for (;
		psNode != NULL;
		psNode = psNode->psNext)
	{
		/* Find the tif or ntf filename */
		if (psNode->eType != CXT_Element
			|| !(EQUAL(psNode->pszValue, "ipdf")))
			continue;

		/* -------------------------------------------------------------------- */
		/*      Fetch ipdf image. COuld be either tif or ntf.                   */
		/*      Replace / by \\                                                 */
		/* -------------------------------------------------------------------- */
		const char *pszBasedFilename = CPLGetXMLValue(psNode, "", "");
		if (*pszBasedFilename == '\0')
			continue;

		/* Count number of image names within ipdf tag*/
		imageBandFileCount++;	

		CPLString pszUpperBasedFilename(CPLString(pszBasedFilename).toupper());

		const bool bEndsWithNTF =
			strlen(pszUpperBasedFilename.c_str()) > 4 &&
			EQUAL(pszUpperBasedFilename.c_str() + strlen(pszUpperBasedFilename.c_str()) - 4, ".NTF");

		if (bEndsWithNTF)
		{
			/* Found it! It would not exist one more */
			bIsNITF = true;
			pszNITF_Filename = pszBasedFilename;
			break;
		}
		else {
			/* Keep adding polarizations filename according to the pole */
			const char *pszPole = CPLGetXMLValue(psNode, "pole", "");
			if (*pszPole == '\0')
				/* Overhere, when no pole has been set but it is NOT a NTF, skip it */
				continue;
			
			imageBandList.AddString(CPLString(pszPole).toupper());
			imageBandFileList.AddString(pszBasedFilename);
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Incidence Angle in a sub-folder                                 */
	/*      called 'calibration' from the 'metadata' folder                 */
	/* -------------------------------------------------------------------- */

	const char *pszIncidenceAngleFileName = CPLGetXMLValue(psImageReferenceAttributes,
		"incidenceAngleFileName", "");

	if (pszIncidenceAngleFileName != NULL) {
		CPLString osIncidenceAnglePath;
		osIncidenceAnglePath.append(CALIBRATION_FOLDER);
		osIncidenceAnglePath.append(szPathSeparator);
		osIncidenceAnglePath.append(pszIncidenceAngleFileName);

		/* Check if the file exist */
		if (IsValidXMLFile(pszPath, osIncidenceAnglePath)) {
			CPLString osIncidenceAngleFilePath = CPLFormFilename(pszPath, osIncidenceAnglePath,
				NULL);

			CPLXMLNode *psIncidenceAngle = CPLParseXMLFile(osIncidenceAngleFilePath);

			int pixelFirstLutValue = atoi(CPLGetXMLValue(psIncidenceAngle, "=incidenceAngles.pixelFirstAnglesValue", "0"));

			int stepSize = atoi(CPLGetXMLValue(psIncidenceAngle, "=incidenceAngles.stepSize", "0"));

			int numberOfValues = atoi(CPLGetXMLValue(psIncidenceAngle, "=incidenceAngles.numberOfValues", "0"));

			/* Get the Pixel Per range */
			int tableSize = abs(stepSize) * abs(numberOfValues);

			CPLString angles;
			// Loop through all nodes with spaces
			CPLXMLNode *psNextNode =
				CPLGetXMLNode(psIncidenceAngle,
					"=incidenceAngles");

			CPLXMLNode *psNodeInc;
			for (psNodeInc = psNextNode->psChild; psNodeInc != NULL;
				psNodeInc = psNodeInc->psNext)
			{
				if (EQUAL(psNodeInc->pszValue, "angles")) {
					if (angles.length() > 0) {
						angles.append(" "); /* separator */
					}
					const char * valAngle = CPLGetXMLValue(psNodeInc, "", "");
					angles.append(valAngle);
				}
			}

			char **papszAngleList = CSLTokenizeString2(angles, " ", CSLT_HONOURSTRINGS);
	
			/* Allocate the right LUT size according to the product range pixel */
			poDS->m_IncidenceAngleTableSize = tableSize;
			poDS->m_nfIncidenceAngleTable = InterpolateValues(papszAngleList, tableSize, stepSize, numberOfValues, pixelFirstLutValue);

			CSLDestroy(papszAngleList);
		}
	}


	for (int iPoleInx=0; iPoleInx<nPolarizationsGridCount; iPoleInx++)
	{
		// Search for a specific band name
		const CPLString pszPole = CPLString(papszPolarizationsGrids[iPoleInx]).toupper();

		double *tableNoiseLevelsBetaNought = NULL;
		double *tableNoiseLevelsSigmaNought = NULL;
		double *tableNoiseLevelsGamma = NULL;

		//Look if the NoiseLevel file xml exist for the 
		CPLXMLNode *psRefNode = psImageReferenceAttributes->psChild;
		for (;
			psRefNode != NULL;
			psRefNode = psRefNode->psNext)
		{
			if (EQUAL(psRefNode->pszValue, "noiseLevelFileName") && bCanCalib) {
				/* Determine which incidence angle correction this is */
				const char *pszPoleToMatch = CPLGetXMLValue(psRefNode,
					"pole", "");
				const char *pszNoiseLevelFile = CPLGetXMLValue(psRefNode, "", "");

				if (*pszPoleToMatch == '\0')
					continue;

				if (EQUAL(pszNoiseLevelFile, ""))
					continue;

				/* -------------------------------------------------------------------- */
				/*      With RCM, LUT file is different per polarizarion (pole)         */
				/*      The following code make sure to loop through all possible       */
				/*      'noiseLevelFileName' and match the <ipdf 'pole' name           */
				/* -------------------------------------------------------------------- */
				if (pszPole.compare(pszPoleToMatch) != 0) {
					continue;
				}

				/* -------------------------------------------------------------------- */
				/*      LUT calibration is unique per pole in a sub-folder              */
				/*      called 'calibration' from the 'metadata' folder                 */
				/* -------------------------------------------------------------------- */
				CPLString oNoiseLevelPath;
				oNoiseLevelPath.append(CALIBRATION_FOLDER);
				oNoiseLevelPath.append(szPathSeparator);
				oNoiseLevelPath.append(pszNoiseLevelFile);

				CPLString osLUTFilePath = CPLFormFilename(pszPath, oNoiseLevelPath,	NULL);

				if (IsValidXMLFile(pszPath, oNoiseLevelPath)) {
				   // File exists for this band
					CPLString osNoiseLevelFilePath = CPLFormFilename(pszPath, oNoiseLevelPath,
						NULL);

					pszNoiseLevelsValues = VSIStrdup(oNoiseLevelPath);
				}
			}
		}

		// Search again with different file
		psRefNode = psImageReferenceAttributes->psChild;
		for (;
			psRefNode != NULL;
			psRefNode = psRefNode->psNext)
		{
			if (EQUAL(psRefNode->pszValue, "lookupTableFileName") && bCanCalib) {
				/* Determine which incidence angle correction this is */
				const char *pszLUTType = CPLGetXMLValue(psRefNode,
					"sarCalibrationType", "");
				const char *pszPoleToMatch = CPLGetXMLValue(psRefNode,
					"pole", "");
				const char *pszLUTFile = CPLGetXMLValue(psRefNode, "", "");

				if (*pszPoleToMatch == '\0')
					continue;

				if (*pszLUTType == '\0')
					continue;

				if (EQUAL(pszLUTType, ""))
					continue;

				/* -------------------------------------------------------------------- */
				/*      With RCM, LUT file is different per polarizarion (pole)         */
				/*      The following code make sure to loop through all possible       */
				/*      'lookupTableFileName' and match the <ipdf 'pole' name           */
				/* -------------------------------------------------------------------- */
				if ( pszPole.compare(pszPoleToMatch) != 0 ) {
					continue;
				}

				/* -------------------------------------------------------------------- */
				/*      LUT calibration is unique per pole in a sub-folder              */
				/*      called 'calibration' from the 'metadata' folder                 */
				/* -------------------------------------------------------------------- */
				CPLString osCalibPath;
				osCalibPath.append(CALIBRATION_FOLDER);
				osCalibPath.append(szPathSeparator);
				osCalibPath.append(pszLUTFile);

				CPLString osLUTFilePath = CPLFormFilename(pszPath, osCalibPath,
					NULL);

				if (EQUAL(pszLUTType, "Beta Nought") &&
					IsValidXMLFile(pszPath, osCalibPath))
				{
					poDS->papszExtraFiles =
						CSLAddString(poDS->papszExtraFiles, osLUTFilePath);

					CPLString pszBuf(FormatCalibration(szBETA0, osMDFilename.c_str()));
					pszBeta0LUT = VSIStrdup(osCalibPath);
					
					const char *oldLut = poDS->GetMetadataItem("BETA_NOUGHT_LUT");
					if (oldLut == NULL) {
						poDS->SetMetadataItem("BETA_NOUGHT_LUT", osCalibPath);
					}
					else {
						/* Keep adding LUT file for all bands, should be planty of space */
						char* ptrConcatLut = static_cast<char*>(VSIMalloc(2048));
						if (ptrConcatLut == NULL) {
							char msgError[256] = "";
							sprintf(msgError, "ERROR: Cannot allocate memory to create BETA LUT metadata for dataType=%s and bitsPerSample=%d", pszDataType, nBitsPerSample);   
							write_to_file_error(msgError, "");

							CPLFree(imageBandList);
							CSLDestroy(papszPolarizationsGrids);
							CPLFree(imageBandFileList);
							CPLFree(pszPath);
							delete poDS;

							CPLError(CE_Failure, CPLE_AppDefined, "%s", msgError);
							return NULL;
						}
						ptrConcatLut[0] = '\0'; /* Just initialize the first byte */
						strcat(ptrConcatLut, oldLut);
						strcat(ptrConcatLut, ",");
						strcat(ptrConcatLut, osCalibPath);
						poDS->SetMetadataItem("BETA_NOUGHT_LUT", ptrConcatLut);
						CPLFree(ptrConcatLut);
					}

					poDS->papszSubDatasets = CSLSetNameValue(
						poDS->papszSubDatasets, "SUBDATASET_3_NAME", pszBuf);
					poDS->papszSubDatasets = CSLSetNameValue(
						poDS->papszSubDatasets, "SUBDATASET_3_DESC",
						"Beta Nought calibrated");

				}
				else if (EQUAL(pszLUTType, "Sigma Nought") &&
					IsValidXMLFile(pszPath, osCalibPath))
				{
					poDS->papszExtraFiles =
						CSLAddString(poDS->papszExtraFiles, osLUTFilePath);

					CPLString pszBuf(FormatCalibration(szSIGMA0, osMDFilename.c_str()));
					pszSigma0LUT = VSIStrdup(osCalibPath);

					const char *oldLut = poDS->GetMetadataItem("SIGMA_NOUGHT_LUT");
					if (oldLut == NULL) {
						poDS->SetMetadataItem("SIGMA_NOUGHT_LUT", osCalibPath);
					}
					else {
						/* Keep adding LUT file for all bands, should be planty of space */
						char* ptrConcatLut = static_cast<char*>(VSIMalloc(2048));
						if (ptrConcatLut == NULL) {
							char msgError[256] = "";
							sprintf(msgError, "ERROR: Cannot allocate memory to create SIGMA0 LUT metadata for dataType=%s and bitsPerSample=%d", pszDataType, nBitsPerSample);
							write_to_file_error(msgError, "");

							CPLFree(imageBandList);
							CSLDestroy(papszPolarizationsGrids);
							CPLFree(imageBandFileList);
							CPLFree(pszPath);
							delete poDS;

							CPLError(CE_Failure, CPLE_AppDefined, "%s", msgError);
							return NULL;
						}
						ptrConcatLut[0] = '\0'; /* Just initialize the first byte */
						strcat(ptrConcatLut, oldLut);
						strcat(ptrConcatLut, ",");
						strcat(ptrConcatLut, osCalibPath);
						poDS->SetMetadataItem("SIGMA_NOUGHT_LUT", ptrConcatLut);
						CPLFree(ptrConcatLut);
					}

					poDS->papszSubDatasets = CSLSetNameValue(
						poDS->papszSubDatasets, "SUBDATASET_2_NAME", pszBuf);
					poDS->papszSubDatasets = CSLSetNameValue(
						poDS->papszSubDatasets, "SUBDATASET_2_DESC",
						"Sigma Nought calibrated");

				}
				else if (EQUAL(pszLUTType, "Gamma") &&
					IsValidXMLFile(pszPath, osCalibPath))
				{
					poDS->papszExtraFiles =
						CSLAddString(poDS->papszExtraFiles, osLUTFilePath);

					CPLString pszBuf(FormatCalibration(szGAMMA, osMDFilename.c_str()));
					pszGammaLUT = VSIStrdup(osCalibPath);
	
					const char *oldLut = poDS->GetMetadataItem("GAMMA_LUT");
					if (oldLut == NULL) {
						poDS->SetMetadataItem("GAMMA_LUT", osCalibPath);
					}
					else {
						/* Keep adding LUT file for all bands, should be planty of space */
						char* ptrConcatLut = static_cast<char*>(VSIMalloc(2048));
						if (ptrConcatLut == NULL) {
							char msgError[256] = "";
							sprintf(msgError, "ERROR: Cannot allocate memory to create GAMMA LUT metadata for dataType=%s and bitsPerSample=%d", pszDataType, nBitsPerSample);
							write_to_file_error(msgError, "");

							CPLFree(imageBandList);
							CSLDestroy(papszPolarizationsGrids);
							CPLFree(imageBandFileList);
							CPLFree(pszPath);
							delete poDS;

							CPLError(CE_Failure, CPLE_AppDefined, "%s", msgError);
							return NULL;
						}
						ptrConcatLut[0] = '\0'; /* Just initialize the first byte */
						strcat(ptrConcatLut, oldLut);
						strcat(ptrConcatLut, ",");
						strcat(ptrConcatLut, osCalibPath);
						poDS->SetMetadataItem("GAMMA_LUT", ptrConcatLut);
						CPLFree(ptrConcatLut);
					}

					poDS->papszSubDatasets = CSLSetNameValue(
						poDS->papszSubDatasets, "SUBDATASET_4_NAME", pszBuf);
					poDS->papszSubDatasets = CSLSetNameValue(
						poDS->papszSubDatasets, "SUBDATASET_4_DESC",
						"Gamma calibrated");

				}
			}
		}


		/* -------------------------------------------------------------------- */
		/*      Fetch ipdf image. Could be either tif or ntf.                   */
		/*      Replace / by \\                                                 */
		/* -------------------------------------------------------------------- */
		const char *pszBasedFilename;
		if (bIsNITF) {
			pszBasedFilename = pszNITF_Filename;
		}
		else {
			const int bandPositionIndex = imageBandList.FindString(pszPole);
			if (bandPositionIndex < 0) {
				char msgError[256] = "";
				sprintf(msgError, "ERROR: RCM cannot find the polarization %s. Please contact your data provider for a corrected dataset", pszPole.c_str());
				write_to_file_error(msgError, "");

				CPLFree(imageBandList);
				CSLDestroy(papszPolarizationsGrids);
				CPLFree(imageBandFileList);
				CPLFree(pszPath);
				delete poDS;

				CPLError(CE_Failure, CPLE_OpenFailed, "%s", msgError);
				return NULL;
			}

			pszBasedFilename = imageBandFileList[bandPositionIndex];
		}

		int nLength = static_cast<int>(strlen(pszBasedFilename));
		char * pszBasename = (char*)CPLCalloc(nLength + 1, sizeof(char));

		/* -------------------------------------------------------------------- */
		/*      Change folder separator if given full path in wrong OS          */
		/* -------------------------------------------------------------------- */
		for (int i = 0; i < nLength; i++) {
			if (pszBasedFilename[i] == cOppositePathSeparator)
				pszBasename[i] = cPathSeparator;
			else
				pszBasename[i] = pszBasedFilename[i];
		}

		/* -------------------------------------------------------------------- */
		/*      Form full filename (path of product.xml + basename).            */
		/* -------------------------------------------------------------------- */
		char *pszFullname =
			CPLStrdup(CPLFormFilename(pszPath, pszBasename, NULL));

		CPLFree(pszBasename);

		/* -------------------------------------------------------------------- */
		/*      Try and open the file.                                          */
		/* -------------------------------------------------------------------- */
		GDALDataset *poBandFile = reinterpret_cast<GDALDataset *>(
			GDALOpen(pszFullname, GA_ReadOnly));
		if (poBandFile == NULL)
		{
			CPLFree(pszFullname);
			continue;
		}
		if (poBandFile->GetRasterCount() == 0)
		{
			GDALClose(reinterpret_cast<GDALRasterBandH>(poBandFile));
			CPLFree(pszFullname);
			continue;
		}

		poDS->papszExtraFiles = CSLAddString(poDS->papszExtraFiles,
			pszFullname);

		/* Some CFloat32 NITF files have nBitsPerSample incorrectly reported    */
		/* as 16, and get misinterpreted as CInt16.  Check the underlying NITF  */
		/* and override if this is the case.                                    */
		if (poBandFile->GetRasterBand(1)->GetRasterDataType() == GDT_CFloat32)
			eDataType = GDT_CFloat32;

		BandMappingRCM b = checkBandFileMappingRCM(eDataType, poBandFile, bIsNITF);
		if (b == BANDERROR)
		{
			GDALClose((GDALRasterBandH)poBandFile);
			CPLFree(pszFullname);
			delete poDS;
			CPLError(CE_Failure, CPLE_AppDefined,
				"The underlying band files do not have an appropriate data type.");
			return NULL;
		}
		bool twoBandComplex = b == TWOBANDCOMPLEX;
		bool isOneFilePerPol = (imageBandCount == imageBandFileCount);	

		/* -------------------------------------------------------------------- */
		/*      Create the band.                                                */
		/* -------------------------------------------------------------------- */     		
		int bandNum = poDS->GetRasterCount()+1;
		if (eCalib == None || eCalib == Uncalib) {
			RCMRasterBand *poBand
				= new RCMRasterBand(poDS, bandNum, eDataType, pszPole, 
									poBandFile, twoBandComplex, isOneFilePerPol, bIsNITF);
			
			poDS->SetBand(poDS->GetRasterCount() + 1, poBand);

		}
		else {
			const char *pszLUT = NULL;
			switch (eCalib) {
			case Sigma0:
				pszLUT = pszSigma0LUT;
				break;
			case Beta0:
				pszLUT = pszBeta0LUT;
				break;
			case Gamma:
				pszLUT = pszGammaLUT;
				break;
			default:
				/* we should bomb gracefully... */
				pszLUT = pszSigma0LUT;
			}

			// The variable 'pszNoiseLevelsValues' is always the same for a ban name except the XML contains different calibration name
			if (poDS->isComplexData) {
				// If Complex, always 32 bits
				RCMCalibRasterBand *poBand
					= new RCMCalibRasterBand(poDS, pszPole, GDT_Float32, poBandFile, eCalib,
						CPLFormFilename(pszPath, pszLUT, NULL),
						CPLFormFilename(pszPath, pszNoiseLevelsValues, NULL), eDataType);
				poDS->SetBand(poDS->GetRasterCount() + 1, poBand);
			}
			else {
				// Whatever the datatype was previoulsy set
				RCMCalibRasterBand *poBand
					= new RCMCalibRasterBand(poDS, pszPole, eDataType, poBandFile, eCalib,
						CPLFormFilename(pszPath, pszLUT, NULL),
						CPLFormFilename(pszPath, pszNoiseLevelsValues, NULL), eDataType);
				poDS->SetBand(poDS->GetRasterCount() + 1, poBand);
			}

		}

		CPLFree(pszFullname);
	}

	if (poDS->papszSubDatasets != NULL && eCalib == None) {
		// must be removed const size_t nBufLen = nFLen + 28;
		CPLString pszBuf = FormatCalibration(szUNCALIB, osMDFilename.c_str());
		poDS->papszSubDatasets = CSLSetNameValue(poDS->papszSubDatasets,
			"SUBDATASET_1_NAME", pszBuf);
		poDS->papszSubDatasets = CSLSetNameValue(poDS->papszSubDatasets,
			"SUBDATASET_1_DESC", "Uncalibrated digital numbers");

	}
	else if (poDS->papszSubDatasets != NULL) {
		CSLDestroy(poDS->papszSubDatasets);
		poDS->papszSubDatasets = NULL;

	}

	/* -------------------------------------------------------------------- */
	/*      Set the appropriate MATRIX_REPRESENTATION.                      */
	/* -------------------------------------------------------------------- */

	if (poDS->GetRasterCount() == 4 && (eDataType == GDT_CInt16 ||
		eDataType == GDT_CFloat32))
	{
		poDS->SetMetadataItem("MATRIX_REPRESENTATION", "SCATTERING");

	}

	/* -------------------------------------------------------------------- */
	/*      Collect a few useful metadata items.                            */
	/* -------------------------------------------------------------------- */


	pszItem = CPLGetXMLValue(psSourceAttrs,
		"satellite", "");
	poDS->SetMetadataItem("SATELLITE_IDENTIFIER", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs,
		"sensor", "");
	poDS->SetMetadataItem("SENSOR_IDENTIFIER", pszItem);

	/* Get beam mode mnemonic */
	pszItem = CPLGetXMLValue(psSourceAttrs, "beamMode", "UNK");
	poDS->SetMetadataItem("BEAM_MODE", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs, "beamModeMnemonic", "UNK");
	poDS->SetMetadataItem("BEAM_MODE_MNEMONIC", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs, "beamModeDefinitionId", "UNK");
	poDS->SetMetadataItem("BEAM_MODE_DEFINITION_ID", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs, "rawDataStartTime", "UNK");
	poDS->SetMetadataItem("ACQUISITION_START_TIME", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs, "inputDatasetFacilityId", "UNK");
	poDS->SetMetadataItem("FACILITY_IDENTIFIER", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs,
		"orbitAndAttitude.orbitInformation.passDirection", "UNK");
	poDS->SetMetadataItem("ORBIT_DIRECTION", pszItem);
	pszItem = CPLGetXMLValue(psSourceAttrs,
		"orbitAndAttitude.orbitInformation.orbitDataSource", "UNK");
	poDS->SetMetadataItem("ORBIT_DATA_SOURCE", pszItem);
	pszItem = CPLGetXMLValue(psSourceAttrs,
		"orbitAndAttitude.orbitInformation.orbitDataFileName", "UNK");
	poDS->SetMetadataItem("ORBIT_DATA_FILE", pszItem);


	/* Get incidence angle information. DONE */
	pszItem = CPLGetXMLValue(psSceneAttributes,
		"imageAttributes.incAngNearRng", "UNK");
	poDS->SetMetadataItem("NEAR_RANGE_INCIDENCE_ANGLE", pszItem);

	pszItem = CPLGetXMLValue(psSceneAttributes,
		"imageAttributes.incAngFarRng", "UNK");
	poDS->SetMetadataItem("FAR_RANGE_INCIDENCE_ANGLE", pszItem);

	pszItem = CPLGetXMLValue(psSceneAttributes,
		"imageAttributes.slantRangeNearEdge", "UNK");
	poDS->SetMetadataItem("SLANT_RANGE_NEAR_EDGE", pszItem);

	pszItem = CPLGetXMLValue(psSceneAttributes,
		"imageAttributes.slantRangeFarEdge", "UNK");
	poDS->SetMetadataItem("SLANT_RANGE_FAR_EDGE", pszItem);

	/*--------------------------------------------------------------------- */
	/*      Collect Map projection/Geotransform information, if present.DONE     */
	/*      In RCM, there is no file that indicates                         */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psMapProjection =
		CPLGetXMLNode(psImageReferenceAttributes,
			"geographicInformation.mapProjection");

	if (psMapProjection != NULL) {
		CPLXMLNode *psPos =
			CPLGetXMLNode(psMapProjection, "positioningInformation");

		pszItem = CPLGetXMLValue(psMapProjection,
			"mapProjectionDescriptor", "UNK");
		poDS->SetMetadataItem("MAP_PROJECTION_DESCRIPTOR", pszItem);

		pszItem = CPLGetXMLValue(psMapProjection,
			"mapProjectionOrientation", "UNK");
		poDS->SetMetadataItem("MAP_PROJECTION_ORIENTATION", pszItem);

		pszItem = CPLGetXMLValue(psMapProjection,
			"resamplingKernel", "UNK");
		poDS->SetMetadataItem("RESAMPLING_KERNEL", pszItem);

		pszItem = CPLGetXMLValue(psMapProjection,
			"satelliteHeading", "UNK");
		poDS->SetMetadataItem("SATELLITE_HEADING", pszItem);

		if (psPos != NULL) {
			const double tl_x = CPLStrtod(CPLGetXMLValue(
				psPos, "upperLeftCorner.mapCoordinate.easting", "0.0"), NULL);
			const double tl_y = CPLStrtod(CPLGetXMLValue(
				psPos, "upperLeftCorner.mapCoordinate.northing", "0.0"), NULL);
			const double bl_x = CPLStrtod(CPLGetXMLValue(
				psPos, "lowerLeftCorner.mapCoordinate.easting", "0.0"), NULL);
			const double bl_y = CPLStrtod(CPLGetXMLValue(
				psPos, "lowerLeftCorner.mapCoordinate.northing", "0.0"), NULL);
			const double tr_x = CPLStrtod(CPLGetXMLValue(
				psPos, "upperRightCorner.mapCoordinate.easting", "0.0"), NULL);
			const double tr_y = CPLStrtod(CPLGetXMLValue(
				psPos, "upperRightCorner.mapCoordinate.northing", "0.0"), NULL);
			poDS->adfGeoTransform[1] = (tr_x - tl_x) / (poDS->nRasterXSize - 1);
			poDS->adfGeoTransform[4] = (tr_y - tl_y) / (poDS->nRasterXSize - 1);
			poDS->adfGeoTransform[2] = (bl_x - tl_x) / (poDS->nRasterYSize - 1);
			poDS->adfGeoTransform[5] = (bl_y - tl_y) / (poDS->nRasterYSize - 1);
			poDS->adfGeoTransform[0] = (tl_x - 0.5*poDS->adfGeoTransform[1]
				- 0.5*poDS->adfGeoTransform[2]);
			poDS->adfGeoTransform[3] = (tl_y - 0.5*poDS->adfGeoTransform[4]
				- 0.5*poDS->adfGeoTransform[5]);

			/* Use bottom right pixel to test geotransform */
			const double br_x = CPLStrtod(CPLGetXMLValue(
				psPos, "lowerRightCorner.mapCoordinate.easting", "0.0"), NULL);
			const double br_y = CPLStrtod(CPLGetXMLValue(
				psPos, "lowerRightCorner.mapCoordinate.northing", "0.0"), NULL);
			const double testx = poDS->adfGeoTransform[0] + poDS->adfGeoTransform[1] *
				(poDS->nRasterXSize - 0.5) + poDS->adfGeoTransform[2] *
				(poDS->nRasterYSize - 0.5);
			const double testy = poDS->adfGeoTransform[3] + poDS->adfGeoTransform[4] *
				(poDS->nRasterXSize - 0.5) + poDS->adfGeoTransform[5] *
				(poDS->nRasterYSize - 0.5);

			/* Give 1/4 pixel numerical error leeway in calculating location
			based on affine transform */
			if ((fabs(testx - br_x) >
				fabs(0.25*(poDS->adfGeoTransform[1] + poDS->adfGeoTransform[2])))
				|| (fabs(testy - br_y) > fabs(0.25*(poDS->adfGeoTransform[4] +
					poDS->adfGeoTransform[5]))))
			{
				const char msgError[] = "WARNING: Unexpected error in calculating affine transform: corner coordinates inconsistent.";
				write_to_file_error(msgError, "");

				CPLError(CE_Warning, CPLE_AppDefined, "%s", msgError);
			}
			else
			{
				poDS->bHaveGeoTransform = TRUE;
			}
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Collect Projection String Information.DONE                      */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psEllipsoid =
		CPLGetXMLNode(psImageReferenceAttributes,
			"geographicInformation.ellipsoidParameters");

	if (psEllipsoid != NULL) {
		OGRSpatialReference oLL, oPrj;

		const char *pszGeodeticTerrainHeight
			= CPLGetXMLValue(psEllipsoid, "geodeticTerrainHeight", "UNK");
		poDS->SetMetadataItem("GEODETIC_TERRAIN_HEIGHT", pszGeodeticTerrainHeight);

		const char *pszEllipsoidName
			= CPLGetXMLValue(psEllipsoid, "ellipsoidName", "");
		double minor_axis
			= CPLAtof(CPLGetXMLValue(psEllipsoid, "semiMinorAxis", "0.0"));
		double major_axis
			= CPLAtof(CPLGetXMLValue(psEllipsoid, "semiMajorAxis", "0.0"));

		if (EQUAL(pszEllipsoidName, "") || (minor_axis == 0.0) ||
			(major_axis == 0.0))
		{
			oLL.SetWellKnownGeogCS("WGS84");
			oPrj.SetWellKnownGeogCS("WGS84");
			const char msgError[] = "WARNING: Incomplete ellipsoid information.  Using wgs-84 parameters.";
			write_to_file_error(msgError, "");

			CPLError(CE_Warning, CPLE_AppDefined, "%s", msgError);
		}
		else if (EQUAL(pszEllipsoidName, "WGS84") || EQUAL(pszEllipsoidName, "WGS 1984")) {
			oLL.SetWellKnownGeogCS("WGS84");
			oPrj.SetWellKnownGeogCS("WGS84");
		}
		else {
			const double inv_flattening = major_axis / (major_axis - minor_axis);
			oLL.SetGeogCS("", "", pszEllipsoidName, major_axis,
				inv_flattening);
			oPrj.SetGeogCS("", "", pszEllipsoidName, major_axis,
				inv_flattening);
		}

		if (psMapProjection != NULL) {
			const char *pszProj = CPLGetXMLValue(
				psMapProjection, "mapProjectionDescriptor", "");
			bool bUseProjInfo = false;

			CPLXMLNode *psUtmParams =
				CPLGetXMLNode(psMapProjection,
					"utmProjectionParameters");

			CPLXMLNode *psNspParams =
				CPLGetXMLNode(psMapProjection,
					"nspProjectionParameters");

			if ((psUtmParams != NULL) && poDS->bHaveGeoTransform) {
				/* double origEasting, origNorthing; */
				bool bNorth = true;

				const int utmZone = atoi(CPLGetXMLValue(psUtmParams, "utmZone", ""));
				const char *pszHemisphere = CPLGetXMLValue(
					psUtmParams, "hemisphere", "");
#if 0
				origEasting = CPLStrtod(CPLGetXMLValue(
					psUtmParams, "mapOriginFalseEasting", "0.0"), NULL);
				origNorthing = CPLStrtod(CPLGetXMLValue(
					psUtmParams, "mapOriginFalseNorthing", "0.0"), NULL);
#endif
				if (STARTS_WITH_CI(pszHemisphere, "southern"))
					bNorth = FALSE;

				if (STARTS_WITH_CI(pszProj, "UTM")) {
					oPrj.SetUTM(utmZone, bNorth);
					bUseProjInfo = true;
				}
			}
			else if ((psNspParams != NULL) && poDS->bHaveGeoTransform) {
				const double origEasting = CPLStrtod(CPLGetXMLValue(
					psNspParams, "mapOriginFalseEasting", "0.0"), NULL);
				const double origNorthing = CPLStrtod(CPLGetXMLValue(
					psNspParams, "mapOriginFalseNorthing", "0.0"), NULL);
				const double copLong = CPLStrtod(CPLGetXMLValue(
					psNspParams, "centerOfProjectionLongitude", "0.0"), NULL);
				const double copLat = CPLStrtod(CPLGetXMLValue(
					psNspParams, "centerOfProjectionLatitude", "0.0"), NULL);
				const double sP1 = CPLStrtod(CPLGetXMLValue(psNspParams,
					"standardParallels1", "0.0"), NULL);
				const double sP2 = CPLStrtod(CPLGetXMLValue(psNspParams,
					"standardParallels2", "0.0"), NULL);

				if (STARTS_WITH_CI(pszProj, "ARC")) {
					/* Albers Conical Equal Area */
					oPrj.SetACEA(sP1, sP2, copLat, copLong, origEasting,
						origNorthing);
					bUseProjInfo = true;
				}
				else if (STARTS_WITH_CI(pszProj, "LCC")) {
					/* Lambert Conformal Conic */
					oPrj.SetLCC(sP1, sP2, copLat, copLong, origEasting,
						origNorthing);
					bUseProjInfo = true;
				}
				else if (STARTS_WITH_CI(pszProj, "STPL")) {
					/* State Plate
					ASSUMPTIONS: "zone" in product.xml matches USGS
					definition as expected by ogr spatial reference; NAD83
					zones (versus NAD27) are assumed. */

					const int nSPZone = atoi(CPLGetXMLValue(psNspParams,
						"zone", "1"));

					oPrj.SetStatePlane(nSPZone, TRUE, NULL, 0.0);
					bUseProjInfo = true;
				}
			}

			if (bUseProjInfo) {
				CPLFree(poDS->pszProjection);
				poDS->pszProjection = NULL;
				oPrj.exportToWkt(&(poDS->pszProjection));
			}
			else {
				const char msgError[] = "WARNING: Unable to interpret projection information; check mapProjection info in product.xml!";
				write_to_file_error(msgError, "");

				CPLError(CE_Warning, CPLE_AppDefined, "%s", msgError);
			}
		}

		CPLFree(poDS->pszGCPProjection);
		poDS->pszGCPProjection = NULL;
		oLL.exportToWkt(&(poDS->pszGCPProjection));

	}

	/* -------------------------------------------------------------------- */
	/*      Collect GCPs.DONE                                                   */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psGeoGrid =
		CPLGetXMLNode(psImageReferenceAttributes,
			"geographicInformation.geolocationGrid");

	if (psGeoGrid != NULL) {
		/* count GCPs */
		poDS->nGCPCount = 0;

		for (psNode = psGeoGrid->psChild; psNode != NULL;
			psNode = psNode->psNext)
		{
			if (EQUAL(psNode->pszValue, "imageTiePoint"))
				poDS->nGCPCount++;
		}

		poDS->pasGCPList = reinterpret_cast<GDAL_GCP *>(
			CPLCalloc(sizeof(GDAL_GCP), poDS->nGCPCount));

		poDS->nGCPCount = 0;

		for (psNode = psGeoGrid->psChild; psNode != NULL;
			psNode = psNode->psNext)
		{
			GDAL_GCP   *psGCP = poDS->pasGCPList + poDS->nGCPCount;

			if (!EQUAL(psNode->pszValue, "imageTiePoint"))
				continue;

			poDS->nGCPCount++;

			char szID[32];
			snprintf(szID, sizeof(szID), "%d", poDS->nGCPCount);
			psGCP->pszId = CPLStrdup(szID);
			psGCP->pszInfo = CPLStrdup("");
			psGCP->dfGCPPixel =
				CPLAtof(CPLGetXMLValue(psNode, "imageCoordinate.pixel", "0"));
			psGCP->dfGCPLine =
				CPLAtof(CPLGetXMLValue(psNode, "imageCoordinate.line", "0"));
			psGCP->dfGCPX =
				CPLAtof(CPLGetXMLValue(psNode, "geodeticCoordinate.longitude", ""));
			psGCP->dfGCPY =
				CPLAtof(CPLGetXMLValue(psNode, "geodeticCoordinate.latitude", ""));
			psGCP->dfGCPZ =
				CPLAtof(CPLGetXMLValue(psNode, "geodeticCoordinate.height", ""));

		}
	}

	CPLFree(pszPath);
	if (pszBeta0LUT) CPLFree(pszBeta0LUT);
	if (pszSigma0LUT) CPLFree(pszSigma0LUT);
	if (pszGammaLUT) CPLFree(pszGammaLUT);

	/* -------------------------------------------------------------------- */
	/*      Collect RPC.DONE                                                   */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psRationalFunctions =
		CPLGetXMLNode(psImageReferenceAttributes,
			"geographicInformation.rationalFunctions");
	if (psRationalFunctions != NULL) {
		char ** papszRPC = NULL;
		static const char* const apszXMLToGDALMapping[] =
		{
			"biasError", "ERR_BIAS",
			"randomError", "ERR_RAND",
			//"lineFitQuality", "????",
			//"pixelFitQuality", "????",
			"lineOffset", "LINE_OFF",
			"pixelOffset", "SAMP_OFF",
			"latitudeOffset", "LAT_OFF",
			"longitudeOffset", "LONG_OFF",
			"heightOffset", "HEIGHT_OFF",
			"lineScale", "LINE_SCALE",
			"pixelScale", "SAMP_SCALE",
			"latitudeScale", "LAT_SCALE",
			"longitudeScale", "LONG_SCALE",
			"heightScale", "HEIGHT_SCALE",
			"lineNumeratorCoefficients", "LINE_NUM_COEFF",
			"lineDenominatorCoefficients", "LINE_DEN_COEFF",
			"pixelNumeratorCoefficients", "SAMP_NUM_COEFF",
			"pixelDenominatorCoefficients", "SAMP_DEN_COEFF",
		};
		for (size_t i = 0; i < CPL_ARRAYSIZE(apszXMLToGDALMapping); i += 2)
		{
			const char* pszValue = CPLGetXMLValue(psRationalFunctions, apszXMLToGDALMapping[i], NULL);
			if (pszValue)
				papszRPC = CSLSetNameValue(papszRPC, apszXMLToGDALMapping[i + 1], pszValue);
		}
		poDS->GDALDataset::SetMetadata(papszRPC, "RPC");
		CSLDestroy(papszRPC);
	}

	/* -------------------------------------------------------------------- */
	/*      Initialize any PAM information.                                 */
	/* -------------------------------------------------------------------- */
	CPLString osDescription;
	CPLString osSubdatasetName;
	bool useSubdatasets = true;

	switch (eCalib) {
	case Sigma0:
	{
		osSubdatasetName = szSIGMA0;
		CPLString pszDescriptionSigma = FormatCalibration(szSIGMA0, osMDFilename.c_str());
		osDescription = pszDescriptionSigma;
	}
	break;
	case Beta0:
	{
		osSubdatasetName = szBETA0;
		CPLString pszDescriptionBeta = FormatCalibration(szBETA0, osMDFilename.c_str());
		osDescription = pszDescriptionBeta;
	}
	break;
	case Gamma:
	{
		osSubdatasetName = szGAMMA;
		CPLString pszDescriptionGamma = FormatCalibration(szGAMMA, osMDFilename.c_str());
		osDescription = pszDescriptionGamma;
	}
	break;
	case Uncalib:
	{
		osSubdatasetName = szUNCALIB;
		CPLString pszDescriptionUncalib = FormatCalibration(szUNCALIB, osMDFilename.c_str());
		osDescription = pszDescriptionUncalib;
	}
	break;
	default:
		osSubdatasetName = szUNCALIB;
		osDescription = osMDFilename;
		useSubdatasets = false;
	}

	if (eCalib != None)
		poDS->papszExtraFiles =
		CSLAddString(poDS->papszExtraFiles, osMDFilename);

	/* -------------------------------------------------------------------- */
	/*      Initialize any PAM information.                                 */
	/* -------------------------------------------------------------------- */
	poDS->SetDescription(osDescription);

	poDS->SetPhysicalFilename(osMDFilename);
	poDS->SetSubdatasetName(osDescription);

	poDS->TryLoadXML();

	/* -------------------------------------------------------------------- */
	/*      Check for overviews.                                            */
	/* -------------------------------------------------------------------- */
	if (useSubdatasets)
		poDS->oOvManager.Initialize(poDS, ":::VIRTUAL:::");
	else
		poDS->oOvManager.Initialize(poDS, osMDFilename);

	return poDS;
}

/************************************************************************/
/*                            GetGCPCount()                             */
/************************************************************************/

int RCMDataset::GetGCPCount()

{
	return nGCPCount;
}

/************************************************************************/
/*                          GetGCPProjection()                          */
/************************************************************************/

const char *RCMDataset::GetGCPProjection()

{
	return pszGCPProjection;
}

/************************************************************************/
/*                               GetGCPs()                              */
/************************************************************************/

const GDAL_GCP *RCMDataset::GetGCPs()

{
	return pasGCPList;
}

/************************************************************************/
/*                          GetProjectionRef()                          */
/************************************************************************/

const char *RCMDataset::GetProjectionRef()

{
	return pszProjection;
}

/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr RCMDataset::GetGeoTransform(double * padfTransform)

{
	memcpy(padfTransform, adfGeoTransform, sizeof(double) * 6);

	if (bHaveGeoTransform)
	{
		return CE_None;
	}
	
	return CE_Failure;
}

/************************************************************************/
/*                      GetMetadataDomainList()                         */
/************************************************************************/

char **RCMDataset::GetMetadataDomainList()
{
	return BuildMetadataDomainList(GDALDataset::GetMetadataDomainList(),
		TRUE,
		"SUBDATASETS", NULL);
}

/************************************************************************/
/*                            GetMetadata()                             */
/************************************************************************/

char **RCMDataset::GetMetadata(const char *pszDomain)

{
	if (pszDomain != NULL && STARTS_WITH_CI(pszDomain, "SUBDATASETS") &&
		papszSubDatasets != NULL)
		return papszSubDatasets;

	return GDALDataset::GetMetadata(pszDomain);
}

/************************************************************************/
/*                         GDALRegister_RCM()                           */
/************************************************************************/

void GDALRegister_RCM()

{
	if (GDALGetDriverByName("RCM") != NULL)
	{
		return;
	}

	GDALDriver *poDriver = new GDALDriver();

	poDriver->SetDescription("RCM");
	poDriver->SetMetadataItem(GDAL_DCAP_RASTER, "YES");
	poDriver->SetMetadataItem(GDAL_DMD_LONGNAME, "Radarsat Constellation Mission XML Product");
	poDriver->SetMetadataItem(GDAL_DMD_HELPTOPIC, "frmt_rcm.html");
	poDriver->SetMetadataItem(GDAL_DMD_SUBDATASETS, "YES");

	poDriver->pfnOpen = RCMDataset::Open;
	poDriver->pfnIdentify = RCMDataset::Identify;

	GetGDALDriverManager()->RegisterDriver(poDriver);
}

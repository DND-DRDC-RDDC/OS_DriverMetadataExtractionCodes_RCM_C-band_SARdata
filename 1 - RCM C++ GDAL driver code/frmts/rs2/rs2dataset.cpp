/******************************************************************************
 *
 * Project:  Polarimetric Workstation
 * Purpose:  Radarsat 2 - XML Products (product.xml) driver
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 2004, Frank Warmerdam <warmerdam@pobox.com>
 * Copyright (c) 2009-2013, Even Rouault <even dot rouault at mines-paris dot org>
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

#include "cpl_minixml.h"
#include "gdal_frmts.h"
#include "gdal_pam.h"
#include "ogr_spatialref.h"
#include "rs2dataset.h"
#include "gdal_io_error.h"
#include "gdal_lut.h"

// MDA Roberto Caron revised based on version 2.2.4 on 2019-02-25
// MDA Shawn Gong ported Roberto's changes to version 2.4.4 on 2020-04-23

CPL_CVSID("$Id: rs2dataset.cpp af7aeb0d4b98257c72c0db4c1bdb1f5880b80af8 2019-06-24 16:49:43 +0300 an-ivanov $")

/*** Function to test for valid LUT files ***/
static bool IsValidXMLFile( const char *pszPath, const char *pszLut)
{
    /* Return true for valid xml file, false otherwise */
    char *pszLutFile
        = VSIStrdup(CPLFormFilename(pszPath, pszLut, nullptr));

    CPLXMLTreeCloser psLut(CPLParseXMLFile(pszLutFile));

    CPLFree(pszLutFile);

    return psLut.get() != nullptr;
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
}BandMappingRS2;
static BandMappingRS2 checkBandFileMappingRS2(GDALDataType dataType, GDALDataset* poBandFile)
{

	GDALRasterBand* band1 = poBandFile->GetRasterBand(1);
	GDALDataType bandfileType = band1->GetRasterDataType();
	//if there is one band and it has the same datatype, the band file gets passed straight through
	if (poBandFile->GetRasterCount() == 1 && dataType == bandfileType)
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
	return BANDERROR;   //don't accept any other combinations
}

/************************************************************************/
/*                            RS2RasterBand                             */
/************************************************************************/

RS2RasterBand::RS2RasterBand( RS2Dataset *poDSIn,
                              GDALDataType eDataTypeIn,
                              const char *pszPole,
                              GDALDataset *poBandFileIn,
							  bool twoBandComplex) :
    poBandFile(poBandFileIn),
	m_eCalib(eCalibration::Uncalib),
	m_nfTable(nullptr),
	m_nTableSize(0),
	m_nfOffset(0),
	m_pszLUTFile(nullptr)
{
    poDS = poDSIn;

    GDALRasterBand *poSrcBand = poBandFile->GetRasterBand( 1 );

    poSrcBand->GetBlockSize( &nBlockXSize, &nBlockYSize );

    eDataType = eDataTypeIn;

	this->twoBandComplex = twoBandComplex;

	if (pszPole != nullptr && strlen(pszPole) != 0) {
		SetMetadataItem("POLARIMETRIC_INTERP", pszPole);
	}
}

/************************************************************************/
/*                            RSRasterBand()                            */
/************************************************************************/

RS2RasterBand::~RS2RasterBand()

{
    if( poBandFile != nullptr )
        GDALClose( reinterpret_cast<GDALRasterBandH>( poBandFile ) );
}


double RS2RasterBand::GetLUT(int pixel)
{
	return NAN;
}

int RS2RasterBand::GetLUTsize()
{
	return 0;
}

const char *RS2RasterBand::GetLUTFilename()
{
	return nullptr;
}

double RS2RasterBand::GetLUTOffset()
{
	return 0.0f;
}

bool RS2RasterBand::IsComplex()
{
	if (this->m_eType == GDT_CInt16 || this->m_eType == GDT_CInt32 || this->m_eType == GDT_CFloat32 || this->m_eType == GDT_CFloat64) {
		return true;
	}
	else {
		return false;
	}
}

bool RS2RasterBand::IsExistLUT()
{
	return false;
}

eCalibration RS2RasterBand::GetCalibration()
{
	return this->m_eCalib;
}

void RS2RasterBand::SetPartialLUT(int pixel_offset, int pixel_width)
{
	// nothing to do
}

double RS2RasterBand::GetNoiseLevels(int pixel)
{
	return 0.0f;
}

int RS2RasterBand::GetNoiseLevelsSize()
{
	return 0;
}

const char *RS2RasterBand::GetNoiseLevelsFilename()
{
	return nullptr;
}

bool RS2RasterBand::IsExistNoiseLevels()
{
	return false;
}

/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr RS2RasterBand::IReadBlock( int nBlockXOff, int nBlockYOff,
                                  void * pImage )

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

	//case: 2 bands representing I+Q -> one complex band
	if (twoBandComplex)
	{
		int dataTypeSize = GDALGetDataTypeSize(eDataType) / 8;
		GDALDataType bandFileType = poBandFile->GetRasterBand(1)->GetRasterDataType();
		int bandFileSize = GDALGetDataTypeSize(bandFileType) / 8;

		//this data type is the complex version of the band file
		CPLAssert(dataTypeSize == bandFileSize * 2);

		return
			//I and Q from each band are pixel-interleaved into this complex band
			poBandFile->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pImage, nRequestXSize, nRequestYSize,
				bandFileType,
				2, nullptr, dataTypeSize, nBlockXSize * dataTypeSize, bandFileSize, nullptr);
	}

	//case: band file == this band
	//NOTE: if the underlying band is opened with the NITF driver, it may combine 2 band I+Q -> complex band
	else if (poBandFile->GetRasterCount() == 1 && poBandFile->GetRasterBand(1)->GetRasterDataType() == eDataType)
	{
		int dataTypeSize = GDALGetDataTypeSize(eDataType) / 8;
		return
			poBandFile->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pImage, nRequestXSize, nRequestYSize,
				eDataType,
				1, nullptr, 0, dataTypeSize*nBlockXSize, 0, nullptr);
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
/************************************************************************/
void RS2CalibRasterBand::ReadLUT() {
    CPLXMLNode *psLUT = CPLParseXMLFile(m_pszLUTFile);

	char bandNumber[6];
	//itoa(poDS->GetRasterCount() + 1, bandNumber, 10);
        sprintf(bandNumber, "%d", poDS->GetRasterCount() + 1);
    this->m_nfOffset = CPLAtof(CPLGetXMLValue( psLUT, "=lut.offset", "0.0" ));

    char **papszLUTList = CSLTokenizeString2( CPLGetXMLValue(psLUT,
        "=lut.gains", ""), " ", CSLT_HONOURSTRINGS);

	this->m_nTableSize = CSLCount(papszLUTList);

	this->m_nfTable = reinterpret_cast<double *>(
        CPLMalloc( sizeof(double) * this->m_nTableSize ) );
	memset(this->m_nfTable, 0, sizeof(double) * this->m_nTableSize);

	const size_t nLen = this->m_nTableSize * max_space_for_string; // 12 max + space + 11 reserved
	char *lut_gains = static_cast<char *>(CPLMalloc(nLen));
	memset(lut_gains, 0, nLen);

    for (int i = 0; i < this->m_nTableSize; i++) {
		this->m_nfTable[i] = CPLAtof(papszLUTList[i]);
		char lut[max_space_for_string];
		// 2.390641e+02 %e Scientific annotation
		sprintf(lut, "%e ", this->m_nfTable[i]);
		strcat(lut_gains, lut);
    }

#ifdef _TRACE_RCM
	write_to_file("RS2 ReadLUT m_pszLUTFile=", m_pszLUTFile);
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
/*                        RS2CalibRasterBand()                          */
/************************************************************************/

RS2CalibRasterBand::RS2CalibRasterBand(
    RS2Dataset *poDataset, const char *pszPolarization, GDALDataType eType,
    GDALDataset *poBandDataset, eCalibration eCalib,
    const char *pszLUT, struct NoiseLevel *noiseLevel, GDALDataType eOriginalType) :
    m_eCalib(eCalib),
    m_poBandDataset(poBandDataset),
    m_eType(eType),
	m_eOriginalType(eOriginalType),
    m_nfTable(nullptr),
    m_nTableSize(0),
    m_nfOffset(0),
	m_nfTableNoiseLevels(nullptr),
    pixelFirstLutValueNoiseLevels(0),
    stepSizeNoiseLevels(0),
    numberOfValuesNoiseLevels(0),
    m_nTableNoiseLevelsSize(0),
    m_pszLUTFile(VSIStrdup(pszLUT))
{
    this->poDS = poDataset;

	if (pszPolarization != nullptr && strlen(pszPolarization) != 0) {
        SetMetadataItem( "POLARIMETRIC_INTERP", pszPolarization );
    }

	if ((eType == GDT_CInt16) || (eType == GDT_CFloat32))
		this->eDataType = GDT_CFloat32;
	else
		this->eDataType = GDT_Float32;

    GDALRasterBand *poRasterBand = poBandDataset->GetRasterBand( 1 );
    poRasterBand->GetBlockSize( &nBlockXSize, &nBlockYSize );

    ReadLUT();

	if (noiseLevel != nullptr) {
		this->pixelFirstLutValueNoiseLevels = noiseLevel->pixelFirstLutValueNoiseLevels;
		this->stepSizeNoiseLevels = noiseLevel->stepSizeNoiseLevels;
		this->numberOfValuesNoiseLevels = noiseLevel->numberOfValuesNoiseLevels;
		this->m_nTableNoiseLevelsSize = noiseLevel->m_nTableNoiseLevelsSize;
		/* Copy values */
		this->m_nfTableNoiseLevels = (double *)malloc(sizeof(double) * m_nTableNoiseLevelsSize);
		memcpy(this->m_nfTableNoiseLevels, noiseLevel->m_nfTableNoiseLevels, sizeof(double) * m_nTableNoiseLevelsSize);
	}
}

double RS2CalibRasterBand::GetLUT(int pixel)
{
	return this->m_nfTable[pixel];
}

int RS2CalibRasterBand::GetLUTsize()
{
	return this->m_nTableSize;
}

const char *RS2CalibRasterBand::GetLUTFilename()
{
	return this->m_pszLUTFile;
}

double RS2CalibRasterBand::GetLUTOffset()
{
	return this->m_nfOffset;
}

bool RS2CalibRasterBand::IsComplex()
{
	if (this->m_eType == GDT_CInt16 || this->m_eType == GDT_CInt32 || this->m_eType == GDT_CFloat32 || this->m_eType == GDT_CFloat64) {
		return true;
	}
	else {
		return false;
	}
}

bool RS2CalibRasterBand::IsExistLUT()
{
	if (this->m_nfTable == nullptr || this->m_pszLUTFile == nullptr || strlen(this->m_pszLUTFile) == 0 || this->m_nTableSize == 0) {
		return false;
	}
	else {
		return true;
	}
}

eCalibration RS2CalibRasterBand::GetCalibration()
{
	return this->m_eCalib;
}

void RS2CalibRasterBand::SetPartialLUT(int pixel_offset, int pixel_width)
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

			const size_t nLen = pixel_width * max_space_for_string; // 32 max + space 
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
			sprintf(bandNumber, "%d", this->GetBand());

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

double * RS2CalibRasterBand::CloneLUT()
{
	double *values = nullptr;

	if (this->m_nfTable != nullptr) {
		values = reinterpret_cast<double *>(CPLMalloc(sizeof(double) * this->m_nTableSize));
		memcpy(values, this->m_nfTable, sizeof(double) * this->m_nTableSize);
	}

	return values;
}

double * RS2CalibRasterBand::CloneNoiseLevels()
{
	double *values = nullptr;

	if (this->m_nfTableNoiseLevels != nullptr) {
		values = (double *)malloc(sizeof(double) * m_nTableNoiseLevelsSize);
		memcpy(values, this->m_nfTableNoiseLevels, sizeof(double) * m_nTableNoiseLevelsSize);
	}

	return values;
}

double RS2CalibRasterBand::GetNoiseLevels(int pixel)
{
	return this->m_nfTableNoiseLevels[pixel];
}

int RS2CalibRasterBand::GetNoiseLevelsSize()
{
	return this->m_nTableNoiseLevelsSize;
}

const char *RS2CalibRasterBand::GetNoiseLevelsFilename()
{
	return "product.xml"; // The file it-self
}

bool RS2CalibRasterBand::IsExistNoiseLevels()
{
	if (this->m_nfTableNoiseLevels == nullptr || this->m_nTableNoiseLevelsSize == 0) {
		return false;
	}
	else {
		return true;
	}
}

/************************************************************************/
/*                       ~RS2CalibRasterBand()                          */
/************************************************************************/

RS2CalibRasterBand::~RS2CalibRasterBand() {
    CPLFree(m_nfTable);
    CPLFree(m_pszLUTFile);
	CPLFree(m_nfTableNoiseLevels);
    GDALClose( m_poBandDataset );
}



/************************************************************************/
/*                        IReadBlock()                                  */
/************************************************************************/

CPLErr RS2CalibRasterBand::IReadBlock( int nBlockXOff, int nBlockYOff,
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

	if (this->m_eOriginalType == GDT_CInt16) { // if (this->m_eType == GDT_CInt16) {
		GInt16 *pnImageTmp;
		/* read the original image complex values in a temporary image space */
		pnImageTmp = (GInt16 *)CPLMalloc(2 * nBlockXSize * nBlockYSize *
			GDALGetDataTypeSize(GDT_Int16) / 8);
		
		// Roberto: Any Complex will be considered as float 32 
		if (m_poBandDataset->GetRasterCount() == 2) {
			eErr = m_poBandDataset->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pnImageTmp, nRequestXSize, nRequestYSize,
				this->m_eOriginalType,
				2, nullptr, 4, nBlockXSize * 4, 4, nullptr);

			/*
			eErr = m_poBandDataset->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pnImageTmp, nRequestXSize, nRequestYSize,
				GDT_Int16,
				2, nullptr, 4, nBlockXSize * 4, 2, nullptr);
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
					1, nullptr, 4, nBlockXSize * 4, 0, nullptr);
			/*
			eErr =
				m_poBandDataset->RasterIO(GF_Read,
					nBlockXOff * nBlockXSize,
					nBlockYOff * nBlockYSize,
					nRequestXSize, nRequestYSize,
					pnImageTmp, nRequestXSize, nRequestYSize,
					GDT_UInt32,
					1, nullptr, 4, nBlockXSize * 4, 0, nullptr);
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

				((float *)pImage)[nTruePixOff] = calibValue;
			}
		}
		
		CPLFree(pnImageTmp);
	}

	//If the underlying file is NITF CFloat32
	else if (this->m_eOriginalType == GDT_CFloat32 && m_poBandDataset->GetRasterCount() == 1) //else if (this->m_eType == GDT_CFloat32 && m_poBandDataset->GetRasterCount() == 1)
	{
		float *pnImageTmp;

		int dataTypeSize = GDALGetDataTypeSize(this->m_eOriginalType) / 8;
		GDALDataType bandFileType = this->m_eOriginalType;
		int bandFileSize = GDALGetDataTypeSize(bandFileType) / 8;

		/* read the original image complex values in a temporary image space */
		pnImageTmp = (float *)CPLMalloc(2 * nBlockXSize * nBlockYSize * bandFileSize);

		eErr =
			m_poBandDataset->RasterIO(GF_Read,
				nBlockXOff * nBlockXSize,
				nBlockYOff * nBlockYSize,
				nRequestXSize, nRequestYSize,
				pnImageTmp, nRequestXSize, nRequestYSize,
				bandFileType,
				1, nullptr, 8, nBlockXSize * 8, 0, nullptr);

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
	
	else if (this->m_eOriginalType == GDT_UInt16) { // else if (this->m_eType == GDT_UInt16) {
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
			1, nullptr, 2, nBlockXSize * 2, 0, nullptr);

		/* iterate over detected values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				int nPixOff = (i * nBlockXSize) + j;

				// Formula for Non-Complex
				float digitalValue = (float)pnImageTmp[nPixOff];
				float A = m_nfTable[nBlockXOff*nBlockXSize + j];
				((float *)pImage)[nPixOff] = ((digitalValue * digitalValue) + this->m_nfOffset) / A;
			}
		}

		CPLFree(pnImageTmp);
	} 
	
	/* Ticket #2104: Support for ScanSAR products */
	else if (this->m_eOriginalType == GDT_Byte) { // else if (this->m_eType == GDT_Byte) {
		GByte *pnImageTmp;
		
		pnImageTmp = (GByte *)CPLMalloc(nBlockXSize * nBlockYSize *
			GDALGetDataTypeSize(GDT_Byte) / 8);

		eErr = m_poBandDataset->RasterIO(GF_Read,
			nBlockXOff * nBlockXSize,
			nBlockYOff * nBlockYSize,
			nRequestXSize, nRequestYSize,
			pnImageTmp, nRequestXSize, nRequestYSize,
			GDT_Byte,
			1, nullptr, 1, nBlockXSize, 0, nullptr);

		/* iterate over detected values */
		for (int i = 0; i < nRequestYSize; i++) {
			for (int j = 0; j < nRequestXSize; j++) {
				int nPixOff = (i * nBlockXSize) + j;

				// Formula for Non-Complex
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
/*                              RS2Dataset                              */
/* ==================================================================== */
/************************************************************************/

/************************************************************************/
/*                             RS2Dataset()                             */
/************************************************************************/

RS2Dataset::RS2Dataset() :
    psProduct(nullptr),
    nGCPCount(0),
    pasGCPList(nullptr),
    pszGCPProjection(CPLStrdup("")),
    papszSubDatasets(nullptr),
    pszProjection(CPLStrdup("")),
    bHaveGeoTransform(FALSE),
    papszExtraFiles(nullptr)
{
    adfGeoTransform[0] = 0.0;
    adfGeoTransform[1] = 1.0;
    adfGeoTransform[2] = 0.0;
    adfGeoTransform[3] = 0.0;
    adfGeoTransform[4] = 0.0;
    adfGeoTransform[5] = 1.0;
}

/************************************************************************/
/*                            ~RS2Dataset()                             */
/************************************************************************/

RS2Dataset::~RS2Dataset()

{
    RS2Dataset::FlushCache();

    CPLDestroyXMLNode( psProduct );
    CPLFree( pszProjection );

    CPLFree( pszGCPProjection );
    if( nGCPCount > 0 )
    {
        GDALDeinitGCPs( nGCPCount, pasGCPList );
        CPLFree( pasGCPList );
    }

    RS2Dataset::CloseDependentDatasets();

	/* Roberto's Fix */
	if (papszSubDatasets != nullptr)
		CSLDestroy( papszSubDatasets );

	if (papszExtraFiles != nullptr)
		CSLDestroy( papszExtraFiles );

	psProduct = nullptr;
	pszProjection = nullptr;
	pszGCPProjection = nullptr;
}

/************************************************************************/
/*                      CloseDependentDatasets()                        */
/************************************************************************/

int RS2Dataset::CloseDependentDatasets()
{
    int bHasDroppedRef = GDALPamDataset::CloseDependentDatasets();

    if (nBands != 0)
        bHasDroppedRef = TRUE;

    for( int iBand = 0; iBand < nBands; iBand++ )
    {
       delete papoBands[iBand];
    }
    nBands = 0;

    return bHasDroppedRef;
}

/************************************************************************/
/*                            GetFileList()                             */
/************************************************************************/

char **RS2Dataset::GetFileList()

{
    char **papszFileList = GDALPamDataset::GetFileList();

    papszFileList = CSLInsertStrings( papszFileList, -1, papszExtraFiles );

    return papszFileList;
}

/************************************************************************/
/*                             Identify()                               */
/************************************************************************/
/* Roberto's Fix */
int RS2Dataset::Identify(GDALOpenInfo *poOpenInfo)
{
	/* Check for the case where we're trying to read the calibrated data: */
	if (STARTS_WITH_CI(poOpenInfo->pszFilename, "RADARSAT_2_CALIB:")) {
		return TRUE;
	}

	/* Check for directory access when there is a product.xml file in the
	directory. */
	if (poOpenInfo->bIsDirectory)
	{
		CPLString osMDFilename =
			CPLFormCIFilename(poOpenInfo->pszFilename, "product.xml", nullptr);

		VSIStatBufL sStat;
		if (VSIStatL(osMDFilename, &sStat) == 0)
		{
			/* Roberto's Fix : Must skip if product.xml is still not what we expect to be a RS2 */
			CPLXMLNode *psProduct = CPLParseXMLFile(osMDFilename);
			if (psProduct == nullptr)
				return FALSE;

			CPLXMLNode *psProductAttributes = CPLGetXMLNode(psProduct, "=product");
			if (psProductAttributes == nullptr)
			{
				CPLDestroyXMLNode(psProduct);
				return FALSE;
			}

			/* Check the namespace only, should be rs2 */
			const char *szNamespace = CPLGetXMLValue(psProductAttributes, "xmlns", "");

			if (strstr(szNamespace, "rs2") == nullptr)
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

	if (strstr((const char *)poOpenInfo->pabyHeader, "/rs2") == nullptr
		|| strstr((const char *)poOpenInfo->pabyHeader, "<product") == nullptr)
		return FALSE;

	return TRUE;
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *RS2Dataset::Open( GDALOpenInfo * poOpenInfo )

{
	/* -------------------------------------------------------------------- */
	/*      Is this a RADARSAT-2 Product.xml definition?                   */
	/* -------------------------------------------------------------------- */
	if (!RS2Dataset::Identify(poOpenInfo)) {
		return nullptr;
	}

	/* -------------------------------------------------------------------- */
	/*        Get subdataset information, if relevant                            */
	/* -------------------------------------------------------------------- */
	CPLString osMDFilename;
	const char *pszFilename = poOpenInfo->pszFilename;
	eCalibration eCalib = None;

	if (EQUALN("RADARSAT_2_CALIB:", pszFilename, 17)) {
		pszFilename += 17;

		if (EQUALN("BETA0", pszFilename, 5))
			eCalib = Beta0;
		else if (EQUALN("SIGMA0", pszFilename, 6))
			eCalib = Sigma0;
		else if (EQUALN("GAMMA", pszFilename, 5))
			eCalib = Gamma;
		else if (EQUALN("UNCALIB", pszFilename, 7))
			eCalib = Uncalib;
		else
			eCalib = None;

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

	if (poOpenInfo->bIsDirectory)
	{
		osMDFilename =
			CPLFormCIFilename(pszFilename, "product.xml", nullptr);
	}
	else
		osMDFilename = pszFilename;

	/* -------------------------------------------------------------------- */
	/*      Ingest the Product.xml file.                                    */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psProduct, *psImageAttributes, *psImageGenerationParameters;

	psProduct = CPLParseXMLFile(osMDFilename);
	if (psProduct == nullptr)
		return nullptr;

	/* -------------------------------------------------------------------- */
	/*      Confirm the requested access is supported.                      */
	/* -------------------------------------------------------------------- */
	if (poOpenInfo->eAccess == GA_Update)
	{
		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_NotSupported,
			"The RS2 driver does not support update access to existing"
			" datasets.\n");
		return nullptr;
	}

	psImageAttributes = CPLGetXMLNode(psProduct, "=product.imageAttributes");
	if (psImageAttributes == nullptr)
	{
		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed,
			"Failed to find <imageAttributes> in document.");
		return nullptr;
	}

	psImageGenerationParameters = CPLGetXMLNode(psProduct,
		"=product.imageGenerationParameters");
	if (psImageGenerationParameters == nullptr) {
		CPLDestroyXMLNode(psProduct);
		CPLError(CE_Failure, CPLE_OpenFailed,
			"Failed to find <imageGenerationParameters> in document.");
		return nullptr;
	}

	/* -------------------------------------------------------------------- */
	/*      Create the dataset.                                             */
	/* -------------------------------------------------------------------- */
	RS2Dataset *poDS = new RS2Dataset();

	poDS->psProduct = psProduct;

	/* -------------------------------------------------------------------- */
	/*      Get overall image information.                                  */
	/* -------------------------------------------------------------------- */
	poDS->nRasterXSize =
		atoi(CPLGetXMLValue(psImageAttributes,
			"rasterAttributes.numberOfSamplesPerLine",
			"-1"));
	poDS->nRasterYSize =
		atoi(CPLGetXMLValue(psImageAttributes,
			"rasterAttributes.numberofLines",
			"-1"));
	if (poDS->nRasterXSize <= 1 || poDS->nRasterYSize <= 1) {
		CPLError(CE_Failure, CPLE_OpenFailed,
			"Non-sane raster dimensions provided in product.xml. If this is "
			"a valid RADARSAT-2 scene, please contact your data provider for "
			"a corrected dataset.");
		delete poDS;
		return nullptr;
	}

	/* -------------------------------------------------------------------- */
	/*      Check product type, as to determine if there are LUTs for       */
	/*      calibration purposes.                                           */
	/* -------------------------------------------------------------------- */

	int bCanCalib = 0;

	const char *pszProductType = CPLGetXMLValue(psImageGenerationParameters,
		"generalProcessingInformation.productType",
		"UNK");

	poDS->SetMetadataItem("PRODUCT_TYPE", pszProductType);

	/* the following cases can be assumed to have no LUTs, as per
	* RN-RP-51-2713, but also common sense
	*/
	if (!(EQUALN(pszProductType, "UNK", 3) ||
		EQUALN(pszProductType, "SSG", 3) ||
		EQUALN(pszProductType, "SPG", 3)))
	{
		bCanCalib = 1;
	}

	/* -------------------------------------------------------------------- */
	/*      Get dataType (so we can recognise complex data), and the        */
	/*      bitsPerSample.                                                  */
	/* -------------------------------------------------------------------- */
	GDALDataType eDataType;

	const char *pszDataType =
		CPLGetXMLValue(psImageAttributes, "rasterAttributes.dataType",
			"");
	poDS->SetMetadataItem("DATA_TYPE", pszDataType);

	const char *pszBitsPerSample = CPLGetXMLValue(psImageAttributes,
		"rasterAttributes.bitsPerSample", "");

	int nBitsPerSample = atoi(pszBitsPerSample);
	poDS->SetMetadataItem("BITS_PER_SAMPLE", pszBitsPerSample);

	bool isComplexType = false;
	if (nBitsPerSample == 16 && EQUAL(pszDataType, "Complex")) {
		eDataType = GDT_CInt16;
		isComplexType = true;
	}
	else if (nBitsPerSample == 32 && EQUAL(pszDataType, "Complex")) { //NITF datasets can come in this configuration
		eDataType = GDT_CFloat32;
		isComplexType = true;
	}
	else if (nBitsPerSample == 16 && EQUALN(pszDataType, "Mag", 3))
		eDataType = GDT_UInt16;
	else if (nBitsPerSample == 8 && EQUALN(pszDataType, "Mag", 3))
		eDataType = GDT_Byte;
	else
	{
		delete poDS;
		CPLError(CE_Failure, CPLE_AppDefined,
			"dataType=%s, bitsPerSample=%d: not a supported configuration.",
			pszDataType, nBitsPerSample);
		return nullptr;
	}

	/* while we're at it, extract the pixel spacing information */
	const char *pszPixelSpacing = CPLGetXMLValue(psImageAttributes,
		"rasterAttributes.sampledPixelSpacing", "UNK");
	poDS->SetMetadataItem("PIXEL_SPACING", pszPixelSpacing);

	const char *pszLineSpacing = CPLGetXMLValue(psImageAttributes,
		"rasterAttributes.sampledLineSpacing", "UNK");
	poDS->SetMetadataItem("LINE_SPACING", pszLineSpacing);

	/* -------------------------------------------------------------------- */
	/*      Load Reference Noise Level                                      */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psRadarParameters = CPLGetXMLNode(psProduct,
		"=product.sourceAttributes.radarParameters");

	struct NoiseLevel noiseLevelBeta0;
	struct NoiseLevel noiseLevelSigma0;
	struct NoiseLevel noiseLevelGamma;

	if (psRadarParameters != nullptr) {
		const char *pszPolarizations = CPLGetXMLValue(psRadarParameters,
			"polarizations", "");
		poDS->SetMetadataItem("POLARIZATIONS", pszPolarizations);

		const char *psAcquisitionType = CPLGetXMLValue(psRadarParameters,
			"acquisitionType", "UNK");
		poDS->SetMetadataItem("ACQUISITION_TYPE", psAcquisitionType);

		const char *psBeams = CPLGetXMLValue(psRadarParameters,
			"beams", "UNK");
		poDS->SetMetadataItem("BEAMS", psBeams);

		// Load Beta Nought, Sigma Nought, Gamma noise levels
		// Loop through all nodes with spaces
		CPLXMLNode *psNodeInc;
		for (psNodeInc = psRadarParameters->psChild; psNodeInc != nullptr;
			psNodeInc = psNodeInc->psNext)
		{
			if (EQUAL(psNodeInc->pszValue, "referenceNoiseLevel")) {
				const char *pszLUTType = CPLGetXMLValue(psNodeInc,
					"incidenceAngleCorrection", nullptr);
				CPLXMLNode *psPixelFirstNoiseValue =
					CPLGetXMLNode(psNodeInc, "pixelFirstNoiseValue");
				CPLXMLNode *psStepSize =
					CPLGetXMLNode(psNodeInc, "stepSize");
				CPLXMLNode *psNumberOfValues =
					CPLGetXMLNode(psNodeInc, "numberOfNoiseLevelValues");
				CPLXMLNode *psNoiseLevelValues =
					CPLGetXMLNode(psNodeInc, "noiseLevelValues");

				if (pszLUTType != nullptr && psPixelFirstNoiseValue != nullptr &&
					psStepSize != nullptr && psNumberOfValues != nullptr &&
					psNoiseLevelValues != nullptr) {
					int pixelFirstLutValueNoiseLevels = atoi(CPLGetXMLValue(psPixelFirstNoiseValue, "", "0"));
					int stepSizeNoiseLevels = atoi(CPLGetXMLValue(psStepSize, "", "0"));
					int numberOfValuesNoiseLevels = atoi(CPLGetXMLValue(psNumberOfValues, "", "0"));
					const char * noiseLevelValues = CPLGetXMLValue(psNoiseLevelValues, "", "");
					char **papszNoiseLevelList = CSLTokenizeString2(noiseLevelValues, " ", CSLT_HONOURSTRINGS);
					/* Get the Pixel Per range */
					int m_nTableNoiseLevelsSize = abs(stepSizeNoiseLevels) * abs(numberOfValuesNoiseLevels);

					if (EQUAL(pszLUTType, "Beta Nought")) {
						noiseLevelBeta0.pixelFirstLutValueNoiseLevels = pixelFirstLutValueNoiseLevels;
						noiseLevelBeta0.stepSizeNoiseLevels = stepSizeNoiseLevels;
						noiseLevelBeta0.numberOfValuesNoiseLevels = numberOfValuesNoiseLevels;
						noiseLevelBeta0.m_nTableNoiseLevelsSize = m_nTableNoiseLevelsSize;
						
						/* Allocate the right Noise Levels size according to the product range pixel */
						noiseLevelBeta0.m_nfTableNoiseLevels = InterpolateValues(papszNoiseLevelList,
																					m_nTableNoiseLevelsSize,
																					stepSizeNoiseLevels,
																					numberOfValuesNoiseLevels,
																					pixelFirstLutValueNoiseLevels);

					}
					else if (EQUAL(pszLUTType, "Sigma Nought")) {
						noiseLevelSigma0.pixelFirstLutValueNoiseLevels = pixelFirstLutValueNoiseLevels;
						noiseLevelSigma0.stepSizeNoiseLevels = stepSizeNoiseLevels;
						noiseLevelSigma0.numberOfValuesNoiseLevels = numberOfValuesNoiseLevels;
						noiseLevelSigma0.m_nTableNoiseLevelsSize = m_nTableNoiseLevelsSize;

						/* Allocate the right Noise Levels size according to the product range pixel */
						noiseLevelSigma0.m_nfTableNoiseLevels = InterpolateValues(papszNoiseLevelList,
							m_nTableNoiseLevelsSize,
							stepSizeNoiseLevels,
							numberOfValuesNoiseLevels,
							pixelFirstLutValueNoiseLevels);
					}
					else if (EQUAL(pszLUTType, "Gamma")) {
						noiseLevelGamma.pixelFirstLutValueNoiseLevels = pixelFirstLutValueNoiseLevels;
						noiseLevelGamma.stepSizeNoiseLevels = stepSizeNoiseLevels;
						noiseLevelGamma.numberOfValuesNoiseLevels = numberOfValuesNoiseLevels;
						noiseLevelGamma.m_nTableNoiseLevelsSize = m_nTableNoiseLevelsSize;

						/* Allocate the right Noise Levels size according to the product range pixel */
						noiseLevelGamma.m_nfTableNoiseLevels = InterpolateValues(papszNoiseLevelList,
							m_nTableNoiseLevelsSize,
							stepSizeNoiseLevels,
							numberOfValuesNoiseLevels,
							pixelFirstLutValueNoiseLevels);
					}

					CSLDestroy(papszNoiseLevelList);
				}
			}
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Open each of the data files as a complex band.                  */
	/* -------------------------------------------------------------------- */
    CPLString osBeta0LUT;
    CPLString osGammaLUT;
    CPLString osSigma0LUT;

	CPLXMLNode *psNode;
	char *pszPath = CPLStrdup(CPLGetPath(osMDFilename));
	char *pszBuf;
	int nFLen = strlen(osMDFilename);

	for (psNode = psImageAttributes->psChild;
		psNode != nullptr;
		psNode = psNode->psNext)
	{
		const char *pszBasename;
		if (psNode->eType != CXT_Element
			|| !(EQUAL(psNode->pszValue, "fullResolutionImageData")
				|| EQUAL(psNode->pszValue, "lookupTable")))
			continue;

		if (EQUAL(psNode->pszValue, "lookupTable") && bCanCalib) {
			/* Determine which incidence angle correction this is */
			const char *pszLUTType = CPLGetXMLValue(psNode,
				"incidenceAngleCorrection", "");
			const char *pszLUTFile = CPLGetXMLValue(psNode, "", "");
			CPLString osLUTFilePath = CPLFormFilename(pszPath, pszLUTFile,
				nullptr);

			if (EQUAL(pszLUTType, ""))
				continue;
			else if (EQUAL(pszLUTType, "Beta Nought") &&
				IsValidXMLFile(pszPath, pszLUTFile))
			{
				poDS->papszExtraFiles =
					CSLAddString(poDS->papszExtraFiles, osLUTFilePath);

				pszBuf = (char *)CPLMalloc(nFLen + 27);
				osBeta0LUT = pszLUTFile;
				poDS->SetMetadataItem("BETA_NOUGHT_LUT", pszLUTFile);

				sprintf(pszBuf, "RADARSAT_2_CALIB:BETA0:%s",
					osMDFilename.c_str());
				poDS->papszSubDatasets = CSLSetNameValue(
					poDS->papszSubDatasets, "SUBDATASET_3_NAME", pszBuf);
				poDS->papszSubDatasets = CSLSetNameValue(
					poDS->papszSubDatasets, "SUBDATASET_3_DESC",
					"Beta Nought calibrated");
				CPLFree(pszBuf);
			}
			else if (EQUAL(pszLUTType, "Sigma Nought") &&
				IsValidXMLFile(pszPath, pszLUTFile))
			{
				poDS->papszExtraFiles =
					CSLAddString(poDS->papszExtraFiles, osLUTFilePath);

				pszBuf = (char *)CPLMalloc(nFLen + 27);
				osSigma0LUT = pszLUTFile;
				poDS->SetMetadataItem("SIGMA_NOUGHT_LUT", pszLUTFile);

				sprintf(pszBuf, "RADARSAT_2_CALIB:SIGMA0:%s",
					osMDFilename.c_str());
				poDS->papszSubDatasets = CSLSetNameValue(
					poDS->papszSubDatasets, "SUBDATASET_2_NAME", pszBuf);
				poDS->papszSubDatasets = CSLSetNameValue(
					poDS->papszSubDatasets, "SUBDATASET_2_DESC",
					"Sigma Nought calibrated");
				CPLFree(pszBuf);
			}
			else if (EQUAL(pszLUTType, "Gamma") &&
				IsValidXMLFile(pszPath, pszLUTFile))
			{
				poDS->papszExtraFiles =
					CSLAddString(poDS->papszExtraFiles, osLUTFilePath);

				pszBuf = (char *)CPLMalloc(nFLen + 27);
				osGammaLUT = pszLUTFile;
				poDS->SetMetadataItem("GAMMA_LUT", pszLUTFile);
				sprintf(pszBuf, "RADARSAT_2_CALIB:GAMMA:%s",
					osMDFilename.c_str());
				poDS->papszSubDatasets = CSLSetNameValue(
					poDS->papszSubDatasets, "SUBDATASET_4_NAME", pszBuf);
				poDS->papszSubDatasets = CSLSetNameValue(
					poDS->papszSubDatasets, "SUBDATASET_4_DESC",
					"Gamma calibrated");
				CPLFree(pszBuf);
			}
			continue;
		}

		/* -------------------------------------------------------------------- */
		/*      Fetch filename.                                                 */
		/* -------------------------------------------------------------------- */
		pszBasename = CPLGetXMLValue(psNode, "", "");
		if (*pszBasename == '\0')
			continue;

		/* -------------------------------------------------------------------- */
		/*      Form full filename (path of product.xml + basename).            */
		/* -------------------------------------------------------------------- */
		char *pszFullname =
			CPLStrdup(CPLFormFilename(pszPath, pszBasename, nullptr));

		/* -------------------------------------------------------------------- */
		/*      Try and open the file.                                          */
		/* -------------------------------------------------------------------- */
		GDALDataset *poBandFile;

		poBandFile = (GDALDataset *)GDALOpen(pszFullname, GA_ReadOnly);
		if (poBandFile == nullptr)
		{
			CPLFree(pszFullname);
			continue;
		}
		if (poBandFile->GetRasterCount() == 0)
		{
			GDALClose((GDALRasterBandH)poBandFile);
			CPLFree(pszFullname);
			continue;
		}

		/* Some CFloat32 NITF files have nBitsPerSample incorrectly reported    */
		/* as 16, and get misinterpreted as CInt16.  Check the underlying NITF  */
		/* and override if this is the case.                                    */
		if (poBandFile->GetRasterBand(1)->GetRasterDataType() == GDT_CFloat32)
			eDataType = GDT_CFloat32;

		BandMappingRS2 b = checkBandFileMappingRS2(eDataType, poBandFile);
		if (b == BANDERROR)
		{
			GDALClose((GDALRasterBandH)poBandFile);
			CPLFree(pszFullname);
			delete poDS;
			CPLError(CE_Failure, CPLE_AppDefined,
				"The underlying band files do not have an appropriate data type.");
			return nullptr;
		}
		bool twoBandComplex = b == TWOBANDCOMPLEX;


		poDS->papszExtraFiles = CSLAddString(poDS->papszExtraFiles,
			pszFullname);

		/* -------------------------------------------------------------------- */
		/*      Create the band.                                                */
		/* -------------------------------------------------------------------- */
		if (eCalib == None || eCalib == Uncalib) {
			RS2RasterBand *poBand;
			// A Complex type remains the same type
			poBand = new RS2RasterBand(poDS, eDataType,
				CPLGetXMLValue(psNode, "pole", ""),
				poBandFile, twoBandComplex);

			poDS->SetBand(poDS->GetRasterCount() + 1, poBand);
		}
		else {
			struct NoiseLevel *noiseLevel  = nullptr;
			const char *pszLUT;
			switch (eCalib) {
			case Sigma0:
				noiseLevel = &noiseLevelSigma0;
				pszLUT = osSigma0LUT;
				break;
			case Beta0:
				noiseLevel = &noiseLevelBeta0;
				pszLUT = osBeta0LUT;
				break;
			case Gamma:
				noiseLevel = &noiseLevelGamma;
				pszLUT = osGammaLUT;
				break;
			default:
				/* we should bomb gracefully... */
				pszLUT = osSigma0LUT;
			}
			RS2CalibRasterBand *poBand;
			
			if (isComplexType) {
				// If Complex, always 32 bits
				poBand = new RS2CalibRasterBand(poDS, CPLGetXMLValue(psNode,
					"pole", ""), GDT_Float32, poBandFile, eCalib,
					CPLFormFilename(pszPath, pszLUT, nullptr),
					noiseLevel, eDataType);
			}
			else {
				// Whatever the datatype was previoulsy set
				poBand = new RS2CalibRasterBand(poDS, CPLGetXMLValue(psNode,
					"pole", ""), eDataType, poBandFile, eCalib,
					CPLFormFilename(pszPath, pszLUT, nullptr),
					noiseLevel, eDataType);
			}

			poDS->SetBand(poDS->GetRasterCount() + 1, poBand);
		}

		CPLFree(pszFullname);
	}

	if (noiseLevelBeta0.m_nfTableNoiseLevels != nullptr) {
		CPLFree(noiseLevelBeta0.m_nfTableNoiseLevels);
	}
	
	if (noiseLevelSigma0.m_nfTableNoiseLevels != nullptr) {
		CPLFree(noiseLevelSigma0.m_nfTableNoiseLevels);
	}
	
	if (noiseLevelGamma.m_nfTableNoiseLevels != nullptr) {
		CPLFree(noiseLevelGamma.m_nfTableNoiseLevels);
	}

	if (poDS->papszSubDatasets != nullptr && eCalib == None) {
		char *pszBuf2;
		pszBuf2 = (char *)CPLMalloc(nFLen + 28);
		sprintf(pszBuf2, "RADARSAT_2_CALIB:UNCALIB:%s",
			osMDFilename.c_str());
		poDS->papszSubDatasets = CSLSetNameValue(poDS->papszSubDatasets,
			"SUBDATASET_1_NAME", pszBuf2);
		poDS->papszSubDatasets = CSLSetNameValue(poDS->papszSubDatasets,
			"SUBDATASET_1_DESC", "Uncalibrated digital numbers");
		CPLFree(pszBuf2);
	}
	else if (poDS->papszSubDatasets != nullptr) {
		CSLDestroy(poDS->papszSubDatasets);
		poDS->papszSubDatasets = nullptr;
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
	/*      Collect a few useful metadata items                             */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psSourceAttrs = CPLGetXMLNode(psProduct,
		"=product.sourceAttributes");
	const char *pszItem;

	pszItem = CPLGetXMLValue(psProduct, "=product.productId", "UNK");
	poDS->SetMetadataItem("PRODUCT_ID", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs,
		"satellite", "");
	poDS->SetMetadataItem("SATELLITE_IDENTIFIER", pszItem);

	pszItem = CPLGetXMLValue(psSourceAttrs,
		"sensor", "");
	poDS->SetMetadataItem("SENSOR_IDENTIFIER", pszItem);

	if (psSourceAttrs != nullptr) {
		/* Get beam mode mnemonic */
		pszItem = CPLGetXMLValue(psSourceAttrs, "beamModeMnemonic", "UNK");
		poDS->SetMetadataItem("BEAM_MODE", pszItem);
		pszItem = CPLGetXMLValue(psSourceAttrs, "rawDataStartTime", "UNK");
		poDS->SetMetadataItem("ACQUISITION_START_TIME", pszItem);
	}

	CPLXMLNode *psSarProcessingInformation =
		CPLGetXMLNode(psProduct, "=product.imageGenerationParameters");

	if (psSarProcessingInformation != nullptr) {
		/* Get incidence angle information */
		pszItem = CPLGetXMLValue(psSarProcessingInformation,
			"sarProcessingInformation.incidenceAngleNearRange", "UNK");
		poDS->SetMetadataItem("NEAR_RANGE_INCIDENCE_ANGLE", pszItem);

		pszItem = CPLGetXMLValue(psSarProcessingInformation,
			"sarProcessingInformation.incidenceAngleFarRange", "UNK");
		poDS->SetMetadataItem("FAR_RANGE_INCIDENCE_ANGLE", pszItem);

		pszItem = CPLGetXMLValue(psSarProcessingInformation,
			"sarProcessingInformation.slantRangeNearEdge", "UNK");
		poDS->SetMetadataItem("SLANT_RANGE_NEAR_EDGE", pszItem);

		pszItem = CPLGetXMLValue(psSarProcessingInformation,
			"sarProcessingInformation.satelliteHeight", "UNK");
		poDS->SetMetadataItem("SATELLITE_HEIGHT", pszItem);
	}

	/*--------------------------------------------------------------------- */
	/*      Collect Map projection/Geotransform information, if present     */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psMapProjection =
		CPLGetXMLNode(psImageAttributes,
			"geographicInformation.mapProjection");

	if (psMapProjection != nullptr) {
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

		if (psPos != nullptr) {
			double testx, testy, br_x, br_y, tl_x, tl_y, tr_x, tr_y,
				bl_x, bl_y;

			tl_x = strtod(CPLGetXMLValue(psPos,
				"upperLeftCorner.mapCoordinate.easting", "0.0"), nullptr);
			tl_y = strtod(CPLGetXMLValue(psPos,
				"upperLeftCorner.mapCoordinate.northing", "0.0"), nullptr);
			bl_x = strtod(CPLGetXMLValue(psPos,
				"lowerLeftCorner.mapCoordinate.easting", "0.0"), nullptr);
			bl_y = strtod(CPLGetXMLValue(psPos,
				"lowerLeftCorner.mapCoordinate.northing", "0.0"), nullptr);
			tr_x = strtod(CPLGetXMLValue(psPos,
				"upperRightCorner.mapCoordinate.easting", "0.0"), nullptr);
			tr_y = strtod(CPLGetXMLValue(psPos,
				"upperRightCorner.mapCoordinate.northing", "0.0"), nullptr);
			poDS->adfGeoTransform[1] = (tr_x - tl_x) / (poDS->nRasterXSize - 1);
			poDS->adfGeoTransform[4] = (tr_y - tl_y) / (poDS->nRasterXSize - 1);
			poDS->adfGeoTransform[2] = (bl_x - tl_x) / (poDS->nRasterYSize - 1);
			poDS->adfGeoTransform[5] = (bl_y - tl_y) / (poDS->nRasterYSize - 1);
			poDS->adfGeoTransform[0] = (tl_x - 0.5*poDS->adfGeoTransform[1]
				- 0.5*poDS->adfGeoTransform[2]);
			poDS->adfGeoTransform[3] = (tl_y - 0.5*poDS->adfGeoTransform[4]
				- 0.5*poDS->adfGeoTransform[5]);

			/* Use bottom right pixel to test geotransform */
			br_x = strtod(CPLGetXMLValue(psPos,
				"lowerRightCorner.mapCoordinate.easting", "0.0"), nullptr);
			br_y = strtod(CPLGetXMLValue(psPos,
				"lowerRightCorner.mapCoordinate.northing", "0.0"), nullptr);
			testx = poDS->adfGeoTransform[0] + poDS->adfGeoTransform[1] *
				(poDS->nRasterXSize - 0.5) + poDS->adfGeoTransform[2] *
				(poDS->nRasterYSize - 0.5);
			testy = poDS->adfGeoTransform[3] + poDS->adfGeoTransform[4] *
				(poDS->nRasterXSize - 0.5) + poDS->adfGeoTransform[5] *
				(poDS->nRasterYSize - 0.5);

			/* Give 1/4 pixel numerical error leeway in calculating location
			based on affine transform */
			if ((fabs(testx - br_x) >
				fabs(0.25*(poDS->adfGeoTransform[1] + poDS->adfGeoTransform[2])))
				|| (fabs(testy - br_y) > fabs(0.25*(poDS->adfGeoTransform[4] +
					poDS->adfGeoTransform[5]))))
			{
				CPLError(CE_Warning, CPLE_AppDefined,
					"Unexpected error in calculating affine transform: "
					"corner coordinates inconsistent.");
			}
			else
				poDS->bHaveGeoTransform = TRUE;

		}

	}


	/* -------------------------------------------------------------------- */
	/*      Collect Projection String Information                           */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psEllipsoid =
		CPLGetXMLNode(psImageAttributes,
			"geographicInformation.referenceEllipsoidParameters");

	if (psEllipsoid != nullptr) {
		const char *pszEllipsoidName;
		double minor_axis, major_axis, inv_flattening;
		OGRSpatialReference oLL, oPrj;

		const char *pszGeodeticTerrainHeight
			= CPLGetXMLValue(psEllipsoid, "geodeticTerrainHeight", "UNK");
		poDS->SetMetadataItem("GEODETIC_TERRAIN_HEIGHT", pszGeodeticTerrainHeight);

		pszEllipsoidName = CPLGetXMLValue(psEllipsoid, "ellipsoidName", "");
		minor_axis = atof(CPLGetXMLValue(psEllipsoid, "semiMinorAxis",
			"0.0"));
		major_axis = atof(CPLGetXMLValue(psEllipsoid, "semiMajorAxis",
			"0.0"));

		if (EQUAL(pszEllipsoidName, "") || (minor_axis == 0.0) ||
			(major_axis == 0.0))
		{
			CPLError(CE_Warning, CPLE_AppDefined, "Warning- incomplete"
				" ellipsoid information.  Using wgs-84 parameters.\n");
			oLL.SetWellKnownGeogCS("WGS84");
			oPrj.SetWellKnownGeogCS("WGS84");
		}
		else if (EQUAL(pszEllipsoidName, "WGS84") || EQUAL(pszEllipsoidName, "WGS 1984")) {
			oLL.SetWellKnownGeogCS("WGS84");
			oPrj.SetWellKnownGeogCS("WGS84");
		}
		else {
			inv_flattening = major_axis / (major_axis - minor_axis);
			oLL.SetGeogCS(nullptr, nullptr, pszEllipsoidName, major_axis,
				inv_flattening);
			oPrj.SetGeogCS(nullptr, nullptr, pszEllipsoidName, major_axis,
				inv_flattening);
		}

		if (psMapProjection != nullptr) {
			const char *pszProj = CPLGetXMLValue(
				psMapProjection, "mapProjectionDescriptor", "");
			bool bUseProjInfo = FALSE;

			CPLXMLNode *psUtmParams =
				CPLGetXMLNode(psMapProjection,
					"utmProjectionParameters");

			CPLXMLNode *psNspParams =
				CPLGetXMLNode(psMapProjection,
					"nspProjectionParameters");

			if ((psUtmParams != nullptr) && poDS->bHaveGeoTransform) {
				const char *pszHemisphere;
				int utmZone;
				double origEasting, origNorthing;
				bool bNorth = TRUE;

				utmZone = atoi(CPLGetXMLValue(psUtmParams, "utmZone", ""));
				pszHemisphere = CPLGetXMLValue(psUtmParams,
					"hemisphere", "");
				origEasting = strtod(CPLGetXMLValue(psUtmParams,
					"mapOriginFalseEasting", "0.0"), nullptr);
				origNorthing = strtod(CPLGetXMLValue(psUtmParams,
					"mapOriginFalseNorthing", "0.0"), nullptr);

				if (EQUALN(pszHemisphere, "southern", 8))
					bNorth = FALSE;

				if (EQUALN(pszProj, "UTM", 3)) {
					oPrj.SetUTM(utmZone, bNorth);
					bUseProjInfo = TRUE;
				}
			}
			else if ((psNspParams != nullptr) && poDS->bHaveGeoTransform) {
				double origEasting, origNorthing, copLong, copLat, sP1, sP2;

				origEasting = strtod(CPLGetXMLValue(psNspParams,
					"mapOriginFalseEasting", "0.0"), nullptr);
				origNorthing = strtod(CPLGetXMLValue(psNspParams,
					"mapOriginFalseNorthing", "0.0"), nullptr);
				copLong = strtod(CPLGetXMLValue(psNspParams,
					"centerOfProjectionLongitude", "0.0"), nullptr);
				copLat = strtod(CPLGetXMLValue(psNspParams,
					"centerOfProjectionLatitude", "0.0"), nullptr);
				sP1 = strtod(CPLGetXMLValue(psNspParams,
					"standardParallels1", "0.0"), nullptr);
				sP2 = strtod(CPLGetXMLValue(psNspParams,
					"standardParallels2", "0.0"), nullptr);

				if (EQUALN(pszProj, "ARC", 3)) {
					/* Albers Conical Equal Area */
					oPrj.SetACEA(sP1, sP2, copLat, copLong, origEasting,
						origNorthing);
					bUseProjInfo = TRUE;
				}
				else if (EQUALN(pszProj, "LCC", 3)) {
					/* Lambert Conformal Conic */
					oPrj.SetLCC(sP1, sP2, copLat, copLong, origEasting,
						origNorthing);
					bUseProjInfo = TRUE;
				}
				else if (EQUALN(pszProj, "STPL", 3)) {
					/* State Plate
					ASSUMPTIONS: "zone" in product.xml matches USGS
					definition as expected by ogr spatial reference; NAD83
					zones (versus NAD27) are assumed. */

					int nSPZone = atoi(CPLGetXMLValue(psNspParams,
						"zone", "1"));

					oPrj.SetStatePlane(nSPZone, TRUE, nullptr, 0.0);
					bUseProjInfo = TRUE;
				}
			}

			if (bUseProjInfo) {
				CPLFree(poDS->pszProjection);
				poDS->pszProjection = nullptr;
				oPrj.exportToWkt(&(poDS->pszProjection));
			}
			else {
				CPLError(CE_Warning, CPLE_AppDefined, "Unable to interpret "
					"projection information; check mapProjection info in "
					"product.xml!");
			}
		}

		CPLFree(poDS->pszGCPProjection);
		poDS->pszGCPProjection = nullptr;
		oLL.exportToWkt(&(poDS->pszGCPProjection));

	}

	/* -------------------------------------------------------------------- */
	/*      Collect GCPs.                                                   */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psGeoGrid =
		CPLGetXMLNode(psImageAttributes,
			"geographicInformation.geolocationGrid");

	if (psGeoGrid != nullptr) {
		/* count GCPs */
		poDS->nGCPCount = 0;

		for (psNode = psGeoGrid->psChild; psNode != nullptr;
			psNode = psNode->psNext)
		{
			if (EQUAL(psNode->pszValue, "imageTiePoint"))
				poDS->nGCPCount++;
		}

		poDS->pasGCPList = (GDAL_GCP *)
			CPLCalloc(sizeof(GDAL_GCP), poDS->nGCPCount);

		poDS->nGCPCount = 0;

		for (psNode = psGeoGrid->psChild; psNode != nullptr;
			psNode = psNode->psNext)
		{
			char    szID[32];
			GDAL_GCP   *psGCP = poDS->pasGCPList + poDS->nGCPCount;

			if (!EQUAL(psNode->pszValue, "imageTiePoint"))
				continue;

			poDS->nGCPCount++;

			sprintf(szID, "%d", poDS->nGCPCount);
			psGCP->pszId = CPLStrdup(szID);
			psGCP->pszInfo = CPLStrdup("");
			psGCP->dfGCPPixel =
				atof(CPLGetXMLValue(psNode, "imageCoordinate.pixel", "0"));
			psGCP->dfGCPLine =
				atof(CPLGetXMLValue(psNode, "imageCoordinate.line", "0"));
			psGCP->dfGCPX =
				atof(CPLGetXMLValue(psNode, "geodeticCoordinate.longitude", ""));
			psGCP->dfGCPY =
				atof(CPLGetXMLValue(psNode, "geodeticCoordinate.latitude", ""));
			psGCP->dfGCPZ =
				atof(CPLGetXMLValue(psNode, "geodeticCoordinate.height", ""));
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Collect RPCs.                                                   */
	/* -------------------------------------------------------------------- */
	CPLXMLNode *psRPCNode =
		CPLGetXMLNode(psImageAttributes,
			"geographicInformation.rationalFunctions");

	if (psRPCNode != nullptr) {
		//the R2 dataset presents the RPC information exactly as needed by GDAL
		//a simple mapping of metadata names will do
		static const char* rpcMap[] = { "ERR_BIAS",        "biasError",
			"ERR_RAND",        "randomError",
			"LINE_OFF",        "lineOffset",
			"SAMP_OFF",        "pixelOffset",
			"LAT_OFF",         "latitudeOffset",
			"LONG_OFF",        "longitudeOffset",
			"HEIGHT_OFF",      "heightOffset",
			"LINE_SCALE",      "lineScale",
			"SAMP_SCALE",      "pixelScale",
			"LAT_SCALE",       "latitudeScale",
			"LONG_SCALE",      "longitudeScale",
			"HEIGHT_SCALE",    "heightScale",
			"LINE_NUM_COEFF",  "lineNumeratorCoefficients",
			"LINE_DEN_COEFF",  "lineDenominatorCoefficients",
			"SAMP_NUM_COEFF",  "pixelNumeratorCoefficients",
			"SAMP_DEN_COEFF",  "pixelDenominatorCoefficients" };

		for (int i = 0; i<32; i += 2)        //16 pairs in above list
		{
			pszItem = CPLGetXMLValue(psRPCNode, rpcMap[i + 1], nullptr);
			if (pszItem)
				poDS->SetMetadataItem(rpcMap[i], pszItem, "RPC");
		}
	}

	CPLFree(pszPath);

	/* -------------------------------------------------------------------- */
	/*      Initialize any PAM information.                                 */
	/* -------------------------------------------------------------------- */
	CPLString osDescription;
	CPLString osSubdatasetName;
	bool useSubdatasets = true;

	switch (eCalib) {
	case Sigma0:
		osSubdatasetName = "SIGMA0";
		osDescription.Printf("RADARSAT_2_CALIB:SIGMA0:%s",
			osMDFilename.c_str());
		break;
	case Beta0:
		osSubdatasetName = "BETA0";
		osDescription.Printf("RADARSAT_2_CALIB:BETA0:%s",
			osMDFilename.c_str());
		break;
	case Gamma:
		osSubdatasetName = "GAMMA0";
		osDescription.Printf("RADARSAT_2_CALIB:GAMMA0:%s",
			osMDFilename.c_str());
		break;
	case Uncalib:
		osSubdatasetName = "UNCALIB";
		osDescription.Printf("RADARSAT_2_CALIB:UNCALIB:%s",
			osMDFilename.c_str());
		break;
	default:
		osDescription = osMDFilename;
		osSubdatasetName = "UNCALIB";
		useSubdatasets = false;
		break;
	}

	if (eCalib != None)
		poDS->papszExtraFiles =
		CSLAddString(poDS->papszExtraFiles, osMDFilename);

	/* -------------------------------------------------------------------- */
	/*      Initialize any PAM information.                                 */
	/* -------------------------------------------------------------------- */
	poDS->SetDescription(osDescription);

	poDS->SetPhysicalFilename(osMDFilename);
	poDS->SetSubdatasetName(osSubdatasetName);

	poDS->TryLoadXML();

	/* -------------------------------------------------------------------- */
	/*      Check for overviews.                                            */
	/* -------------------------------------------------------------------- */
	if (useSubdatasets)
		poDS->oOvManager.Initialize(poDS, ":::VIRTUAL:::");
	else
		poDS->oOvManager.Initialize(poDS, osMDFilename);

	return(poDS);
}

/************************************************************************/
/*                            GetGCPCount()                             */
/************************************************************************/

int RS2Dataset::GetGCPCount()

{
    return nGCPCount;
}

/************************************************************************/
/*                          GetGCPProjection()                          */
/************************************************************************/

const char *RS2Dataset::GetGCPProjection()

{
    return pszGCPProjection;
}

/************************************************************************/
/*                               GetGCPs()                              */
/************************************************************************/

const GDAL_GCP *RS2Dataset::GetGCPs()

{
    return pasGCPList;
}

/************************************************************************/
/*                          GetProjectionRef()                          */
/************************************************************************/

const char *RS2Dataset::GetProjectionRef()

{
    return pszProjection;
}

/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr RS2Dataset::GetGeoTransform( double * padfTransform )

{
    memcpy( padfTransform,  adfGeoTransform, sizeof(double) * 6 );

    if (bHaveGeoTransform)
        return CE_None;

    return CE_Failure;
}

/************************************************************************/
/*                      GetMetadataDomainList()                         */
/************************************************************************/

char **RS2Dataset::GetMetadataDomainList()
{
    return BuildMetadataDomainList(GDALDataset::GetMetadataDomainList(),
                                   TRUE,
                                   "SUBDATASETS", nullptr);
}

/************************************************************************/
/*                            GetMetadata()                             */
/************************************************************************/

char **RS2Dataset::GetMetadata( const char *pszDomain )

{
    if( pszDomain != nullptr && STARTS_WITH_CI(pszDomain, "SUBDATASETS") &&
        papszSubDatasets != nullptr)
        return papszSubDatasets;

    return GDALDataset::GetMetadata( pszDomain );
}

/************************************************************************/
/*                         GDALRegister_RS2()                          */
/************************************************************************/

void GDALRegister_RS2()

{
    if( GDALGetDriverByName( "RS2" ) != nullptr )
        return;

    GDALDriver *poDriver = new GDALDriver();

    poDriver->SetDescription( "RS2" );
    poDriver->SetMetadataItem( GDAL_DCAP_RASTER, "YES" );
    poDriver->SetMetadataItem( GDAL_DMD_LONGNAME, "RadarSat 2 XML Product" );
    poDriver->SetMetadataItem( GDAL_DMD_HELPTOPIC, "frmt_rs2.html" );
    poDriver->SetMetadataItem( GDAL_DMD_SUBDATASETS, "YES" );
    poDriver->SetMetadataItem( GDAL_DCAP_VIRTUALIO, "YES" );

    poDriver->pfnOpen = RS2Dataset::Open;
    poDriver->pfnIdentify = RS2Dataset::Identify;

    GetGDALDriverManager()->RegisterDriver( poDriver );
}

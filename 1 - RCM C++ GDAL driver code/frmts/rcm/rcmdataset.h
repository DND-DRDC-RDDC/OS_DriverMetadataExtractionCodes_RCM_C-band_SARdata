#ifndef GDAL_RCM_H_INCLUDED
#define GDAL_RCM_H_INCLUDED

/* Roberto Fix */
#include "gdal_pam.h"
#include "gdal_lut.h"


// Should be size of larged possible filename.
static const int CPL_PATH_BUF_SIZE = 2048;
static const char szLayerCalibration[] = "RCM_CALIB";
static const char szLayerSeparator[] = ":";
static const char szSIGMA0[] = "SIGMA0";
static const char szGAMMA[] = "GAMMA";
static const char szBETA0[] = "BETA0";
static const char szUNCALIB[] = "UNCALIB";
static const char szPathSeparator[] =
#ifdef _WIN32 /* Defined if Win32 and Win64 */
"\\";
#else
"/";
#endif
static const char cPathSeparator =
#ifdef _WIN32 /* Defined if Win32 and Win64 */
'\\';
#else
'/';
#endif
static const char cOppositePathSeparator =
#ifdef _WIN32 /* Defined if Win32 and Win64 */
'/';
#else
'\\';
#endif

/************************************************************************/
/* ==================================================================== */
/*                               RCMDataset                             */
/* ==================================================================== */
/************************************************************************/

class RCMDataset : public GDALPamDataset
{
	CPLXMLNode *psProduct;

	int         nGCPCount;
	GDAL_GCP   *pasGCPList;
	char       *pszGCPProjection;
	char      **papszSubDatasets;
	char       *pszProjection;
	char       *pszLutApplied;
	double      adfGeoTransform[6];
	bool        bHaveGeoTransform;
	bool		bPerPolarizationScaling;
	bool        isComplexData;
	int         magnitudeBits;
	int         realBitsComplexData;
	int         imaginaryBitsComplexData;
	char      **papszExtraFiles;
	double     *m_nfIncidenceAngleTable;
	int         m_IncidenceAngleTableSize;

protected:
	virtual int         CloseDependentDatasets() override;

public:
	RCMDataset();
	virtual ~RCMDataset();

	virtual int    GetGCPCount() override;
	virtual const char *GetGCPProjection() override;
	virtual const GDAL_GCP *GetGCPs() override;

	virtual const char *GetProjectionRef(void) override;
	virtual CPLErr GetGeoTransform(double *) override;

	virtual char      **GetMetadataDomainList() override;
	virtual char **GetMetadata(const char * pszDomain = "") override;
	virtual char **GetFileList(void) override;

	static GDALDataset *Open(GDALOpenInfo *);
	static int Identify(GDALOpenInfo *);

	CPLXMLNode *GetProduct() { return psProduct; }

	/* If False, this is Magnitude,   True, Complex data with Real and Imaginary*/
	bool IsComplexData() { return isComplexData; }

	/* These 2 variables are used in case of Complex Data */
	int GetRealBitsComplexData() { return realBitsComplexData; }
	int GetImaginaryBitsComplexData() { return imaginaryBitsComplexData; }

	/* This variable is used in case of Magnitude */
	int GetMagnitudeBits() { return magnitudeBits; }

	/* This variable is used to hold the Incidence Angle */
	double * GetIncidenceAngle() { return m_nfIncidenceAngleTable; }

	/* This variable is used to hold the Incidence Angle Table Size */
	int GetIncidenceAngleSize() { return m_IncidenceAngleTableSize; }
};

/************************************************************************/
/* ==================================================================== */
/*                    RCMRasterBand                           */
/* ==================================================================== */
/************************************************************************/

class RCMRasterBand : public GDALPamRasterBand
{
	eCalibration m_eCalib;
	GDALDataset     *poBandFile;
	RCMDataset      *poRCMDataset;
	GDALDataset *m_poBandDataset;
	GDALDataType m_eType; /* data type of data being ingested */
	

	double *m_nfTable;
	int m_nTableSize;
	double m_nfOffset;
	char *m_pszLUTFile;
	int pixelFirstLutValue;
	int stepSize;
	int numberOfValues;
	GDALRasterBand *poBand;

	//2 bands representing I+Q -> one complex band
	//otherwise poBandFile is passed straight through
	bool            twoBandComplex;

	bool 		isOneFilePerPol;
	bool		isNITF;

public:
	RCMRasterBand(RCMDataset *poDSIn,
		int nBandIn,
		GDALDataType eDataTypeIn,
		const char *pszPole,
		GDALDataset *poBandFile,
		bool bTwoBandComplex, bool isOneFilePerPol, bool isNITF);

	virtual     ~RCMRasterBand();

	virtual CPLErr IReadBlock(int, int, void *) override;

	bool IsExistLUT();

	double GetLUT(int pixel);

	const char *GetLUTFilename();

	int GetLUTsize();

	double GetLUTOffset();

	void SetPartialLUT(int pixel_offset, int pixel_width);

	bool IsComplex();

	bool IsExistNoiseLevels();

	double GetNoiseLevels(int pixel);

	const char *GetNoiseLevelsFilename();

	int GetNoiseLevelsSize();

	eCalibration GetCalibration();

	static GDALDataset *Open(GDALOpenInfo *);
};



/************************************************************************/
/* ==================================================================== */
/*                         RCMCalibRasterBand                           */
/* ==================================================================== */
/************************************************************************/
/* Returns data that has been calibrated to sigma nought, gamma         */
/* or beta nought.                                                      */
/************************************************************************/

class RCMCalibRasterBand : public GDALPamRasterBand {
private:
	eCalibration m_eCalib;
	RCMDataset *m_poRCMDataset;
	GDALDataset *m_poBandDataset;
	GDALDataType m_eType; /* data type of data being ingested */
	GDALDataType m_eOriginalType; /* data type that used to be before transformation */

	double *m_nfTable;
	int m_nTableSize;
	double m_nfOffset;
	char *m_pszLUTFile;
	int pixelFirstLutValue;
	int stepSize;
	int numberOfValues;

	char *m_pszNoiseLevelsFile;
	double *m_nfTableNoiseLevels;
	int pixelFirstLutValueNoiseLevels;
	int stepSizeNoiseLevels;
	int numberOfValuesNoiseLevels;
	int m_nTableNoiseLevelsSize;

	void ReadLUT();
	void ReadNoiseLevels();
public:
	RCMCalibRasterBand(
		RCMDataset *poDataset, const char *pszPolarization,
		GDALDataType eType, GDALDataset *poBandDataset, eCalibration eCalib,
		const char *pszLUT, const char *pszNoiseLevels, 
		GDALDataType eOriginalType);
	~RCMCalibRasterBand();

	CPLErr IReadBlock(int nBlockXOff, int nBlockYOff, void *pImage) override;

	bool IsExistLUT();

	double GetLUT(int pixel);

	const char *GetLUTFilename();

	int GetLUTsize();

	double GetLUTOffset();

	void SetPartialLUT(int pixel_offset, int pixel_width);

	bool IsExistNoiseLevels();

	double GetNoiseLevels(int pixel);

	const char *GetNoiseLevelsFilename();

	int GetNoiseLevelsSize();

	bool IsComplex();

	eCalibration GetCalibration();

	double * CloneLUT();

	double * CloneNoiseLevels();

};


#endif /* ndef GDAL_RCM_H_INCLUDED */

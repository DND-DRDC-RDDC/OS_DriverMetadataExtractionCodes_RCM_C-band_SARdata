#ifndef GDAL_RS2_H_INCLUDED
#define GDAL_RS2_H_INCLUDED


/* Roberto's Fixed */

#include "gdal_pam.h"

struct NoiseLevel {
	double *m_nfTableNoiseLevels = NULL;
	int pixelFirstLutValueNoiseLevels = 0;
	int stepSizeNoiseLevels = 0;
	int numberOfValuesNoiseLevels = 0;
	int m_nTableNoiseLevelsSize = 0;
};

/************************************************************************/
/* ==================================================================== */
/*                         RS2CalibRasterBand                           */
/* ==================================================================== */
/************************************************************************/
/* Returns data that has been calibrated to sigma nought, gamma         */
/* or beta nought.                                                      */
/************************************************************************/

/************************************************************************/
/* ==================================================================== */
/*                               RS2Dataset                             */
/* ==================================================================== */
/************************************************************************/

class RS2Dataset final: public GDALPamDataset
{
	CPLXMLNode *psProduct;

	int           nGCPCount;
	GDAL_GCP      *pasGCPList;
	char          *pszGCPProjection;
	char        **papszSubDatasets;
	char          *pszProjection;
	double      adfGeoTransform[6];
	bool        bHaveGeoTransform;

	char        **papszExtraFiles;

protected:
	virtual int         CloseDependentDatasets() override;

public:
	RS2Dataset();
	virtual ~RS2Dataset();

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

	/* This variable is used to hold the Incidence Angle, RS2 has no information */
	double * GetIncidenceAmgle() { return NULL; }

};

/************************************************************************/
/* ==================================================================== */
/*                    RS2RasterBand                           */
/* ==================================================================== */
/************************************************************************/

class RS2RasterBand final: public GDALPamRasterBand
{
private:
	eCalibration     m_eCalib;
	GDALDataset     *poBandFile;
	GDALDataType m_eType; /* data type of data being ingested */
	double *m_nfTable;
	int m_nTableSize;
	double m_nfOffset;
	char *m_pszLUTFile;

	//2 bands representing I+Q -> one complex band
	//otherwise poBandFile is passed straight through
	bool            twoBandComplex;

public:
	RS2RasterBand(RS2Dataset *poDSIn,
		GDALDataType eDataTypeIn,
		const char *pszPole,
		GDALDataset *poBandFile,
		bool twoBandComplex);
	virtual     ~RS2RasterBand();

	virtual CPLErr IReadBlock(int, int, void *) override;

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

	static GDALDataset *Open(GDALOpenInfo *);
};

class RS2CalibRasterBand : public GDALPamRasterBand {
private:
	eCalibration m_eCalib;
	GDALDataset *m_poBandDataset;
	GDALDataType m_eType; /* data type of data being ingested */
	GDALDataType m_eOriginalType; /* data type that used to be before transformation */
	double *m_nfTable;
	int m_nTableSize;
	double m_nfOffset;
	char *m_pszLUTFile;

	double *m_nfTableNoiseLevels;
	int pixelFirstLutValueNoiseLevels;
	int stepSizeNoiseLevels;
	int numberOfValuesNoiseLevels;
	int m_nTableNoiseLevelsSize;

	void ReadLUT();
public:
	RS2CalibRasterBand(
		RS2Dataset *poDataset, const char *pszPolarization,
		GDALDataType eType, GDALDataset *poBandDataset, eCalibration eCalib,
		const char *pszLUT, struct NoiseLevel *noiseLevel, GDALDataType eOriginalType);
	~RS2CalibRasterBand();

	CPLErr IReadBlock(int nBlockXOff, int nBlockYOff, void *pImage) override;

	bool IsExistLUT();

	double GetLUT(int pixel);

	const char *GetLUTFilename();

	int GetLUTsize();

	double GetLUTOffset();

	bool IsExistNoiseLevels();

	double GetNoiseLevels(int pixel);

	const char *GetNoiseLevelsFilename();

	int GetNoiseLevelsSize();

	bool IsComplex();

	eCalibration GetCalibration();

	void SetPartialLUT(int pixel_offset, int pixel_width);

	double * CloneLUT();

	double * CloneNoiseLevels();
};

#endif /* ndef GDAL_RS2_H_INCLUDED */

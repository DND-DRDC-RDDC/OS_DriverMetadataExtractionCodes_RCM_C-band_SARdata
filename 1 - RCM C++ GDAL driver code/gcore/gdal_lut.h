#ifndef GDAL_LUT_H_INCLUDED
#define GDAL_LUT_H_INCLUDED

#include "gdal_pam.h"

/* Start: Roberto July, 2018 */
int CPL_DLL CPL_STDCALL GetMetadataLutValues(GDALDataset *ds, double **values, char *bandNumber);
void CPL_DLL CPL_STDCALL CalculateComplexSigmaLutDB(GDALDataset *dst, float pix_real, float pix_imaginary, int pixel, double *lut_value, double *lut_valueDB, double *magnitude, double *sigma0, char *bandNumber);
void CPL_DLL CPL_STDCALL CalculateMagnitudeLutDB(GDALDataset *dst, float pix, int pixel, double *lut_value, double *lut_valueDB, double *magnitude, char *bandNumber);
void CPL_DLL CPL_STDCALL SetRasterDataLUTPartial(GDALDataset *dst, int pixel_offset, int pixel_width, char *bandNumber);
double CPL_DLL * CPL_STDCALL InterpolateValues(char **papszList, int tableSize, int stepSize, int numberOfValues, int pixelFirstLutValue);
/* End: Roberto July, 2018 */

#endif /* ndef GDAL_LUT_H_INCLUDED */
This package was developed by Defence Research and Development Canada (DRDC) Ottawa Research Centre, in 2020.
It contains RCM GDAL Raster driver, and adds RCM to GDAL supported Raster Formats.
  
This driver was compiled as part of the GDAL-2.4.4
The cpp and header files are stored in gdal-2.4.4\frmts\rcm\, which also contains 'makefile.vc' and 'GNUmakefile' for Windows and Linux build, respectively.

The gdal-2.4.4\frmts\rs2\ and gdal-2.4.4\frmts\nitf\  cpp and header files are also provided as RS2 SWIG functions are included in the gdalrasterband.cpp.

Files in this package:
  gdal-2.4.4\frmts\gdalallregister.cpp  add '#ifdef  FRMT_rcm GDALRegister_RCM();  #endif'
  gdal-2.4.4\frmts\makefile.vc          Windows makefile: add -DFRMT_rcm in EXTRAFLAGS
  gdal-2.4.4\frmts\formats_list.html    add an entry for RCM
  
  gdal-2.4.4\frmts\rcm\rcmdataset.cpp   RCM C++ driver
  gdal-2.4.4\frmts\rcm\rcmdataset.h     RCM header
  gdal-2.4.4\frmts\rcm\makefile.vc      Windows makefile
  gdal-2.4.4\frmts\rcm\GNUmakefile      Linux makefile
  gdal-2.4.4\frmts\rcm\frmt_rcm.html    RCM format HTML
  
  gdal-2.4.4\frmts\rs2\rs2dataset.cpp   (modified) RS2 C++ driver
  gdal-2.4.4\frmts\rs2\rs2dataset.h     (modified) RS2 header
  gdal-2.4.4\frmts\rs2\makefile.vc      Windows makefile
  gdal-2.4.4\frmts\rs2\GNUmakefile      Linux makefile
  gdal-2.4.4\frmts\rs2\frmt_rs2.html    RS2 format HTML
  
  gdal-2.4.4\frmts\nitf\nitfdataset.cpp     NITF C++ driver
  gdal-2.4.4\frmts\nitf\nitfdataset.h       NITF header
  gdal-2.4.4\frmts\nitf\nitfrasterband.cpp  NITF C++
  
  gdal-2.4.4\gcore\makefile.vc          Windows makefile: add gdal_io_error.obj
  gdal-2.4.4\gcore\GNUmakefile          Linux makefile: add gdal_io_error.o
  gdal-2.4.4\gcore\gdal.h               add SWIG functions for RS2 and RCM
  gdal-2.4.4\gcore\gdal_frmts.h         add void CPL_DLL GDALRegister_RCM(void);
  gdal-2.4.4\gcore\gdal_io_error.h      new header file for dubugging RS2 and RCM
  gdal-2.4.4\gcore\gdal_lut.h           new header file for RS2 and RCM LUTs
  gdal-2.4.4\gcore\gdal_pam.h           moved eCalibration definition from rs2dataset.cpp and rcmdataset.cpp
  gdal-2.4.4\gcore\gdal_priv.h          add functions: DeleteOneBand(), DeleteAllBands()
  gdal-2.4.4\gcore\gdal_io_error.cpp    new C++ file for dubugging RS2 and RCM
  gdal-2.4.4\gcore\gdaldataset.cpp      add SWIG functions for RS2 and RCM
  gdal-2.4.4\gcore\gdalrasterband.cpp   add SWIG functions for RS2 and RCM
  gdal-2.4.4\gcore\gdalarraybandblockcache.cpp  one change made in AdoptBlock()
  
The Linux user is required to edit the following file which comes with GDAL:
gdal-2.4.4\GDALmake.opt                 after running the command ./configure, add a line 'GDAL_FORMATS += rcm' towards the end of the file, following all the other lines of GDAL_FORMATS statements

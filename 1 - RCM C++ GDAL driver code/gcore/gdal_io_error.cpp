#include <time.h>
#include <stdio.h>
//#include <conio.h>
#include <string.h>
#include "gdal_pam.h"
#include "gdal_io_error.h"

char *get_current_time()
{
	time_t current_time;
	char* c_time_string;

	/* Obtain current time. */
	current_time = time(NULL);

	/* Convert to local time format. */
	c_time_string = ctime(&current_time);

	for (int i = 0; i < strlen(c_time_string); i++)
	{
		if (c_time_string[i] == '\n' || c_time_string[i] == '\r') c_time_string[i] = ' ';
	}

	return c_time_string;
}

void write_to_file_error(const char *header, const char *value)
{
	/* Convert to local time format. */
	//char *c_time_string = get_current_time();

	//FILE *pFile2 = fopen("gdal_error.log", "a");
	//if (value == NULL)
	//{
	//	fprintf(pFile2, "%s: %s\n", c_time_string, header);
	//}
	//else
	//{
	//	fprintf(pFile2, "%s: %s %s\n", c_time_string, header, value);
	//}
	//fclose(pFile2);
}



void write_to_file(const char *header, const char *value)
{
	/* Convert to local time format. */
	char *c_time_string = get_current_time();
	//CPLDebug(header, value);

	FILE *pFile2 = fopen("gdal_trace.log", "a");
	fprintf(pFile2, "%s: %s %f\n", c_time_string, header, value);
	fclose(pFile2);
}

#ifdef _TRACE_RCM

void write_to_file_dbl(const char *header, const double value)
{
	///* Convert to local time format. */
	//char *c_time_string = get_current_time();

	//char str_value[5000];
	//sprintf(str_value, "%f", value);
	//CPLDebug(header, str_value);

	//FILE *pFile2 = fopen("gdal_trace.log", "a");
	//fprintf(pFile2, "%s: %s %f\n", c_time_string, header, value);
	//fclose(pFile2);
}

void write_to_file_flt(const char *header, const float value)
{
	///* Convert to local time format. */
	//char *c_time_string = get_current_time();

	//char str_value[5000];
	//sprintf(str_value, "%f", value);
	//CPLDebug(header, str_value);

	//FILE *pFile2 = fopen("gdal_trace.log", "a");
	//fprintf(pFile2, "%s: %s %f\n", c_time_string, header, value);
	//fclose(pFile2);
}


void write_to_file_int(const char *header, const int value)
{
	///* Convert to local time format. */
	//char *c_time_string = get_current_time();

	//char str_value[5000];
	//sprintf(str_value, "%d", value);
	//CPLDebug(header, str_value);

	//FILE *pFile2 = fopen("gdal_trace.log", "a");
	//fprintf(pFile2, "%s: %s %d\n", c_time_string, header, value);
	//fclose(pFile2);
}
#endif

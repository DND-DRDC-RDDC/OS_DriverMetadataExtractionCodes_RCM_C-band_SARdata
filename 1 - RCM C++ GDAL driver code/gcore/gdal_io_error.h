#define _NO_TRACE_RCM

char *get_current_time();

void write_to_file_error(const char *header, const char *value);

void write_to_file(const char *header, const char *value);


#ifdef _TRACE_RCM

void write_to_file_dbl(const char *header, const double value);

void write_to_file_flt(const char *header, const float value);

void write_to_file_int(const char *header, const int value);

#endif

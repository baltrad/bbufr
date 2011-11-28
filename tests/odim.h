/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

/* Interface for the zlib library */
#include <zlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Interface for the HDF5 library */
#include <hdf5.h>

#include "desc.h"

extern int _bufr_edition;

unsigned char * my_compress(varfl * buf, int n, int * ncomp);
varfl * my_decompress(unsigned char * buf, unsigned long n, int * ndecomp);
char * set_fuseau(char * nouveau);
time_t build_time(char * date, char * time);
int split(char * args, char * delim, char* subdelim);
char * getsplit_typ(char * args, char * delim, char * subdelim, int i0);
unsigned char * getsplit_val(char * args, char * delim, char * subdelim, int i0);
herr_t H5Acreatestring(hid_t root_id, char * name, char * s);
char * H5Areadstring(hid_t root_id, char * name);
void header_dump(sect_1_t *s1);

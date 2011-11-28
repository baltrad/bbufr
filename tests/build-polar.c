/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Interface for the PROJ library of geographic projections */
#include <proj_api.h>

/* Interface definition of the OPERA coder/decoder */
#include "desc.h"
#include "bufr.h"
#include "bitio.h"
#include "rlenc.h"

#define MAX_DESC 500               /* maximum number of descriptors in a message */

static varfl *global_val;          /* array of values corresponding to a given BUFR descriptor;
                                      must be big enouhg to hold the pixels of an image */
static int global_size;            /* current size of 'global_val' */
static int mindex = 0;             /* index to next uninitialized item of 'global_val' array */

static char *table_dir = NULL;     /* directory for BUFR tables */

static int get_lens (char *buf, long len, size_t secl[6])
/* This function reads from a bufr-message the length of data- and
   data-descriptor-section. Therefore the buffer is opened as a bitstream
   and data is read.

   parameters:
   BUF:   Memory-area containing the BUFR-message.
   LEN:   Number of bytes of the complete BUFR-message determined from
          the length ob the input-file.
   SECL:  Array containing the lengths of the BUFR-sections.

   The return-value is 1 on success, 0 on a fault.
*/

{
  int h, co, i, totlen, lens0;
  unsigned long l;
  long sum;

/******* The length of section 0 is constant, but get the length of the
         whole BUFR message */

  h = bitio_i_open (buf, 8);
  bitio_i_input (h, &l, 32);        /* skip that 'BUFR' */
  bitio_i_input (h, &l, 24);        /* length of whole message */
  bitio_i_close (h);
  lens0 = l;

  secl[0] = 8;
  co = 8;
  sum = 8;

/******* length of section 1 */

  h = bitio_i_open (buf + co, 20);
  if (h == -1) return 0;
  bitio_i_input (h, &l, 24);
  secl[1] = (size_t) l;
  co += secl[1];
  bitio_i_close (h);
  sum += l;
  if (sum > len) goto err;

/******* there is no section 2 */

  secl[2] = 0;

/******* length of section 3 */

  h = bitio_i_open (buf + co, 20);
  if (h == -1) return 0;
  bitio_i_input (h, &l, 24);
  secl[3] = (size_t) l;
  co += secl[3];
  bitio_i_close (h);
  sum += l;
  if (sum > len) goto err;

/******* length of section 4 */

  h = bitio_i_open (buf + co, 20);
  if (h == -1) return 0;
  bitio_i_input (h, &l, 24);
  secl[4] = (size_t) l;
  co += secl[4];
  bitio_i_close (h);
  sum += l;
  if (sum > len) goto err;

/******* length of section 5 is constant */

  secl[5] = 4;
  sum += 4;
  if (sum > len) goto err;

/******* Check the total length of the message against the sum of the lengths 
         of the sections. */

  totlen = 0;
  for (i = 0; i < 6; i ++) {
    totlen += secl[i];
  }
  if (totlen != lens0) {
    fprintf (stderr, "WARNING: Total length of message doesn't match with the lengths\n"
                     "of the individual sections !\n");
  }

  return 1;

/******* Lengths of BUFR-sections not correct */

err:
  fprintf (stderr, "FATAL: Lengths of BUFR-sections > size of input-file !\n");
  return 0;

}

static int ascii_out (varfl val, int ind)
/* This function outputs one character of an ASCII string

   parameters:
   VAL:    Data-value to be output.
   IND:    Index to the global array DES[] holding the description of
           known data-descriptors.

   The function returns 1 on success, 0 on a fault.
*/

{
  fprintf (stderr, "%c", (int) val);
  return 1;
}

static void array_check(int index, int *global_size, varfl **global_val, char * info) 
/* This function checks if index is valid for array 'global_val' */
{
  if (index >= *global_size) {
    *global_size += (200000 > (index - *global_size) ? 200000 : index - *global_size);
    *global_val = realloc(*global_val, *global_size*sizeof(varfl));
    if (!*global_val) {
      fprintf (stderr, "FATAL: Unable to reallocate bytes for %s to hold the output BUFR message !\n", info);
      exit (EXIT_FAILURE);
    }
    fprintf(stderr, "Reallocated memory for %s to size %d\n", info, *global_size);
  } 
  return;
}

static int value_out (varfl val, int ind)
/* This function outputs one descriptor

   parameters:
   VAL:    Data-value to be output.
   IND:    Index to the global array DES[] holding the description of
           known data-descriptors.

   The function returns 1 on success, 0 on a fault.
*/

{
  array_check(mindex, &global_size, &global_val, "global_val");
  global_val[mindex++] = val;
  return 1;
}

static int null_out (varfl val, int ind)
/* This function outputs nothing

   parameters:
   VAL:    Data-value to be output.
   IND:    Index to the global array DES[] holding the description of
           known data-descriptors.

   The function returns 1 on success, 0 on a fault.
*/

{
  return 1;
}

/*===========================================================================*/
static void decode_sec_1 (char *data, size_t datalen, sect_1_t *s1)
/* Fills the section size structure s1 of the message pointed to by data */
{
  unsigned long l;
  int h = bitio_i_open (data, datalen);

  if (h == -1) {
    fprintf (stderr, "FATAL: Unable to read section 1 !\n");
    exit (EXIT_FAILURE);
  }

  bitio_i_input (h, &l, 24);                 /* length of section */

  bitio_i_input (h, &l, 8);  s1->mtab = l;    /* master table used */
  bitio_i_input (h, &l, 8);  s1->subcent = l;
  bitio_i_input (h, &l, 8);  s1->gencent = l;
  bitio_i_input (h, &l, 8);  s1->updsequ = l; /* original BUFR message */
  bitio_i_input (h, &l, 8);  s1->opsec = l;   /* no optional section */
  bitio_i_input (h, &l, 8);  s1->dcat = l;    /* message type */
  bitio_i_input (h, &l, 8);  s1->dcatst = l;  /* message subtype */
  bitio_i_input (h, &l, 8);  s1->vmtab = l;   /* version number of master table used */
  bitio_i_input (h, &l, 8);  s1->vltab = l;   /* version number of local table used */
  bitio_i_input (h, &l, 8);  s1->year = l;    /* year */
  bitio_i_input (h, &l, 8);  s1->mon = l;     /* month */
  bitio_i_input (h, &l, 8);  s1->day = l;     /* day */
  bitio_i_input (h, &l, 8);  s1->hour = l;    /* hour */
  bitio_i_input (h, &l, 8);  s1->min = l;     /* minute */
  bitio_i_close (h);

  _bufr_edition = 2;

  /******* Write all that stuff to an ASCII file */

  fprintf (stderr, "%5d    master table used                  \n", s1->mtab);
  fprintf (stderr, "%5d    subcenter                          \n", s1->subcent);
  fprintf (stderr, "%5d    generating center                  \n", s1->gencent);
  fprintf (stderr, "%5d    original BUFR message              \n", s1->updsequ);
  fprintf (stderr, "%5d    no optional section                \n", s1->opsec);
  fprintf (stderr, "%5d    message type                       \n", s1->dcat);
  fprintf (stderr, "%5d    message subtype                    \n", s1->dcatst);
  fprintf (stderr, "%5d    version number of master table used\n", s1->vmtab);
  fprintf (stderr, "%5d    version number of local table used \n", s1->vltab);
  fprintf (stderr, "%5d    year                               \n", s1->year);
  fprintf (stderr, "%5d    month                              \n", s1->mon);
  fprintf (stderr, "%5d    day                                \n", s1->day);
  fprintf (stderr, "%5d    hour                               \n", s1->hour);
  fprintf (stderr, "%5d    minute                             \n", s1->min);

  return;
}

static int test_fxy(dd *d, int ff, int xx, int yy)
/* Tests equality of descriptor d with (ff,xx,yy) */
{
  return (d->f == ff) && (d->x == xx) && (d->y == yy);
}

static int transcode_it(sect_1_t *s1_in, dd *descr, int iout, varfl *v, char *outfile)
/* Create a new OPERA BUFR message */
{
  int i, ret;
  char *sec[6];      /* 6 sections in BUFR-message */
  size_t secl[6];    /* length of 6 sections */
  sect_1_t s1;       /* Here we store section 1 of BUFR message */


  /****** Initialize basic data for output */

  for (i = 0; i < 6; i ++) sec[i] = NULL;

  /******* Prepare data for Section 1 (output) */

  s1.year = s1_in->year;
  s1.mon  = s1_in->mon;
  s1.day = s1_in->day;
  s1.hour = s1_in->hour;
  s1.min  = s1_in->min;
  s1.mtab = s1_in->mtab;                      /* master table used */
  s1.subcent = 255;
  s1.gencent = 255;
  s1.updsequ = s1_in->updsequ;                   /* original BUFR message */
  s1.opsec = s1_in->opsec;                     /* no optional section */
  s1.dcat = s1_in->dcat;                      /* message type */
  s1.dcatst = s1_in->dcatst;                    /* message subtype */
  s1.vmtab = 13;                     /* version number of master table used */
  s1.vltab = 6;                     /* version number of local table used */

  /******* Write all that stuff to an ASCII file */

  fprintf (stderr, "%5d    master table used                  \n", s1.mtab);
  fprintf (stderr, "%5d    subcenter                          \n", s1.subcent);
  fprintf (stderr, "%5d    generating center                  \n", s1.gencent);
  fprintf (stderr, "%5d    original BUFR message              \n", s1.updsequ);
  fprintf (stderr, "%5d    no optional section                \n", s1.opsec);
  fprintf (stderr, "%5d    message type                       \n", s1.dcat);
  fprintf (stderr, "%5d    message subtype                    \n", s1.dcatst);
  fprintf (stderr, "%5d    version number of master table used\n", s1.vmtab);
  fprintf (stderr, "%5d    version number of local table used \n", s1.vltab);
  fprintf (stderr, "%5d    year                               \n", s1.year);
  fprintf (stderr, "%5d    month                              \n", s1.mon);
  fprintf (stderr, "%5d    day                                \n", s1.day);
  fprintf (stderr, "%5d    hour                               \n", s1.hour);
  fprintf (stderr, "%5d    minute                             \n", s1.min);

  /******* read supported data descriptors */
  /******* parameters are the table-directory (NULL to search in current dirctory), 
	   version of mtab, version of ltab, orcenter */

  if (read_tables(table_dir, s1.vmtab, s1.vltab, s1.subcent, s1.gencent)) {
    fprintf (stderr, "FATAL: Unable to read tables\n");
    exit (EXIT_FAILURE);
  }

  /******* Code the data (section 3 and 4) */

  ret = bufr_create_msg (descr, iout, v, (void *)&sec[4], (void *)&sec[3], &secl[4], &secl[3]);

  /******* Setup section 0, 1, 2, 5 */

  if (!setup_sec0125 (sec, secl, s1)) {
    fprintf (stderr, "WARNING: Unable to create section 0, 1, 2 and/or 5\n");
    return 0;
  }
    
  /******* Save coded data */

  if (!save_sections (sec, secl, outfile)) return 0;

  free(v);
  free(global_val);
  return 1;
}

#define fill_desc(ff,xx,yy) if (iout > MAX_DESC) { \
    fprintf (stderr, "FATAL: Increase the MAX_DESC value\n");\
  } else {\
     descr[iout].f = ff; descr[iout].x = xx; descr[iout++].y = yy;\
     fprintf (stderr, "F X Y output: %d %d %d\n", ff, xx, yy);\
  }

#define fill_v(val) (array_check(jout, &v_size, &v, "v"), v[jout++] = (val))

static int decode_it (char *dad, size_t dadlen, char *da, 
		      size_t dalen, char *outfile, sect_1_t *s1_in)
/* This function decodes a BUFR-message consisting of a data-descriptor-section
   and a data-section and saves the data we look for

   parameters:
   DAD:      Is where the data-descriptor-section is stored.
   DADLEN:   Number of bytes of the data-descriptor-section.
   DA:       Is where the data-section is stored.
   DALEN:    Number of bytes of the data-section.

   The function returns 1 on success, 0 on a fault.
*/

{
  dd *dds, *d;
  int ndds, i, ii, iout, nrows = MISSVAL, ncols = MISSVAL, v_size;
  unsigned int vali, jout;
  size_t nvals;
  varfl *vals;
  char *unit;

  varfl    region = MISSVAL, station = MISSVAL;
  varfl    lat_rad = MISSVAL, lon_rad = MISSVAL, hauteur = MISSVAL;

  varfl    an_3_1_11 = MISSVAL;
  varfl    mois_3_1_11 = MISSVAL;
  varfl    jour_3_1_11 = MISSVAL;
  varfl    heure_3_1_12 = MISSVAL;
  varfl    min_3_1_12 = MISSVAL;
  varfl    dx = MISSVAL, dy = MISSVAL, rotation = MISSVAL;
  varfl    picture_type = MISSVAL, coord_type = MISSVAL, proj_type = MISSVAL;
  varfl    lat_NW = MISSVAL, lon_NW = MISSVAL, lat_orig = MISSVAL, lon_ref = MISSVAL;
  varfl    gate_length = MISSVAL, ngates = MISSVAL, npulses = MISSVAL, echos = MISSVAL;

  int nseuils;
  varfl * seuils = NULL;

  int nombre_seuils = 0;
  varfl * pp = NULL, seuil0 = MISSVAL;
  unsigned char *pixmap = NULL;

  int nscans = 0;
  varfl * scans_info = NULL;

  dd descr[MAX_DESC]; /* This array must be huge enough to hold all required descriptors */
  varfl *v;           /* This array must be huge enough to hold all corresponding data values */

  v_size = 1000;
  v = calloc(v_size, sizeof(varfl));
  global_size = 500000;
  global_val = calloc(global_size, sizeof(varfl));
  if (!v || !global_val) {
    fprintf (stderr, "FATAL: Unable to allocate bytes to hold the output BUFR message !\n");
    return 0;
  }

  i = bufr_read_msg (da, dad, dalen, dadlen, &dds, &ndds, &vals, &nvals);
  if (!i) {
    fprintf (stderr, "FATAL: Unable to read the input BUFR message !\n");
    return 0;
  }

  /* The data we are interested in is read to array 'v' 
     and the corresponding descriptors are in array 'descr' */

  vali = 0; i = 0;

  /* Iteration through the input data */
  while (i < ndds) {
    d = &dds[i];
    unit = get_unit (d);
    if (test_fxy(d, 3, 1, 1)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      region = global_val[0];
      station = global_val[1];
      mindex = 0;
    } else if (test_fxy(d, 3, 1, 11)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      an_3_1_11 = global_val[0];
      mois_3_1_11 = global_val[1];
      jour_3_1_11 = global_val[2];
      mindex = 0;
    } else if (test_fxy(d, 3, 1, 13)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      heure_3_1_12 = global_val[0];
      min_3_1_12 = global_val[1];
      mindex = 0;
    } else if (test_fxy(d, 3, 1, 21)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      lat_rad = global_val[0];
      lon_rad = global_val[1];
      mindex = 0;
    } else if (test_fxy(d, 0, 7, 2)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      hauteur = global_val[0]; 
      mindex = 0;
    } else if (test_fxy(d, 0, 5, 33)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      dx = global_val[0]; 
      mindex = 0;
    } else if (test_fxy(d, 0, 6, 33)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      dy = global_val[0]; 
      mindex = 0;
    } else if (test_fxy(d, 0, 30, 21)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      ncols = global_val[0]; 
      mindex = 0;
    } else if (test_fxy(d, 0, 30, 22)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      nrows = global_val[0]; 
      mindex = 0;
    } else if (test_fxy(d, 3, 21, 10)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      hauteur = global_val[1]; /* 0 7 2 */
      mindex = 0;
    } else if (test_fxy(d, 3, 21, 11)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      picture_type = global_val[0]; /* 0 30 31 */
      coord_type = global_val[2]; /* 0 29 2 */
      mindex = 0;
    } else if (test_fxy(d, 0, 2, 202)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      rotation = global_val[0];
      mindex = 0;
    } else if (test_fxy(d, 0, 30, 31)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      picture_type = global_val[0];
      mindex = 0;
    } else if (test_fxy(d, 0, 29, 2)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      coord_type = global_val[0];
      mindex = 0;
    } else if (test_fxy(d, 0, 29, 1)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      proj_type = global_val[0];
      mindex = 0;
    } else if (test_fxy(d, 3, 21, 6)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      gate_length = global_val[0];
      ngates = global_val[1];
      npulses = global_val[2];
      echos = global_val[3];
      mindex = 0;
    } else if (test_fxy(d, 3, 29, 192)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      lat_NW = global_val[0];
      lon_NW = global_val[1];
      lat_orig = global_val[2];
      lon_ref = global_val[5];
      mindex = 0;
    } else if (test_fxy(d, 1, 6, 0)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      nseuils = global_val[0];
      seuils = calloc(11*nseuils, sizeof(varfl));
      memcpy(seuils, global_val+1, 11*nseuils*sizeof(varfl));
      mindex = 0;
    } else if (test_fxy(d, 1, 4, 0)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      nscans = global_val[0];
      scans_info = calloc(4*nscans, sizeof(varfl));
      memcpy(scans_info, global_val+1, 4*nscans*sizeof(varfl));
      mindex = 0;
    } else if (test_fxy(d, 3, 21, 193)) {
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      seuil0 = global_val[3];
      nombre_seuils = global_val[0] - 1;;
      pp = calloc(3*nombre_seuils, sizeof(varfl));
      memcpy(pp, global_val+4, 3*nombre_seuils*sizeof(varfl));
      mindex = 0;
    } else if (test_fxy(d, 1, 3, 0)) {
      int ii, jj, valiold = vali;
      bufr_parse (d, 0, 0, vals, &vali, value_out);
      if (vali - valiold - 1 > 1000) {
	/* horrible heuristic test used to detect the pixels 
	   array position in the French radar data */
	fprintf(stderr, "ncols= %d nrows= %d\n", ncols, nrows);
	fprintf(stderr, "taille vraie: %d taille calculee: %d\n", (vali - valiold - 1), ncols*nrows);
	fprintf(stderr, "taille: %d\n", 720*256);
	nrows = 256;
	ncols = 720;
	pixmap = calloc(ncols*nrows, sizeof(unsigned char));
	for (ii = 0; ii < nrows; ii++) {
	  for (jj = 0; jj < ncols; jj++) {
	    pixmap[ii*ncols + jj] = global_val[ii*ncols + jj + 1];
	  }
	}
      }
      mindex = 0;
    } else {

/******* special treatment for ASCII data */

      if (unit != NULL && strcmp (unit, "CCITT IA5") == 0) {
	if (!bufr_parse (d, 0, 0, vals, &vali, ascii_out)) return 0;
      }

/******* "add associated field" descriptor */

      else if (d->f == 2 && d->x == 4) {
	if (!bufr_parse (d, 0, 0, vals, &vali, ascii_out)) return 0;
      }

/******* "ordinary" data */

      else {
	if (!bufr_parse (d, 0, 0, vals, &vali, null_out)) return 0;
      }

    }
/******* If there is a replication factor -> increase i by the number of 
         descriptors to be replicated */

    if      (d->f == 1 && d->y == 0) i += d->x + 2;
    else if (d->f == 1 && d->y != 0) i += d->x + 1;
    else                             i ++;

  }

  /* Output phase */

  iout = 0; jout = 0;

  {
    fill_desc(3,1,31);
      fill_v(region); 
      fill_v(station);

      fill_v(0); /* automatic station */
  
      fill_v(an_3_1_11); 
      fill_v(mois_3_1_11);
      fill_v(jour_3_1_11);

      fill_v(heure_3_1_12);
      fill_v(min_3_1_12);
  
      /* fill_desc(0,5,1); */
      fill_v(lat_rad); 

      /* fill_desc(0,6,1); */
      fill_v(lon_rad); 

      /* fill_desc(0,7,1); */
      fill_v(hauteur); 

    fill_desc(0,30,31);
    fill_v(picture_type);

    fill_desc(0,29,2);
    fill_v(coord_type); 

    fill_desc(0,30,196);
    fill_v(4);

    fill_desc(0,33,3);
    fill_v(0); 
  }

  if (pp) {
    fill_desc(2,1,129);
    fill_desc(3,13,9);
    fill_v(seuil0);
    fill_v(nombre_seuils);
    for (ii = 0; ii < nombre_seuils; ++ii) {
      fill_v(pp[3*ii+2]);
    }
    fill_desc(2,1,0);
    free(pp);
  }

  fill_desc(3,21,199);

  fill_v(1); /* one scan */

  fill_v(heure_3_1_12); /* may be wrong by one hour */
  if (scans_info && (nscans >= 1)) {
    fill_v(scans_info[2]); 
    fill_v(scans_info[3]);
    fill_v(scans_info[1]);
    fill_v(scans_info[0]);
    free(scans_info);
  } else {
    fprintf (stderr, "FATAL: no polar scan\n");
    exit(1);
  } 

  fill_v(gate_length);
  fill_v(ngates);
  fill_v(npulses);
  fill_v(echos);     

  fill_v(gate_length*ngates);
  fill_v(0.5);

  fill_v(rotation);

  fprintf(stderr, "ngates= %g npulses= %g\n", ngates, npulses);
  nrows = 256;
  fill_v(nrows);
  ncols = 720;
  fill_v(ncols);

  {
    int jj;
    unsigned char *p = calloc(ncols, sizeof(unsigned char));
    varfl *ptr = NULL;
    size_t count = 0;
    for (ii = 0; ii < nrows; ii++) {
      for (jj = 0; jj < ncols; jj++) {
	p[jj] = pixmap[ii*ncols + jj];
      }
      if (!rlenc_compress_line (ii, p, ncols, &ptr, &count)) {
	fprintf (stderr, "WARNING: Error coding RLL line %d\n", ii);
      }
      /*fprintf (stderr, "line %d %d ok\n", ii, count);*/
    }
    fill_v(nrows);
    for (ii = 0; ii < count; ii++) {
      fill_v(ptr[ii]);  
    }
    free(p);
    free(ptr);
    free(pixmap);
  }

  free (dds);
  free (vals);
  free_descs();

  /* Now building new transcoded BUFR message */

  fprintf (stderr, "Output file header:\n");
  return transcode_it(s1_in, descr, iout, v, outfile);

}

/*===========================================================================
                                  Main program
  ===========================================================================*/
int main (int argc, char *argv[])
{
  char *usage = "usage: transcode [-v] [-d tabdir] input_file output_file\n";
  char *version = "transcode V2.3, 02-09-2005\n";

  char *buf,*b1 = NULL, *b7777;
  FILE *fp;                        /* File-Pointer for input file */
  long len, co, l;
  int i;
  size_t secl[6];
  sect_1_t s1;
  char *sec[6];

/******* check command line parameter */

  while (argc > 1 && *argv[1] == '-')
  {
    if (*(argv[1] + 1) == 'v')
      fprintf (stderr, "%s", version);
    else if (*(argv[1] + 1) == 'd')
    {
      if (argc < 2)
      {
        fprintf (stderr, "Missing parameter for -d\n\n%s", usage);
        exit (EXIT_FAILURE);
      }
      table_dir = argv[2];
      argc--;
      argv++;
    }
    else
    {
        fprintf (stderr, "Invalid parameter %s\n\n%s", argv[1], usage);
        exit (EXIT_FAILURE);
    }
    argc--;
    argv++;
  }

/******* Get input- and output-filenames from the command-line */

  if (argc != 3)
  {
    fprintf (stderr, "%s", usage);
    exit (EXIT_FAILURE);
  }

/****** read source-file. Therefore allocate memory to hold the complete
        BUFR-message */

  fp = fopen (argv[1], "rb");
  if (fp == NULL) {
    fprintf (stderr, "FATAL: Unable to open input file '%s'\n", argv[1]);
    exit (EXIT_FAILURE);
  }
  fseek (fp, 0L, SEEK_END);
  len = ftell (fp);
  fseek (fp, 0L, SEEK_SET);
  b1 = (char *) malloc ((size_t) len);
  if (b1 == NULL) {
    fprintf (stderr, "FATAL: Unable to allocate %ld bytes to hold BUFR-message !\n", len);
    goto err;
  }
  if (fread (b1, 1, (size_t) len, fp) != (size_t) len) {
    fprintf (stderr, "FATAL: Unable to read the BUFR-message !\n");
    goto err;
  }

/****** Search for "BUFR" */

  buf = NULL;
  for (l = 0; l < len - 4 && buf == NULL; l ++) {
    if (*(b1 + l)     == 'B' && 
        *(b1 + l + 1) == 'U' &&
        *(b1 + l + 2) == 'F' &&
        *(b1 + l + 3) == 'R') buf = b1 + l;
  }
  if (buf == NULL) {
    fprintf (stderr, "FATAL: 'BUFR' not found in BUFR-message !\n");
    goto err;
  }

/****** Check for the ending "7777" */

  b7777 = NULL;
  for (l = 0; l < len - 3 && b7777 == NULL; l ++) {
    if (*(b1 + l)     == '7' && 
        *(b1 + l + 1) == '7' &&
        *(b1 + l + 2) == '7' &&
        *(b1 + l + 3) == '7') b7777 = b1 + l;
  }
  if (b7777 == NULL) {
    fprintf (stderr, "FATAL: '7777' not found in BUFR-message !\n");
    goto err;
  }

/****** Get length of all 6 sections */

  if (!get_lens (buf, len, secl)) {
    fprintf (stderr, "FATAL: Unable to read lengths of BUFR-sections !\n");
    goto err;
  }

/******* Calculate a pointer to the beginning of each section */

  co = 0L;
  for (i = 0; i < 6; i ++) {
    sec[i] = buf + co;
    co += secl[i];
  }

/******* decode section 1 */

  fprintf (stderr, "Input file header:\n");
  decode_sec_1 (sec[1], secl[1], &s1);

/* read descriptor tables */

  if (read_tables (table_dir, s1.vmtab, s1.vltab, s1.subcent, s1.gencent)) {
    fprintf (stderr, "FATAL: Unable to read tables\n");
    goto err;
  }

/******* decode data descriptor- and data-section now */

  if (!decode_it (sec[3], secl[3], sec[4], secl[4], argv[2], &s1)) {
    fprintf (stderr, "FATAL: Unable to decode BUFR-message !\n");
    goto err;
  }
  fclose (fp);
  free (b1);
  free_descs();
  exit (EXIT_SUCCESS);

/******* An error occoured: */

err:
  fclose (fp);
  free (b1);
  free_descs();
  exit (EXIT_FAILURE);
}


/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <float.h>

/* Interface for the PROJ library of geographic projections */
#include <proj_api.h>

/* Interface for the PNG graphic format */
#include <png.h>

/* Interface definition of the OPERA coder/decoder */
#include "desc.h"
#include "bufr.h"
#include "bitio.h"
#include "rlenc.h"

#include "radar.h"
#include "polar.h"
#include "odim.h"

FILE* initPngFile(char *filename, png_structpp png_ptr_ptr, png_infopp info_ptr_ptr){
  FILE *fp = NULL;
  if (filename != NULL) {
    fp = fopen(filename , "wb");
    if (!fp) return fp;
  }
  else fp = stdout;  
  *png_ptr_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr_ptr){
    fclose(fp);
    return NULL;
  }
  *info_ptr_ptr = png_create_info_struct(*png_ptr_ptr);
  if (!info_ptr_ptr){
    png_destroy_write_struct(png_ptr_ptr, (png_infopp)NULL);
    fclose(fp);
    return NULL;
  }
 
  if (setjmp(png_jmpbuf(*png_ptr_ptr))){  /*set longjump on error */
    png_destroy_write_struct(png_ptr_ptr, info_ptr_ptr);
    fclose(fp);
    return NULL;
  }
  png_init_io(*png_ptr_ptr, fp);
 
  return fp;
}

void deinit(FILE** fp, png_structpp png_ptr_ptr, png_infopp info_ptr_ptr){
  png_write_end(*png_ptr_ptr, *info_ptr_ptr);
  png_destroy_write_struct(png_ptr_ptr, info_ptr_ptr);
  fclose(*fp);
}

void writeData(png_structp png_ptr,png_infop info_ptr,odim_polar_t* od){
  varfl* p = od->datasets[0].dataset_data[0].data;
  
  int height = od->datasets[0].dataset_where.nrays;
  int width = od->datasets[0].dataset_where.nbins;

  unsigned char * row = calloc(height*width, sizeof(unsigned char));
  int i;

  for (i = 0; i < height*width; i++) {
    if (p[i] == od->datasets[0].dataset_data[0].dataset_what.nodata) {
      row[i] = 255;
    } else if (p[i] < 0) {
      row[i] = 0;
    } else {
      row[i] = p[i];
    }
  }
  
  png_bytep row_pointer = row;
  int count = 0;

  for (count = 0; count < height; count++){
    png_write_row(png_ptr, row_pointer);
    row_pointer += width;
  }
  free(row);
}

void writeInfo(png_structp png_ptr, png_infop info_ptr, odim_polar_t* od){
  int count;
  int size_palette;
  png_color palette[256] = {{255, 255, 255}, {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, /* "White" */
			    {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, 
			    {72,  61, 139}, {72,  61, 139}, {72,  61, 139}, {72,  61, 139},     /* "DarkSlateBlue" */
			    {72,  61, 139}, {72,  61, 139}, {72,  61, 139}, {72,  61, 139}, 
			    {0, 0, 205}, {0, 0, 205}, {0, 0, 205}, {0, 0, 205},                 /* "MediumBlue" */
			    {30, 144, 255}, {30, 144, 255}, {30, 144, 255}, {30, 144, 255},     /* "DodgerBlue" */
			    {173, 216, 230}, {173, 216, 230}, {173, 216, 230}, {173, 216, 230}, /* "LightBlue" */
			    {85, 107,  47}, {85, 107,  47}, {85, 107,  47}, {85, 107,  47},     /* "DarkOliveGreen" */
			    {60, 179, 113}, {60, 179, 113}, {60, 179, 113}, {60, 179, 113},     /* "MediumSeaGreen" */
			    {127, 255, 212}, {127, 255, 212}, {127, 255, 212}, {127, 255, 212}, /* "Aquamarine" */
			    {127, 255,   0}, {127, 255,   0}, {127, 255,   0}, {127, 255,   0}, /* "Chartreuse" */
			    {255, 255,   0}, {255, 255,   0}, {255, 255,   0}, {255, 255,   0}, /* "Yellow" */
			    {240, 230, 140}, {240, 230, 140}, {240, 230, 140}, {240, 230, 140}, /* "khaki" */
			    {222, 184, 135}, {222, 184, 135}, {222, 184, 135}, {222, 184, 135}, /* "burlywood" */
			    {255, 165,   0}, {255, 165,   0}, {255, 165,   0}, {255, 165,   0}, /* "Orange" */
			    {160,  82,  45}, {160,  82,  45}, {160,  82,  45}, {160,  82,  45}, /* "sienna" */
			    {255,   0,   0}};                                                   /* "Red" */
  png_byte trans[256];
  png_color no_value = {255, 192, 203}; /* Pink */

  png_set_invert_alpha(png_ptr);
  png_set_IHDR(png_ptr, info_ptr,od->datasets[0].dataset_where.nbins, od->datasets[0].dataset_where.nrays, 8,
      PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  size_palette = 256; /*(image->depth == PIXEL4BIT) ?  16 : 256;*/
    
  for (count = 64; count < size_palette; count++){
    palette[count].red = 255;
    palette[count].green = 0;
    palette[count].blue = 0;
/*     palette[count].blue = 256*count/size_palette; */
/*     palette[count].red = 256*(size_palette-count-1)/size_palette; */
/*     palette[count].green = 256*count/size_palette; */
  }
  for (count = 0; count < size_palette; count++) 
    trans[count] = 0;
  palette[size_palette-1] = no_value;
  trans[size_palette-1] = 255;
  png_set_tRNS(png_ptr, info_ptr, trans, size_palette, NULL);
  png_set_PLTE(png_ptr, info_ptr, palette, size_palette);
  png_write_info(png_ptr, info_ptr);
}

odim_polar_t odim_data; /* structure holding ODIM data */

static int write_hdf5(odim_polar_t * od, char * file)
{
  hid_t       file_id, dataset_id, dataspace_id, root_id, group_id, attr_id, space_id;  /* identifiers */
  herr_t      status;
  int i, j;

  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;

  FILE* fp=initPngFile(file, &png_ptr, &info_ptr); 
  if (!fp) return 1;

  writeInfo(png_ptr, info_ptr, od);
  writeData(png_ptr, info_ptr, od);
  deinit(&fp, &png_ptr, &info_ptr);

  return 1;
}

static int our_callback (varfl val, int ind) {

    odim_polar_t * od = &odim_data;   /* our global data structure */

    /* do nothing if data modification descriptor or replication descriptor */
    if (ind == _desc_special) {
      if (des[ind]->el->d.f != 2) fprintf (stderr,
	       " fxy: %d %d %d\n", des[ind]->el->d.f, des[ind]->el->d.x, des[ind]->el->d.y);
      return 1;
    }

    /* sequence descriptor */

    if (des[ind]->id == SEQDESC) {
      /* get descriptor */
      dd *d = &(des[ind]->seq->d);
      varfl *vv;

      /* open array for values */
      bufrval_t* v = bufr_open_val_array ();
      if (v == NULL) return 0;
      
      /* decode sequence to global array */
      if (!bufr_parse_out (des[ind]->seq->del, 0, des[ind]->seq->nel - 1,
			   bufr_val_to_global, 0)) {
	bufr_close_val_array ();
	return 0;
      }
      vv = v->vals;

      /* WMO block and station number */
        if (bufr_check_fxy (d, 3,1,31)) { 
	  radar_data_t b;
	  int i = 0;
	  {
	    varfl tempo = vv[i++];
	    /* fprintf (stderr, "OMM: %g %d %d\n", tempo, (tempo == MISSVAL), od->what.nstations);*/
	    if (tempo == MISSVAL) {
	      if (od->what.nstations) od->what.nstations--; /* we supposed WMO station was not missing */
	      assert(vv[i++] == MISSVAL);
	    } else {
	      b.wmoblock = tempo;
	      b.wmostat = vv[i++];
	      assert(b.wmostat != MISSVAL);
	      od->what.source[od->what.nstations-1].identifier = "WMO";
	      od->what.source[od->what.nstations-1].value = calloc(6, sizeof(char));
	      sprintf(od->what.source[od->what.nstations-1].value, "%02d%03d", b.wmoblock, b.wmostat);
	    }
	  }
	    i++;
            b.meta.year = vv[i++];       /* Date */
            b.meta.month = vv[i++];
            b.meta.day = vv[i++];
            b.meta.hour = vv[i++];       /* Time */
            b.meta.min = vv[i++];
	    b.meta.radar.lat = vv[i++];
            b.meta.radar.lon = vv[i++];
            b.meta.radar_height = vv[i++];
 	  /* what */
	  od->what.object = "PVOL";
	  od->what.version = "H5rad 2.0";

	  /* where */
	  od->where.lat = b.meta.radar.lat;
	  od->where.lon = b.meta.radar.lon;
	  od->where.height = b.meta.radar_height;
	}
	/* ODIM */
        else if (bufr_check_fxy (d, 3,21,203)) {
	  int ii, jj, i = 0;
	  od->nscans = vv[i++];
	  od->datasets = calloc(od->nscans, sizeof(odim_polar_dataset_t));
	  for (ii = 0; ii < od->nscans; ++ii) {
	    odim_polar_dataset_t * ds = &od->datasets[ii];
    
	    /* what */
	    ds->dataset_what.product = "SCAN";
	    {
	      char * tz = set_fuseau("TZ=UTC");
	      struct tm local;
	      local.tm_year = vv[i++] - 1900;
	      local.tm_mon = vv[i++] - 1;
	      local.tm_mday = vv[i++];
	      local.tm_hour = vv[i++];
	      local.tm_min = vv[i++];
	      local.tm_sec = vv[i++];
	      local.tm_wday = 0;
	      local.tm_yday = 0;    
	      local.tm_isdst = 0;
	      ds->dataset_what.start_time = mktime(&local);
	      
	      local.tm_year = vv[i++] - 1900;
	      local.tm_mon = vv[i++] - 1;
	      local.tm_mday = vv[i++];
	      local.tm_hour = vv[i++];
	      local.tm_min = vv[i++];
	      local.tm_sec = vv[i++];
	      local.tm_wday = 0;
	      local.tm_yday = 0;    
	      local.tm_isdst = 0;
	      ds->dataset_what.end_time = mktime(&local);

	      set_fuseau(tz);
	    }
	    {
	      varfl val = vv[i++];
	      assert(val == 90);
	    }

	    /* where */
	    ds->dataset_where.elangle = vv[i++];
	    ds->dataset_where.nbins = vv[i++];
	    ds->dataset_where.rscale = vv[i++];
	    ds->dataset_where.rstart = vv[i++];
	    ds->dataset_where.nrays = vv[i++];
	    ds->dataset_where.a1gate = vv[i++];
    
	    /* parameters */
	    ds->nparams = vv[i++];
	    ds->dataset_data = calloc(ds->nparams, sizeof(odim_polar_dataset_data_t));
    
	    for (jj = 0; jj < ds->nparams; ++jj) {
	      odim_polar_dataset_data_t * p = &ds->dataset_data[jj];
	      int param = vv[i++];

	      /* what */
	      if (param == 0) {
		p->dataset_what.quantity = "DBZH";
	      } else if (param == 40) {
		p->dataset_what.quantity = "VRAD";
	      } else if (param == 91) {
		p->dataset_what.quantity = "TH";
	      } else if(param == 92) {
		p->dataset_what.quantity = "WRAD";
	      } else {
		p->dataset_what.quantity = "XXXX";
		fprintf(stderr, "unknown quantity %d\n", param);
	      }
	      fprintf(stderr, "Quantity %s\n" , p->dataset_what.quantity);
	      p->dataset_what.gain = 1;
	      p->dataset_what.offset = 0;
	      p->dataset_what.nodata = DBL_MAX;
	      p->dataset_what.undetect = -DBL_MAX;
      
	      param = vv[i++];
	      assert(param == 0); /* zlib compresssion */
	      
	      /* data */
	      {
		int ndecomp = ds->dataset_where.nbins*ds->dataset_where.nrays*sizeof(varfl);
		int k, n, j = 0, niter = vv[i++];
		unsigned char * tempo = calloc(niter*65534, sizeof(unsigned char));
		for (n = 0; n < niter; ++n) {
		  int ncomp = vv[i++];
		  for (k = 0; k < ncomp; ++k) {
		    tempo[j++] = (vv[i+k] == MISSVAL ? 255 : vv[i+k]);
		  }
		  i += ncomp;
		}
		p->data = my_decompress(tempo, j, &ndecomp);
		free(tempo);
		fprintf (stderr,"--> %g %g %g\n", p->data[0], p->data[1], p->data[2]);
	      }
	    }
	  }
	}
	/* ODIM stations */
        else if (bufr_check_fxy (d, 3,21,204)) {
	  int i = 0, j;
	  od->what.nstations = vv[i++] + 1; /* WMO station comes after */
	  od->what.source = calloc(od->what.nstations, sizeof(odim_station_t));
	  for (j = 0; j < od->what.nstations-1; ++j) {
	    int k; 
	    od->what.source[j].identifier = calloc(4, sizeof(char)); 
	    for (k = 0; k < 3; ++k) {
	      od->what.source[j].identifier[k] = vv[i++];
	    }
	    od->what.source[j].value = calloc(17, sizeof(char));
	    for (k = 0; k < 16; ++k) {
	      od->what.source[j].value[k] = vv[i++];
	    }
	  }
        } else {
            fprintf (stderr,
                     "Unknown sequence descriptor %d %d %d\n", d->f, d->x, d->y);
        }
        /* close the global value array */
        bufr_close_val_array ();
    }
    /* element descriptor */
    else if (des[ind]->id == ELDESC) {
      dd *d = &(des[ind]->el->d);
      fprintf (stderr,
	       "Unknown element descriptor %d %d %d\n", d->f, d->x, d->y);
    }
    return 1;
}

/*===========================================================================
                                  Main program
  ===========================================================================*/
int main (int argc, char *argv[])
{
  char *usage = "usage: transcode [-v] [-d tabdir] input_file output_file\n";
  char *version = "transcode V2.3, 02-10-2008\n";

  char *table_dir = NULL;     /* directory for BUFR tables */
  sect_1_t s1;
  bufr_t msg;

/******* check command line parameter */

  while (argc > 1 && *argv[1] == '-') {
    if (*(argv[1] + 1) == 'v')
      fprintf (stderr, "%s", version);
    else if (*(argv[1] + 1) == 'd') {
      if (argc < 2) {
        fprintf (stderr, "Missing parameter for -d\n\n%s", usage);
        exit (EXIT_FAILURE);
      }
      table_dir = argv[2];
      argc--;
      argv++;
    } else {
        fprintf (stderr, "Invalid parameter %s\n\n%s", argv[1], usage);
        exit (EXIT_FAILURE);
    }
    argc--;
    argv++;
  }

  /******* Get input- and output-filenames from the command-line */
  if (argc != 3) {
    fprintf (stderr, "%s", usage);
    exit (EXIT_FAILURE);
  }

  if (!bufr_read_file(&msg, argv[1])) {
        fprintf (stderr, "FATAL: Unable to read the BUFR-message in %s!\n", argv[1]);
        exit (EXIT_FAILURE);
  }

  /******* decode section 1 */
  fprintf (stderr, "Input file header:\n");
  if (!bufr_decode_sections01(&s1, &msg)) {
    fprintf (stderr, "FATAL: Unable to decode section 1\n");
    exit (EXIT_FAILURE);
  }
  header_dump(&s1);

  {
    char * tz = set_fuseau("TZ=UTC");
    struct tm local;
    local.tm_sec = s1.sec;
    local.tm_min = s1.min;
    local.tm_hour = s1.hour;
    local.tm_mday = s1.day;
    local.tm_mon = s1.mon - 1;
    local.tm_year = s1.year - 1900;
    local.tm_wday = 0;
    local.tm_yday = 0;    
    local.tm_isdst = 0;
    
    odim_data.what.nominal_time = mktime(&local);
    set_fuseau(tz);
  }

  /* read descriptor tables */
  if (read_tables (table_dir, s1.vmtab, s1.vltab, s1.subcent, s1.gencent)) {
    fprintf (stderr, "FATAL: Unable to read tables\n");
    exit (EXIT_FAILURE);
  }

  /* decode data descriptor and data-section now */
  
  {
    int ok, desch, ndescs;
    dd* dds = NULL;

    /* open bitstreams for section 3 and 4 */
    desch = bufr_open_descsec_r(&msg, NULL);
    ok = (desch >= 0);
    if (ok) ok = (bufr_open_datasect_r(&msg) >= 0);

    /* calculate number of data descriptors  */
    ndescs = bufr_get_ndescs (&msg);

    /* allocate memory and read data descriptors from bitstream */
    if (ok) ok = bufr_in_descsec (&dds, ndescs, desch);

    /* output data to our global data structure */
    if (ok) ok = bufr_parse_out (dds, 0, ndescs - 1, our_callback, 1);

    /* fill the HDF5 ODIM data structure */
    if (ok) write_hdf5(&odim_data, argv[2]);

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }

  free_descs();
  exit (EXIT_SUCCESS);
}


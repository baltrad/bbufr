/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <float.h>

/* Interface for the PROJ library of geographic projections */
#include <proj_api.h>

/* Interface definition of the OPERA coder/decoder */
#include "desc.h"
#include "bufr.h"
#include "bitio.h"
#include "rlenc.h"

#include "radar.h"
#include "compo.h"
#include "odim.h"

odim_comp_t odim_data; /* structure holding ODIM data */

static int write_hdf5(odim_comp_t * od, char * file)
{
  hid_t       file_id, dataset_id, dataspace_id, root_id, group_id, attr_id, space_id;  /* identifiers */
  herr_t      status;
  int i, j;

  /* Create a new file using default properties. */
  file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/");

  status = H5Acreatestring(root_id, "Conventions", Conventions);

  group_id = H5Gcreate(root_id, "what", 0);
  status = H5Acreatestring(group_id, "object", od->what.object);
  status = H5Acreatestring(group_id, "version", od->what.version);
  {
    char str[100];
    int i;
    char * p;
    struct tm * local = gmtime(&od->what.nominal_time);
    sprintf(str, "%04d%02d%02d", local->tm_year+1900, local->tm_mon+1, local->tm_mday);
    status = H5Acreatestring(group_id, "date", str);
    sprintf(str, "%02d%02d%02d", local->tm_hour, local->tm_min, local->tm_sec);
    status = H5Acreatestring(group_id, "time", str);

    p = str;
    fprintf(stderr, "stations: %d\n", od->what.nstations);
    for (i = 0; i < od->what.nstations; ++i) {
      sprintf(p, "%s:%s", od->what.source[i].identifier, od->what.source[i].value);
      p += strlen(p); 
      if (i < od->what.nstations-1) {
	sprintf(p, ",");
	++p;
      }
    }
    if (od->what.nstations) {
      status = H5Acreatestring(group_id, "source", str);
    } else {
      fprintf(stderr, "Warning: source missing\n");
    }
  }
  status = H5Gclose(group_id);

  group_id = H5Gcreate(root_id, "where", 0);

  attr_id = H5Acreate(group_id, "xsize", H5T_STD_I32BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_INT, &od->where.ncols);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "ysize", H5T_STD_I32BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_INT, &od->where.nrows);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "xscale", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.psizex);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "yscale", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.psizey);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "LL_lon", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.sw.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "LL_lat", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.sw.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "UL_lon", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.nw.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "UL_lat", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.nw.lat);
  status = H5Aclose(attr_id);

  fprintf(stderr, "NW: %g %g\n", od->where.nw.lon, od->where.nw.lat);

  attr_id = H5Acreate(group_id, "UR_lon", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.ne.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "UR_lat", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.ne.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "LR_lon", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.se.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "LR_lat", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.se.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Acreatestring(group_id, "projdef", od->where.projdef);

  status = H5Gclose(group_id);

  for (i = 0; i < od->nscans; ++i) {
    hid_t g_id; 
    hsize_t dims[2];
    char str[20];

    sprintf(str, "dataset%d", i+1);
    group_id = H5Gcreate(root_id, str, 0);

    g_id = H5Gcreate(group_id, "what", 0);
    status = H5Acreatestring(g_id, "product", od->datasets[i].dataset_what.product);
    status = H5Acreatestring(g_id, "prodpar", "CAPPI");
    status = H5Acreatestring(g_id, "quantity", od->datasets[i].dataset_what.quantity);
    {
      char str[100];
      struct tm * local = gmtime(&od->datasets[i].dataset_what.start_time);
      sprintf(str, "%04d%02d%02d", local->tm_year+1900, local->tm_mon+1, local->tm_mday);
      status = H5Acreatestring(g_id, "startdate", str);
      sprintf(str, "%02d%02d%02d", local->tm_hour, local->tm_min, local->tm_sec);
      status = H5Acreatestring(g_id, "starttime", str);
      local = gmtime(&od->datasets[i].dataset_what.end_time);
      sprintf(str, "%04d%02d%02d", local->tm_year+1900, local->tm_mon+1, local->tm_mday);
      status = H5Acreatestring(g_id, "enddate", str);
      sprintf(str, "%02d%02d%02d", local->tm_hour, local->tm_min, local->tm_sec);
      status = H5Acreatestring(g_id, "endtime", str);
    }
    attr_id = H5Acreate(g_id, "gain", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.gain);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "offset", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.offset);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "nodata", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.nodata);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "undetect", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.undetect);
    status = H5Aclose(attr_id);
 
    status = H5Gclose(g_id);
    dims[0] = od->where.ncols; 
    dims[1] = od->where.nrows;

    {
      int j = 0; /* always for composites */
      hid_t plist, data_id; 
      char str[20];
      hsize_t cdims[2];
      cdims[0] = dims[0]/2;
      cdims[1] = dims[1]/2;
      sprintf(str, "data%d", j+1);
      
      if (od->datasets[i].data) {
	hid_t g_id;
	data_id = H5Gcreate(group_id, str, 0);
	dataspace_id = H5Screate_simple(2, dims, NULL);

	plist  = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(plist, 2, cdims);
	status = H5Pset_deflate(plist, 6); 

	dataset_id = H5Dcreate(data_id, "data", H5T_IEEE_F64BE, dataspace_id, plist);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, 
			  H5S_ALL, H5P_DEFAULT, od->datasets[i].data);
	status = H5Acreatestring(dataset_id, "CLASS", "IMAGE");
	status = H5Acreatestring(dataset_id, "IMAGE_VERSION", "1.2");
	status = H5Dclose(dataset_id);
	status = H5Pclose(plist);
	status = H5Sclose(dataspace_id);

	status = H5Gclose(data_id);
	}
    }

    status = H5Gclose(group_id);
  }
  
  status = H5Gclose(root_id);
  /* Close the file. */
  status = H5Fclose(file_id);

  return 1;
}

static int count = 0;
static char bufr[100];
static bufr_char_0_29_205 (varfl val, int ind) {
  bufr[count++] = val;
  return 1;
}

static int our_callback (varfl val, int ind) {

    odim_comp_t * od = &odim_data;   /* our global data structure */
    static radar_data_t b;
    static int ncorners = 0, nscans = 0;

    /* do nothing if data modification descriptor or replication descriptor */
    /*if (ind == _desc_special) {
      if (des[ind]->el->d.f != 2) fprintf (stderr,
	       " fxy: %d %d %d\n", des[ind]->el->d.f, des[ind]->el->d.x, des[ind]->el->d.y);
      if (strcmp (des[ind]->el->unit, "CCITT IA5") == 0) {
	fprintf (stderr,"vv CCITT: %g\n", val);
	return 1;
      }
      fprintf (stderr,"vv: %g\n", val);
      return 1;
      }*/

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
        if (bufr_check_fxy (d, 3,1,11)) { 
	  int i = 0;
	  b.meta.year = vv[i++];       /* Date */
	  b.meta.month = vv[i++];
	  b.meta.day = vv[i++];
	} 
	else if (bufr_check_fxy (d, 3,1,13)) { 
	  int i = 0;
	  b.meta.hour = vv[i++];       /* Time */
	  b.meta.min = vv[i++];
	  b.meta.sec = vv[i++];

	}
	/* ODIM */
        else if (bufr_check_fxy (d, 3,1,21) && (ncorners == 0)) {
	  ++ncorners;
	  od->where.nw.lat = vv[0];
	  od->where.nw.lon = vv[1];
	}
	else if (bufr_check_fxy (d, 3,1,21) && (ncorners == 1)) {
	  ++ncorners;
	  od->where.ne.lat = vv[0];
	  od->where.ne.lon = vv[1];
	}
	else if (bufr_check_fxy (d, 3,1,21) && (ncorners == 2)) {
	  ++ncorners;
	  od->where.se.lat = vv[0];
	  od->where.se.lon = vv[1];
	}
	else if (bufr_check_fxy (d, 3,1,21) && (ncorners == 3)) {
	  ++ncorners;
	  od->where.sw.lat = vv[0];
	  od->where.sw.lon = vv[1];
	}
	else if (bufr_check_fxy (d, 3,21,206)) {
	  int param, i = 0;
	      od->datasets[nscans].dataset_what.gain = 1;
	      od->datasets[nscans].dataset_what.offset = 0;
	      od->datasets[nscans].dataset_what.nodata = DBL_MAX;
	      od->datasets[nscans].dataset_what.undetect = -DBL_MAX;
      
	      param = vv[i++];
	      assert(param == 0); /* zlib compresssion */
	      
	      /* data */
	      {
		int ndecomp = od->where.ncols*od->where.nrows*sizeof(varfl);
		int k, n, j = 0, niter = vv[i++];
		unsigned char * tempo = calloc(niter*65534, sizeof(unsigned char));
		fprintf (stderr,"%d %d\n", od->where.nrows, od->where.ncols);
		fprintf (stderr,"%d %d\n", ndecomp, niter);
		for (n = 0; n < niter; ++n) {
		  int ncomp = vv[i++];
		  for (k = 0; k < ncomp; ++k) {
		    tempo[j++] = (vv[i+k] == MISSVAL ? 255 : vv[i+k]);
		  }
		  i += ncomp;
		}
		fprintf (stderr,"%d\n", i);
		od->datasets[nscans].data = my_decompress(tempo, j, &ndecomp);
		free(tempo);
		fprintf (stderr,"--> %g %g %g\n", od->datasets[nscans].data[0], 
			 od->datasets[nscans].data[1], od->datasets[nscans].data[2]);
	      }
	      ++nscans;
	    }
	/* ODIM stations */
        else if (bufr_check_fxy (d, 3,21,204)) {
	  int i = 0, j;

 	  /* what */
	  od->what.object = "COMP";
	  od->what.version = "H5rad 2.0";
	    {
	      char * tz = set_fuseau("TZ=UTC");
	      struct tm local;
	      local.tm_year = b.meta.year - 1900;
	      local.tm_mon = b.meta.month - 1;
	      local.tm_mday = b.meta.day;
	      local.tm_hour = b.meta.hour;
	      local.tm_min = b.meta.min;
	      local.tm_sec = b.meta.sec;
	      local.tm_wday = 0;
	      local.tm_yday = 0;    
	      local.tm_isdst = 0;
	      od->what.nominal_time = mktime(&local);
	    }

	  od->what.nstations = vv[i++];
	  od->what.source = calloc(od->what.nstations, sizeof(odim_station_t));
	  for (j = 0; j < od->what.nstations; ++j) {
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
	}
	/* ODIM stations */
        else if (bufr_check_fxy (d, 3,21,205)) {
	  int i = 0;
	      char * tz = set_fuseau("TZ=UTC");
	      struct tm local;

	      local.tm_year = vv[i++] - 1900;
	      fprintf (stderr, "Date start an %g\n", vv[0]);
	      local.tm_mon = vv[i++] - 1;
	      fprintf (stderr, "Date start mois %d\n", local.tm_mon);
	      local.tm_mday = vv[i++];
	      local.tm_hour = vv[i++];
	      local.tm_min = vv[i++];
	      local.tm_sec = vv[i++];
	      local.tm_wday = 0;
	      local.tm_yday = 0;    
	      local.tm_isdst = 0;
	      od->datasets[nscans].dataset_what.start_time = mktime(&local);

	      local.tm_year = vv[i++] - 1900;
	      local.tm_mon = vv[i++] - 1;
	      local.tm_mday = vv[i++];
	      local.tm_hour = vv[i++];
	      local.tm_min = vv[i++];
	      local.tm_sec = vv[i++];
	      local.tm_wday = 0;
	      local.tm_yday = 0;    
	      local.tm_isdst = 0;
	      od->datasets[nscans].dataset_what.end_time = mktime(&local);
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

      if (bufr_check_fxy (d, 0,29,205)) {
	int k; 
	od->where.projdef = calloc(101, sizeof(char));
	count = 0;
	bufr_parse_out (d, 0, 0, bufr_char_0_29_205, 0);
	fprintf (stderr, "count: %d\n", count);
	for (k = 0; k < 100; ++k) {
	  od->where.projdef[k] = bufr[k];
	}
	fprintf (stderr, "Projection: %s\n", od->where.projdef);
	return 1;
      }
      else if (bufr_check_fxy (d, 0,5,33)) {
	od->where.psizex = val;
      }
      else if (bufr_check_fxy (d, 0,6,33)) {
	od->where.psizey = val;
      }
      else if (bufr_check_fxy (d, 0,30,21)) {
	od->where.ncols = val;
      }
      else if (bufr_check_fxy (d, 0,30,22)) {
	od->where.nrows = val;
      }
      else if (bufr_check_fxy (d, 0,31,1)) {
	  od->nscans = val;
	  fprintf (stderr, "Nparams:% d\n", od->nscans);
	  od->datasets = calloc(od->nscans,sizeof(odim_comp_dataset_t));
      }
      else if (bufr_check_fxy (d, 0,30,31)) {
	int param = val;
	if (param == 1) {
	  od->datasets[nscans].dataset_what.product = "COMP";
	} else {
	  od->datasets[nscans].dataset_what.product = "XXXX";
	  fprintf(stderr, "unknown product %d\n", param);
	}
	fprintf(stderr, "Product %s\n" , od->datasets[nscans].dataset_what.product);
      }
      else if (bufr_check_fxy (d, 0,30,196)) {
	int param  = val;
	if (param == 4) {
	  od->datasets[nscans].dataset_what.quantity = "DBZH";
	} else if (param == 20) {
	  od->datasets[nscans].dataset_what.quantity = "RATE";
	} else if (param == 21) {
	  od->datasets[nscans].dataset_what.quantity = "ACCR";
	} else if (param == 99999) {
	  od->datasets[nscans].dataset_what.quantity = "QIND";
	} else {
	  od->datasets[nscans].dataset_what.quantity = "XXXX";
	  fprintf(stderr, "unknown quantity %d\n", param);
	}
	fprintf(stderr, "Quantity %s\n" , od->datasets[nscans].dataset_what.quantity);
      }
      else 
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


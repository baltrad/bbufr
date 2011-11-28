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
#include "polar.h"
#include "odim.h"

odim_polar_t odim_data; /* structure holding ODIM data */

static int write_hdf5(odim_polar_t * od, char * file)
{
   hid_t       file_id, dataset_id, dataspace_id, root_id, group_id, attr_id, space_id, prop_id, fap_id;  /* identifiers */
  herr_t      status;
  int i, j;

  /* Create a new file using optimized file-creation properties. */
  prop_id = H5Pcreate(H5P_FILE_CREATE);
  fap_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_userblock(prop_id, (hsize_t)0);
  H5Pset_sizes(prop_id, (size_t)4, (size_t)4);
  H5Pset_sym_k(prop_id, (int)1, (int)1);
  H5Pset_istore_k(prop_id, (long)1);
  H5Pset_meta_block_size(fap_id, (long)0);

  file_id = H5Fcreate(file, H5F_ACC_TRUNC, prop_id, fap_id);
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
    /* fprintf(stderr, "stations: %d\n", od->what.nstations);*/
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

  attr_id = H5Acreate(group_id, "lon", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "lat", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Acreate(group_id, "height", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id,  H5T_NATIVE_DOUBLE, &od->where.height);
  status = H5Aclose(attr_id);

  status = H5Gclose(group_id);

  for (i = 0; i < od->nscans; ++i) {
    hid_t g_id; 
    hsize_t dims[2];
    char str[20];

    sprintf(str, "dataset%d", i+1);
    group_id = H5Gcreate(root_id, str, 0);

    g_id = H5Gcreate(group_id, "what", 0);
    status = H5Acreatestring(g_id, "product", "SCAN");
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
    status = H5Gclose(g_id);

    g_id = H5Gcreate(group_id, "where", 0);

    attr_id = H5Acreate(g_id, "elangle", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_where.elangle);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "nbins", H5T_STD_I32BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &od->datasets[i].dataset_where.nbins);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "rstart", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_where.rstart);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "rscale", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_where.rscale);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "nrays", H5T_STD_I32BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &od->datasets[i].dataset_where.nrays);
    status = H5Aclose(attr_id);

    attr_id = H5Acreate(g_id, "a1gate", H5T_STD_I32BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &od->datasets[i].dataset_where.a1gate);
    status = H5Aclose(attr_id);
      
    status = H5Gclose(g_id);
    dims[0] = od->datasets[i].dataset_where.nrays; 
    dims[1] = od->datasets[i].dataset_where.nbins; 

    for (j = 0; j < od->datasets[i].nparams; ++j) {
      hid_t plist, data_id; 
      char str[20];
      hsize_t cdims[2];
      cdims[0] = dims[0];
      cdims[1] = dims[1];
      sprintf(str, "data%d", j+1);
      
      if (od->datasets[i].dataset_data[j].data) {
	hid_t g_id;
	data_id = H5Gcreate(group_id, str, 0);
	g_id = H5Gcreate(data_id, "what", 0);

	status = H5Acreatestring(g_id, "quantity", od->datasets[i].dataset_data[j].dataset_what.quantity);

	attr_id = H5Acreate(g_id, "gain", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.gain);
	status = H5Aclose(attr_id);

	attr_id = H5Acreate(g_id, "offset", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.offset);
	status = H5Aclose(attr_id);

	attr_id = H5Acreate(g_id, "nodata", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.nodata);
	status = H5Aclose(attr_id);

	attr_id = H5Acreate(g_id, "undetect", H5T_IEEE_F64BE, H5Screate(H5S_SCALAR), H5P_DEFAULT);
	status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.undetect);
	status = H5Aclose(attr_id);
      
	status = H5Gclose(g_id);
	dataspace_id = H5Screate_simple(2, dims, NULL);

	plist  = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(plist, 2, cdims);
	status = H5Pset_deflate(plist, 6); 

	dataset_id = H5Dcreate(data_id, "data", H5T_IEEE_F64BE, dataspace_id, plist);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, 
			  H5S_ALL, H5P_DEFAULT, od->datasets[i].dataset_data[j].data);
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
  status = H5Pclose(prop_id);
  status = H5Pclose(fap_id);
  /* Close the file. */
  status = H5Fclose(file_id);

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
    local.tm_year = (s1.year < 50 ? s1.year+2000 : s1.year-1900);
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


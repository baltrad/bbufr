/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <float.h>

/* Interface definition of the OPERA coder/decoder */
#include "desc.h"
#include "bufr.h"
#include "bitio.h"
#include "rlenc.h"

#include "radar.h"
#include "polar.h"
#include "odim.h"

int loc = 8; /* default local BUFR table */
odim_polar_t odim_data; /* structure holding ODIM data */

static char *table_dir = NULL;     /* directory for BUFR tables */

#define START  { herr_t (*old_func)(void*);\
      void *old_client_data;\
      H5Eget_auto(&old_func, &old_client_data);\
      H5Eset_auto(NULL, NULL);

#define END H5Eset_auto(old_func, old_client_data);}

static herr_t aiter_cb(hid_t location_id, const char *attr_name, void *op_data) 
{
  char * s;
  double val;
  hid_t attr_id;
  fprintf(stderr, "%s:", attr_name);
  {
    herr_t status;
    hid_t attr_id = H5Aopen_name(location_id, attr_name);
    status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &val);
    if (status == 0) {
      fprintf(stderr, " %g\n", val);
      if (op_data) {
	odim_polar_t * od = (odim_polar_t *)op_data;
	odim_polar_how_double_t * p = od->how_double;
	odim_polar_how_double_t * nouveau = malloc(sizeof(odim_polar_how_double_t));
	nouveau->next = p;
	nouveau->value = val;
	nouveau->id = strdup(attr_name);
	od->how_double = nouveau;
      }
      H5Aclose(attr_id);
      return 0;
    } else {
      H5Aclose(attr_id);
      s = H5Areadstring(location_id, attr_name);
      if (s) {
	fprintf(stderr, " %s\n", s);
	if (op_data) {
	  odim_polar_t * od = (odim_polar_t *)op_data;
	  odim_polar_how_string_t * p = od->how_string;
	  odim_polar_how_string_t * nouveau = malloc(sizeof(odim_polar_how_string_t));
	  nouveau->next = p;
	  nouveau->value = s;
	  nouveau->id = strdup(attr_name);
	  od->how_string = nouveau;
	} else {
	  free(s);
	}
	return 0;
      }
    }
  }
  fprintf(stderr, "Unknown type\n");
  H5Aclose(attr_id);
  return 0;
}

static time_t startepochs = 0;

static int read_hdf5(odim_polar_t * od, char * file)
{
  hid_t       file_id, dataset_id, dataspace_id, root_id, how_id, group_id, attr_id, space_id, ftype;  /* identifiers */
  herr_t      status;
  int i, j;
  char * result;

  /* Create a new file using default properties. */
  file_id = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);
  root_id = H5Gopen(file_id, "/");

  START
  result = H5Areadstring(root_id, "Conventions");
  END
  if (!result) {
    fprintf(stderr, "Warning: Conventions missing\n");
  } else {
    if (strcmp(result, "ODIM_H5/V2_1") == 0) {
      loc = 9;
    }
  }

  START
  how_id = H5Gopen(root_id, "how");
  if (how_id >= 0) {
    fprintf(stderr, "===============\n");
    fprintf(stderr, "Top how:\n");
    {
      int ret;
      unsigned idx = 0;            /* Index in the attribute list */
      od->how_string = NULL;
      od->how_double = NULL;
      while((ret = H5Aiterate(how_id, &idx, aiter_cb, od)) > 0) {
	fprintf(stderr, "%d %d\n", ret, idx);
      }
      
    }
    fprintf(stderr, "===============\n");

    {
      double begins, ends;
      attr_id = H5Aopen_name(how_id, "startepochs");
      status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &begins);
      if (!status) {
	fprintf(stderr, "startepochs: %.15lg\n", begins);
	startepochs = (time_t)begins;
      }
      status = H5Aclose(attr_id);
      attr_id = H5Aopen_name(how_id, "stopepochs");
      status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &ends);
      if (!status) {
	fprintf(stderr, "stopepochs: %.15lg\n", ends);
	od->what.nominal_time = (time_t)ends;
      }
      status = H5Aclose(attr_id);
    }
	  
    status = H5Gclose(how_id);
  }
  END

  group_id = H5Gopen(root_id, "what");
  od->what.object = H5Areadstring(group_id, "object");
  START
  od->what.version = H5Areadstring(group_id, "version");
  END
  if (!od->what.version) { 
    od->what.version = "H5rad 2.0"; 
  }
  {
    char * date, * time;
    START
    date = H5Areadstring(group_id, "date");
    time = H5Areadstring(group_id, "time");
    END
    if (date && time) {
      od->what.nominal_time = build_time(date, time);
    } else {
      fprintf(stderr, "Warning: date or time missing, filling with stopepochs\n");
    }
  }
  {
    char * s;
    START
    s = H5Areadstring(group_id, "source");
    END
    if (!s) {
      fprintf(stderr, "Warning: source missing\n");
    } else {
      int i, n = split(s, ",", ":");
      od->what.nstations = n;
      if (n == 0) {
	fprintf(stderr, "Warning: source missing\n");
      } else {
	od->what.source = calloc(od->what.nstations, sizeof(odim_station_t));
	for (i = 0; i < n; ++i) {
	  od->what.source[i].identifier = getsplit_typ(H5Areadstring(group_id, "source"), ",", ":", i+1);
	  od->what.source[i].value = getsplit_val(H5Areadstring(group_id, "source"), ",", ":", i+1);
	  fprintf(stderr, "source: %s:%s\n", od->what.source[i].identifier, od->what.source[i].value);
	}
      }
    }
  }

  status = H5Gclose(group_id);

  group_id = H5Gopen(root_id, "where");

  attr_id = H5Aopen_name(group_id, "lon");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "lat");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "height");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.height);
  status = H5Aclose(attr_id);
  
  fprintf(stderr, "lon: %g lat: %g height: %g\n", 
	  od->where.lon, od->where.lat, od->where.height);

  status = H5Gclose(group_id);

  /* Number of scans determination */
  {
    int i = 0;
    char str[20];

    START
    while (1) {
      sprintf(str, "dataset%d", i+1);
      if ((group_id = H5Gopen(root_id, str)) < 0) {
	od->nscans = i;
	break;
      }
      H5Gclose(group_id);
      ++i;
    }
    fprintf(stderr, "nscans: %d\n", od->nscans);
    END
  }

  od->datasets = calloc(od->nscans, sizeof(odim_polar_dataset_t));
  for (i = 0; i < od->nscans; ++i) {
    hid_t g_id;
    hsize_t dims[2];

    char str[20];
    sprintf(str, "dataset%d", i+1);

    group_id = H5Gopen(root_id, str);

    START
    how_id = H5Gopen(group_id, "how");
    if (how_id >= 0) {
      fprintf(stderr, "===============\n");
      fprintf(stderr, "Dataset%d how:\n", i+1);
      {  
	int ret;
	unsigned idx = 0;            /* Index in the attribute list */
	od->datasets[i].how_string = NULL;
	od->datasets[i].how_double = NULL;
	while((ret = H5Aiterate(how_id, &idx, aiter_cb, &od->datasets[i])) > 0) {
	  fprintf(stderr, "%d %d\n", ret, idx);
	}
	
      }
      fprintf(stderr, "===============\n");
      status = H5Gclose(how_id);
    } else  {
      fprintf(stderr, "Warning: no dataset%d how\n", i+1);
    }
    END

    START
    g_id = H5Gopen(group_id, "what");
    END
    if (g_id < 0) {
      fprintf(stderr, "Warning: no dataset%d what\n", i+1);
      od->datasets[i].dataset_what.product = "SCAN  ";
      od->datasets[i].dataset_what.start_time = startepochs;
      od->datasets[i].dataset_what.end_time = od->what.nominal_time;
    } else {
      od->datasets[i].dataset_what.product = result = H5Areadstring(g_id, "product");
      assert((strcmp(result, "SCAN") == 0) || (strcmp(result, "SCAN  ") == 0));

      START
	{
	  char * date, * time;
	  START
	  date = H5Areadstring(g_id, "startdate");
	  time = H5Areadstring(g_id, "starttime");
	  END
	    if (date && time) {
	      od->datasets[i].dataset_what.start_time = build_time(date, time);
	    } else {
	      fprintf(stderr, "Warning: startdate or starttime missing, filling with bogus\n");
	      od->datasets[i].dataset_what.start_time = build_time("19990101", "120000");
	    }
	}
    
      {
	char * date = H5Areadstring(g_id, "enddate");
	char * time = H5Areadstring(g_id, "endtime");
	if (!date || !time) {
	  fprintf(stderr, "Warning: no enddate or endtime, using startdate and starttime\n");
	  od->datasets[i].dataset_what.end_time = od->datasets[i].dataset_what.start_time;
	} else {
	  od->datasets[i].dataset_what.end_time = build_time(date, time);
	}
      }
      END

      status = H5Gclose(g_id);
    }

    START
    g_id = H5Gopen(group_id, "where");
    END
    if (g_id < 0) {
      fprintf(stderr, "Warning: no dataset%d where, trying global where...\n", i+1);
      g_id = H5Gopen(root_id, "where");
    }
    START
    attr_id = H5Aopen_name(g_id, "elangle");
    END
    if (attr_id < 0) {
      fprintf(stderr, "Warning: e1angle missing, trying angle\n");
      attr_id = H5Aopen_name(g_id, "angle");
    }
    status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_where.elangle);
    status = H5Aclose(attr_id);
      
    START
    attr_id = H5Aopen_name(g_id, "nbins");
    END
    if (attr_id < 0) {
      fprintf(stderr, "Warning: nbins missing, trying rsize\n");
      attr_id = H5Aopen_name(g_id, "rsize");
    }
    status = H5Aread(attr_id, H5T_NATIVE_INT, &od->datasets[i].dataset_where.nbins);
    status = H5Aclose(attr_id);
    
    attr_id = H5Aopen_name(g_id, "rstart");
    status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_where.rstart);
    status = H5Aclose(attr_id);

    attr_id = H5Aopen_name(g_id, "rscale");
    status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_where.rscale);
    status = H5Aclose(attr_id);
    if (od->datasets[i].dataset_where.rscale < 10) {
      fprintf(stderr, "Warning: rscale too small, multiplying it by 1000\n");
      od->datasets[i].dataset_where.rscale *= 1000;
    }

    START
    attr_id = H5Aopen_name(g_id, "nrays");
    END
    if (attr_id < 0) {
      fprintf(stderr, "Warning: nrays missing, trying asize\n");
      attr_id = H5Aopen_name(g_id, "asize");
    }

    status = H5Aread(attr_id, H5T_NATIVE_INT, &od->datasets[i].dataset_where.nrays);
    status = H5Aclose(attr_id);

    attr_id = H5Aopen_name(g_id, "a1gate");
    status = H5Aread(attr_id, H5T_NATIVE_INT, &od->datasets[i].dataset_where.a1gate);
    status = H5Aclose(attr_id);
      
    fprintf(stderr, "e1angle: %g nbins: %d rstart: %g\n", od->datasets[i].dataset_where.elangle,
	    od->datasets[i].dataset_where.nbins, od->datasets[i].dataset_where.rstart);
    fprintf(stderr, "rscale: %g nrays: %d a1gate: %d\n", od->datasets[i].dataset_where.rscale,
	    od->datasets[i].dataset_where.nrays, od->datasets[i].dataset_where.a1gate);

    status = H5Gclose(g_id);
    dims[0] = od->datasets[i].dataset_where.nrays;
    dims[1] = od->datasets[i].dataset_where.nbins;

    /* Number of params determination */
    {
      int j = 0;
      char str[20];

      START
      while (1) {
	sprintf(str, "data%d", j+1);
	if ((space_id = H5Gopen(group_id, str)) < 0) {
	  od->datasets[i].nparams = j;
	break;
	}
	H5Gclose(space_id);
	++j;
      }
      END
    }
    fprintf(stderr, "nparams for scan %d: %d (nrays: %d nbins: %d)\n", i+1, od->datasets[i].nparams, 
	    od->datasets[i].dataset_where.nrays,
	    od->datasets[i].dataset_where.nbins);
    
    od->datasets[i].dataset_data = calloc(od->datasets[i].nparams, sizeof(odim_polar_dataset_data_t));

    for (j = 0; j < od->datasets[i].nparams; ++j) {
      hid_t plist, data_id;
      char str[20];
      hsize_t cdims[2];
      cdims[0] = 200;
      cdims[1] = 200;
      sprintf(str, "data%d", j+1);
      
      od->datasets[i].dataset_data[j].data = 
	calloc(od->datasets[i].dataset_where.nbins*od->datasets[i].dataset_where.nrays, sizeof(varfl));
      if (od->datasets[i].dataset_data[j].data) {
	hid_t g_id;
	data_id = H5Gopen(group_id, str);

	START
	how_id = H5Gopen(data_id, "how");
	if (how_id >= 0) {
	  fprintf(stderr, "===============\n");
	  fprintf(stderr, "Data%d how:\n", j+1);
	  {  
	    int ret;
	    unsigned idx = 0;            /* Index in the attribute list */
	    od->datasets[i].dataset_data[j].how_string = NULL;
	    od->datasets[i].dataset_data[j].how_double = NULL;
	    while((ret = H5Aiterate(how_id, &idx, aiter_cb, &od->datasets[i].dataset_data[j])) > 0) {
	      fprintf(stderr, "%d %d\n", ret, idx);
	    }
	  }
	  fprintf(stderr, "===============\n");
	  status = H5Gclose(how_id);
	}
	END

	g_id = H5Gopen(data_id, "what");

	od->datasets[i].dataset_data[j].dataset_what.quantity = H5Areadstring(g_id, "quantity");
	if (strcmp(od->datasets[i].dataset_data[j].dataset_what.quantity, "TH") &&
	    strcmp(od->datasets[i].dataset_data[j].dataset_what.quantity, "DBZH") &&
	    strcmp(od->datasets[i].dataset_data[j].dataset_what.quantity, "WRAD") &&
	    strcmp(od->datasets[i].dataset_data[j].dataset_what.quantity, "WRAD") &&
	    strcmp(od->datasets[i].dataset_data[j].dataset_what.quantity, "RATE")) {
	  fprintf(stderr, "Warning: Unknown quantity %s\n", od->datasets[i].dataset_data[j].dataset_what.quantity);
	}
	fprintf(stderr, "Quantity %s\n", od->datasets[i].dataset_data[j].dataset_what.quantity);

	attr_id = H5Aopen_name(g_id, "gain");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.gain);
	status = H5Aclose(attr_id);

	attr_id = H5Aopen_name(g_id, "offset");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.offset);
	status = H5Aclose(attr_id);

	attr_id = H5Aopen_name(g_id, "nodata");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.nodata);
	status = H5Aclose(attr_id);

	attr_id = H5Aopen_name(g_id, "undetect");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_data[j].dataset_what.undetect);
	status = H5Aclose(attr_id);
	fprintf(stderr, "gain: %g offset: %g nodata: %g undetect: %g\n", 
		od->datasets[i].dataset_data[j].dataset_what.gain,
		od->datasets[i].dataset_data[j].dataset_what.offset,
		od->datasets[i].dataset_data[j].dataset_what.nodata,
		od->datasets[i].dataset_data[j].dataset_what.undetect);
	status = H5Gclose(g_id);
	dataset_id = H5Dopen(data_id, "data");
	{
	  START
	  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, od->datasets[i].dataset_data[j].data);
	  if (status < 0) {
	    unsigned char * p = calloc(od->datasets[i].dataset_where.nrays*od->datasets[i].dataset_where.nbins, 
				       sizeof(unsigned char));
	    /*int kk; 
	    for (kk = 0; kk < od->datasets[i].dataset_where.nrays; kk++) {
	      p[kk] = calloc(od->datasets[i].dataset_where.nbins, sizeof(unsigned char));
	      }*/
	   
	    fprintf(stderr, "Warning: trying unsigned char for data %d %d\n", 
		    od->datasets[i].dataset_where.nrays,
		    od->datasets[i].dataset_where.nbins);
	    status = H5Dread(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);
	    if (status == 0) {
	      int k; 
	      fprintf(stderr, "Ok unsigned char for data\n");
	      for (k = 0; k < od->datasets[i].dataset_where.nbins*od->datasets[i].dataset_where.nrays; k++) {
		od->datasets[i].dataset_data[j].data[k] = p[k];
	      }
	    } else {
	      unsigned short * p = calloc(od->datasets[i].dataset_where.nrays*od->datasets[i].dataset_where.nbins, 
				       sizeof(unsigned short));
	      fprintf(stderr, "Warning: trying unsigned short (%d) for data %d %d\n", 
		      sizeof(unsigned char),
		      od->datasets[i].dataset_where.nrays,
		      od->datasets[i].dataset_where.nbins);
	      status = H5Dread(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);
	      if (status == 0) {
		int k; 
		fprintf(stderr, "Ok unsigned short for data\n");
		for (k = 0; k < od->datasets[i].dataset_where.nbins*od->datasets[i].dataset_where.nrays; k++) {
		  od->datasets[i].dataset_data[j].data[k] = p[k];
		}
	      }
	      free(p);
	    }
	    free(p);
	  }
	  END
	}
    
	{
	  int k, n_nodata = 0, n_undetect = 0;
	  for (k = 0; k < od->datasets[i].dataset_where.nbins*od->datasets[i].dataset_where.nrays; k++) {
	    if (od->datasets[i].dataset_data[j].data[k] == od->datasets[i].dataset_data[j].dataset_what.nodata) {
	      od->datasets[i].dataset_data[j].data[k] = DBL_MAX;
	      ++n_nodata;
	    } else if (od->datasets[i].dataset_data[j].data[k] == od->datasets[i].dataset_data[j].dataset_what.undetect) {
	      od->datasets[i].dataset_data[j].data[k] = -DBL_MAX;
	      ++n_undetect;
	    } else {
	      od->datasets[i].dataset_data[j].data[k] =
		od->datasets[i].dataset_data[j].data[k]*od->datasets[i].dataset_data[j].dataset_what.gain +
		od->datasets[i].dataset_data[j].dataset_what.offset;
	    }
	  }
	  fprintf(stderr, "nodata: %d undetect: %d\n", n_nodata, n_undetect);
	}

	START
	result = H5Areadstring(dataset_id, "CLASS");
	if (result) assert(strcmp(result, "IMAGE") == 0);
	
	result = H5Areadstring(dataset_id, "IMAGE_VERSION");
	if (!result) {
	  fprintf(stderr, "Trying VERSION instead...\n");
	  result = H5Areadstring(dataset_id, "VERSION");
	}
	if (result) assert(strcmp(result, "1.2") == 0);
	END
	
	status = H5Dclose(dataset_id);
	fprintf(stderr, "--> %g %g %g\n", od->datasets[i].dataset_data[j].data[0],
		od->datasets[i].dataset_data[j].data[1],
		od->datasets[i].dataset_data[j].data[2]);

	status = H5Gclose(data_id);
      }
    end: ;
    }

    status = H5Gclose(group_id);
  }
  
  status = H5Gclose(root_id);
  /* Close the file. */
  status = H5Fclose(file_id);

  return 1;
}

static int write_opera(sect_1_t *s1_in, dd *descr, int iout, varfl *v, char *outfile)
/* Create a new OPERA BUFR message */
{
  int i, ret;
  sect_1_t s1;       /* Here we store section 1 of BUFR message */
  bufr_t msg;

  /******* Prepare data for Section 1 (output) */
  s1.year = s1_in->year;
  s1.mon  = s1_in->mon;
  s1.day = s1_in->day;
  s1.hour = s1_in->hour;
  s1.min  = s1_in->min;
  s1.sec  = s1_in->sec;
  s1.mtab = s1_in->mtab;            /* master table used */
  s1.subcent = 0;
  s1.gencent = 247;
  s1.updsequ = s1_in->updsequ;      /* original BUFR message */
  s1.opsec = s1_in->opsec;          /* no optional section */
  s1.dcat = s1_in->dcat;            /* message type */
  s1.dcatst = s1_in->dcatst;        /* message subtype */
  s1.idcatst = s1_in->idcatst;
  s1.vmtab = (loc == 9 ? 16 : 13);    /* version number of master table used */
  s1.vltab = loc;                     /* version number of local table used */

  _bufr_edition = 4;
  header_dump(&s1);

  /******* read supported data descriptors */
  if (read_tables(table_dir, s1.vmtab, s1.vltab, s1.subcent, s1.gencent)) {
    fprintf (stderr, "FATAL: Unable to read tables\n");
    exit (EXIT_FAILURE);
  }

  /******* Code the data (section 3 and 4) */
  if (!bufr_encode_sections34(descr, iout, v, &msg)) return 0;

  /******* Setup section 0, 1, 2, 5 */
  if (!bufr_encode_sections0125(&s1, &msg)) {
    fprintf (stderr, "WARNING: Unable to create section 0, 1, 2 and/or 5\n");
    return 0;
  }
    
  /******* Save coded data */
  if (!bufr_write_file(&msg, outfile)) return 0;

  free(v);
  return 1;
}

#define fill_desc(ff,xx,yy) {dd dds; dds.f=ff; dds.x=xx; dds.y=yy; bufr_desc_to_array(descr, dds, &iout);}

#define fill_v(val) bufr_val_to_array(&v, val, &jout)

static int recode_opera_local (odim_polar_t *od, char *outfile, sect_1_t *s1_in)
/* This function recodes a BUFR-message; returns 1 on success, 0 on a fault.*/
{
  int iout = 0;
  unsigned int jout = 0;

  dd descr[MAX_DESCS]; /* This array must be huge enough to hold all required descriptors */
  varfl *v = 0;        /* This array must be huge enough to hold all corresponding data values */

  /* The data we are interested in is read to array 'v' 
     and the corresponding descriptors are in array 'descr' */

  int ii, flag_omm = 0;
  
  for (ii = 0; ii < od->what.nstations; ++ii) {
    if (strcmp(od->what.source[ii].identifier, "WMO") == 0) {
      flag_omm = 1; 
      break;
    }
  }

  /* Output phase */
  {
    int ii, jj; 

    fill_desc(3,21,204);
    /*fprintf(stderr, "nstations %d %d\n", od->what.nstations, flag_omm);*/
    if (flag_omm) {
      fill_v(od->what.nstations-1);
    } else {
      fill_v(od->what.nstations);
    }
    for (ii = 0; ii < od->what.nstations; ++ii) {
      if (strcmp(od->what.source[ii].identifier, "WMO")) {
	for (jj = 0; jj < 3; ++jj) {      
	  fill_v(od->what.source[ii].identifier[jj]);
	}
	for (jj = 0; jj < strlen(od->what.source[ii].value); ++jj) {      
	  fill_v(od->what.source[ii].value[jj]);
	}
	for (jj = strlen(od->what.source[ii].value); jj < 16; ++jj) {      
	  fill_v(0);
	}
      }
    }

    fill_desc(3,1,31);
    if (flag_omm) {
	int j;
	for (j = 0; j < od->what.nstations; ++j) {
	  if (strcmp(od->what.source[j].identifier, "WMO") == 0) {
	    break;
	  }
	}
	{
	  long omm = atol(od->what.source[j].value);
	  if (omm > 99999) {
	    fprintf(stderr, "bogus OMM: %ld\n", omm);
	    /* in case we get a bogus WMO number */
	    fill_v(MISSVAL); 
	    fill_v(MISSVAL);
	  } else {
	    fill_v(omm / 1000); 
	    fill_v(omm % 1000);
	  }
	}
    } else {
      /*fprintf(stderr, "OMM missing\n");*/
      fill_v(MISSVAL); 
      fill_v(MISSVAL);
    }

    fill_v(0); /* automatic station */
  
    {
      struct tm * local = gmtime(&od->what.nominal_time);
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
    }

    /* fill_desc(0,5,1); */
    fill_v(od->where.lat); 
    
    /* fill_desc(0,6,1); */
    fill_v(od->where.lon); 
    
    /* fill_desc(0,7,1); */
    fill_v(od->where.height); 

   fill_desc(3,21,203);
    fill_v(od->nscans);
    for (ii = 0; ii < od->nscans; ++ii) {
      struct tm * local = gmtime(&od->datasets[ii].dataset_what.start_time);
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
      fill_v(local->tm_sec);
      local = gmtime(&od->datasets[ii].dataset_what.end_time);
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
      fill_v(local->tm_sec);
      fill_v(90); /* full polar volume */
      fill_v(od->datasets[ii].dataset_where.elangle);
      fill_v(od->datasets[ii].dataset_where.nbins);
      fill_v(od->datasets[ii].dataset_where.rscale);
      fill_v(od->datasets[ii].dataset_where.rstart);
      fill_v(od->datasets[ii].dataset_where.nrays);
      fill_v(od->datasets[ii].dataset_where.a1gate);
      fill_v(od->datasets[ii].nparams); fprintf(stderr, "nparms = %d\n", od->datasets[ii].nparams);
      for (jj = 0; jj < od->datasets[ii].nparams; ++jj) {
	if ((strcmp(od->datasets[ii].dataset_data[jj].dataset_what.quantity, "DBZH") == 0)) {
	  fill_v(0);
	} else if ((strcmp(od->datasets[ii].dataset_data[jj].dataset_what.quantity, "VRAD") == 0)) {
	  fill_v(40);
	} else if ((strcmp(od->datasets[ii].dataset_data[jj].dataset_what.quantity, "TH") == 0)) {
	  fill_v(91);
	} else if ((strcmp(od->datasets[ii].dataset_data[jj].dataset_what.quantity, "RATE") == 0)) {
	  fill_v(93);
	} else if ((strcmp(od->datasets[ii].dataset_data[jj].dataset_what.quantity, "WRAD") == 0)) {
	  fill_v(92);
	} else { /* polar volume reflectivity or radial wind */
	  fprintf(stderr, "Warning: Unknown quantity %s\n", od->datasets[ii].dataset_data[jj].dataset_what.quantity);
	}
	fill_v(0); /* zlib compression */
	{
	int i, ncomp, n = od->datasets[ii].dataset_where.nbins*od->datasets[ii].dataset_where.nrays;
	unsigned char * result;
	/*assert(od->datasets[ii].dataset_data[jj].dataset_what.gain == 1);
	  assert(od->datasets[ii].dataset_data[jj].dataset_what.offset == 0);*/
	/*assert(od->datasets[ii].dataset_data[jj].dataset_what.nodata == DBL_MAX);
	  assert(od->datasets[ii].dataset_data[jj].dataset_what.undetect == -DBL_MAX);*/
	assert(result = my_compress(od->datasets[ii].dataset_data[jj].data, n, &ncomp));
	fprintf(stderr, "%d (%d packets, rest %d)\n", ncomp, ncomp/65534 + 1, ncomp%65534);
	fill_v(ncomp/65534 + 1);
	for (i = 0; i < ncomp/65534; ++i) {
	  fill_v(65534); 
	  /*fprintf(stderr, "%d\n", 65534);*/
	  for (n = 0; n < 65534; ++n) {
	    fill_v(result[n + i*65534]); 
	  }
	}
	fill_v(ncomp%65534);
	for (n = 0; n < ncomp%65534; ++n) {
	  fill_v(result[n + (ncomp/65534)*65534]); 
	}
	free(result);
	}
      }
    }
  }

  free_descs();

  /* Now building new transcoded BUFR message */
  fprintf (stderr, "Output file header:\n");
  return write_opera(s1_in, descr, iout, v, outfile);

}

#define fill_how(p0, pp0) \
    { \
      int len = 0; \
      odim_polar_how_string_t * p = p0; \
      while (p) { \
	len++; \
	p = p->next; \
      } \
      fill_v(len); \
 \
      p = p0; \
      while (p) { \
	{ \
	  int jj, i = 0;  \
	  const char * s = p->id; \
	  for (jj = 0; (i < 16) && (jj < strlen(s)); ++jj) {       \
	    fill_v(s[jj]); ++i; \
	  } \
	  while (i < 16) { \
	    fill_v(' '); ++i; \
	  } \
	  s = p->value; i = 0; \
	  for (jj = 0; (i < 16) && (jj < strlen(s)); ++jj) {       \
	    fill_v(s[jj]); ++i; \
	  } \
	  while (i < 16) { \
	    fill_v(' '); ++i; \
	  } \
	} \
	p = p->next; \
      } \
    } \
    { \
      int len = 0; \
      odim_polar_how_double_t * p = pp0; \
      while (p) { \
	len++; \
	p = p->next; \
      } \
      fill_v(len); \
 \
      p = pp0; \
      while (p) { \
	{ \
	  int jj, i = 0;  \
	  const char * s = p->id; \
	  const unsigned char * ss; \
	  for (jj = 0; (i < 16) && (jj < strlen(s)); ++jj) {       \
	    fill_v(s[jj]); ++i; \
	  } \
	  while (i < 16) { \
	    fill_v(' '); ++i; \
	  } 	\
	  ss = (unsigned char *)&p->value;		\
	  for (jj = 0; jj < 8; ++jj) {       \
	    fill_v(ss[jj]); \
	    }	\
	} \
	p = p->next; \
      } \
    }

static int recode_opera_local_v2_1 (odim_polar_t *od, char *outfile, sect_1_t *s1_in)
/* This function recodes a BUFR-message; returns 1 on success, 0 on a fault.*/
{
  int iout = 0;
  unsigned int jout = 0;

  dd descr[MAX_DESCS]; /* This array must be huge enough to hold all required descriptors */
  varfl *v = 0;        /* This array must be huge enough to hold all corresponding data values */

  /* The data we are interested in is read to array 'v' 
     and the corresponding descriptors are in array 'descr' */

  int ii, flag_omm = 0;
  
  for (ii = 0; ii < od->what.nstations; ++ii) {
    if (strcmp(od->what.source[ii].identifier, "WMO") == 0) {
      flag_omm = 1; 
      break;
    }
  }

  /* Output phase */
  {
    int ii, jj; 

    fill_desc(3,21,204);
    /*fprintf(stderr, "nstations %d %d\n", od->what.nstations, flag_omm);*/
    if (flag_omm) {
      fill_v(od->what.nstations-1);
    } else {
      fill_v(od->what.nstations);
    }
    for (ii = 0; ii < od->what.nstations; ++ii) {
      if (strcmp(od->what.source[ii].identifier, "WMO")) {
	for (jj = 0; jj < 3; ++jj) {      
	  fill_v(od->what.source[ii].identifier[jj]);
	}
	for (jj = 0; jj < strlen(od->what.source[ii].value); ++jj) {      
	  fill_v(od->what.source[ii].value[jj]);
	}
	for (jj = strlen(od->what.source[ii].value); jj < 16; ++jj) {      
	  fill_v(0);
	}
      }
    }

    fill_desc(3,1,31);
    if (flag_omm) {
	int j;
	for (j = 0; j < od->what.nstations; ++j) {
	  if (strcmp(od->what.source[j].identifier, "WMO") == 0) {
	    break;
	  }
	}
	{
	  long omm = atol(od->what.source[j].value);
	  if (omm > 99999) {
	    fprintf(stderr, "bogus OMM: %ld\n", omm);
	    /* in case we get a bogus WMO number */
	    fill_v(MISSVAL); 
	    fill_v(MISSVAL);
	  } else {
	    fill_v(omm / 1000); 
	    fill_v(omm % 1000);
	  }
	}
    } else {
      /*fprintf(stderr, "OMM missing\n");*/
      fill_v(MISSVAL); 
      fill_v(MISSVAL);
    }

    fill_v(0); /* automatic station */
  
    {
      struct tm * local = gmtime(&od->what.nominal_time);
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
    }

    /* fill_desc(0,5,1); */
    fill_v(od->where.lat); 
    
    /* fill_desc(0,6,1); */
    fill_v(od->where.lon); 
    
    /* fill_desc(0,7,1); */
    fill_v(od->where.height); 

    fill_desc(3,21,207);
    fill_how(od->how_string, od->how_double)
    fill_v(od->nscans);
    for (ii = 0; ii < od->nscans; ++ii) {
      struct tm * local = gmtime(&od->datasets[ii].dataset_what.start_time);
      fill_how(od->datasets[ii].how_string, od->datasets[ii].how_double)
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
      fill_v(local->tm_sec);
      local = gmtime(&od->datasets[ii].dataset_what.end_time);
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
      fill_v(local->tm_sec);
      {
	char * s = od->datasets[ii].dataset_what.product;
	assert(strlen(s) <= 6);
	for (jj = 0; jj < strlen(s); ++jj) {      
	  fill_v(s[jj]);
	}
	for (jj = strlen(s); jj < 6; ++jj) {      
	  fill_v(' ');
	}
      }
      fill_v(od->datasets[ii].dataset_where.elangle);
      fill_v(od->datasets[ii].dataset_where.nbins);
      fill_v(od->datasets[ii].dataset_where.rscale);
      fill_v(od->datasets[ii].dataset_where.rstart);
      fill_v(od->datasets[ii].dataset_where.nrays);
      fill_v(od->datasets[ii].dataset_where.a1gate);
      fill_v(od->datasets[ii].nparams); fprintf(stderr, "nparms = %d\n", od->datasets[ii].nparams);
      for (jj = 0; jj < od->datasets[ii].nparams; ++jj) {
        fill_how(od->datasets[ii].dataset_data[jj].how_string, od->datasets[ii].dataset_data[jj].how_double)
	{
	  int k;
	  char * s = od->datasets[ii].dataset_data[jj].dataset_what.quantity;
	  for (k = 0; k < strlen(s); ++k) {      
	    fill_v(s[k]);
	  }
	  for (k = strlen(s); k < 6; ++k) {      
	    fill_v(' ');
	  }
	}
	fill_v(0); /* zlib compression */
	{
	int i, ncomp, n = od->datasets[ii].dataset_where.nbins*od->datasets[ii].dataset_where.nrays;
	unsigned char * result;
	/*assert(od->datasets[ii].dataset_data[jj].dataset_what.gain == 1);
	  assert(od->datasets[ii].dataset_data[jj].dataset_what.offset == 0);*/
	/*assert(od->datasets[ii].dataset_data[jj].dataset_what.nodata == DBL_MAX);
	  assert(od->datasets[ii].dataset_data[jj].dataset_what.undetect == -DBL_MAX);*/
	assert(result = my_compress(od->datasets[ii].dataset_data[jj].data, n, &ncomp));
	fprintf(stderr, "%d (%d packets, rest %d)\n", ncomp, ncomp/65534 + 1, ncomp%65534);
	fill_v(ncomp/65534 + 1);
	for (i = 0; i < ncomp/65534; ++i) {
	  fill_v(65534); 
	  /*fprintf(stderr, "%d\n", 65534);*/
	  for (n = 0; n < 65534; ++n) {
	    fill_v(result[n + i*65534]); 
	  }
	}
	fill_v(ncomp%65534);
	for (n = 0; n < ncomp%65534; ++n) {
	  fill_v(result[n + (ncomp/65534)*65534]); 
	}
	free(result);
	}
      }
    }
  }

  free_descs();

  /* Now building new transcoded BUFR message */
  fprintf (stderr, "Output file header:\n");
  return write_opera(s1_in, descr, iout, v, outfile);

}

/*===========================================================================
                                  Main program
  ===========================================================================*/
int main (int argc, char *argv[])
{
  char *usage = "usage: transcode [-v] [-d tabdir] input_file output_file\n";
  char *version = "transcode V2.3, 02-10-2008\n";

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
  
  {
    int ok;
    /* fill the HDF5 ODIM data structure */
    ok = read_hdf5(&odim_data, argv[1]);
  }

  /******* recode data descriptor- and data-section now */
      
  {
    struct tm * local = gmtime(&odim_data.what.nominal_time);
    s1.year = local->tm_year+1900; /* we use BUFR edition 4 */
    s1.mon = local->tm_mon+1;
    s1.day = local->tm_mday;
    s1.hour = local->tm_hour;
    s1.min = local->tm_min;
    s1.sec = local->tm_sec;
  }
  s1.mtab = 0;            /* master table used */
  s1.updsequ = 0;      /* original BUFR message */
  s1.opsec = 0;          /* no optional section */
  s1.dcat = 6;            /* message type */
  s1.dcatst = 0;        /* message subtype */
  s1.idcatst = 0;        /* message subtype */
  {
    int status;
    if (loc == 8) {
      status = recode_opera_local (&odim_data, argv[2], &s1);
    } else { 
      status = recode_opera_local_v2_1 (&odim_data, argv[2], &s1);
    }
    if (!status) {
      fprintf (stderr, "FATAL: Unable to recode BUFR-message !\n");
      exit (EXIT_FAILURE);
    }
  }
  free_descs();
  exit (EXIT_SUCCESS);
}


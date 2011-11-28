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
#include "compo.h"
#include "odim.h"

odim_comp_t odim_data; /* structure holding ODIM data */

static char *table_dir = NULL;     /* directory for BUFR tables */

#define START  { herr_t (*old_func)(void*);\
      void *old_client_data;\
      H5Eget_auto(&old_func, &old_client_data);\
      H5Eset_auto(NULL, NULL);

#define END H5Eset_auto(old_func, old_client_data);}

static int read_hdf5(odim_comp_t * od, char * file)
{
  hid_t       file_id, dataset_id, dataspace_id, root_id, group_id, attr_id, space_id, ftype;  /* identifiers */
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
  }

  group_id = H5Gopen(root_id, "what");
  od->what.object = H5Areadstring(group_id, "object");
  od->what.version = H5Areadstring(group_id, "version");

  {
    char * date = H5Areadstring(group_id, "date");
    char * time = H5Areadstring(group_id, "time");
    od->what.nominal_time = build_time(date, time);
  }
  {
    char * s;
    START
    s = H5Areadstring(group_id, "source");
    END
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

  status = H5Gclose(group_id);

  group_id = H5Gopen(root_id, "where");

  od->where.projdef = H5Areadstring(group_id, "projdef");

  attr_id = H5Aopen_name(group_id, "xscale");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.psizex);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "yscale");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.psizey);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "xsize");
  status = H5Aread(attr_id, H5T_NATIVE_INT, &od->where.ncols);
  status = H5Aclose(attr_id);
  
  attr_id = H5Aopen_name(group_id, "ysize");
  status = H5Aread(attr_id, H5T_NATIVE_INT, &od->where.nrows);
  status = H5Aclose(attr_id);
  
  attr_id = H5Aopen_name(group_id, "LL_lon");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.sw.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "LL_lat");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.sw.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "UL_lon");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.nw.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "UL_lat");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.nw.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "UR_lon");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.ne.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "UR_lat");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.ne.lat);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "LR_lon");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.se.lon);
  status = H5Aclose(attr_id);

  attr_id = H5Aopen_name(group_id, "LR_lat");
  status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->where.se.lat);
  status = H5Aclose(attr_id);

  fprintf(stderr, "ncols: %d nrows: %d xscale: %g yscale: %g\n", 
	  od->where.ncols, od->where.nrows, od->where.psizex, od->where.psizey);

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

  od->datasets = calloc(od->nscans, sizeof(odim_comp_dataset_t));
  for (i = 0; i < od->nscans; ++i) {
    hid_t g_id;
    hsize_t dims[2];

    char str[20];
    sprintf(str, "dataset%d", i+1);
    group_id = H5Gopen(root_id, str);

    g_id = H5Gopen(group_id, "what");
    od->datasets[i].dataset_what.product = H5Areadstring(g_id, "product");
    assert(strcmp(od->datasets[i].dataset_what.product, "COMP") == 0);

    {
      char * date = H5Areadstring(g_id, "startdate");
      char * time = H5Areadstring(g_id, "starttime");
      od->datasets[i].dataset_what.start_time = build_time(date, time);
    }
    
    {
      char * date = H5Areadstring(g_id, "enddate");
      char * time = H5Areadstring(g_id, "endtime");
      od->datasets[i].dataset_what.end_time = build_time(date, time);
    }

	od->datasets[i].dataset_what.quantity = H5Areadstring(g_id, "quantity");
	if (strcmp(od->datasets[i].dataset_what.quantity, "DBZH") &&
	    strcmp(od->datasets[i].dataset_what.quantity, "RATE") &&
	    strcmp(od->datasets[i].dataset_what.quantity, "ACCR") &&
	    strcmp(od->datasets[i].dataset_what.quantity, "QIND")) {
	  fprintf(stderr, "Warning: Unknown quantity %s\n", od->datasets[i].dataset_what.quantity);
	}
	fprintf(stderr, "Quantity %s\n", od->datasets[i].dataset_what.quantity);

	attr_id = H5Aopen_name(g_id, "gain");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.gain);
	status = H5Aclose(attr_id);

	attr_id = H5Aopen_name(g_id, "offset");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.offset);
	status = H5Aclose(attr_id);

	attr_id = H5Aopen_name(g_id, "nodata");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.nodata);
	status = H5Aclose(attr_id);

	attr_id = H5Aopen_name(g_id, "undetect");
	status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, &od->datasets[i].dataset_what.undetect);
	status = H5Aclose(attr_id);
	fprintf(stderr, "gain: %g offset: %g nodata: %g undetect: %g\n", 
		od->datasets[i].dataset_what.gain,
		od->datasets[i].dataset_what.offset,
		od->datasets[i].dataset_what.nodata,
		od->datasets[i].dataset_what.undetect);

    status = H5Gclose(g_id);

    dims[0] = od->where.nrows;
    dims[1] = od->where.ncols;

    /* Number of params determination */
    {
      int j = 0;
      char str[20];

      START
      while (1) {
	sprintf(str, "data%d", j+1);
	if ((space_id = H5Gopen(group_id, str)) < 0) {
	  break;
	}
	H5Gclose(space_id);
	++j;
	assert(j == 1);
      }
      END
    }

    {
      hid_t plist, data_id;
      char str[20];
      hsize_t cdims[2];
      cdims[0] = 200;
      cdims[1] = 200;
      sprintf(str, "data%d", 1);
      
      od->datasets[i].data = 
	calloc(od->where.nrows*od->where.ncols, sizeof(varfl));
      if (od->datasets[i].data) {
	hid_t g_id;
	data_id = H5Gopen(group_id, str);
	dataset_id = H5Dopen(data_id, "data");
	{
	  START
	  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, od->datasets[i].data);
	  if (status < 0) {
	    unsigned char * p = calloc(od->where.nrows*od->where.ncols, 
				       sizeof(unsigned char));
	    fprintf(stderr, "Warning: trying unsigned char for data %d %d\n", od->where.nrows,
		    od->where.ncols);
	    status = H5Dread(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);
	    if (status == 0) {
	      int k; 
	      fprintf(stderr, "Ok unsigned char for data\n");
	      for (k = 0; k < od->where.ncols*od->where.nrows; k++) {
		od->datasets[i].data[k] = p[k];
	      }
	      free(p);
	    }
	  }
	  END
	}
    
	{
	  int k, n_nodata = 0, n_undetect = 0;
	  for (k = 0; k < od->where.ncols*od->where.nrows; k++) {
	    if (od->datasets[i].data[k] == od->datasets[i].dataset_what.nodata) {
	      od->datasets[i].data[k] = DBL_MAX;
	      ++n_nodata;
	    } else if (od->datasets[i].data[k] == od->datasets[i].dataset_what.undetect) {
	      od->datasets[i].data[k] = -DBL_MAX;
	      ++n_undetect;
	    } else {
	      od->datasets[i].data[k] =
		od->datasets[i].data[k]*od->datasets[i].dataset_what.gain +
		od->datasets[i].dataset_what.offset;
	    }
	  }
	  fprintf(stderr, "nodata: %d undetect: %d\n", n_nodata, n_undetect);
	}

	result = H5Areadstring(dataset_id, "CLASS");
	assert(strcmp(result, "IMAGE") == 0);
	
	START
	result = H5Areadstring(dataset_id, "IMAGE_VERSION");
	END
	if (!result) {
	  fprintf(stderr, "Warning: IMAGE_VERSION missing\n");
	  result = H5Areadstring(dataset_id, "VERSION");
	  if (!result) {
	    fprintf(stderr, "Warning: VERSION missing\n");
	  }
	}
	assert(strcmp(result, "1.2") == 0);
	
	status = H5Dclose(dataset_id);
	fprintf(stderr, "--> %g %g %g\n", od->datasets[i].data[0],
		od->datasets[i].data[1],
		od->datasets[i].data[2]);

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
  s1.vmtab = 14;                    /* version number of master table used */
  s1.vltab = 8;                     /* version number of local table used */

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

static int recode_opera_local (odim_comp_t *od, char *outfile, sect_1_t *s1_in)
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

    {
      struct tm * local = gmtime(&od->what.nominal_time);
      fill_desc(3,1,11);
      fill_v(local->tm_year+1900); 
      fill_v(local->tm_mon+1);
      fill_v(local->tm_mday);
      fill_desc(3,1,13);
      fill_v(local->tm_hour);
      fill_v(local->tm_min);
      fill_v(local->tm_sec);
    }

    fill_desc(3,21,204);
    fprintf(stderr, "nstations %d %d\n", od->what.nstations, flag_omm);
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

    fill_desc(0,29,205);
    for (ii = 0; ii < 100; ++ii) {
      fill_v(od->where.projdef[ii]);
    }

    fill_desc(0,5,33);
    fill_v(od->where.psizex); 
    fill_desc(0,6,33);
    fill_v(od->where.psizey);
    fill_desc(2,1,129);
    fill_desc(0,30,21);
    fill_v(od->where.ncols); 
    fill_desc(0,30,22);
    fill_v(od->where.nrows); 
    fill_desc(2,1,0);
  
    /* Corners ccordinates */
    fill_desc(3,1,21);
    fill_v(od->where.nw.lat);
    fill_v(od->where.nw.lon);
    fill_desc(3,1,21);
    fill_v(od->where.ne.lat);
    fill_v(od->where.ne.lon);
    fill_desc(3,1,21);
    fill_v(od->where.se.lat);
    fill_v(od->where.se.lon);
    fill_desc(3,1,21);
    fill_v(od->where.sw.lat);
    fill_v(od->where.sw.lon);

    fill_desc(1,4,0);
    fill_desc(0,31,1);
    fill_v(od->nscans);

    fill_desc(3,21,205);
    fill_desc(0,30,31);
    fill_desc(0,30,196);
    fill_desc(3,21,206);

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

      if ((strcmp(od->datasets[ii].dataset_what.product, "COMP") == 0)) {
	fill_v(1);
      } else {
	fprintf(stderr, "Warning: Unknown quantity %s\n", od->datasets[ii].dataset_what.product);
	fill_v(1); /* COMP */
      }
      if ((strcmp(od->datasets[ii].dataset_what.quantity, "DBZH") == 0)) {
	fill_v(4);
      } else if ((strcmp(od->datasets[ii].dataset_what.quantity, "RATE") == 0)) {
	fill_v(20);
      } else if ((strcmp(od->datasets[ii].dataset_what.quantity, "ACCR") == 0)) {
	fill_v(21);
      } else if ((strcmp(od->datasets[ii].dataset_what.quantity, "QIND") == 0)) {
	fill_v(255);
      } else {
	fprintf(stderr, "Warning: Unknown quantity %s\n", od->datasets[ii].dataset_what.quantity);
	fill_v(255);
      }

      {
	int i;
	fill_v(0); /* zlib compression */
	{
	  int i, ncomp, n = od->where.nrows*od->where.ncols;
	  unsigned char * result;
	  fprintf(stderr, "rows: %d cols: %d\n", od->where.nrows, od->where.ncols);
	  assert(result = my_compress(od->datasets[ii].data, n, &ncomp));
	  fprintf(stderr, "%d\n", ncomp);
	  fill_v(ncomp/65534 + 1);
	  for (i = 0; i < ncomp/65534; ++i) {
	    fill_v(65534); 
	    fprintf(stderr, "%d\n", 65534);
	    for (n = 0; n < 65534; ++n) {
	      fill_v(result[n + i*65534]); 
	    }
	  }
	  fill_v(ncomp%65534);
	  fprintf(stderr, "reste %d\n", ncomp%65534);
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
  if (!recode_opera_local (&odim_data, argv[2], &s1)) {
    fprintf (stderr, "FATAL: Unable to recode BUFR-message !\n");
    exit (EXIT_FAILURE);
  }
  free_descs();
  exit (EXIT_SUCCESS);
}


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

radar_data_t our_data1, our_data2, our_data3; /* structure holding our decoded data */
odim_polar_t odim_data; /* structure holding ODIM data */

static char *table_dir = NULL;     /* directory for BUFR tables */

static int fill_odim(radar_data_t * b, radar_data_t * b2, radar_data_t * b3, odim_polar_t * d, sect_1_t *s1_in)
{
  /* what */
  d->what.object = "PVOL";
  d->what.version = "H5rad 2.0";
  {
    char * tz = set_fuseau("TZ=UTC");
    struct tm local;
    local.tm_sec = b->meta.sec;
    local.tm_min = b->meta.min;
    local.tm_hour = b->meta.hour;
    local.tm_mday = b->meta.day;
    local.tm_mon = b->meta.month - 1;
    local.tm_year = b->meta.year - 1900;
    local.tm_wday = 0;
    local.tm_yday = 0;    
    local.tm_isdst = 0;
    
    d->what.nominal_time = mktime(&local);
    set_fuseau(tz);
    /* ensures coherence with section 1 */ 
    s1_in->year = b->meta.year - 2000;
    s1_in->mon = b->meta.month;
    s1_in->day = b->meta.day;
    s1_in->hour = b->meta.hour;
    s1_in->min = b->meta.min;
    s1_in->sec = b->meta.sec;
  }

  d->what.nstations = 1;
  d->what.source = calloc(d->what.nstations, sizeof(odim_station_t));
  d->what.source[0].identifier = "WMO";
  d->what.source[0].value = calloc(6, sizeof(char));
  sprintf(d->what.source[0].value, "%02d%03d", b->wmoblock, b->wmostat);

  /* where */
  d->where.lat = b->meta.radar.lat;
  d->where.lon = b->meta.radar.lon;
  d->where.height = b->meta.radar_height;

  /* datasets */
  d->nscans = 3;
  d->datasets = calloc(d->nscans, sizeof(odim_polar_dataset_t));

  /* datasets[0] */
  {
    odim_polar_dataset_t * ds = &d->datasets[0];
    
    /* what */
    ds->dataset_what.product = "SCAN";
    ds->dataset_what.start_time = d->what.nominal_time - 300;
    ds->dataset_what.end_time = d->what.nominal_time;
    
    /* where */
    ds->dataset_where.elangle = b->polar.elangle;
    ds->dataset_where.nbins = 256;
    ds->dataset_where.rstart = 0;
    ds->dataset_where.rscale = b->polar.gate_length*b->polar.ngates;
    ds->dataset_where.nrays = 720;
    ds->dataset_where.a1gate = 0;
    
    /* parameters */
    ds->nparams = 1;
    ds->dataset_data = calloc(ds->nparams, sizeof(odim_polar_dataset_data_t));
    
    /* parameters[0] */
    {
      odim_polar_dataset_data_t * p = &ds->dataset_data[0];
      
      /* what */
      p->dataset_what.quantity = "DBZH";
      p->dataset_what.gain = 1;
      p->dataset_what.offset = 0;
      p->dataset_what.nodata = DBL_MAX;
      p->dataset_what.undetect = -DBL_MAX;
      
      /* data */
      p->data = calloc(ds->dataset_where.nbins*ds->dataset_where.nrays, sizeof(varfl));
      {
	int i;
	srand48(0);
	for (i = 0; i < ds->dataset_where.nbins*ds->dataset_where.nrays; ++i) {
	  if (b->img.data[i] == 255) {
	    p->data[i] = p->dataset_what.nodata;
	  } else if (b->img.data[i] < 0) {
	    p->data[i] = p->dataset_what.undetect;
	  } else {
	    varfl low = ((b->img.data[i] == 0) ? 0 : b->img.scale.vals[b->img.data[i] - 1]);
	    varfl high = b->img.scale.vals[b->img.data[i]];
	    varfl ran = drand48();
	    /*fprintf(stderr, "%g %d %g\n", b->img.scale.vals[b->img.data[i]], b->img.data[i], 
	      ((int)((low + ran*(high - low))*2))/2.0);*/
	    p->data[i] = ((int)((low + ran*(high - low))*2))/2.0; /* 0.5 dBZ resolution */
	  }
	}
      }
    }
  }

  /* datasets[1] */
  {
    odim_polar_dataset_t * ds = &d->datasets[1];
    
    /* what */
    ds->dataset_what.product = "SCAN";
    ds->dataset_what.start_time = d->what.nominal_time - 300;
    ds->dataset_what.end_time = d->what.nominal_time;
    
    /* where */
    ds->dataset_where.elangle = b2->polar.elangle;
    ds->dataset_where.nbins = 256;
    ds->dataset_where.rstart = 0;
    ds->dataset_where.rscale = b2->polar.gate_length*b2->polar.ngates;
    ds->dataset_where.nrays = 720;
    ds->dataset_where.a1gate = 0;
    
    /* parameters */
    ds->nparams = 1;
    ds->dataset_data = calloc(ds->nparams, sizeof(odim_polar_dataset_data_t));
    
    /* parameters[0] */
    {
      odim_polar_dataset_data_t * p = &ds->dataset_data[0];
      
      /* what */
      p->dataset_what.quantity = "DBZH";
      p->dataset_what.gain = 1;
      p->dataset_what.offset = 0;
      p->dataset_what.nodata = DBL_MAX;
      p->dataset_what.undetect = -DBL_MAX;
      
      /* data */
      p->data = calloc(ds->dataset_where.nbins*ds->dataset_where.nrays, sizeof(varfl));
      {
	int i;
	srand48(0);
	for (i = 0; i < ds->dataset_where.nbins*ds->dataset_where.nrays; ++i) {
	  if (b2->img.data[i] == 255) {
	    p->data[i] = p->dataset_what.nodata;
	  } else if (b2->img.data[i] < 0) {
	    p->data[i] = p->dataset_what.undetect;
	  } else {
	    varfl low = ((b2->img.data[i] == 0) ? 0 : b2->img.scale.vals[b2->img.data[i] - 1]);
	    varfl high = b2->img.scale.vals[b2->img.data[i]];
	    varfl ran = drand48();
	    /*fprintf(stderr, "%g %d %g\n", b->img.scale.vals[b->img.data[i]], b->img.data[i], 
	      ((int)((low + ran*(high - low))*2))/2.0);*/
	    p->data[i] = ((int)((low + ran*(high - low))*2))/2.0; /* 0.5 dBZ resolution */
	  }
	}
      }
    }
  }

  /* datasets[2] */
  {
    odim_polar_dataset_t * ds = &d->datasets[2];
    
    /* what */
    ds->dataset_what.product = "SCAN";
    ds->dataset_what.start_time = d->what.nominal_time - 300;
    ds->dataset_what.end_time = d->what.nominal_time;
    
    /* where */
    ds->dataset_where.elangle = b3->polar.elangle;
    ds->dataset_where.nbins = 256;
    ds->dataset_where.rstart = 0;
    ds->dataset_where.rscale = b3->polar.gate_length*b3->polar.ngates;
    ds->dataset_where.nrays = 720;
    ds->dataset_where.a1gate = 0;
    
    /* parameters */
    ds->nparams = 1;
    ds->dataset_data = calloc(ds->nparams, sizeof(odim_polar_dataset_data_t));
    
    /* parameters[0] */
    {
      odim_polar_dataset_data_t * p = &ds->dataset_data[0];
      
      /* what */
      p->dataset_what.quantity = "DBZH";
      p->dataset_what.gain = 1;
      p->dataset_what.offset = 0;
      p->dataset_what.nodata = DBL_MAX;
      p->dataset_what.undetect = -DBL_MAX;
      
      /* data */
      p->data = calloc(ds->dataset_where.nbins*ds->dataset_where.nrays, sizeof(varfl));
      {
	int i;
	srand48(0);
	for (i = 0; i < ds->dataset_where.nbins*ds->dataset_where.nrays; ++i) {
	  if (b3->img.data[i] == 255) {
	    p->data[i] = p->dataset_what.nodata;
	  } else if (b3->img.data[i] < 0) {
	    p->data[i] = p->dataset_what.undetect;
	  } else {
	    varfl low = ((b3->img.data[i] == 0) ? 0 : b3->img.scale.vals[b3->img.data[i] - 1]);
	    varfl high = b3->img.scale.vals[b3->img.data[i]];
	    varfl ran = drand48();
	    /*fprintf(stderr, "%g %d %g\n", b->img.scale.vals[b->img.data[i]], b->img.data[i], 
	      ((int)((low + ran*(high - low))*2))/2.0);*/
	    p->data[i] = ((int)((low + ran*(high - low))*2))/2.0; /* 0.5 dBZ resolution */
	  }
	}
      }
    }
  }

  return 1;
}

static int write_opera(sect_1_t *s1_in, dd *descr, int iout, varfl *v, char *outfile)
/* Create a new OPERA BUFR message */
{
  int i, ret;
  sect_1_t s1;       /* Here we store section 1 of BUFR message */
  bufr_t msg;

  /******* Prepare data for Section 1 (output) */
  s1.year = s1_in->year + 2000; /* we use bufr edition 4 */
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
  s1.idcatst = 0;
  s1.vmtab = 13;                    /* version number of master table used */
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

#define fill_string(str,nsize) \
{\
  int i;\
  for (i = 0; i < strlen(str); ++i) {\
    fill_v(str[i]);\
  }\
  for (i = strlen(str); i < nsize; ++i) {\
    fill_v(' ');\
  }\
}

static int recode_opera_local (bufr_t* msg, char *outfile, sect_1_t *s1_in)
/* This function recodes a BUFR-message; returns 1 on success, 0 on a fault.*/
{
  int iout = 0;
  unsigned int jout = 0;

  odim_polar_t *od = &odim_data;

  dd descr[MAX_DESCS]; /* This array must be huge enough to hold all required descriptors */
  varfl *v = 0;        /* This array must be huge enough to hold all corresponding data values */

  /* The data we are interested in is read to array 'v' 
     and the corresponding descriptors are in array 'descr' */

  /* Output phase */
  {
    int ii, jj; 

    /* Station identification */
    fill_desc(3,21,204);
    fill_v(1);
    fill_string("PLC",3);
    fill_string("Abbeville",16);

    fill_desc(3,1,31);
      {
	int omm = atoi(od->what.source[0].value);
	fill_v(omm / 1000); 
	fill_v(omm % 1000);
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
      fill_v(1);
     for (jj = 0; jj < 1; ++jj) {
      int i;
      fill_v(jj*40); /* polar volume reflectivity or radial wind */
      fill_v(0); /* zlib compression */
      if (jj == 0) {
	int i, ncomp, n = od->datasets[ii].dataset_where.nbins*od->datasets[ii].dataset_where.nrays;
	unsigned char * result;
	assert(result = my_compress(od->datasets[ii].dataset_data[jj].data, n, &ncomp));
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
      } else {
	fill_v(0);
      }
     }
    }
  }

  free_descs();

  /* Now building new transcoded BUFR message */
  fprintf (stderr, "Output file header:\n");
  return write_opera(s1_in, descr, iout, v, outfile);

}

static int our_callback1 (varfl val, int ind) {

    radar_data_t* b = &our_data1;   /* our global data structure */
    static int nrepet = 0, nmax;

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
        if  (bufr_check_fxy (d, 3,1,1)) {   
            b->wmoblock = vv[0];
            b->wmostat = vv[1];
        }
        else if  (bufr_check_fxy (d, 3,1,11)) {   
            b->meta.year = vv[0];       /* Date */
            b->meta.month = vv[1];
            b->meta.day = vv[2];
        }
        else if  (bufr_check_fxy (d, 3,1,13)) {   
            b->meta.hour = vv[0];       /* Time */
            b->meta.min = vv[1];
            b->meta.sec = vv[2];
        }
        /* Latitude, longitude of station */
        else if (bufr_check_fxy (d, 3,1,21)) { 
            b->meta.radar.lat = vv[0];
            b->meta.radar.lon = vv[1];
        }
        else if (bufr_check_fxy (d, 3,21,6)) { 
            b->polar.gate_length = vv[0];
            b->polar.ngates = vv[1];
	}
        else if (bufr_check_fxy (d, 3,21,10)) { 
            b->meta.radar_height = vv[1];
	}
        /* slicing table */
        else if (bufr_check_fxy (d, 3,21,193)) {
            int j;
            b->img.scale.nvals = vv[0] - 1;  /* number of scale values */ 
            assert(b->img.scale.nvals < 256);
            b->img.scale.vals[0] = vv[3];
            for (j = 0; j < b->img.scale.nvals; j++) {
                b->img.scale.vals[j+1] = vv[6+3*j];
            }
        }
        /* elevations */
        else if (bufr_check_fxy (d, 3,21,196)) {
            /* vv[0] number of elevations */ 
            assert(vv[0] == 1);
	    b->polar.elangle = vv[1];
        }
        else {
            fprintf (stderr,
                     "Unknown sequence descriptor %d %d %d\n", d->f, d->x, d->y);
        }
        /* close the global value array */
        bufr_close_val_array ();
    }
    /* element descriptor */
    else if (des[ind]->id == ELDESC) {
        dd *d = &(des[ind]->el->d);

	if (bufr_check_fxy (d, 0,30,1)) {
	  if (nrepet) {
	    --nrepet;
	    b->img.data[nmax-1-nrepet] = val;
	  }
	} else if (bufr_check_fxy (d, 0,31,192)) {
            nrepet = nmax = val;
	    b->img.data = calloc(nmax, sizeof(unsigned short));
        }
	else {
            fprintf (stderr,
                     "Unknown element descriptor %d %d %d\n", d->f, d->x, d->y);
        }
    }
    return 1;
}

static int our_callback2 (varfl val, int ind) {

    radar_data_t* b = &our_data2;   /* our global data structure */
    static int nrepet = 0, nmax;

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
        if  (bufr_check_fxy (d, 3,1,1)) {   
            b->wmoblock = vv[0];
            b->wmostat = vv[1];
        }
        else if  (bufr_check_fxy (d, 3,1,11)) {   
            b->meta.year = vv[0];       /* Date */
            b->meta.month = vv[1];
            b->meta.day = vv[2];
        }
        else if  (bufr_check_fxy (d, 3,1,13)) {   
            b->meta.hour = vv[0];       /* Time */
            b->meta.min = vv[1];
            b->meta.sec = vv[2];
        }
        /* Latitude, longitude of station */
        else if (bufr_check_fxy (d, 3,1,21)) { 
            b->meta.radar.lat = vv[0];
            b->meta.radar.lon = vv[1];
        }
        else if (bufr_check_fxy (d, 3,21,6)) { 
            b->polar.gate_length = vv[0];
            b->polar.ngates = vv[1];
	}
        else if (bufr_check_fxy (d, 3,21,10)) { 
            b->meta.radar_height = vv[1];
	}
        /* slicing table */
        else if (bufr_check_fxy (d, 3,21,193)) {
            int j;
            b->img.scale.nvals = vv[0] - 1;  /* number of scale values */ 
            assert(b->img.scale.nvals < 256);
            b->img.scale.vals[0] = vv[3];
            for (j = 0; j < b->img.scale.nvals; j++) {
                b->img.scale.vals[j+1] = vv[6+3*j];
            }
        }
        /* elevations */
        else if (bufr_check_fxy (d, 3,21,196)) {
            /* vv[0] number of elevations */ 
            assert(vv[0] == 1);
	    b->polar.elangle = vv[1];
        }
        else {
            fprintf (stderr,
                     "Unknown sequence descriptor %d %d %d\n", d->f, d->x, d->y);
        }
        /* close the global value array */
        bufr_close_val_array ();
    }
    /* element descriptor */
    else if (des[ind]->id == ELDESC) {
        dd *d = &(des[ind]->el->d);

	if (bufr_check_fxy (d, 0,30,1)) {
	  if (nrepet) {
	    --nrepet;
	    b->img.data[nmax-1-nrepet] = val;
	  }
	} else if (bufr_check_fxy (d, 0,31,192)) {
            nrepet = nmax = val;
	    b->img.data = calloc(nmax, sizeof(unsigned short));
        }
	else {
            fprintf (stderr,
                     "Unknown element descriptor %d %d %d\n", d->f, d->x, d->y);
        }
    }
    return 1;
}

static int our_callback3 (varfl val, int ind) {

    radar_data_t* b = &our_data3;   /* our global data structure */
    static int nrepet = 0, nmax;

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
        if  (bufr_check_fxy (d, 3,1,1)) {   
            b->wmoblock = vv[0];
            b->wmostat = vv[1];
        }
        else if  (bufr_check_fxy (d, 3,1,11)) {   
            b->meta.year = vv[0];       /* Date */
            b->meta.month = vv[1];
            b->meta.day = vv[2];
        }
        else if  (bufr_check_fxy (d, 3,1,13)) {   
            b->meta.hour = vv[0];       /* Time */
            b->meta.min = vv[1];
            b->meta.sec = vv[2];
        }
        /* Latitude, longitude of station */
        else if (bufr_check_fxy (d, 3,1,21)) { 
            b->meta.radar.lat = vv[0];
            b->meta.radar.lon = vv[1];
        }
        else if (bufr_check_fxy (d, 3,21,6)) { 
            b->polar.gate_length = vv[0];
            b->polar.ngates = vv[1];
	}
        else if (bufr_check_fxy (d, 3,21,10)) { 
            b->meta.radar_height = vv[1];
	}
        /* slicing table */
        else if (bufr_check_fxy (d, 3,21,193)) {
            int j;
            b->img.scale.nvals = vv[0] - 1;  /* number of scale values */ 
            assert(b->img.scale.nvals < 256);
            b->img.scale.vals[0] = vv[3];
            for (j = 0; j < b->img.scale.nvals; j++) {
                b->img.scale.vals[j+1] = vv[6+3*j];
            }
        }
        /* elevations */
        else if (bufr_check_fxy (d, 3,21,196)) {
            /* vv[0] number of elevations */ 
            assert(vv[0] == 1);
	    b->polar.elangle = vv[1];
        }
        else {
            fprintf (stderr,
                     "Unknown sequence descriptor %d %d %d\n", d->f, d->x, d->y);
        }
        /* close the global value array */
        bufr_close_val_array ();
    }
    /* element descriptor */
    else if (des[ind]->id == ELDESC) {
        dd *d = &(des[ind]->el->d);

	if (bufr_check_fxy (d, 0,30,1)) {
	  if (nrepet) {
	    --nrepet;
	    b->img.data[nmax-1-nrepet] = val;
	  }
	} else if (bufr_check_fxy (d, 0,31,192)) {
            nrepet = nmax = val;
	    b->img.data = calloc(nmax, sizeof(unsigned short));
        }
	else {
            fprintf (stderr,
                     "Unknown element descriptor %d %d %d\n", d->f, d->x, d->y);
        }
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

  sect_1_t s1;
  bufr_t msg, msg_t1, msg_t2, msg_t3;

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
  if (argc != 5) {
    fprintf (stderr, "%s", usage);
    exit (EXIT_FAILURE);
  }

  /* scan t1 */
  
  if (!bufr_read_file(&msg_t1, argv[1])) {
        fprintf (stderr, "FATAL: Unable to read the BUFR-message in %s!\n", argv[1]);
        exit (EXIT_FAILURE);
  }

  /******* decode section 1 */
  fprintf (stderr, "Input file header:\n");
  if (!bufr_decode_sections01(&s1, &msg_t1)) {
    fprintf (stderr, "FATAL: Unable to decode section 1\n");
    exit (EXIT_FAILURE);
  }
  header_dump(&s1);

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
    desch = bufr_open_descsec_r(&msg_t1, NULL);
    ok = (desch >= 0);
    if (ok) ok = (bufr_open_datasect_r(&msg_t1) >= 0);

    /* calculate number of data descriptors  */
    ndescs = bufr_get_ndescs (&msg_t1);

    /* allocate memory and read data descriptors from bitstream */
    if (ok) ok = bufr_in_descsec (&dds, ndescs, desch);

    /* output data to our global data structure */
    if (ok) ok = bufr_parse_out (dds, 0, ndescs - 1, our_callback1, 1);

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }
  free_descs();

  /* scan t2 */
  
  if (!bufr_read_file(&msg_t2, argv[2])) {
        fprintf (stderr, "FATAL: Unable to read the BUFR-message in %s!\n", argv[2]);
        exit (EXIT_FAILURE);
  }

  /******* decode section 1 */
  fprintf (stderr, "Input file header:\n");
  if (!bufr_decode_sections01(&s1, &msg_t2)) {
    fprintf (stderr, "FATAL: Unable to decode section 1\n");
    exit (EXIT_FAILURE);
  }
  header_dump(&s1);

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
    desch = bufr_open_descsec_r(&msg_t2, NULL);
    ok = (desch >= 0);
    if (ok) ok = (bufr_open_datasect_r(&msg_t2) >= 0);

    /* calculate number of data descriptors  */
    ndescs = bufr_get_ndescs (&msg_t2);

    /* allocate memory and read data descriptors from bitstream */
    if (ok) ok = bufr_in_descsec (&dds, ndescs, desch);

    /* output data to our global data structure */
    if (ok) ok = bufr_parse_out (dds, 0, ndescs - 1, our_callback2, 1);

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }
  free_descs();

  /* scan t3 */
  
  if (!bufr_read_file(&msg_t3, argv[3])) {
        fprintf (stderr, "FATAL: Unable to read the BUFR-message in %s!\n", argv[3]);
        exit (EXIT_FAILURE);
  }

  /******* decode section 1 */
  fprintf (stderr, "Input file header:\n");
  if (!bufr_decode_sections01(&s1, &msg_t3)) {
    fprintf (stderr, "FATAL: Unable to decode section 1\n");
    exit (EXIT_FAILURE);
  }
  header_dump(&s1);

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
    desch = bufr_open_descsec_r(&msg_t3, NULL);
    ok = (desch >= 0);
    if (ok) ok = (bufr_open_datasect_r(&msg_t3) >= 0);

    /* calculate number of data descriptors  */
    ndescs = bufr_get_ndescs (&msg_t3);

    /* allocate memory and read data descriptors from bitstream */
    if (ok) ok = bufr_in_descsec (&dds, ndescs, desch);

    /* output data to our global data structure */
    if (ok) ok = bufr_parse_out (dds, 0, ndescs - 1, our_callback3, 1);

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }
  free_descs();

  /* fill the polar ODIM data structure */
  fill_odim(&our_data1, &our_data2, &our_data3, &odim_data, &s1);

  /******* recode data descriptor- and data-section now */
  if (!recode_opera_local (&msg, argv[4], &s1)) {
    fprintf (stderr, "FATAL: Unable to recode BUFR-message !\n");
    exit (EXIT_FAILURE);
  }
  exit (EXIT_SUCCESS);
}


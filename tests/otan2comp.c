/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <math.h>
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

static char *table_dir = NULL;     /* directory for BUFR tables */

static int count = 0;
static char bufr[100];
static bufr_char_0_29_205 (varfl val, int ind) {
  bufr[count++] = val;
  return 1;
}

static int our_callback (varfl val, int ind) {

    odim_comp_t * od = &odim_data;   /* our global data structure */
    static radar_data_t b;
    static int nscans = 0;
    static varfl Za, Zb;

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
	else if (bufr_check_fxy (d, 3,1,12)) { 
	  int i = 0;
	  b.meta.hour = vv[i++];       /* Time */
	  b.meta.min = vv[i++];
	  b.meta.sec = 0;
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
	      od->datasets[nscans].dataset_what.start_time = od->what.nominal_time;
	      od->datasets[nscans].dataset_what.end_time = od->what.nominal_time;
	    }
	}
        else if (bufr_check_fxy (d, 3,21,8)) {
	  Za = vv[1];
	  Zb = vv[2];
	}
        else if (bufr_check_fxy (d, 3,13,10)) {
            int j;
            int i = 0;
            b.img.scale.vals[0] = vv[i++];
            b.img.scale.nvals = vv[i++] + 1;  /* number of scale values */ 
            assert(b.img.scale.nvals < 256);
            for (j = 1; j < b.img.scale.nvals; j++) {
	      fprintf(stderr, "%d: %g mm/h  ", j, vv[i]);
	      if (vv[i] == 0) {
		b.img.scale.vals[j] = 0; i++;
		fprintf(stderr, "%g dBZ\n", b.img.scale.vals[j]);
	      } else {
		b.img.scale.vals[j] = 10*Zb*log10(vv[i++]) + 10*log10(Za);
		fprintf(stderr, "%g dBZ\n", b.img.scale.vals[j]);
	      }
            }
	}
        else if (bufr_check_fxy (d, 3,21,193)) {
	  int i;
	  unsigned short * img = NULL;
	  int nvals, ncols, nrows, status = rldec_to_mem (vv, &img, &nvals, &nrows, &ncols);
	  if (status == 0) {
	    fprintf(stderr,"Error rldec_to_mem\n");
	  }
	  fprintf(stderr,"ncols= %d nrows=%d nvals = %d\n", ncols, nrows, nvals);
	  if ((ncols != od->where.ncols) || (nrows != od->where.nrows)) {
	    fprintf(stderr,"Error consistency for ncols or nrows\n");
	  }
	  od->datasets[nscans].data = calloc(ncols*nrows, sizeof(varfl));
	  for (i = 0; i < ncols*nrows; ++i) { 
	    if (img[i] == 255) {
	      od->datasets[nscans].data[i] = od->datasets[nscans].dataset_what.nodata;
	    } else if (img[i] < 0) {
	      od->datasets[nscans].data[i] = od->datasets[nscans].dataset_what.undetect;
	    } else {
	      varfl low = ((img[i] == 0) ? 0 : b.img.scale.vals[img[i] - 1]);
	      varfl high = b.img.scale.vals[img[i]];
	      varfl ran = drand48();
	      /*fprintf(stderr, "%g %d %g\n", b->img.scale.vals[b->img.data[i]], b->img.data[i], 
		((int)((low + ran*(high - low))*2))/2.0);*/
	      varfl res = 0.5; /* 0.5 dBZ resolution */
	      od->datasets[nscans].data[i] = low + ((int)(ran*(high - low)/res))*res; 
	    }
	  }
	  free(img);
	}
        else if (bufr_check_fxy (d, 3,1,1)) {
	  int j;

 	  /* what */
	  od->what.object = "COMP";
	  od->what.version = "H5rad 2.0";
	  od->what.nstations = 1;
	  od->what.source = calloc(od->what.nstations, sizeof(odim_station_t));
	  for (j = 0; j < od->what.nstations; ++j) {
	    int k;
	    char * typ = "ORG";
	    char * name = "Pilot Data Hub  ";
	    od->what.source[j].identifier = calloc(4, sizeof(char)); 
	    for (k = 0; k < 3; ++k) {
	      od->what.source[j].identifier[k] = typ[k];
	    }
	    od->what.source[j].value = calloc(17, sizeof(char));
	    for (k = 0; k < 16; ++k) {
	      od->what.source[j].value[k] = name[k];
	    }
	  }
	  od->nscans = 1;
	  fprintf (stderr, "Nparams:% d\n", od->nscans);
	  od->datasets = calloc(od->nscans,sizeof(odim_comp_dataset_t));
	  od->datasets[nscans].dataset_what.product = "COMP";
	  fprintf(stderr, "Product %s\n" , od->datasets[nscans].dataset_what.product);
	  od->datasets[nscans].dataset_what.quantity = "DBZH";
	  fprintf(stderr, "Quantity %s\n" , od->datasets[nscans].dataset_what.quantity);
	  od->datasets[nscans].dataset_what.gain = 1;
	  od->datasets[nscans].dataset_what.offset = 0;
	  od->datasets[nscans].dataset_what.nodata = DBL_MAX;
	  od->datasets[nscans].dataset_what.undetect = -DBL_MAX;
	  od->datasets[nscans].data = NULL;
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

      if (bufr_check_fxy (d, 0,5,33)) {
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
      else if (bufr_check_fxy (d, 0,29,199)) {
	/* Semi-major axis or rotation ellipsoid */
	b.proj.majax = val;
	fprintf (stderr, " Semi-major axis %g\n", val);
      }
      else if (bufr_check_fxy (d, 0,29,200)) {
	/* Semi-minor axis or rotation ellipsoid */
	b.proj.minax = val;
	fprintf (stderr, " Semi-minor axis %g\n", val);
      }
      else if (bufr_check_fxy (d, 0,29,193)) {
	/* Longitude Origin */
	b.proj.orig.lon = val;
	fprintf (stderr, " Lon 0 %g\n", val);
      }
      else if (bufr_check_fxy (d, 0,29,194)) {
	/* Latitude Origin */
	b.proj.orig.lat = val;
	fprintf (stderr, " Lat 0 %g\n", val);
      }
      else if (bufr_check_fxy (d, 0,29,195)) {
	/* False Easting */
	b.proj.xoff = val;
	fprintf (stderr, " Easting %g\n", val);
      }
      else if (bufr_check_fxy (d, 0,29,196)) {
	/* False Northing */
	b.proj.yoff = val;
	fprintf (stderr, " Northing %g\n", val);
     }
      else if (bufr_check_fxy (d, 0,29,197)) {
	b.proj.stdpar1 = val;
	fprintf (stderr, " Reference lat. %g\n", val);
	  {
	    int i; 
	    {
	      double Requatorial = b.proj.majax;
	      double Rpolaire = b.proj.minax;
	      double excentricite = 1 - (Rpolaire/Requatorial)*(Rpolaire/Requatorial); 
	  
	      projPJ pj;
	      projUV p, pout;

	      od->where.projdef = calloc(100, sizeof(char));
	      snprintf(od->where.projdef, 100, "+proj=stere +lat_0=%g +lon_0=%g +lat_ts=%g +x_0=%ld +y_0=%ld +a=%.8g +es=%.8g",
		       b.proj.orig.lat, b.proj.orig.lon, b.proj.stdpar1, -b.proj.xoff, -b.proj.yoff, Requatorial, excentricite);
	      for (i = strlen(od->where.projdef); i < 100; ++i) {
		od->where.projdef[i] = ' ';
	      }

	      if (!(pj = pj_init_plus(od->where.projdef))) {
		fprintf(stderr, "FATAL: error init proj %s\n", pj_strerrno(pj_errno));
		return 0;
	      }

	      p.v = 0;
	      p.u = 0;
	      pout = pj_inv(p, pj);
	      od->where.nw.lat = pout.v / DEG_TO_RAD;
	      od->where.nw.lon = pout.u / DEG_TO_RAD;

	      fprintf(stderr, "nw %g %g\n", od->where.nw.lat, od->where.nw.lon);
	      pout.v = od->where.nw.lat*DEG_TO_RAD;
	      pout.u = od->where.nw.lon*DEG_TO_RAD;
	      p = pj_fwd(pout,pj);

	      p.u += od->where.ncols*od->where.psizex; /* ne */
	      pout = pj_inv(p,pj);
	      od->where.ne.lat = pout.v / DEG_TO_RAD;
	      od->where.ne.lon = pout.u / DEG_TO_RAD;
	      fprintf(stderr, "ne %g %g\n", od->where.ne.lat, od->where.ne.lon);

	      p.v -= od->where.nrows*od->where.psizey; /* se */
	      pout = pj_inv(p,pj);
	      od->where.se.lat = pout.v / DEG_TO_RAD;
	      od->where.se.lon = pout.u / DEG_TO_RAD;
	      fprintf(stderr, "se %g %g\n", od->where.se.lat, od->where.se.lon);

	      p.u -= od->where.ncols*od->where.psizex; /* sw */
	      pout = pj_inv(p,pj);
	      od->where.sw.lat = pout.v / DEG_TO_RAD;
	      od->where.sw.lon = pout.u / DEG_TO_RAD;
	      fprintf(stderr, "sw %g %g\n", od->where.sw.lat, od->where.sw.lon);

	      p.v += od->where.nrows*od->where.psizey; /* nw verification */
	      pout = pj_inv(p,pj);
	      fprintf(stderr, "nw again %g %g\n", pout.v / DEG_TO_RAD, pout.u / DEG_TO_RAD);
	    }
	  }
      }  else if (bufr_check_fxy (d, 0,29,201)) {
	/* Projection type */
	b.proj.type = val;
	assert(val == 1); /* stereo */
      }
      else 
	fprintf (stderr,
		 "Unknown element descriptor %d %d %d\n", d->f, d->x, d->y);
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
    /*if (ok) write_hdf5(&odim_data, argv[2]);*/

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }

  free_descs();

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


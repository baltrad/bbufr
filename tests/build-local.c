/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* Interface for the PROJ library of geographic projections */
#include <proj_api.h>

/* Interface definition of the OPERA coder/decoder */
#include "desc.h"
#include "bufr.h"
#include "bitio.h"
#include "rlenc.h"

#include "radar.h"

radar_data_t our_data; /* structure holding our decoded data */

static char *table_dir = NULL;     /* directory for BUFR tables */

static void header_dump(sect_1_t *s1)
{
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
  s1.mtab = s1_in->mtab;            /* master table used */
  s1.subcent = 255;
  s1.gencent = 255;
  s1.updsequ = s1_in->updsequ;      /* original BUFR message */
  s1.opsec = s1_in->opsec;          /* no optional section */
  s1.dcat = s1_in->dcat;            /* message type */
  s1.dcatst = s1_in->dcatst;        /* message subtype */
  s1.vmtab = 13;                    /* version number of master table used */
  s1.vltab = 6;                     /* version number of local table used */

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

static int recode_opera_local (bufr_t* msg, char *outfile, sect_1_t *s1_in)
/* This function recodes a BUFR-message; returns 1 on success, 0 on a fault.*/
{
  int iout = 0;
  unsigned int jout = 0;

  radar_data_t *d = &our_data;

  dd descr[MAX_DESCS]; /* This array must be huge enough to hold all required descriptors */
  varfl *v = 0;        /* This array must be huge enough to hold all corresponding data values */

  /* The data we are interested in is read to array 'v' 
     and the corresponding descriptors are in array 'descr' */

  /* Output phase */
  {
    double rterre = 6370997; /* radius of the spheric earth used */
    char str[100];
    projPJ pj;
    projUV p, pout;

    sprintf(str, "+proj=aeqd +lat_0=%g +lon_0=%g +a=%g +es=0", 
	    d->meta.radar.lat, d->meta.radar.lon, rterre);
    if (!(pj = pj_init_plus(str))) {
      fprintf(stderr, "FATAL: error init proj %s\n", pj_strerrno(pj_errno));
      return 0;
    }
    p.v = d->meta.radar.lat*DEG_TO_RAD;
    p.u = d->meta.radar.lon*DEG_TO_RAD;
    p = pj_fwd(p,pj);

    fill_desc(3,1,1);
      fill_v(d->wmoblock); 
      fill_v(d->wmostat);
  
    fill_desc(3,1,194);
      fill_v(d->meta.year);       /* Date */
      fill_v(d->meta.month);
      fill_v(d->meta.day);

      fill_v(d->meta.hour);       /* Time */
      fill_v(d->meta.min);

      p.u -= d->img.ncols*d->img.psizex/2;
      p.v += d->img.nrows*d->img.psizey/2;
      pout = pj_inv(p,pj);
      /* (3,1,21); */
      fill_v(pout.v / DEG_TO_RAD); /* latitude NW corner */
      fill_v(pout.u / DEG_TO_RAD); /* longitude NW corner */
    
      p.u += d->img.ncols*d->img.psizex;
      pout = pj_inv(p,pj);
      fill_v(pout.v / DEG_TO_RAD); /* latitude NE corner */
      fill_v(pout.u / DEG_TO_RAD); /* longitude NE corner */
    
      p.v -= d->img.nrows*d->img.psizey;
      pout = pj_inv(p,pj);
      fill_v(pout.v / DEG_TO_RAD); /* latitude SE corner */
      fill_v(pout.u / DEG_TO_RAD); /* longitude SE corner */
    
      p.u -= d->img.ncols*d->img.psizex;
      pout = pj_inv(p,pj);
      fill_v(pout.v / DEG_TO_RAD); /* latitude SW corner */
      fill_v(pout.u / DEG_TO_RAD); /* longitude SW corner */
      p.v += d->img.nrows*d->img.psizey;

      /* (0,29,201); */
      fill_v(d->proj.type);
      if (d->proj.type != 4) { /* proj. radar */
	fprintf (stderr, "WARNING: not radar projection, output will be inconsistent\n");
      }
  
      /* (0,5,1); */
      fill_v(d->meta.radar.lat);        /* Latitude of radar */

      /* (0,6,1); */
      fill_v(d->meta.radar.lon);        /* Longitude of radar */

      /* (0,7,1); */
      fill_v(d->meta.radar_height); 

      /* (0,5,33); */
      fill_v(d->img.psizex);            /* Pixel size along x coordinate */

      /* (0,6,33); */
      fill_v(d->img.psizey);            /* Pixel size along y coordinate */

      /* (0,30,21); */
      fill_v(d->img.ncols);             /* Number of pixels per column */

      /* (0,30,22); */
      fill_v(d->img.nrows);             /* Number of pixels per row */
      
    fill_desc(0,30,196);
    fill_v(0);
    fill_desc(3,1,193);
    fill_v(rterre);
    fill_v(rterre);
    fill_v(d->meta.radar.lon);        /* Longitude origin */
    fill_v(d->meta.radar.lat);        /* Latitude origin */
    fill_v(p.u); 
    fill_v(p.v); 
    fill_v(MISSVAL); /* latitude reference unique */
    fill_v(MISSVAL);

    fill_desc(0,30,31);
    fill_v(d->img.type);

    fill_desc(0,29,2);
    fill_v(d->img.grid); 

    fill_desc(0,33,3);
    fill_v(0); 
  }

  if (d->img.scale.nvals) {
    int ii;
    fill_desc(2,1,129);
    fill_desc(3,13,9);
    fill_v(d->img.scale.vals[0]);
    fill_v(d->img.scale.nvals);
    for (ii = 0; ii < d->img.scale.nvals; ++ii) {
      fill_v(d->img.scale.vals[ii+1]);
    }
    fill_desc(2,1,0);
  }

  if (d->img.data) {
    fill_desc(3,21,193);
    if (!rlenc_from_mem (d->img.data, d->img.nrows, d->img.ncols, &v, 
			 &jout)) {
      fprintf (stderr, "WARNING: Error RLL coding\n");
    }
    free(d->img.data);
  }

  free_descs();

  /* Now building new transcoded BUFR message */
  fprintf (stderr, "Output file header:\n");
  return write_opera(s1_in, descr, iout, v, outfile);

}

static int our_callback (varfl val, int ind) {

    radar_data_t* b = &our_data;   /* our global data structure */
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
        }
        /* Meta information */
        else if (bufr_check_fxy (d, 3,1,192)) { 
            int i = 0;
            b->meta.year = vv[i++];       /* Date */
            b->meta.month = vv[i++];
            b->meta.day = vv[i++];
            b->meta.hour = vv[i++];       /* Time */
            b->meta.min = vv[i++];
            b->img.nw.lat = vv[i++];      /* Lat. / lon. of NW corner */
            b->img.nw.lon = vv[i++];
            b->img.ne.lat = vv[i++];      /* Lat. / lon. of NE corner */
            b->img.ne.lon = vv[i++];
            b->img.se.lat = vv[i++];      /* Lat. / lon. of SE corner */
            b->img.se.lon = vv[i++];
            b->img.sw.lat = vv[i++];      /* Lat. / lon. of SW corner */
            b->img.sw.lon = vv[i++];
            b->proj.type = vv[i++];       /* Projection type */
            b->meta.radar.lat = vv[i++];        /* Latitude of radar */
            b->meta.radar.lon = vv[i++];        /* Longitude of radar */
            b->img.psizex = vv[i++];      /* Pixel size along x coordinate */
            b->img.psizey = vv[i++];      /* Pixel size along y coordinate */
            b->img.nrows = vv[i++];     /* Number of pixels per row */
            b->img.ncols = vv[i++];     /* Number of pixels per column */
        }
        /* Latitude, longitude of station */
        else if (bufr_check_fxy (d, 3,1,21)) { 
            b->meta.radar.lat = vv[0];
            b->meta.radar.lon = vv[1];
        }
        /* Latitude, longitude and height of station */
        else if (bufr_check_fxy (d, 3,1,22)) { 
            b->meta.radar.lat = vv[0];
            b->meta.radar.lon = vv[1];
            b->meta.radar_height = vv[2];
        }
        /* Reflectivity scale */
        else if (bufr_check_fxy (d, 3,13,9)) { 
            int j;
            int i = 0;
            b->img.scale.vals[0] = vv[i++];
            b->img.scale.nvals = vv[i++] + 1;  /* number of scale values */ 
            assert(b->img.scale.nvals < 256);
            for (j = 1; j < b->img.scale.nvals; j++) {
                b->img.scale.vals[j] = vv[i++];
            }
        }
        else if (bufr_check_fxy (d, 3,21,10)) { 
            b->meta.radar_height = vv[1];
	}
        else if (bufr_check_fxy (d, 3,21,11)) { 
            b->img.type = vv[0];
            b->img.grid = vv[2];
	}
        /* slicing table */
        else if (bufr_check_fxy (d, 3,21,193)) {
            int j;
            int i = 0;
            b->img.scale.vals[0] = vv[i+3];
            b->img.scale.nvals = vv[i] - 1;  /* number of scale values */ 
            assert(b->img.scale.nvals < 256);
            for (j = 0; j < b->img.scale.nvals; j++) {
                b->img.scale.vals[j+1] = vv[i+6+3*j];
            }
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

        if (bufr_check_fxy (d, 0,5,33))
            b->img.psizex = val;
        else if (bufr_check_fxy (d, 0,6,33))
            b->img.psizey = val;
	else if (bufr_check_fxy (d, 0,30,1)) {
	  if (nrepet) {
	    --nrepet;
	    b->img.data[nmax-1-nrepet] = val;
	  }
	} else if (bufr_check_fxy (d, 0,30,21))
            b->img.ncols = val;
	else if (bufr_check_fxy (d, 0,30,22))
            b->img.nrows = val;
	else if (bufr_check_fxy (d, 0,29,199))
            /* Semi-major axis or rotation ellipsoid */
            b->proj.majax = val;
        else if (bufr_check_fxy (d, 0,29,200))
            /* Semi-minor axis or rotation ellipsoid */
            b->proj.minax = val;
        else if (bufr_check_fxy (d, 0,29,193))
            /* Longitude Origin */
            b->proj.orig.lon = val;
        else if (bufr_check_fxy (d, 0,29,194))
            /* Latitude Origin */
            b->proj.orig.lat = val;
        else if (bufr_check_fxy (d, 0,29,195))
            /* False Easting */
            b->proj.xoff = val;
        else if (bufr_check_fxy (d, 0,29,196))
            /* False Northing */
            b->proj.yoff = val;
        else if (bufr_check_fxy (d, 0,29,197))
            /* 1st Standard Parallel */
            b->proj.stdpar1 = val;
        else if (bufr_check_fxy (d, 0,29,198))
            /* 2nd Standard Parallel */
            b->proj.stdpar2 = val;
        else if (bufr_check_fxy (d, 0,30,31))
            /* Image type */
            b->img.type = val;
        else if (bufr_check_fxy (d, 0,29,1))
            /* Co-ordinate grid */
            b->proj.type = val;
        else if (bufr_check_fxy (d, 0,29,2))
            /* Co-ordinate grid */
            b->img.grid = val;
        else if (bufr_check_fxy (d, 0,33,3))
            /* Quality information */
            b->img.qual = val;
        else if (bufr_check_fxy (d, 0,21,198))
            /* dBZ Value offset */
            b->img.scale.offset = val;
        else if (bufr_check_fxy (d, 0,21,199))
            /* dBZ Value increment */
            b->img.scale.increment = val;
	else if (bufr_check_fxy (d, 0,31,192)) {
            nrepet = nmax = val;
	    b->img.data = calloc(b->img.ncols*b->img.nrows, sizeof(unsigned short));
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

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }

  /******* recode data descriptor- and data-section now */
  if (!recode_opera_local (&msg, argv[2], &s1)) {
    fprintf (stderr, "FATAL: Unable to recode BUFR-message !\n");
    exit (EXIT_FAILURE);
  }
  free_descs();
  exit (EXIT_SUCCESS);
}


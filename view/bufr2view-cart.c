/* Some standard headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <math.h>

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

#include <proj_api.h>

static int kml = 0;
static FILE *fp = NULL;

char * create_filename(char * filename, char * suffix)
{
  char * complete_filename = malloc(strlen(filename) + strlen(suffix) + 1);
  strcpy(complete_filename, filename);
  strcat(complete_filename, suffix);
  return complete_filename;
}

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

#define size_palette 256

void writeInfo(png_structp png_ptr, png_infop info_ptr, odim_polar_t* od, int width, int height){
  int count;
  png_color palette[size_palette] = {{255, 255, 255}, {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, /* White */
			    {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, 
			    {72,  61, 139}, {72,  61, 139}, {72,  61, 139}, {72,  61, 139},     /* DarkSlateBlue */
			    {72,  61, 139}, {72,  61, 139}, {72,  61, 139}, {72,  61, 139}, 
			    {0, 0, 205}, {0, 0, 205}, {0, 0, 205}, {0, 0, 205},                 /* MediumBlue */
			    {30, 144, 255}, {30, 144, 255}, {30, 144, 255}, {30, 144, 255},     /* DodgerBlue */
			    {173, 216, 230}, {173, 216, 230}, {173, 216, 230}, {173, 216, 230}, /* LightBlue */
			    {85, 107,  47}, {85, 107,  47}, {85, 107,  47}, {85, 107,  47},     /* DarkOliveGreen */
			    {60, 179, 113}, {60, 179, 113}, {60, 179, 113}, {60, 179, 113},     /* MediumSeaGreen */
			    {127, 255, 212}, {127, 255, 212}, {127, 255, 212}, {127, 255, 212}, /* Aquamarine */
			    {127, 255,   0}, {127, 255,   0}, {127, 255,   0}, {127, 255,   0}, /* Chartreuse */
			    {255, 255,   0}, {255, 255,   0}, {255, 255,   0}, {255, 255,   0}, /* Yellow */
			    {240, 230, 140}, {240, 230, 140}, {240, 230, 140}, {240, 230, 140}, /* khaki */
			    {222, 184, 135}, {222, 184, 135}, {222, 184, 135}, {222, 184, 135}, /* burlywood */
			    {255, 165,   0}, {255, 165,   0}, {255, 165,   0}, {255, 165,   0}, /* Orange */
			    {160,  82,  45}, {160,  82,  45}, {160,  82,  45}, {160,  82,  45}, /* sienna */
			    {255,   0,   0}};                                                   /* Red */
  png_byte trans[size_palette];
  png_color no_value = {255, 192, 203}; /* Pink */

  png_set_invert_alpha(png_ptr);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8,
      PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  for (count = 64; count < size_palette; count++){
    palette[count].red = 255;
    palette[count].green = 0;
    palette[count].blue = 0;
  }
  for (count = 0; count < size_palette; count++) {
    /* KML files are half transparent */
    trans[count] = (kml ? 127 : 0);
  }
  palette[size_palette-1] = no_value;
  trans[size_palette-1] = 127;
  png_set_tRNS(png_ptr, info_ptr, trans, size_palette, NULL);
  png_set_PLTE(png_ptr, info_ptr, palette, size_palette);
  png_write_info(png_ptr, info_ptr);
}

void writeData(png_structp png_ptr,png_infop info_ptr,odim_polar_t* od, char * filename, int n, int m){
  varfl* p = od->datasets[n].dataset_data[m].data;
  varfl range = od->datasets[n].dataset_where.rstart +
    od->datasets[n].dataset_where.nbins*od->datasets[n].dataset_where.rscale;
  varfl pixsizex, pixsizey;
  int width, height;

  varfl north, south, east, west;

  char str[100];
  projPJ pj, pj2, pj3, pj4;
  projUV pp, pp2, ppN, ppW, ppE, ppS;
 
  sprintf(str, "+proj=polar +lat_0=%g +lon_0=%g +ellps=WGS84", od->where.lat, od->where.lon);
  if (!(pj = pj_init_plus(str))) {
    fprintf(stderr, "FATAL: error init proj %s (%s)\n", str, pj_strerrno(pj_errno));
    exit(1);
  }
  fprintf(stderr, "polar projection: %s rstart: %g a1gate: %d rscale: %g range: %g\n", str, 
	  od->datasets[n].dataset_where.rstart, od->datasets[n].dataset_where.a1gate,
	  od->datasets[n].dataset_where.rscale, range);
  pj4 = pj_latlong_from_proj(pj); 
  fprintf(stderr, "Latlon PROJ equivalent definition of polar projection: %s\n", pj_get_def(pj4, 0));

  ppW.u = -90;
  ppW.v = range;
  pj_transform(pj, pj4, 1, 1, &ppW.u, &ppW.v, NULL);
  fprintf(stderr, "West border: %g %g\n", west = ppW.u/DEG_TO_RAD, ppW.v/DEG_TO_RAD);
  ppN.u = 0;
  ppN.v = range;
  pj_transform(pj, pj4, 1, 1, &ppN.u, &ppN.v, NULL);
  fprintf(stderr, "North border: %g %g\n", ppN.u/DEG_TO_RAD, north = ppN.v/DEG_TO_RAD);
  ppE.u = 90;
  ppE.v = range;
  pj_transform(pj, pj4, 1, 1, &ppE.u, &ppE.v, NULL);
  fprintf(stderr, "East border: %g %g\n", east = ppE.u/DEG_TO_RAD, ppE.v/DEG_TO_RAD);
  ppS.u = 180;
  ppS.v = range;
  pj_transform(pj, pj4, 1, 1, &ppS.u, &ppS.v, NULL);
  fprintf(stderr, "South border: %g %g\n", ppS.u/DEG_TO_RAD, south = ppS.v/DEG_TO_RAD);

  if (kml) { 
    fprintf(fp, "  <GroundOverlay>\n");
    fprintf(fp, "    <name>Reflectivity at elevation %d (%g degree)</name>\n", 
	    n, od->datasets[n].dataset_where.elangle);
    fprintf(fp, "    <description>Reflectivity at elevation %d (%g degree)</description>\n", 
	    n, od->datasets[n].dataset_where.elangle);
    fprintf(fp, "    <Icon>\n");
    fprintf(fp, "      <href>%s.png</href>\n", (strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/')+1));
    fprintf(fp, "    </Icon>\n");
    fprintf(fp, "    <LatLonBox>\n");
    fprintf(fp, "      <north>%g</north>\n", north);
    fprintf(fp, "      <south>%g</south>\n", south);
    fprintf(fp, "      <east>%g</east>\n", east);
    fprintf(fp, "      <west>%g</west>\n", west);
    fprintf(fp, "    </LatLonBox>\n");
    fprintf(fp, "  </GroundOverlay>\n");
  }

  if (kml) {
    pj2 = pj4;
    pj3 = pj_latlong_from_proj(pj2);

    /* SW corner in lat/lon */
    pp2.u = west*DEG_TO_RAD;
    pp2.v = south*DEG_TO_RAD;
    pp.u = east*DEG_TO_RAD;
    pp.v = north*DEG_TO_RAD;
    width = 2*od->datasets[n].dataset_where.nbins;
    height = 2*od->datasets[n].dataset_where.nbins;
    pixsizex = (pp.u - pp2.u)/width;
    pixsizey = (pp.v - pp2.v)/height;
  } else {
    sprintf(str, "+proj=aeqd +lat_0=%g +lon_0=%g +ellps=WGS84", od->where.lat, od->where.lon);
    if (!(pj2 = pj_init_plus(str))) {
      fprintf(stderr, "FATAL: error init proj %s (%s)\n", str, pj_strerrno(pj_errno));
      exit(1);
    }
    fprintf(stderr, "aeqd: %s\n", str);
    pj3 = pj_latlong_from_proj(pj2);
    fprintf(stderr, "aeqd PROJ equivalent definition: %s\n", pj_get_def(pj3, 0));

    /* SW corner in aeqd */
    pp2.u = - range;
    pp2.v = - range;
    width = 2*od->datasets[n].dataset_where.nbins;
    height = 2*od->datasets[n].dataset_where.nbins;
    pixsizex = od->datasets[n].dataset_where.rscale;
    pixsizey = od->datasets[n].dataset_where.rscale;
  }
  
  fprintf(stderr, "x resolution for %d pixels: %g\n", width, pixsizex);
  fprintf(stderr, "y resolution for %d pixels: %g\n", height, pixsizey);

  {
  unsigned char * row = calloc(width*height, sizeof(unsigned char));
  int i, j;

  fprintf(stderr, "Image size: %dx%d name: %s\n", width, width, filename);

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      int k, l, jj; 

      pp.u = (i+0.5)*pixsizex + pp2.u;
      pp.v = (j+0.5)*pixsizey + pp2.v;
      pj_transform(pj2, pj, 1, 1, &pp.u, &pp.v, NULL);
      k = (pp.u < 0 ? 360 + pp.u : pp.u)*od->datasets[n].dataset_where.nrays/360.0 - 0.5;
      l = (pp.v - od->datasets[n].dataset_where.rstart)/od->datasets[n].dataset_where.rscale - 0.5;
      jj = height - j - 1;
 
      if ((k < 0) || (k >= od->datasets[n].dataset_where.nrays)) {
	fprintf(stderr, " (h) radar: %d %d (%d %d) %g %g\n", k, l, i, jj, pp.u, pp.v);
	row[i + jj*width] = size_palette-1;
	continue;
      }
      if ((l < 0) || (l >= od->datasets[n].dataset_where.nbins)) {
	row[i + jj*width] = size_palette-1;
	continue;
      }
      if (p[l + k*od->datasets[n].dataset_where.nbins] == od->datasets[n].dataset_data[m].dataset_what.nodata) {
	row[i + jj*width] = size_palette-1;
      } else if (p[l + k*od->datasets[n].dataset_where.nbins] == od->datasets[n].dataset_data[m].dataset_what.undetect) {
	row[i + jj*width] = size_palette-2;
      } else if (p[l + k*od->datasets[n].dataset_where.nbins] <= 0) {
	row[i + jj*width] = 0;
      } else if (strcmp(od->datasets[n].dataset_data[m].dataset_what.quantity,"RATE") == 0) {
	double Za = 200, Zb = 1.6;
	row[i + jj*width] = 10*Zb*log10(p[l + k*od->datasets[n].dataset_where.nbins]) + 10*log10(Za);
      } else {
	row[i + jj*width] = p[l + k*od->datasets[n].dataset_where.nbins];
      }
    }
  }

  writeInfo(png_ptr, info_ptr, od, width, height);

  png_bytep row_pointer = row;
  int count = 0;
  
  for (count = 0; count < height; count++){
    png_write_row(png_ptr, row_pointer);
    row_pointer += width;
  }
  free(row);
  }
}

odim_polar_t odim_data; /* structure holding ODIM data */

static int write_png(odim_polar_t * od, char * file)
{
  hid_t       file_id, dataset_id, dataspace_id, root_id, group_id, attr_id, space_id;  /* identifiers */
  herr_t      status;
  int i, j, width;

  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;

  char str[20];
  int n,m;
  if (kml) { 
    char * kmlfile = create_filename(file, ".kml");
    fp = fopen(kmlfile, "wb");
    if (!fp) return;
    free(kmlfile);
    fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fp, "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n");
    fprintf(fp, "<Folder>\n");
    fprintf(fp, "  <name>Radar data at station %s (coding type: %s)</name>\n", 
	    od->what.source[0].value, od->what.source[0].identifier);
    {
      char str[100];
      struct tm * local = gmtime(&od->what.nominal_time);
      sprintf(str, "%04d%02d%02d %02d:%02d:%02d UTC", local->tm_year+1900, local->tm_mon+1, local->tm_mday,
	      local->tm_hour, local->tm_min, local->tm_sec);
      fprintf(fp, "  <description>Radar reflectivity at date: %s</description>\n", str);
    }
  }
  for (n = 0; n < od->nscans; ++n) {
    sprintf(str, "-%d-elev-%g", n, od->datasets[n].dataset_where.elangle);
    fprintf(stderr, "Elevation %d %g (over %d)\n", n, od->datasets[n].dataset_where.elangle, od->nscans);
    for (m=0;m<od->datasets[n].nparams;m++)
    {
        printf("ici\n");
      char nameqtty[6];
      sprintf(nameqtty,"-%s",od->datasets[n].dataset_data[m].dataset_what.quantity);
      char * filetmp = create_filename(file,nameqtty);
      char * elevfile = create_filename(filetmp, str);
      char * pngfile = create_filename(elevfile, ".png");
      FILE* fp=initPngFile(pngfile, &png_ptr, &info_ptr); 
      free(pngfile);
      if (!fp) return 1;
      
      writeData(png_ptr, info_ptr, od, elevfile, n, m);
      deinit(&fp, &png_ptr, &info_ptr);
      free(elevfile);
    }
  }
  if (kml) { 
    fprintf(fp, "</Folder>\n");
    fprintf(fp, "</kml>\n");
    fclose(fp);
  }
  
  return 1;
}

static void how(odim_polar_how_string_t ** un, odim_polar_how_double_t ** deux, varfl * vv, int * ii)
{
  int k, n, m, i = *ii;
  n = vv[i++];
  *un = NULL;
  for (m = 0; m < n; ++m) {
    odim_polar_how_string_t * old = *un;
    odim_polar_how_string_t * p = malloc(sizeof(odim_polar_how_string_t));
    char * s1 = calloc(16, sizeof(char));
    char * s2 = calloc(16, sizeof(char));
    for (k = 0; k < 16; ++k) {
      s1[k] = vv[i++];
    }
    for (k = 0; k < 16; ++k) {
      s2[k] = vv[i++];
    }
    p->id = s1;
    p->value = s2;
    p->next = old;
    *un = p;
  }
  n = vv[i++];
  *deux = NULL;
  for (m = 0; m < n; ++m) {
    odim_polar_how_double_t * old = *deux;
    odim_polar_how_double_t * p = malloc(sizeof(odim_polar_how_double_t));
    char * s1 = calloc(16, sizeof(char));
    char * s2 = (char *)&p->value;
    for (k = 0; k < 16; ++k) {
      s1[k] = vv[i++];
    }
    for (k = 0; k < 8; ++k) {
      s2[k] = vv[i++];
    }
    p->id = s1;
    p->next = old;
    *deux = p;
  }
  *ii = i;
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
	      } else if(param == 93) {
		p->dataset_what.quantity = "RATE";
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
	/* ODIM */
        else if (bufr_check_fxy (d, 3,21,207)) {
	  int ii, jj, i = 0;
	  how(&od->how_string, &od->how_double, vv, &i);
	  od->nscans = vv[i++];
	  od->datasets = calloc(od->nscans, sizeof(odim_polar_dataset_t));
	  for (ii = 0; ii < od->nscans; ++ii) {
	    odim_polar_dataset_t * ds = &od->datasets[ii];
	    how(&ds->how_string, &ds->how_double, vv, &i);
	    /* what */
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
	      int k; 
	      char * s = calloc(6, sizeof(char)); 
	      for (k = 0; k < 6; ++k) {
		s[k] = vv[i++];
	      }
	      for (k = 5; k >= 0; --k) {
		if (s[k] ==  ' ') {
		  s[k] = 0;
		} else {
		  break;
		}
	      }		
	      ds->dataset_what.product = s;
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
	      int param;
	      how(&p->how_string, &p->how_double, vv, &i);
	      /* what */
	      {
		int k; 
		char * s = calloc(6, sizeof(char)); 
		for (k = 0; k < 6; ++k) {
		  s[k] = vv[i++];
		}
		for (k = 5; k >= 0; --k) {
		  if (s[k] ==  ' ') {
		    s[k] = 0;
		  } else {
		    break;
		  }
		}		
		p->dataset_what.quantity = s;
	      }
	      fprintf(stderr, "Quantity %s\n" , p->dataset_what.quantity);
	      p->dataset_what.gain = 1;
	      p->dataset_what.offset = 0;
	      p->dataset_what.nodata = DBL_MAX;
	      p->dataset_what.undetect = -DBL_MAX;
      
	      param = vv[i++];
	      fprintf(stderr, "%d\n", param);
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
  char *usage = "usage: bufrview-cart [-v] [-k] [-d tabdir] input_file output_file\n";
  char *version = "bufrview-cart V2.3, 02-10-2008\n";

  char *table_dir = NULL;     /* directory for BUFR tables */
  sect_1_t s1;
  bufr_t msg;

/******* check command line parameter */

  while (argc > 1 && *argv[1] == '-') {
    if (*(argv[1] + 1) == 'v')
      fprintf (stderr, "%s", version);
    else if (*(argv[1] + 1) == 'k') {
      kml = 1;
    } else if (*(argv[1] + 1) == 'd') {
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
    if (ok) write_png(&odim_data, argv[2]);

    /* close bitstreams and free descriptor array */
    if (dds != (dd*) NULL)
        free (dds);
    bufr_close_descsec_r (desch);
    bufr_close_datasect_r ();
  }

  free_descs();
  exit (EXIT_SUCCESS);
}


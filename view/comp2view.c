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
#include "compo.h"
#include "odim.h"

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

#define size_palette 67

void writeInfo(png_structp png_ptr, png_infop info_ptr, odim_comp_t* od, int width, int height){
  int count;
  png_color palette[size_palette] = {{255, 255, 255}, {255, 255, 255}, {255, 255, 255}, {255, 255, 255}, /* "White" */
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
			    {255,   0,   0},                                                    /* "Red" */
			    {255, 192, 203},                                                    /* "Pink" */
			    {0, 0, 0}};                                                         /* "Black" */
  png_byte trans[size_palette];

  png_set_invert_alpha(png_ptr);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8,
      PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  for (count = 0; count < size_palette; count++) {
    /* KML files are half transparent */
    trans[count] = (kml ? 127 : 0);
  }
  trans[size_palette-1] = 255;

  png_set_tRNS(png_ptr, info_ptr, trans, size_palette, NULL);
  png_set_PLTE(png_ptr, info_ptr, palette, size_palette);
  png_write_info(png_ptr, info_ptr);
}

void writeData(png_structp png_ptr, png_infop info_ptr, odim_comp_t* od, char * filename){
  varfl* p = od->datasets[0].data;
  
  int height = od->where.nrows;
  int width = od->where.ncols;
  int width0;
  int height0;

  unsigned char * row;

  fprintf(stderr, "projection: %s\n", od->where.projdef);

  if (kml) {
    varfl north = -DBL_MAX, south = DBL_MAX, east = -DBL_MAX, west = DBL_MAX;
    projPJ pj, pj4;
    projUV pp, pp2, ppNW, ppSW, ppNE, ppSE; 

    if (!(pj = pj_init_plus(od->where.projdef))) {
      fprintf(stderr, "FATAL: error init proj %s (%s)\n", od->where.projdef, pj_strerrno(pj_errno));
      exit(1);
    }
    pj4 = pj_latlong_from_proj(pj); 

    ppNW.u = od->where.nw.lon*DEG_TO_RAD;
    ppNW.v = od->where.nw.lat*DEG_TO_RAD;
    fprintf(stderr, "NW corner coord: %g %g\n", od->where.nw.lon, od->where.nw.lat);
    ppNE.u = od->where.ne.lon*DEG_TO_RAD;
    ppNE.v = od->where.ne.lat*DEG_TO_RAD;
    fprintf(stderr, "NE corner coord: %g %g\n", od->where.ne.lon, od->where.ne.lat);
    ppSE.u = od->where.se.lon*DEG_TO_RAD;
    ppSE.v = od->where.se.lat*DEG_TO_RAD;
    fprintf(stderr, "SE corner coord: %g %g\n", od->where.se.lon, od->where.se.lat);
    ppSW.u = od->where.sw.lon*DEG_TO_RAD;
    ppSW.v = od->where.sw.lat*DEG_TO_RAD;
    fprintf(stderr, "SW corner coord: %g %g\n", od->where.sw.lon, od->where.sw.lat);

    if (north < ppNW.v) north = ppNW.v;
    if (north < ppNE.v) north = ppNE.v;
    if (north < ppSE.v) north = ppSE.v;
    if (north < ppSW.v) north = ppSW.v;

    if (south > ppNW.v) south = ppNW.v;
    if (south > ppNE.v) south = ppNE.v;
    if (south > ppSE.v) south = ppSE.v;
    if (south > ppSW.v) south = ppSW.v;
    ppNW.v = ppNE.v = north;
    ppSE.v = ppSW.v = south;

    if (east < ppNW.u) east = ppNW.u;
    if (east < ppNE.u) east = ppNE.u;
    if (east < ppSE.u) east = ppSE.u;
    if (east < ppSW.u) east = ppSW.u;

    if (west > ppNW.u) west = ppNW.u;
    if (west > ppNE.u) west = ppNE.u;
    if (west > ppSE.u) west = ppSE.u;
    if (west > ppSW.u) west = ppSW.u;
    ppSE.u = ppNE.u = east;
    ppNW.u = ppSW.u = west;
  
    fprintf(fp, "  <GroundOverlay>\n");
    fprintf(fp, "    <name>Composite reflectivity</name>\n");
    fprintf(fp, "    <description>Composite reflectivity</description>\n");
    fprintf(fp, "    <Icon>\n");
    fprintf(fp, "      <href>%s.png</href>\n", (strrchr(filename, '/') == NULL ? filename : strrchr(filename, '/')+1));
    fprintf(fp, "    </Icon>\n");
    fprintf(fp, "    <LatLonBox>\n");
    fprintf(fp, "      <north>%g</north>\n", north/DEG_TO_RAD);
    fprintf(fp, "      <south>%g</south>\n", south/DEG_TO_RAD);
    fprintf(fp, "      <east>%g</east>\n", east/DEG_TO_RAD);
    fprintf(fp, "      <west>%g</west>\n", west/DEG_TO_RAD);
    fprintf(fp, "    </LatLonBox>\n");
    fprintf(fp, "  </GroundOverlay>\n");
    fprintf(fp, "  <Placemark><name>NW</name>\n");
    fprintf(fp, "  <description>Corner of original domain</description>\n");
    fprintf(fp, "  <Point>\n");
    fprintf(fp, "    <coordinates>%g,%g,0</coordinates>\n", od->where.nw.lon, od->where.nw.lat);
    fprintf(fp, "  </Point>\n");
    fprintf(fp, "  </Placemark>\n");
    fprintf(fp, "  <Placemark><name>NE</name>\n");
    fprintf(fp, "  <description>Corner of original domain</description>\n");
    fprintf(fp, "  <Point>\n");
    fprintf(fp, "    <coordinates>%g,%g,0</coordinates>\n", od->where.ne.lon, od->where.ne.lat);
    fprintf(fp, "  </Point>\n");
    fprintf(fp, "  </Placemark>\n");
    fprintf(fp, "  <Placemark><name>SE</name>\n");
    fprintf(fp, "  <description>Corner of original domain</description>\n");
    fprintf(fp, "  <Point>\n");
    fprintf(fp, "    <coordinates>%g,%g,0</coordinates>\n", od->where.se.lon, od->where.se.lat);
    fprintf(fp, "  </Point>\n");
    fprintf(fp, "  </Placemark>\n");
    fprintf(fp, "  <Placemark><name>SW</name>\n");
    fprintf(fp, "  <description>Corner of original domain</description>\n");
    fprintf(fp, "  <Point>\n");
    fprintf(fp, "    <coordinates>%g,%g,0</coordinates>\n", od->where.sw.lon, od->where.sw.lat);
    fprintf(fp, "  </Point>\n");
    fprintf(fp, "  </Placemark>\n");

    fprintf(stderr, "x resolution for %d pixels: %g\n", width, (ppNE.u - ppNW.u)/width);
    fprintf(stderr, "y resolution for %d pixels: %g\n", height, (ppNW.v - ppSW.v)/height);
    
    width0 = width*1.3; /*(ppNE.u - ppNW.u)/od->where.psizex;*/
    height0 = height*1.3; /*(ppNW.v - ppSW.v)/od->where.psizey;*/

    fprintf(stderr, "x resolution for %d pixels: %g\n", width0, (ppNE.u - ppNW.u)/width0);
    fprintf(stderr, "y resolution for %d pixels: %g\n", height0, (ppNW.v - ppSW.v)/height0);

    pp2.u = od->where.sw.lon*DEG_TO_RAD;
    pp2.v = od->where.sw.lat*DEG_TO_RAD;
    pj_transform(pj4, pj, 1, 1, &pp2.u, &pp2.v, NULL );
    fprintf(stderr, "False easting/northing: +x_0=%g +y_0=%g\n", pp2.u, pp2.v);

    fprintf(stderr, "Image projection: %s\n", pj_get_def(pj, 0));
    fprintf(stderr, "Latlon PROJ equivalent definition: %s\n", pj_get_def(pj4, 0));

    {
    int i, j;

    row = calloc(height0*width0, sizeof(unsigned char));
    for (i = 0; i < width0; i++) {
      for (j = 0; j < height0; j++) {
	int k, l, jj, ll;

	pp.u = (i+0.5)*(ppNE.u - ppNW.u)/width0 + ppSW.u;
	pp.v = (j+0.5)*(ppNW.v - ppSW.v)/height0 + ppSW.v;
	pj_transform(pj4, pj, 1, 1, &pp.u, &pp.v, NULL);

	k = (pp.u - pp2.u)/od->where.psizex - 0.5;
	l = (pp.v - pp2.v)/od->where.psizey - 0.5;
 
	jj = height0 - j - 1;
	ll = height - l - 1;
	/* not missing value 255 here; so it
	   gives the border of the original domain on the final image */
	if ((k < 0) || (k >= width)) {
	  row[i + jj*width0] = size_palette-1;
	  continue;
	}
	if ((l < 0) || (l >= height)) {
	  row[i + jj*width0] = size_palette-1;
	  continue;
	}
	
	if (p[k + ll*width] == od->datasets[0].dataset_what.nodata) {
	  row[i + jj*width0] = size_palette-2;
	} else if (p[k + ll*width] == od->datasets[0].dataset_what.undetect) {
	  row[i + jj*width0] = size_palette-1;
	} else if (p[k + ll*width] > 255) {
	  /* on French tests files, high dBZ values are used to locate radars positions */
	  row[i + jj*width0] = p[k + ll*width];
	  fprintf(stderr, " radar: %d %d (%d %d) ", k, ll, i, jj);
	  pp.u = (i+0.5)*(ppNE.u - ppNW.u)/width0 + ppSW.u;
	  pp.v = (j+0.5)*(ppNW.v - ppSW.v)/height0 + ppSW.v;
	  fprintf(stderr, "lon/lat coord: %g %g\n", pp.u/DEG_TO_RAD, pp.v/DEG_TO_RAD);
	  fprintf(fp, "  <Placemark><name>Radar</name>\n");
	  fprintf(fp, "  <description>Center of composite pixel nearest to the radar</description>\n");
	  fprintf(fp, "  <Point>\n");
	  fprintf(fp, "    <coordinates>%g,%g,0</coordinates>\n", pp.u/DEG_TO_RAD, pp.v/DEG_TO_RAD);
	  fprintf(fp, "  </Point>\n");
	  fprintf(fp, "  </Placemark>\n");
	} else if (p[k + ll*width] <= 0) {
	  row[i + jj*width0] = 0;
	} else if (strcmp(od->datasets[0].dataset_what.quantity,"RATE") == 0) {
	  double Za = 200, Zb = 1.6;
	  row[i + jj*width0] = 10*Zb*log10(p[k + ll*width]) + 10*log10(Za);
	} else if (strcmp(od->datasets[0].dataset_what.quantity,"ACCR") == 0) {
	  row[i + jj*width0] = 65;
	  if (p[k + ll*width] < 300) row[i + jj*width0] = 63;
	  if (p[k + ll*width] < 200) row[i + jj*width0] = 59;
	  if (p[k + ll*width] < 150) row[i + jj*width0] = 55;
	  if (p[k + ll*width] < 100) row[i + jj*width0] = 51;
	  if (p[k + ll*width] < 75) row[i + jj*width0] = 47;
	  if (p[k + ll*width] < 50) row[i + jj*width0] = 43;
	  if (p[k + ll*width] < 30) row[i + jj*width0] = 39;
	  if (p[k + ll*width] < 20) row[i + jj*width0] = 35;
	  if (p[k + ll*width] < 10) row[i + jj*width0] = 31;
	  if (p[k + ll*width] < 5) row[i + jj*width0] = 27;
	  if (p[k + ll*width] < 2) row[i + jj*width0] = 23;
	  if (p[k + ll*width] < 1) row[i + jj*width0] = 19;
	  if (p[k + ll*width] < 0.5) row[i + jj*width0] = 15;
	  if (p[k + ll*width] < 0.1) row[i + jj*width0] = 7;
	} else {
	  row[i + jj*width0] = p[k + ll*width];
	}
      }
    }
    }
    width = width0;
    height = height0;
  } else {
    int i;

    fprintf(stderr, "Image size: %dx%d name: %s\n", width, height, filename);
    row = calloc(height*width, sizeof(unsigned char));

    for (i = 0; i < height*width; i++) {
      if (p[i] == od->datasets[0].dataset_what.nodata) {
	row[i] = size_palette-1;
      } else if (p[i] == od->datasets[0].dataset_what.undetect) {
	row[i] = size_palette-2;
      } else if (p[i] <= 0) {
	row[i] = 0;
      } else if (strcmp(od->datasets[0].dataset_what.quantity,"RATE") == 0) {
	double Za = 200, Zb = 1.6;
	row[i] = 10*Zb*log10(p[i]) + 10*log10(Za);
      } else if (strcmp(od->datasets[0].dataset_what.quantity,"ACCR") == 0) {
	row[i] = 65;
	if (p[i] < 300) row[i] = 63;
	if (p[i] < 200) row[i] = 59;
	if (p[i] < 150) row[i] = 55;
	if (p[i] < 100) row[i] = 51;
	if (p[i] < 75) row[i] = 47;
	if (p[i] < 50) row[i] = 43;
	if (p[i] < 30) row[i] = 39;
	if (p[i] < 20) row[i] = 35;
	if (p[i] < 10) row[i] = 31;
	if (p[i] < 5) row[i] = 27;
	if (p[i] < 2) row[i] = 23;
	if (p[i] < 1) row[i] = 19;
	if (p[i] < 0.5) row[i] = 15;
	if (p[i] < 0.1) row[i] = 7;
      } else {
	row[i] = p[i];
      }
    }
  }
  
  writeInfo(png_ptr, info_ptr, od, width, height);

  {
    png_bytep row_pointer = row;
    int count = 0;

    for (count = 0; count < height; count++){
      png_write_row(png_ptr, row_pointer);
      row_pointer += width;
    }
  }
  free(row);
}

odim_comp_t odim_data; /* structure holding ODIM data */

static int write_hdf5(odim_comp_t * od, char * file)
{
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;

  if (kml) { 
    char * kmlfile = create_filename(file, ".kml");
    fp = fopen(kmlfile, "wb");
    if (!fp) return 0;
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

  {
    char * pngfile = create_filename(file, (kml ? ".png" : ""));
    FILE* fp=initPngFile(pngfile, &png_ptr, &info_ptr); 
    if (!fp) {
      free(pngfile);
      return 1;
    }
    
    writeData(png_ptr, info_ptr, od, file);
    deinit(&fp, &png_ptr, &info_ptr);
    free(pngfile);
  }
  if (kml) { 
    fprintf(fp, "</Folder>\n");
    fprintf(fp, "</kml>\n");
    fclose(fp);
  }

  return 1;
}

static int count = 0;
static char bufr[100];
static int bufr_char_0_29_205 (varfl val, int ind) {
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
	} else if (val == MISSVAL) {
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


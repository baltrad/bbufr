#ifndef __ODC_ODIM_COMP_H_
#define __ODC_ODIM_COMP_H_

#include "polar.h"

/* A geographic coordinates pair in WGS84 system */
typedef struct odim_point_s {
  varfl lat;      /* Latitude */
  varfl lon;      /* Longitude */
} odim_point_t;

/* ODC data structure for composite data */

/* All units are S.I., except angles in decimal degrees, reflectivity in dB,
   accumulated rainfall in mm and rainfall rate in mm/h */

typedef struct odim_comp_what_s {
  char *object;       /* Type of product (here COMP), only an ODIM format
			 reminder, not used */
  char *version;      /* Value "H5rad 2.0", another not used ODIM reminder */
  time_t nominal_time; /* Nominal time of composite data, for now H, H+15,
			  H+30 or H+45 */
  int nstations;       /* Number of product identifiers */
  odim_station_t *source;  /* An array of length nstations */
} odim_comp_what_t;

/* Global what structure for a composite */

typedef struct odim_comp_dataset_what_s {
  char *product;        /* Value COMP; an ODIM content reminder, not used */
  char *quantity;       /* We expect DBZH (Zh), RATE (mm/h), ACRR (mm) or
			   QIND */
  varfl gain,            /* Values are coded as gain*code + offset */
    offset, 
    nodata,              /* Special code for not radiated (beyond range) */
    undetect;            /* Special code for below detection theshhold */
  time_t start_time;     /* Start time of the earliest scan used */
  time_t end_time;       /* End time of the latest scan used  */
} odim_comp_dataset_what_t;

/* Individual composite structure */

typedef struct odim_comp_dataset_s {
  odim_comp_dataset_what_t dataset_what;    /* Time information about
						    scans used in composite */ 
  varfl *data;          /* Data decoded with gain*code + offset formula; 
			   dimension is nrows*ncols, see where_cart_t
			   structure for their value;
			   data[i + ncols*j] is value for pixel (i,j) */
} odim_comp_dataset_t;

typedef struct odim_comp_where_s {
  char *projdef; /* Initialization string for PROJ4 */
  odim_point_t nw;     /* Northwest corner of the image */
  odim_point_t ne;     /* NE corner */
  odim_point_t se;     /* SE corner */
  odim_point_t sw;     /* SW corner */
  int nrows;      /* Number of pixels per column */
  int ncols;      /* Number of pixels per row */
  varfl psizex;   /* Pixel size along x coordinate */
  varfl psizey;   /* Pixel size along y coordinate */
} odim_comp_where_t;

/* This is our internal composite data structure */
typedef struct odim_comp_composite_s {
  odim_comp_what_t what;
  odim_comp_where_t where;
  int nscans;           /* Number of parameters in a composite
			   (normally 2 to include quality) */
  /* Array of composite structures of dimension nscans */
  odim_comp_dataset_t *datasets;
} odim_comp_t;

/* If we have a structure composite_t called arr for a composite, 
   accessing pixel (i,j) is referencing:
   arr.datasets[0].dataset_data.data[i + ncols*j], 
   where it's assumed index 0 references a radar param; 
   corresponding quality information map would then be:
   arr.datasets[1].dataset_data.data[i + ncols*j]
*/

#endif /* __ODC_ODIM_COMP_H_ */

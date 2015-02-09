#ifndef __ODC_ODIM_POLAR_H_
#define __ODC_ODIM_POLAR_H_

/* ODC data structure for polar data */

/* All units are S.I., except angles in decimal degrees, and reflectivty in dB */

char * Conventions = "ODIM_H5/V2_0"; /* not used for compositing */
char * Conventions_v2_1 = "ODIM_H5/V2_1"; /* not used for compositing */

/* The ODIM data format allows a radar to be identified by various means, as
   WMO identifier is not available for all radars
 */

typedef struct odim_station_s {
  char *identifier; /* Type of product identifier, see table 3 of the ODIM specs */
  unsigned char *value;      /* Value of the identifier: */
} odim_station_t;

typedef struct odim_polar_what_s {
  char * object;           /* Type of product (here PVOL), only an ODIM format
			      reminder, not used */
  char * version;          /* Value "H5rad 2.0", another not used ODIM reminder */
  time_t nominal_time;     /* Nominal time of polar data, for now H, H+15, H+30 or
			      H+45 */
  int nstations;           /* Number of radar identifiers, at least one
                              type of WMO or RAD must be present for polar data */
  odim_station_t *source;  /* An array of length nstations */
} odim_polar_what_t;

typedef struct odim_polar_how_string_s {
  const char * id;
  char * value;
  struct odim_polar_how_string_s * next;
} odim_polar_how_string_t;

typedef struct odim_polar_how_double_s {
  const char * id;
  varfl value;
  struct odim_polar_how_double_s * next;
} odim_polar_how_double_t;

/* A geographic coordinate + height in WGS84 system */
typedef struct odim_geo_point_3d_s {
  varfl lat;      /* Latitude */
  varfl lon;      /* Longitude */
  varfl height;   /* Height above ellipsoid */
} odim_geo_point_3d_t;

/* Some parameters for a given scan, related to how input data is coded */

typedef struct odim_polar_dataset_data_what_s {
  char * quantity;       /* We expect DBZH here (clutter free horizontally polarized Z) */
  varfl gain,            /* Reflectivity values are coded as gain*code + offset */
    offset, 
    nodata,              /* Special code for not radiated (beyond range, or for incomplete scan) */
    undetect;            /* Special code for below detection threshhold or clutter */
} odim_polar_dataset_data_what_t;

typedef struct odim_polar_dataset_data_s {
  odim_polar_how_string_t * how_string;
  odim_polar_how_double_t * how_double;
  odim_polar_dataset_data_what_t dataset_what;
  varfl *data;          /* Scan data decoded with gain*code + offset formula; 
			   dimension is nrays*nbins, see odim_polar_dataset_where_t
			   structure for their value;
			   data[i + nbins*j] is value for bin i and ray j */
} odim_polar_dataset_data_t;

/* Helper structures for a scan */
typedef struct odim_polar_dataset_what_s {
  char *product;        /* Value SCAN normally; an ODIM content reminder, not
			   used */
  time_t start_time;     /* Start time of the scan */
  time_t end_time;       /* End time of the scan */
} odim_polar_dataset_what_t;

typedef struct odim_polar_dataset_where_s {
  varfl elangle,         /* Elevation of scan */
    rstart,              /* Distance of start of first bin */
    rscale;              /* Length of a bin */ 
  int nbins,             /* Number of bins, so range is rstart +
			    nbins*rscale */ 
    nrays,               /* Number of rays MUST cover 360 degrees */
    a1gate;              /* Azimuth of first ray; recall that first bin in scan data array starts 
			    always at North, and moves clockwise */
} odim_polar_dataset_where_t;

/* Individual scan structure */

typedef struct odim_polar_dataset_s {
  odim_polar_how_string_t * how_string;
  odim_polar_how_double_t * how_double;
  odim_polar_dataset_what_t dataset_what;    /* Time information about scan */ 
  odim_polar_dataset_where_t dataset_where;  /* Scan parameters */
  int nparams;                               /* Number of parameters for a scan (currently one) */
  odim_polar_dataset_data_t *dataset_data;   /* Array of size nparams */
} odim_polar_dataset_t;

/* This is our internal polar data structure for one radar */

typedef struct odim_polar_s {
  odim_polar_how_string_t * how_string;
  odim_polar_how_double_t * how_double;
  odim_polar_what_t what;
  odim_geo_point_3d_t where;
  int nscans;                       /* Number of scans */
  odim_polar_dataset_t * datasets;  /* Array of size nscans */
} odim_polar_t;

/* If we have an odim_polar_t array arr of all the available radars, 
   accessing value for radar n, elevation p, bin i and ray j is referencing:
   arr[n].datasets[p].dataset_data[0].data[i + nbins*j], 
   where it's assumed dataset_data[0] references the DBZH param; corresponding
   clutter map would then be:
   arr[n].datasets[p].dataset_data[1].data[i + nbins*j]
*/

#endif /* __ODC_ODIM_POLAR_H_ */

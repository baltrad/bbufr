/* A coordinate pair */

typedef struct point_s {
  varfl lat;      /* latitude */
  varfl lon;      /* longitude */
} point_t;


/* Meta information about image */

typedef struct meta_s {
  int year;
  int month;
  int day;
  int hour;
  int min;
  int sec;
  point_t radar;  /* Radar position */
  varfl radar_height;
} meta_t;

/* Level slicing table */

typedef struct scale_s {
  /* one method: */
  int nvals;       /* number of values in level slicing table */
  varfl vals[255]; /* scale values */
  
  /* another method: */
  varfl offset;    /* offset */
  varfl increment; /* increment */
} scale_t;

/* Radar image */

typedef struct img_s {
  int type;       /* Image type */
  varfl qual;     /* quality indicator */
  int grid;       /* Co-ordinate grid type */
  point_t nw;     /* Northwest corner of the image */
  point_t ne;     /* NE corner */
  point_t se;     /* SE corner */
  point_t sw;     /* SW corner */
  int nrows;      /* Number of pixels per row */
  int ncols;      /* Number of pixels per column */
  varfl psizex;   /* Pixel size along x coordinate */
  varfl psizey;   /* Pixel size along y coordinate */
  scale_t scale;  /* Level slicing table */
  unsigned short* data; /* Image data */
} img_t;

/* Polar parameters */

typedef struct polar_s {
  double gate_length, elangle;
  int ngates;
} polar_t;

/* Projection information */

typedef struct proj_s {
  int type;       /* Projection type */
  varfl majax;    /* Semi-major axis or rotation ellipsoid */
  varfl minax;    /* Semi-minor axis or rotation ellipsoid */
  point_t orig;   /* Projection origin */
  int xoff;       /* False easting */
  int yoff;       /* False northing */
  varfl stdpar1;  /* 1st standard parallel */
  varfl stdpar2;  /* 2nd standard parallel */
} proj_t;


/* This is our internal data structure */

typedef struct radar_data_s {
  int wmoblock;           /* WMO block number */
  int wmostat;            /* WMO station number */
  meta_t meta;            /* Meta information about the product */
  img_t img;              /* Radar reflectivity image */
  proj_t proj;            /* Projection information */
  polar_t polar;          /* Polar parameters */
} radar_data_t;


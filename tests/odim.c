#include <endian.h>
#include "odim.h"

unsigned char * my_compress(varfl * buf, int n, int * ncomp)
{
  int nchar = n*sizeof(varfl);
  unsigned char * CompDataBuff = malloc(compressBound(nchar));
  unsigned long DestBuffSize = compressBound(nchar);
  int i, j;
  unsigned char * tempo = calloc(nchar,sizeof(unsigned char));
  for (i = 0; i < nchar/sizeof(varfl); ++i) {
    unsigned char * p = (unsigned char *)(&buf[i]);
    for (j = 0; j < sizeof(varfl); ++j) {
#if __BYTE_ORDER == __BIG_ENDIAN
      tempo[i*sizeof(varfl)+j] = p[sizeof(varfl) - 1 - j];
      
#else
      tempo[i*sizeof(varfl)+j] = p[j];
#endif
    }
  }
  if (compress(CompDataBuff, &DestBuffSize,
		tempo, nchar) != Z_OK) {
    free(tempo);
    return NULL;
  }
  free(tempo);
  *ncomp = DestBuffSize;
  return CompDataBuff;
}

varfl * my_decompress(unsigned char * buf, unsigned long n, int * ndecomp)
{
  int i, j;
  unsigned char str[sizeof(varfl)];
  unsigned char * UnCompDataBuff = malloc(*ndecomp);
  unsigned long DestBuffSize = *ndecomp;
  if (uncompress(UnCompDataBuff, &DestBuffSize, buf, n) != Z_OK) {
    return NULL;
  }
#if __BYTE_ORDER == __BIG_ENDIAN
  for (i = 0; i < *ndecomp/sizeof(varfl); ++i) {
    for (j = 0; j < sizeof(varfl); ++j) {
      str[j] = UnCompDataBuff[i*sizeof(varfl)+j];
    }
    for (j = 0; j < sizeof(varfl); ++j) {
      UnCompDataBuff[i*sizeof(varfl)+j] = str[sizeof(varfl) - 1 - j];
    }
  }
#endif
  return (varfl *)UnCompDataBuff;
}

char * set_fuseau(char * nouveau)
{
  if (!nouveau) return nouveau;
  {
    char * ancien = getenv("TZ");
    if (ancien) ancien -= strlen("TZ=");
    putenv(nouveau);
    tzset();
    return ancien;
  }
}

time_t build_time(char * date, char * time)
  {
    char * tz = set_fuseau("TZ=UTC");
    struct tm local;
    time_t tt;

    int d;
    sscanf(time+4, "%2d", &d);
    local.tm_sec = d;
    sscanf(time+2, "%2d", &d);
    local.tm_min = d;
    sscanf(time, "%2d", &d);
    local.tm_hour = d;
    sscanf(date+6, "%2d", &d);
    local.tm_mday = d;
    sscanf(date+4, "%2d", &d);
    local.tm_mon = d - 1;
    sscanf(date, "%4d", &d);
    local.tm_year = d - 1900;
    local.tm_wday = 0;
    local.tm_yday = 0;    
    local.tm_isdst = 0;
    
    tt = mktime(&local);
    set_fuseau(tz);
    return tt;
  }

int split(char * args, char * delim, char* subdelim)
{
  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  int j;

  for (j = 1, str1 = args; ; j++, str1 = NULL) {
    token = strtok_r(str1, delim, &saveptr1);
    if (token == NULL)
      break;
    /*printf("%d: %s\n", j, token);*/
    
    for (str2 = token; ; str2 = NULL) {
      subtoken = strtok_r(str2, subdelim, &saveptr2);
      if (subtoken == NULL)
	break;
      /*printf(" --> %s\n", subtoken);*/
    }
  }
  return j-1;
}

char * getsplit_typ(char * args, char * delim, char * subdelim, int i0)
{
  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  int j;

  for (j = 1, str1 = args; j <= i0 ; j++, str1 = NULL) {
    token = strtok_r(str1, delim, &saveptr1);
    if (token == NULL)
      break;
    /*printf("%d: %s\n", j, token);*/
    
     subtoken = strtok_r(token, subdelim, &saveptr2);
     /*printf(" --> %s\n", subtoken);*/
  }
  return subtoken;
}

unsigned char * getsplit_val(char * args, char * delim, char * subdelim, int i0)
{
  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  int j;

  for (j = 1, str1 = args; j <= i0 ; j++, str1 = NULL) {
    token = strtok_r(str1, delim, &saveptr1);
    if (token == NULL)
      break;
    /*printf("%d: %s\n", j, token);*/
    
    str2 = token;
    subtoken = strtok_r(str2, subdelim, &saveptr2);
    str2 = NULL;
    subtoken = strtok_r(str2, subdelim, &saveptr2);
    /*printf(" --> %s\n", subtoken);*/
  }
  return subtoken;
}

herr_t H5Acreatestring(hid_t root_id, char * name, char * s)
{
  hid_t strtype, attr_id;
  herr_t status;

  if ((strtype=H5Tcopy(H5T_C_S1))<0) {
    return -1;
  }
  /* the length */
  if ((H5Tset_size(strtype,strlen(s)+1)) < 0) {
    return -1;
  }
  attr_id = H5Acreate(root_id, name, strtype, H5Screate(H5S_SCALAR), H5P_DEFAULT);
  status = H5Awrite(attr_id, strtype, s);
  status = H5Aclose(attr_id);
  
  return status;
}

char * H5Areadstring(hid_t root_id, const char * name)
{
  hid_t strtype, attr_id;
  herr_t status;
  char * s;
  size_t size;

  if ((strtype=H5Tcopy(H5T_C_S1))<0) {
    return 0;
  }
  /* the length */
  if ((H5Tset_size(strtype,H5T_VARIABLE) < 0)) {
    return 0;
  }
  attr_id = H5Aopen_name(root_id, name);

  if (attr_id < 0) {
    fprintf(stderr,"Warning: %s inexistant\n", name);
    return 0;
  }
  {
    hid_t type, ftype;
    H5T_class_t type_class;
    htri_t size_var;

    ftype = H5Aget_type(attr_id);

    type_class = H5Tget_class (ftype);   

    /*if (type_class == H5T_STRING) printf ("File datatype has class H5T_STRING\n");*/
    size = H5Tget_size(ftype);

    /*printf(" Size is of the file datatype returned by H5Tget_size %d \n This is a size of char pointer\n Use H5Tis_variable_str call instead \n", size);*/

    if((size_var = H5Tis_variable_str(ftype)) == 1)
      fprintf(stderr, "Warning: %s has variable size\n", name);

    type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);

    s = malloc(size+1);
    status = H5Aread(attr_id, type, s);
    status = H5Tclose(ftype);
    status = H5Tclose(type);
    /*fprintf(stderr, "--%.*s--\n", size, s);*/
  }

  status = H5Aclose(attr_id);
  if (!size) {
    free(s);
    return 0;
  }

  s[size] = '\0';
  return s;
  
}

void header_dump(sect_1_t *s1)
{
  fprintf (stderr, "%5d    BUFR edition                       \n", _bufr_edition);
  fprintf (stderr, "%5d    master table used                  \n", s1->mtab);
  fprintf (stderr, "%5d    subcenter                          \n", s1->subcent);
  fprintf (stderr, "%5d    generating center                  \n", s1->gencent);
  fprintf (stderr, "%5d    original BUFR message              \n", s1->updsequ);
  fprintf (stderr, "%5d    no optional section                \n", s1->opsec);
  fprintf (stderr, "%5d    message type                       \n", s1->dcat);
  fprintf (stderr, "%5d    message subtype                    \n", s1->dcatst);
  fprintf (stderr, "%5d    international message subtype      \n", s1->idcatst);
  fprintf (stderr, "%5d    version number of master table used\n", s1->vmtab);
  fprintf (stderr, "%5d    version number of local table used \n", s1->vltab);
  fprintf (stderr, "%5d    year                               \n", s1->year);
  fprintf (stderr, "%5d    month                              \n", s1->mon);
  fprintf (stderr, "%5d    day                                \n", s1->day);
  fprintf (stderr, "%5d    hour                               \n", s1->hour);
  fprintf (stderr, "%5d    minute                             \n", s1->min);
  fprintf (stderr, "%5d    sec                                \n", s1->sec);

  return;
}


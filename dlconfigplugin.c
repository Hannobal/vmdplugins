/***************************************************************************
 *cr
 *cr       adapted from John Stone's dlpolyhist plugin revision 1.22
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: dlconfigplugin.c,v $
 *      $Author: hannod $       $Locker:  $             $State: Exp $
 *      $Revision: 0.20 $       $Date: 2016/01/28 16:47:49 $
 *
 ***************************************************************************/

#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "molfile_plugin.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

typedef struct {
  FILE *file;
  int numatoms;
  int cellwarnflag;
  char *file_name;
  molfile_atom_t *atomlist;
  #if vmdplugin_ABIVERSION > 10
  molfile_timestep_metadata_t ts_meta;
  #endif
} dlcfgdata;

static void *open_dlcfg_read(const char *filename, const char *filetype, 
                           int *natoms) {
  FILE *fd;
  dlcfgdata *data;
  char fbuffer[4096];
  int scancount, imcon, keytrj;

  fd = fopen(filename, "rb");
  if (!fd) return NULL;

  /* the first line is the title */ 
  if (NULL == fgets(fbuffer, 1024, fd))  
    return NULL;
  /* then follows the line with number of atoms etc */ 
  if (NULL == fgets(fbuffer, 1024, fd))  
    return NULL;
  scancount = sscanf(fbuffer, "%d %d", &keytrj, &imcon);
  if (scancount != 2) {
    fprintf(stderr,"dlcfg structure) invalid header\n");
    return NULL;
  }
  /*read all the atoms to determine number of atoms*/
  *natoms = 0;
  /*skip the cell data*/
  if (imcon > 0) {
    if (NULL == fgets(fbuffer, 1024, fd)) return NULL;
    if (NULL == fgets(fbuffer, 1024, fd)) return NULL;
    if (NULL == fgets(fbuffer, 1024, fd)) return NULL;
  }
  while(NULL != fgets(fbuffer, 1024, fd)) {
    if (0==strlen(fbuffer)) break;
    if (NULL == fgets(fbuffer, 1024, fd)) return NULL;
    /*skip velocities*/
    if (keytrj > 0) {
      if (NULL == fgets(fbuffer, 1024, fd)) return NULL;
    }
    /*skip forces*/
    if (keytrj > 1) {
      if (NULL == fgets(fbuffer, 1024, fd)) return NULL;
    }
    (*natoms)++;
  }
    
  data = (dlcfgdata *) malloc(sizeof(dlcfgdata));
  memset(data, 0, sizeof(dlcfgdata));
  data->file = fd;
  data->numatoms= *natoms;
  data->cellwarnflag = 0;
  
#if vmdplugin_ABIVERSION > 10
  data->ts_meta.count = -1;
  if(keytrj>0) {
    data->ts_meta.has_velocities = 1;
  } else {
    data->ts_meta.has_velocities = 0;
  }
#endif

  rewind(data->file); /* prepare for first read_timestep call */
  return data;
}

static int read_dlcfg_structure(void *mydata, int *optflags,
                              molfile_atom_t *atoms) {
  char fbuffer[4096], buf[4096];
  float x, y, z;
  int i, atomcount, keytrj, imcon, scancount;
  dlcfgdata *data = (dlcfgdata *)mydata;
  
  /* we don't read any optional data */
  *optflags = MOLFILE_NOOPTIONS;
  
  /* the first line is the title */ 
  if (NULL == fgets(fbuffer, 1024, data->file))
    return MOLFILE_EOF;
  /* then follows the line with number of atoms etc */ 
  if (NULL == fgets(fbuffer, 1024, data->file))
    return MOLFILE_EOF;
  atomcount = 0;
  scancount = sscanf(fbuffer, "%d %d %d", &keytrj, &imcon, &atomcount);
  if (scancount < 2) {
    fprintf(stderr,"dlcfg structure) unrecognized header record\n");
    return MOLFILE_ERROR;
  }
  /* check atom count */
  if (atomcount != data->numatoms && atomcount != 0) {
    fprintf(stderr,"dlcfg structure) WARNING: mismatched atom count in CONFIG file\n");
  }
  
  /* read periodic cell vectors */
  if (imcon > 0) {
    float vec1[3];
    float vec2[3];
    float vec3[3];
    
    /* eat the data but don't use it for anything */
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f\n", &vec1[0], &vec1[1], &vec1[2]);
    if (3 != scancount) {
      fprintf(stderr,"dlcfg structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f\n", &vec2[0], &vec2[1], &vec2[2]);
    if (3 != scancount) {
      fprintf(stderr,"dlcfg structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f\n", &vec3[0], &vec3[1], &vec3[2]);
    if (3 != scancount) {
      fprintf(stderr,"dlcfg structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
  }
  
  for (i=0; i<data->numatoms; i++) {
    molfile_atom_t *atom = atoms + i;
    /* read the coordinates */
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%s", buf);
    if (1 != scancount) {
	fprintf(stderr,"dlcfg structure) failed parsing atom name\n");
	return MOLFILE_ERROR;
    }
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f", &x, &y, &z);
    if (3 != scancount) {
	fprintf(stderr,"dlcfg structure) failed parsing atom coordinates\n");
	return MOLFILE_ERROR;
    }
    /* read the velocities */
    if (keytrj > 0) {
      float vx, vy, vz;
      if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
      scancount = sscanf(fbuffer, "%f %f %f", &vx, &vy, &vz);
      if (3 != scancount) {
	  fprintf(stderr,"dlcfg structure) failed parsing atom velocities\n");
	  return MOLFILE_ERROR;
      }
      /* read the forces */
      if (keytrj > 1) {
	float fx, fy, fz;
	if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
	scancount = sscanf(fbuffer, "%f %f %f", &fx, &fy, &fz);
	if (3 != scancount) {
	    fprintf(stderr,"dlcfg structure) failed parsing atom forces\n");
	    return MOLFILE_ERROR;
	}
      }
    }
    strncpy(atom->name, buf, sizeof(atom->name));
    strncpy(atom->type, atom->name, sizeof(atom->type));
    atom->resname[0] = '\0';
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
  }
  rewind(data->file);
  return MOLFILE_SUCCESS;
}

static int read_dlcfg_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
  char fbuffer[4096], buf[4096];
  float x, y, z, vx, vy, vz, fx, fy, fz;
  int i, keytrj, imcon, scancount;
  dlcfgdata *data = (dlcfgdata *)mydata;
  
  /* the first line is the title */ 
  if (NULL == fgets(fbuffer, 1024, data->file))
    return MOLFILE_EOF;
  /* then follows the line with number of atoms etc */ 
  if (NULL == fgets(fbuffer, 1024, data->file))
    return MOLFILE_EOF;
  scancount = sscanf(fbuffer, "%d %d", &keytrj, &imcon);
  if (scancount != 2) {
    fprintf(stderr,"dlcfg structure) invalid header\n");
    return MOLFILE_ERROR;
  }
  
  /* read periodic cell vectors */
  if (imcon > 0) {
    float vec1[3];
    float vec2[3];
    float vec3[3];
    
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f\n", &vec1[0], &vec1[1], &vec1[2]);
    if (3 != scancount) {
      fprintf(stderr,"dlcfg structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f\n", &vec2[0], &vec2[1], &vec2[2]);
    if (3 != scancount) {
      fprintf(stderr,"dlcfg structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f\n", &vec3[0], &vec3[1], &vec3[2]);
    if (3 != scancount) {
      fprintf(stderr,"dlcfg structure) failed reading unit cell basis vectors\n");
      return MOLFILE_ERROR;
    }
    /* check and copy in periodic cell info */
    
    if (ts != NULL) {
      if(imcon<=2) {
	if (fabs(vec1[1]) > 1.0e-10 || fabs(vec1[2]) > 1.0e-10  || 
	    fabs(vec2[0]) > 1.0e-10 || fabs(vec2[2]) > 1.0e-10  || 
	    fabs(vec3[0]) > 1.0e-10 || fabs(vec3[1]) > 1.0e-10 ) {
	  if (data->cellwarnflag != 1)
	    fprintf(stderr,"dlcfg timestep) non-orthogonal cell found when orthogonal cell was expected\n");
	  data->cellwarnflag = 1;
	} else if (vec1[0] < -1.0e-10 || vec2[1] < -1.0e-10 || vec3[2] < -1.0e-10) {
	  if (data->cellwarnflag != 1)
	    fprintf(stderr,"dlcfg timestep) cell vector components smaller than zero\n");
	  data->cellwarnflag = 1;
	} else {
	  ts->A = vec1[0];
	  ts->B = vec2[1];
	  ts->C = vec3[2];
	  if (data->cellwarnflag != 2)
	    fprintf(stderr,"dlcfg timestep) converting orthogonal DLPOLY periodic cell data\n");
	  data->cellwarnflag = 2;
	}
      } else {
	if (fabs(vec1[1]) > 1.0e-10 || fabs(vec1[2]) > 1.0e-10) {
	  if (data->cellwarnflag != 1)
	    fprintf(stderr,"dlcfg timestep) first vector not parallel to x-axis\n");
	  data->cellwarnflag = 1;
	} else if (fabs(vec2[2]) > 1.0e-10) {
	    fprintf(stderr,"dlcfg timestep) second vector not in xy-plane\n");
	  data->cellwarnflag = 1;
	} else if (vec1[0] < -1.0e-10 || vec2[1] < -1.0e-10 || vec3[2] < -1.0e-10) {
	  if (data->cellwarnflag != 1)
	    fprintf(stderr,"dlcfg timestep) cell vector components smaller than zero\n");
	  data->cellwarnflag = 1;
	} else {
	  ts->A = vec1[0];
	  ts->B = sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
	  ts->C = sqrt(vec3[0]*vec3[0]+vec3[1]*vec3[1]+vec3[2]*vec3[2]);
	  ts->alpha = 90.0/M_PI_2 * acos((vec2[0]*vec3[0]+vec2[1]*vec3[1]+vec2[2]*vec3[2])/(ts->B * ts->C));
	  ts->beta  = 90.0/M_PI_2 * acos((vec1[0]*vec3[0]+vec1[1]*vec3[1]+vec1[2]*vec3[2])/(ts->A * ts->C));
	  ts->gamma = 90.0/M_PI_2 * acos((vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2])/(ts->A * ts->B));
	  if (data->cellwarnflag != 2)
	    fprintf(stderr,"dlcfg timestep) converting non-orthogonal DLPOLY periodic cell data\n");
	  data->cellwarnflag = 2;
	}
      }
    }
  }
  
  for (i=0; i<data->numatoms; i++) {
    /* read the coordinates */
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%s", buf);
    if (1 != scancount) {
	fprintf(stderr,"dlcfg timestep) failed parsing atom name\n");
	return MOLFILE_ERROR;
    }
    if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
    scancount = sscanf(fbuffer, "%f %f %f", &x, &y, &z);
    if (3 != scancount) {
	fprintf(stderr,"dlcfg timestep) failed parsing atom coordinates\n");
	return MOLFILE_ERROR;
    }
    /* read the velocities */
    if (keytrj > 0) {
      if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
      scancount = sscanf(fbuffer, "%f %f %f", &vx, &vy, &vz);
      if (3 != scancount) {
	  fprintf(stderr,"dlcfg timestep) failed parsing atom velocities\n");
	  return MOLFILE_ERROR;
      }
      /* read the forces */
      if (keytrj > 1) {
	if (NULL == fgets(fbuffer, 1024, data->file)) return MOLFILE_EOF;
	scancount = sscanf(fbuffer, "%f %f %f", &fx, &fy, &fz);
	if (3 != scancount) {
	    fprintf(stderr,"dlcfg timestep) failed parsing atom forces\n");
	    return MOLFILE_ERROR;
	}
      }
    }
    /* only save coords if we're given a timestep pointer, */
    /* otherwise assume that VMD wants us to skip past it. */
    if (ts != NULL) { 
      int addr = 3 * i;
      ts->coords[addr    ] = x;
      ts->coords[addr + 1] = y;
      ts->coords[addr + 2] = z;
#if vmdplugin_ABIVERSION > 10
      if (ts->velocities != NULL) {
	ts->velocities[3*i] = vx;
	ts->velocities[3*i+1] = vy;
	ts->velocities[3*i+2] = vz;
      }
#endif
    }
  }
  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 10
/***********************************************************/
static int read_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  dlcfgdata *data = (dlcfgdata *)mydata;
  
  meta->count = -1;
  meta->has_velocities = data->ts_meta.has_velocities;
  if (meta->has_velocities) {
    fprintf(stderr,"dlcfg timestep) Importing velocities from "
                      "CONFIG file.\n");
  }
  return MOLFILE_SUCCESS;
}
#endif

static void close_dlcfg_read(void *mydata) {
  dlcfgdata *data = (dlcfgdata *)mydata;
  fclose(data->file);
  free(data);
}

static void *open_dlcfg_write(const char *filename, const char *filetype, 
                           int natoms) {
  FILE *fp;
  dlcfgdata *data;

  fp = fopen(filename, "w");
  if (!fp) { 
    fprintf(stderr,"dlcfg timestep) Unable to open CONFIG file %s for writing\n",
            filename);
    return NULL;
  }
  
  data = (dlcfgdata *)malloc(sizeof(dlcfgdata));
  data->numatoms = natoms;
  data->file = fp;
  data->file_name = strdup(filename);
  return data;
}


static int write_dlcfg_structure(void *mydata, int optflags, 
                               const molfile_atom_t *atoms) {
  dlcfgdata *data = (dlcfgdata *)mydata;
  data->atomlist = (molfile_atom_t *)malloc(data->numatoms*sizeof(molfile_atom_t));
  memcpy(data->atomlist, atoms, data->numatoms*sizeof(molfile_atom_t));
  return MOLFILE_SUCCESS;
}

static int write_dlcfg_timestep(void *mydata, const molfile_timestep_t *ts) {
  dlcfgdata *data = (dlcfgdata *)mydata; 
  const molfile_atom_t *atom;
  const float *pos;
  int i;

  fprintf(data->file, " generated by VMD\n"); // title
  fprintf(data->file, " %9d", 0); // CONFIG key
  if(ts->A!=0.0 && ts->B!=0.0 && ts->C!=0.0) {
    if(ts->alpha==90.0 && ts->beta==90.0 && ts->gamma==90.0) {
      fprintf(data->file, " %9d\n", 2); // periodic key
      fprintf(data->file, " %19.12e %19.12e %19.12e\n", ts->A,0.0,0.0);
      fprintf(data->file, " %19.12e %19.12e %19.12e\n", 0.0,ts->B,0.0);
      fprintf(data->file, " %19.12e %19.12e %19.12e\n", 0.0,0.0,ts->C);
    } else {
      fprintf(data->file, " %9d\n", 3); // periodic key
      double bx = ts->B*cos(ts->alpha/90.0*M_PI_2);
      double by = ts->B*sin(ts->alpha/90.0*M_PI_2);
      double cx = ts->C*cos(ts->beta/90.0*M_PI_2);
      double cy = (ts->B*ts->C*cos(ts->gamma/90.0*M_PI_2))/by;
      double cz = sqrt(ts->C*ts->C - cx*cx - cy*cy);
      fprintf(data->file, " %19.12e %19.12e %19.12e\n", ts->A,0.0,0.0);
      fprintf(data->file, " %19.12e %19.12e %19.12e\n", bx,by,0.0);
      fprintf(data->file, " %19.12e %19.12e %19.12e\n", cx,cy,cz);
    }
  } else {
    fprintf(data->file, " %9d\n", 0); // periodic key
  }
  
  atom = data->atomlist;
  pos = ts->coords;
  
  for (i = 0; i < data->numatoms; ++i) {
    fprintf(data->file, " %-8s %10d\n", atom->name, i+1);
    fprintf(data->file, "%19.12e %19.12e %19.12e\n", pos[0], pos[1], pos[2]);
    ++atom; 
    pos += 3;
  }
  return MOLFILE_SUCCESS;
}

static void close_dlcfg_write(void *mydata) {
  dlcfgdata *data = (dlcfgdata *)mydata;
  fclose(data->file);
  free(data->atomlist);
  free(data->file_name);
  free(data);
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "dlpolycfg";
  plugin.prettyname = "DLPOLY CONFIG file";
  plugin.author = "Hanno Dietrich";
  plugin.majorv = 0;
  plugin.minorv = 5;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "dlcfg,config,cfg";
  plugin.open_file_read = open_dlcfg_read;
  plugin.read_structure = read_dlcfg_structure;
  plugin.read_next_timestep = read_dlcfg_timestep;
#if vmdplugin_ABIVERSION > 10
  plugin.read_timestep_metadata = read_timestep_metadata;
#endif
  plugin.close_file_read = close_dlcfg_read;
  plugin.open_file_write = open_dlcfg_write;
  plugin.write_structure = write_dlcfg_structure;
  plugin.write_timestep = write_dlcfg_timestep;
  plugin.close_file_write = close_dlcfg_write;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}


#ifdef TEST_PLUGIN

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  void *v;
  int natoms;
  int i, nsets, set;
  
  while (--argc) {
    ++argv;
    v = open_dlcfg_read(*argv, "dlcfg", &natoms);
    if (!v) {
      fprintf(stderr, "open_dlcfg_read failed for file %s\n", *argv);
      return 1;
    }
    fprintf(stderr, "open_dlcfg_read succeeded for file %s\n", *argv);
    fprintf(stderr, "number of atoms: %d\n", natoms);

    i = 0;
    timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
    while (!read_dlcfg_timestep(v, natoms, &timestep)) {
      i++;
    }
    fprintf(stderr, "ended read_next_timestep on frame %d\n", i);

    close_dlcfg_read(v);
#if vmdplugin_ABIVERSION > 10
    read_timestep_metadata(v,&ts_meta);
    if (ts_meta.has_velocities) {
      fprintf(stderr, "found timestep velocities metadata.\n");
    }
    timestep.velocities = (float *) malloc(3*natoms*sizeof(float));
#endif
  }
  return 0;
}

#endif
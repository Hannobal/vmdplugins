
#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "molfile_plugin.h"
#include "periodic_table.h"
#include "hash.h"
#include "inthash.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define LINE_LEN 2048

class lammpsdatadata {
public:
  FILE *file;
  char *file_name;
  int natoms, nattypes;
  int nbonds, nbdtypes;
  int nangles, nangtypes;
  int ndihedrals, ndihtypes;
  int nimpropers, nimptypes;
  int nxtrabonds;
  float xlo, xhi, ylo, yhi, zlo, zhi;
  float xy, xz, yz;
  bool cgcmm, readtimestep, has_velocities;
  int *idlist, *bdat1, *bdat2, *attypeint, *resid;
  inthash_t *idmap;
  
  char (*attypestr)[16], (*resname)[8], (*atname)[8];
  float *pos, *vel, *mass, *charge;
  #if vmdplugin_ABIVERSION > 10
  molfile_timestep_metadata_t ts_meta;
  #endif
  
  lammpsdatadata(void) {
//     fprintf(stderr,"constructor\n");
    natoms = 0;
    nbonds = 0;
    nangles = 0;
    ndihedrals = 0;
    nimpropers = 0;
    nxtrabonds = 0;
    nattypes = 0;
    nbdtypes = 0;
    nangtypes = 0;
    nimptypes = 0;
    has_velocities = false;
    readtimestep = false;
    xlo = 0.0; xhi = 0.0;
    ylo = 0.0; yhi = 0.0;
    zlo = 0.0; zhi = 0.0;
    xy = 0.0; xz = 0.0; yz = 0.0;
    cgcmm = false;
    attypeint = NULL;
    resid = NULL;
    idlist = NULL;
    bdat1 = NULL;
    bdat2 = NULL;
    pos = NULL;
    vel = NULL;
    mass = NULL;
    charge = NULL;
    attypestr=NULL;
    resname=NULL;
    atname=NULL;
    idmap=NULL;
  }
  
  ~lammpsdatadata(void) {
//     fprintf(stderr,"destructor\n");
    free(idlist);
    free(bdat1);
    free(bdat2);
    free(attypeint);
    free(resid);
    free(pos);
    free(vel);
    free(mass);
    free(charge);
    free(attypestr);
    free(resname);
    free(atname);
    if (idmap != NULL) {
      inthash_destroy(idmap);
      free(idmap);
    }
  }
};

bool get_comment(char* in, char* out, int len) {
  char *ptr = strstr(in,"#");
  if(ptr == NULL)
    return false;
  strncpy(out, ++ptr, len-(ptr-in));
}

int is_empty(const char *s) {
  while (*s != '\0') {
    if (!isspace((unsigned char)*s))
      return 0;
    s++;
  }
  return 1;
}

int get_pte_idx_from_mass(float mass) {
  for(int i=0; i<nr_pte_entries; i++) {
    if(fabsf(mass-pte_mass[i])<0.5)
      return i;
  }
  return 0;
}

int find_atomindex(int id, int *arr, int arrlen) {
  for(int i=0; i<arrlen; i++) {
    if(arr[i] == id)
      return i;
  }
  return -1;
}

/* sort for integer map. initial call  id_sort(idmap, 0, natoms - 1);
from lammpsplugin */
static void id_sort(int *idmap, int left, int right)
{
  int pivot, l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = idmap[left];
  
  while (left < right) {
    while ((idmap[right] >= pivot) && (left < right))
      right--;
    if (left != right) {
      idmap[left] = idmap[right];
      left++;
    }
    while ((idmap[left] <= pivot) && (left < right))
      left++;
    if (left != right) {
      idmap[right] = idmap[left];
      right--;
    }
  }
  idmap[left] = pivot;
  pivot = left;
  left = l_hold;
  right = r_hold;

  if (left < pivot)
    id_sort(idmap, left, pivot-1);
  if (right > pivot)
    id_sort(idmap, pivot+1, right);
}

static bool read_lammpsdata_header(void *mydata)
{
//   fprintf(stderr,"lammpsdatadata) read header\n");
  lammpsdatadata *data = (lammpsdatadata *) mydata;
  char fbuffer[LINE_LEN], str1[128], str2[128], str3[128];
  int i;
  float f1, f2, f3;
  
  // the first line is the title
  if (NULL == fgets(fbuffer, 1024, data->file)) return false;
  if (strstr(fbuffer, "CGCMM") != NULL)
    data->cgcmm=true;
  
  // here should be an empty line
  if (NULL == fgets(fbuffer, 1024, data->file)) return false;
  if(is_empty(fbuffer)==0) {
    fprintf(stderr,"lammpsdatadata) non-empty line:\n");
    fprintf(stderr,"lammpsdatadata) %s",fbuffer);
    return false;
  }
  long fpos = ftell(data->file);
  
  while(NULL != fgets(fbuffer, 1024, data->file)) {
    if(is_empty(fbuffer)==1) {
      continue;
    }
    if(sscanf(fbuffer, "%f %f %s %s", &f1, &f2, str1, str2)==4) {
      if(strcmp(str1,"xlo")==0 && strcmp(str2,"xhi")==0) {
        data->xlo=f1;
        data->xhi=f2;
      } else if(strcmp(str1,"ylo")==0 && strcmp(str2,"yhi")==0) {
        data->ylo=f1;
        data->yhi=f2;
      } else if(strcmp(str1,"zlo")==0 && strcmp(str2,"zhi")==0) {
        data->zlo=f1;
        data->zhi=f2;
      } else if(sscanf(fbuffer, "%f %f %f %s %s %s", &f1, &f2, &f3, str1, str2, str3)==6
        && strcmp(str1,"xy")==0 && strcmp(str2,"xz")==0 && strcmp(str3,"yz")==0) {
          data->xy=f1;
          data->xz=f2;
          data->yz=f3;
      } else {
        fprintf(stderr,"lammpsdatadata) unknown cell line in header:\n");
        fprintf(stderr,"lammpsdatadata) %s",fbuffer);
        return false;
      }
    } else if(sscanf(fbuffer, "%d %s", &i, str1)==2) {
      if(strcmp(str1,"atoms")==0)
        data->natoms = i;
      else if(strcmp(str1,"bonds")==0)
        data->nbonds = i;
      else if(strcmp(str1,"angles")==0)
        data->nangles = i;
      else if(strcmp(str1,"dihedrals")==0)
        data->ndihedrals = i;
      else if(strcmp(str1,"impropers")==0)
        data->nimpropers = i;
      else if(strcmp(str1,"atom")==0)
        data->nattypes = i;
      else if(strcmp(str1,"bond")==0)
        data->nbdtypes = i;
      else if(strcmp(str1,"angle")==0)
        data->nangtypes = i;
      else if(strcmp(str1,"dihedral")==0)
        data->ndihtypes = i;
      else if(strcmp(str1,"improper")==0)
        data->nimptypes = i;
      else if(strcmp(str1,"extra")==0)
        data->nxtrabonds = i;
      else {
        fprintf(stderr,"lammpsdatadata) unknown keyword in header:\n");
        fprintf(stderr,"lammpsdatadata) %s\n",str1);
        return false;
      }
    } else {
      // assume the header is done
      break;
    }
    fpos = ftell(data->file);
  }
  fseek(data->file, fpos, SEEK_SET);
  return true;
}

static bool read_mylammpsdata_atoms(void *mydata)
{
//   fprintf(stderr,"lammpsdatadata) read atoms\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  char fbuffer[LINE_LEN], subbuffer[LINE_LEN], str1[128], str2[128];
  int i, index, scancount, tpe, res, addr;
  float chg, px, py, pz;
  
  for(i=0; i<data->natoms; i++) {
    if(NULL == fgets(fbuffer, 1024, data->file)) {
      fprintf(stderr,"lammpsdatadata) unexpected EOF in section Atoms\n");
      return false;
    }
    if(is_empty(fbuffer)==1) {
      fprintf(stderr,"lammpsdatadata) unexpected end of section Atoms\n");
      return false;
    }
    scancount = sscanf(fbuffer, "%d %d %d %f %f %f %f",
      &index, &res, &tpe, &chg, &px, &py, &pz);
    if(scancount != 7) {
      fprintf(stderr,"lammpsdatadata) error parsing atom:\n");
      fprintf(stderr,"%s",fbuffer);
      return false;
    }
    if(index < 1) {
      fprintf(stderr,"lammpsdatadata) invalid atom id in section Atoms\n");
      fprintf(stderr,"%s",fbuffer);
      return false;
    }
    if(tpe > data->nattypes || index < 1) {
      fprintf(stderr,"lammpsdatadata) invalid atom type index in section Atoms\n");
      fprintf(stderr,"%s",fbuffer);
      return false;
    }
    if (inthash_insert(data->idmap, index, i) != HASH_FAIL) {
        fprintf(stderr, "lammpsdatadata) Duplicate atomid %d\n", index);
        return false;
    }
    data->idlist[i] = index;
    data->charge[i] = chg;
    data->resid[i] = res;
    data->attypeint[i] = tpe-1;
    addr = 3 * i;
    data->pos[addr    ] = px;
    data->pos[addr + 1] = py;
    data->pos[addr + 2] = pz;
    if(data->cgcmm==true) {
      if(get_comment(fbuffer,subbuffer,LINE_LEN)==true) {
        scancount = sscanf(subbuffer, "%s %s", &str1, &str2);
        if(scancount==2) {
          strncpy(data->atname[i], str1, 8);
          strncpy(data->resname[i], str2, 8);
        } else if(scancount==1) {
          strncpy(data->atname[i], str1, 8);
          data->resname[i][0]='\0';
        } else {
          data->atname[i][0]='\0';
          data->resname[i][0]='\0';
        }
      } else {
        data->atname[i][0]='\0';
        data->resname[i][0]='\0';
      }
    }
  }
  id_sort(data->idlist, 0, data->natoms-1);
  return true;
}

static bool read_mylammpsdata_velocities(void *mydata)
{
//   fprintf(stderr,"lammpsdatadata) read Velocities\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  char fbuffer[LINE_LEN];
  int i, index, scancount, addr;
  float vx, vy, vz;
  
  for(i=0; i<data->natoms; i++) {
    if(NULL == fgets(fbuffer, 1024, data->file)) {
      fprintf(stderr,"lammpsdatadata) unexpected EOF in section Velocities\n");
      return false;
    }
    if(is_empty(fbuffer)==1) {
      fprintf(stderr,"lammpsdatadata) unexpected end of section Velocities\n");
      return false;
    }
    scancount = sscanf(fbuffer, "%d %f %f %f", &index, &vx, &vy, &vz);
    if(scancount != 4) {
      fprintf(stderr,"lammpsdatadata) error parsing velocities:\n");
      fprintf(stderr,"%s",fbuffer);
      return false;
    }
    addr = 3 * i;
    data->vel[addr    ] = vx;
    data->vel[addr + 1] = vy;
    data->vel[addr + 2] = vz;
  }
  data->has_velocities = true;
  return true;
}

static bool read_mylammpsdata_masses(void *mydata)
{
//   fprintf(stderr,"lammpsdatadata) read_masses\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  char fbuffer[LINE_LEN], subbuffer[LINE_LEN], str[128];
  float mass;
  int i, index;
  for(i=0; i<data->nattypes; i++) {
    if(NULL == fgets(fbuffer, 1024, data->file)) {
      fprintf(stderr,"lammpsdatadata) unexpected EOF in section Masses\n");
      return false;
    }
    if(sscanf(fbuffer, "%d %f", &index, &mass) != 2) {
      fprintf(stderr,"lammpsdatadata) error parsing section Masses:\n");
      fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      fprintf(stderr,"%s%s",fbuffer,"\n");
      return false;
    }
    if(index > data->nattypes || index < 1) {
      fprintf(stderr,"lammpsdatadata) invalid bond index in section Masses:\n");
      fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      return false;
    }
    data->mass[--index] = mass;
    if(data->cgcmm==true) {
      if(get_comment(fbuffer,subbuffer,LINE_LEN)
      && sscanf(subbuffer, "%s", &str) == 1) {
        strncpy(data->attypestr[index], str, 16);
      } else {
        data->attypestr[index][0] = '\0';
      }
    }
  }
  return true;
}


static bool read_mylammpsdata_bonds(void *mydata)
{
//   fprintf(stderr,"lammpsdatadata) read bonds\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  char fbuffer[LINE_LEN];
  int i, a1, a2, b1, b2, tpe, index;
  for(i=0; i<data->nbonds; i++) {
    if(NULL == fgets(fbuffer, 1024, data->file)) {
      fprintf(stderr,"lammpsdatadata) unexpected EOF in section Bonds\n");
      return false;
    }
    if(sscanf(fbuffer, "%d %d %d %d", &index, &tpe, &a1, &a2) != 4) {
      fprintf(stderr,"lammpsdatadata) error parsing section Bonds:\n");
        fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      fprintf(stderr,"%s\n",fbuffer);
      return false;
    }
    b1 = find_atomindex(a1, data->idlist, data->natoms);
    if(b1 == -1) {
      fprintf(stderr,"lammpsdatadata) invalid atom id %d in section Bonds.\n",a1);
      fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      return false;
    }
    b2 = find_atomindex(a2, data->idlist, data->natoms);
    if(b2 == -1) {
      fprintf(stderr,"lammpsdatadata) invalid atom id %d in section Bonds.\n",a2);
      fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      return false;
    }
    data->bdat1[i] = b1+1;
    data->bdat2[i] = b2+1;
  }
  return true;
}

static bool skip_lammpsdata_section(void *mydata, int nentries, char *section)
{
//   fprintf(stderr,"%s%s\n","lammpsdatadata) skip section ",section);
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  char fbuffer[LINE_LEN];
  for(int i=0; i<nentries; i++) {
    if(NULL == fgets(fbuffer, 1024, data->file)) {
      fprintf(stderr, "%s%s%s",
               "lammpsdatadata) unexpected EOF in section ",section,"\n");
      return false;
    }
  }
  return true;
}

static void *open_lammpsdata_read(const char *filename, const char *filetype, 
                           int *natoms)
{
//   fprintf(stderr,"lammpsdatadata) open read\n");
  FILE *fp;
  *natoms = 0;

  fp = fopen(filename, "rb");
  if (!fp) return NULL;
  bool success;
  
  lammpsdatadata *data = new lammpsdatadata;
  data->file = fp;
  data->file_name = strdup(filename);
  
  if(read_lammpsdata_header(data) == false || data->natoms == 0) {
    delete data;
    return NULL;
  }
  
  data->idlist    = new int[data->natoms];
  data->bdat1     = new int[data->nbonds];
  data->bdat2     = new int[data->nbonds];
  data->attypeint = new int[data->natoms];
  data->resid     = new int[data->natoms];
  data->pos       = new float[3*data->natoms];
  data->vel       = new float[3*data->natoms];
  data->mass      = new float[data->nattypes];
  data->charge    = new float[data->natoms];
  data->attypestr = new char[data->nattypes][16];
  data->resname   = new char[data->natoms][8];
  data->atname    = new char[data->natoms][8];
  data->idmap     = (inthash_t *)calloc(1, sizeof(inthash_t));
  inthash_init(data->idmap, data->natoms);
  
  char fbuffer[LINE_LEN], str[128];
  while(NULL != fgets(fbuffer, 1024, data->file)) {
//     ignore empty lines between blocks
    if(is_empty(fbuffer)==1)
      continue;
    if(sscanf(fbuffer,"%s",str) != 1) {
      fprintf(stderr,"lammpsdatadata) weird line:\n");
      fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      delete data;
      return NULL;
    }
    if(NULL == fgets(fbuffer, 1024, data->file)) {
      fprintf(stderr,"lammpsdatadata) unexpected end of file\n");
      delete data;
      return NULL;
    }
    if(is_empty(fbuffer)==0) {
      fprintf(stderr,"lammpsdatadata) expected empty line:\n");
      fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      delete data;
      return NULL;
    }
    if(strcmp(str,"Atoms")==0)
      success = read_mylammpsdata_atoms(data);
    else if(strcmp(str,"Bonds")==0)
      success = read_mylammpsdata_bonds(data);
    else if(strcmp(str,"Masses")==0)
      success = read_mylammpsdata_masses(data);
    else if(strcmp(str,"Velocities")==0)
      success = read_mylammpsdata_velocities(data);
    else if(strcmp(str,"Angles")==0)
      success = skip_lammpsdata_section(data,data->nangles,str);
    else if(strcmp(str,"Dihedrals")==0)
      success = skip_lammpsdata_section(data,data->ndihedrals,str);
    else if(strcmp(str,"Impropers")==0)
      success = skip_lammpsdata_section(data,data->nimpropers,str);
    else if(strcmp(str,"Nonbond")==0)
      success = skip_lammpsdata_section(data,data->nattypes,str);
    else if(strcmp(str,"Bond")==0)
      success = skip_lammpsdata_section(data,data->nbdtypes,str);
    else if(strcmp(str,"Angle")==0)
      success = skip_lammpsdata_section(data,data->nangtypes,str);
    else if(strcmp(str,"Dihedral")==0)
      success = skip_lammpsdata_section(data,data->ndihtypes,str);
    else if(strcmp(str,"Improper")==0)
      success = skip_lammpsdata_section(data,data->nimptypes,str);
    else if(strcmp(str,"BondBond")==0)
      success = skip_lammpsdata_section(data,data->nangtypes,str);
    else if(strcmp(str,"BondAngle")==0)
      success = skip_lammpsdata_section(data,data->nangtypes,str);
    else if(strcmp(str,"MiddleBondTorsion")==0)
      success = skip_lammpsdata_section(data,data->ndihtypes,str);
    else if(strcmp(str,"EndBondTorsion")==0)
      success = skip_lammpsdata_section(data,data->ndihtypes,str);
    else if(strcmp(str,"AngleTorsion")==0)
      success = skip_lammpsdata_section(data,data->ndihtypes,str);
    else if(strcmp(str,"AngleAngleTorsion")==0)
      success = skip_lammpsdata_section(data,data->ndihtypes,str);
    else if(strcmp(str,"BondBond13")==0)
      success = skip_lammpsdata_section(data,data->ndihtypes,str);
    else if(strcmp(str,"AngleAngle")==0)
      success = skip_lammpsdata_section(data,data->nimptypes,str);
    else if(strcmp(str,"Pair")==0)
      success = skip_lammpsdata_section(data,data->nattypes,str);
    else if(strcmp(str,"PairIJ")==0)
      success = skip_lammpsdata_section(data,data->nattypes*(data->nattypes+1)/2,str);
    else {
      fprintf(stderr, "%s%s%s",
        "lammpsdatadata) unknown section keyword ",str,":\n");
        fprintf(stderr,"lammpsdatadata) %s",fbuffer);
      delete data;
      return NULL;
    }
    if(success == false) {
      delete data;
      return NULL;
    }
  }
  *natoms = data->natoms;
  return data;
}

static int read_lammpsdata_structure(void *mydata, int *optflags,
molfile_atom_t *atoms) {
//   fprintf(stderr,"lammpsdatadata) read structure\n");
  lammpsdatadata *data = (lammpsdatadata *) mydata;
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_CHARGE | MOLFILE_MASS;
  for (int i=0; i<data->natoms; i++) {
    molfile_atom_t *atom = atoms+i;
    int j = inthash_lookup(data->idmap, data->idlist[i]);
    if(j == HASH_FAIL) {
      fprintf(stderr,"lammpsdatadata) atom id %d was not found\n",data->idlist[i]);
    }
    atom->mass = data->mass[data->attypeint[j]];
    atom->charge = data->charge[j];
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    atom->resid = data->resid[j];
    int idx = get_pte_idx_from_mass(atom->mass);
    atom->radius = get_pte_vdw_radius(idx);
    atom->atomicnumber = idx;
    if(data->cgcmm==true) {
      strncpy(atom->type,data->attypestr[data->attypeint[j]],16);
      strncpy(atom->name, data->atname[j],8);
      strncpy(atom->resname, data->resname[j],8);
    } else {
      strncpy(atom->name, data->atname[j],8);
      snprintf(atom->type, sizeof(atom->type), "%d", data->attypeint[j]+1);
      snprintf(atom->name, sizeof(atom->name), "%d", data->attypeint[j]+1);
      snprintf(atom->resname, sizeof(atom->resname), "%d", data->resid[j]);
    }
  }
  return MOLFILE_SUCCESS;
}

static int read_lammpsdata_timestep(void *mydata, int natoms, molfile_timestep_t *ts) {
//   fprintf(stderr,"lammpsdatadata) read timestep\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  if(natoms != data->natoms || ts == NULL)
    return MOLFILE_ERROR;
  
  for (int i=0; i < natoms; i++) {
    int j = inthash_lookup(data->idmap, data->idlist[i]);
    if(j == HASH_FAIL) {
      fprintf(stderr,"lammpsdatadata) atom id %d was not found\n",data->idlist[i]);
    }
    memcpy(ts->coords+3*i, data->pos+3*j, 3*sizeof(float));
#if vmdplugin_ABIVERSION > 10
    if(data->has_velocities == true)
      memcpy(ts->velocities+3*i, data->vel+3*j, 3*sizeof(float));
#endif
  }
  if(abs(data->xy) <= 1.0e-10)
    data->xy = 0.0;
  if(abs(data->xz) <= 1.0e-10)
    data->xz = 0.0;
  if(abs(data->yz) <= 1.0e-10)
    data->yz = 0.0;
  float lx = data->xhi - data->xlo;
  float ly = data->yhi - data->ylo;
  float lz = data->zhi - data->zlo;
  ts->A = lx;
  ts->B = sqrt(ly*ly + data->xy * data->xy);
  ts->C = sqrt(lz*lz + data->xz * data->xz + data->yz * data->yz);
  ts->alpha = (90.0/M_PI_2) * acos((data->xy * data->xz + ly * data->yz)/(ts->B * ts->C));
  ts->beta  = (90.0/M_PI_2) * acos(data->xz / ts->C);
  ts->gamma = (90.0/M_PI_2) * acos(data->xy / ts->B);
  if(data->readtimestep==true) {
    return MOLFILE_EOF;
  } else {
    data->readtimestep=true;
  }
  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 10
/***********************************************************/
static int read_lammpsdata_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  
  meta->count = -1;
  if(data->has_velocities==true)
    meta->has_velocities = 1;
  else
    meta->has_velocities = 0;
  return MOLFILE_SUCCESS;
}
#endif

#if vmdplugin_ABIVERSION > 14
static int read_lammpsdata_bonds(void *mydata, int *nbonds, int **from, int **to, 
                           float **bondorder,  int **bondtype, 
                           int *nbondtypes, char ***bondtypename) {
#else
static int read_lammpsdata_bonds(void *mydata, int *nbonds, int **from, int **to, float **bondorder) {
#endif
  
//   fprintf(stderr,"lammpsdatadata) read bonds\n");
  lammpsdatadata *data = (lammpsdatadata *) mydata;
  *nbonds  = data->nbonds;
  *from = data->bdat1;
  *to = data->bdat2;
  *bondorder = NULL;
#if vmdplugin_ABIVERSION > 14
  *bondtype     = NULL;
  *nbondtypes   = 0;
  *bondtypename = NULL;
#endif
  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int read_lammpsdata_angles(void *mydata, int *numangles, int **angles, int **angletypes,
                      int *numangletypes, char ***angletypenames, int *numdihedrals,
                      int **dihedrals, int **dihedraltypes, int *numdihedraltypes,
                      char ***dihedraltypenames, int *numimpropers, int **impropers,        
                      int **impropertypes, int *numimpropertypes, char ***impropertypenames,
                      int *numcterms, int **cterms, int *ctermcols, int *ctermrows) {
#else
static int read_lammpsdata_angles(void *mydata,
                int *numangles,    int **angles,    double **angleforces,
                int *numdihedrals, int **dihedrals, double **dihedralforces,
                int *numimpropers, int **impropers, double **improperforces,
                int *numcterms,    int **cterms, 
                int *ctermcols,    int *ctermrows,  double **ctermforces);
#endif

//   fprintf(stderr,"lammpsdatadata) read angles\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  
  *numangles        = data->nangles;
  *numdihedrals     = data->ndihedrals;
  *numimpropers     = data->nimpropers;
  *angles           = NULL;
  *dihedrals        = NULL;
  *impropers        = NULL;
  *numcterms        = 0;
  *ctermcols        = 0;
  *ctermrows        = 0;
  *cterms            = NULL;
  *numangletypes    = 0;
  *numdihedraltypes = 0;
  *numimpropertypes = 0;
  
#if vmdplugin_ABIVERSION > 14
  *angletypenames    = NULL;
  *dihedraltypenames = NULL;
  *impropertypenames = NULL;
  *angletypes        = NULL;
  *dihedraltypes     = NULL;
  *impropertypes     = NULL;
#else
  *angleforces    = NULL;
  *dihedralforces = NULL;
  *improperforces = NULL;
  *ctermforces    = NULL;
#endif
  
  return VMDPLUGIN_SUCCESS;
}

static void close_lammpsdata_read(void *mydata) {
//   fprintf(stderr,"lammpsdatadata) close\n");
  lammpsdatadata *data = (lammpsdatadata *)mydata;
  fclose(data->file);
  delete data;
}

/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "lammpsdata";
  plugin.prettyname = "LAMMPS data file";
  plugin.author = "Hanno Dietrich";
  plugin.majorv = 0;
  plugin.minorv = 1;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "data,lmpdat,lmpdata,lammpsdata";
  plugin.open_file_read = open_lammpsdata_read;
  plugin.read_structure = read_lammpsdata_structure;
  plugin.read_next_timestep = read_lammpsdata_timestep;
#if vmdplugin_ABIVERSION > 10
  plugin.read_timestep_metadata = read_lammpsdata_timestep_metadata;
#endif
  plugin.read_bonds = read_lammpsdata_bonds;
  //plugin.read_angles = read_lammpsdata_angles;
  plugin.close_file_read = close_lammpsdata_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

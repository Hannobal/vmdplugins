/***************************************************************************
 *cr
 *cr                 (C) Copyright 2015-2020 Hanno Dietrich
 *cr                   University of Erlangen-Nuernberg
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: dlfieldplugin.c,v $
 *      $Author: hannod $       $Locker:  $             $State: Exp $
 *      $Revision: 0.10 $       $Date: 2015/11/10 14:12:44 $
 *
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "molfile_plugin.h"
#include "periodic_table.h"

#define LENMOLNAME 128
#define LENATNAME  8

class dlflddata {
public:
  FILE *file;
  int optflags;
  // net values for entire field file
  int nmols, nvdw, ntbp, nfbp, nmetal, ntersoff, nexternal;
  int natoms, nshells, nconstraints, nbonds, nangles, ndiheds;
  int nrigid, ntether, ninv;
  // bonds, angles, dihedrals, etc
  int *bdat1, *bdat2;
  int *angleatoms, *dihedralatoms, *inversionatoms;
  // values for molecules
  int *molnents, *molnatoms, *molnshells, *molnconstraints, *molnbonds;
  int *molnangles, *molndiheds, *molnrigid, *molntether, *molninv;
  // values for atoms
  int *molid, *resid;
  float *atmass, *atcharge;
  char (*atname)[LENATNAME];
  char (*molname)[LENMOLNAME];
  
  dlflddata(void) {
    file = NULL;
    nmols = -1;
    nvdw      = -1;
    ntbp      = -1;
    nfbp      = -1;
    nmetal    = -1;
    ntersoff  = -1;
    nexternal = -1;
    natoms       = 0;
    nshells      = 0;
    nconstraints = 0;
    nbonds       = 0;
    nangles      = 0;
    ndiheds      = 0;
    nrigid       = 0;
    ntether      = 0;
    ninv         = 0;
    bdat1          = NULL;
    bdat2          = NULL;
    angleatoms     = NULL;
    dihedralatoms  = NULL;
    inversionatoms = NULL;
    atmass   = NULL;
    atcharge = NULL;
    atname   = NULL;
    molname  = NULL;
    molid    = NULL;
    resid    = NULL;
    molnents        = NULL;
    molnatoms       = NULL;
    molnshells      = NULL;
    molnconstraints = NULL;
    molnbonds       = NULL;
    molnangles      = NULL;
    molndiheds      = NULL;
    molninv         = NULL;
    molnrigid       = NULL;
    molntether      = NULL;
  }
  
  ~dlflddata(void) {
//     printf("**** called destructor\n");
    if(NULL != bdat1) delete [] bdat1;
    if(NULL != bdat2) delete [] bdat2;
    if(NULL != angleatoms) delete [] angleatoms;
    if(NULL != dihedralatoms) delete [] dihedralatoms;
    if(NULL != inversionatoms) delete [] inversionatoms;
    if(NULL != atmass) delete [] atmass;
    if(NULL != atcharge) delete [] atcharge;
    if(NULL != atname) delete [] atname;
    if(NULL != molname) delete [] molname;
    if(NULL != molid) delete [] molid;
    if(NULL != resid) delete [] resid;
    if(NULL != molnents) delete [] molnents;
    if(NULL != molnatoms) delete [] molnatoms;
    if(NULL != molnshells) delete [] molnshells;
    if(NULL != molnconstraints) delete [] molnconstraints;
    if(NULL != molnbonds) delete [] molnbonds;
    if(NULL != molnangles) delete [] molnangles;
    if(NULL != molndiheds) delete [] molndiheds;
    if(NULL != molninv) delete [] molninv;
    if(NULL != molnrigid) delete [] molnrigid;
    if(NULL != molntether) delete [] molntether;
  }
};

int get_pte_idx_from_mass(float mass) {
  for(int i=0; i<nr_pte_entries; i++) {
    if(fabsf(mass-pte_mass[i])<0.5)
      return i;
  }
  return 0;
}

static void *open_dlfld_read(const char *filename, const char *, 
    int *natoms) {
  FILE *fd;
  dlflddata *data;
  char fbuffer[4096], buf1[16], buf2[16], buf3[16], buf4[16], buf5[16];
  char molname[LENMOLNAME];
  int scancount, t, i, j, numentries, nents, i1, i2, i3, i4, section, molnatoms;
  float r1, r2, r3, r4, r5, r6, r7;
  
  fd = fopen(filename, "rb");
  if (!fd) return NULL;
  
  data = new dlflddata;
  data->file = fd;
  
  while(NULL != fgets(fbuffer, 1024, fd)) {
    for(i = 0; fbuffer[i]; i++){
      fbuffer[i] = toupper(fbuffer[i]);
    }
//     printf("**** read: %s",fbuffer);
    if(strstr(fbuffer,"MOLECULES")) {
//       printf("**** found MOLECULES entry\n");
      scancount = sscanf(fbuffer, "%s %d", buf1, &i1);
      if(scancount != 2) {
	printf("open_dlfld_read) error reading MOLECULES entry in FIELD file\n");
	delete data;
	return NULL;
      }
      data->nmols=i1;
//       printf("**** found %d molecules\n", data->nmols);
      break;
    }
  }
  
  if(data->nmols<0) {
    printf("open_dlfld_read) no MOLECULES entry found in FIELD file\n");
    delete data;
    return NULL;
  }
  
  for(t=0;t<data->nmols;t++) {
    // molecule name
    if(NULL == fgets(molname, LENMOLNAME, fd)) {
      printf("open_dlfld_read) unexpected end of FIELD file 1\n");
      delete data;
      return NULL;
    }
    // number of molecule entities
    if(NULL == fgets(fbuffer, 1024, fd)) {
      printf("open_dlfld_read) unexpected end of FIELD file 2\n");
      delete data;
      return NULL;
    }
    for(i = 0; fbuffer[i]; i++){
      fbuffer[i] = toupper(fbuffer[i]);
    }
    if(! strstr(fbuffer,"NUMMOLS")) {
      printf("open_dlfld_read) expected NUMMOLS but got %s\n",fbuffer);
      delete data;
      return NULL;
    }
    scancount = sscanf(fbuffer, "%s %d", buf1, &nents);
    if(scancount != 2) {
      printf("open_dlfld_read) error reading NUMMOLS entry for molecule %s in FIELD file\n",molname);
      delete data;
      return NULL;
    }
    molnatoms = 0;
    while(NULL != fgets(fbuffer, 1024, fd)) {
      for(i = 0; fbuffer[i]; i++){
	fbuffer[i] = toupper(fbuffer[i]);
      }
//       printf("**** read: %s",fbuffer);
      if(strstr(fbuffer,"FINISH")) {
	break;
      } else if(strstr(fbuffer,"ATOMS")) {
	section=1;
      } else if(strstr(fbuffer,"SHELL")) {
	section=2;
      } else if(strstr(fbuffer,"BONDS")) {
	section=3;
      } else if(strstr(fbuffer,"CONSTRAINTS")) {
	section=4;
      } else if(strstr(fbuffer,"PMF")) {
	section=5;
      } else if(strstr(fbuffer,"ANGLES")) {
	section=6;
      } else if(strstr(fbuffer,"DIHEDRALS")) {
	section=7;
      } else if(strstr(fbuffer,"INVERSIONS")) {
	section=8;
      } else if(strstr(fbuffer,"RIGID")) {
	section=9;
      } else if(strstr(fbuffer,"TETH")) {
	section=10;
      } else {
	printf("open_dlfld_read) unknown keyword in FIELD file:\nopen_dlfld_read) %s\n",fbuffer);
	delete data;
	return NULL;
      }
      if(section!=5) {
	scancount = sscanf(fbuffer, "%s %d", buf1, &numentries);
	if(scancount != 2) {
	  printf("open_dlfld_read) error reading number of entries  in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	  delete data;
	  return NULL;
	}
	if(section==1) {
	  molnatoms = numentries;
	}
      } else {
	numentries=2;
      }
      for(i=0;i<numentries;i++) {
	if(NULL == fgets(fbuffer, 1024, fd)) {
	  printf("open_dlfld_read) unexpected end of FIELD file 3\n");
	  delete data;
	  return NULL;
	}
// 	printf("**** read: %s",fbuffer);
	if(section==1) { // ATOMS
	  scancount = sscanf(fbuffer, "%s %f %f %d %d %d", buf1, &r1, &r2, &i1, &i2, &i3);
	  if(scancount != 6) {
	    printf("open_dlfld_read) error reading ATOM entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>1) numentries -= i1-1; // check nreps
	  data->natoms += i1*nents;
	} else if(section==2) { // SHELL
	  scancount = sscanf(fbuffer, "%d %d %f %f", &i1, &i2, &r1, &r2);
	  if(scancount != 4) {
	    printf("open_dlfld_read) error reading SHELL entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms) {
	    printf("open_dlfld_read) error reading SHELL entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->nshells += nents;
	} else if(section==3) { // BONDS
	  scancount = sscanf(fbuffer, "%s %d %d %f %f %f %f", buf1, &i1, &i2, &r1, &r2, &r3, &r4);
	  if(scancount < 4) {
	    printf("open_dlfld_read) error reading BOND entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms or i2>molnatoms) {
	    printf("open_dlfld_read) error reading BOND entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->nbonds += nents;
	} else if(section==4) { // CONSTRAINTS
	  scancount = sscanf(fbuffer, "%d %d %f", &i1, &i2, &r1);
	  if(scancount != 3) {
	    printf("open_dlfld_read) error reading CONSTRAINT entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms or i2>molnatoms) {
	    printf("open_dlfld_read) error reading CONSTRAINT entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->nbonds += nents;
	} else if(section==5) { // PMF
	    printf("open_dlfld_read) keyword PMF is not implemented yet!\n");
	    delete data;
	    return NULL;
	} else if(section==6) { // ANGLES
	  scancount = sscanf(fbuffer, "%s %d %d %d %f %f %f %f", buf1, &i1, &i2, &i3, &r1, &r2, &r3, &r4);
	  if(scancount < 6) {
	    printf("open_dlfld_read) error reading ANGLE entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms or i2>molnatoms or i3>molnatoms) {
	    printf("open_dlfld_read) error reading ANGLE entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->nangles += nents;
	} else if(section==7) { // DIHEDRALS
	  scancount = sscanf(fbuffer, "%s %d %d %d %d %f %f %f %f %f", buf1, &i1, &i2, &i3, &i4, &r1, &r2, &r3, &r4, &r5);
	  if(scancount < 6) {
	    printf("open_dlfld_read) error reading DIHEDRAL entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms or i2>molnatoms or i3>molnatoms or i4>molnatoms) {
	    printf("open_dlfld_read) error reading DIHEDRAL entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->ndiheds += nents;
	} else if(section==8) { // INVERSIONS
	  scancount = sscanf(fbuffer, "%s %d %d %d %d %f %f", buf1, &i1, &i2, &i3, &i4, &r1, &r2);
	  if(scancount < 6) {
	    printf("open_dlfld_read) error reading INVERSION entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms or i2>molnatoms or i3>molnatoms or i4>molnatoms) {
	    printf("open_dlfld_read) error reading INVERSION entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->ninv += nents;
	} else if(section==9) { // RIGID
	    printf("open_dlfld_read) keyword RIGID is not implemented yet!\n");
	    delete data;
	    return NULL;
	} else if(section==10) { // TETH
	  scancount = sscanf(fbuffer, "%s %d %f %f %f %f", buf1, &i1, &r1, &r2, &r3, &r4);
	  if(scancount < 3) {
	    printf("open_dlfld_read) error reading TETH entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	    delete data;
	    return NULL;
	  }
	  if(i1>molnatoms) {
	    printf("open_dlfld_read) error reading TETH entry: molecule %s has only %d atoms:\nopen_dlfld_read) %s\n",
		   molname, molnatoms, fbuffer);
	    delete data;
	    return NULL;
	  }
	  data->ntether += nents;
	}
      }
    }
  }
//   printf("**** reading nonbonded interactions\n");
  // read non-bonded interactions
  while(NULL != fgets(fbuffer, 1024, fd)) {
    for(i = 0; fbuffer[i]; i++){
      fbuffer[i] = toupper(fbuffer[i]);
    }
//     printf("**** read: %s",fbuffer);
    if(strstr(fbuffer,"CLOSE")) {
      break;
    } else if(strstr(fbuffer,"VDW")) {
      section=11;
    } else if(strstr(fbuffer,"TBP")) {
      section=12;
    } else if(strstr(fbuffer,"FBP")) {
      section=13;
    } else if(strstr(fbuffer,"METAL")) {
      section=14;
    } else if(strstr(fbuffer,"TERSOFF")) {
      section=15;
    } else if(strstr(fbuffer,"EXTERN")) {
      section=16;
    } else {
      printf("open_dlfld_read) unknown keyword in FIELD file:\nopen_dlfld_read) %s\n",fbuffer);
      delete data;
      return NULL;
    }
    if(section!=16) {
      scancount = sscanf(fbuffer, "%s %d", buf1, &numentries);
      if(scancount != 2) {
	printf("open_dlfld_read) error reading number of entries  in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	delete data;
	return NULL;
      }
    } else {
      numentries=1;
    }
    for(i=0;i<numentries;i++) {
      if(NULL == fgets(fbuffer, 1024, fd)) {
	printf("open_dlfld_read) unexpected end of FIELD file 4\n");
	delete data;
	return NULL;
      }
//       printf("**** read: %s",fbuffer);
      if(section==11) { // VDW
	scancount = sscanf(fbuffer, "%s %s %s %f %f %f %f %f", buf1, buf2, buf3, &r1, &r2, &r3, &r4, &r5);
	if(scancount < 4) {
	  printf("open_dlfld_read) error reading VDW entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	  delete data;
	  return NULL;
	}
	data->nvdw++;
      } else if(section==12) { // TPB
	scancount = sscanf(fbuffer, "%s %s %s %s %f %f %f %f %f", buf1, buf2, buf3, buf4, &r1, &r2, &r3, &r4, &r5);
	if(scancount < 5) {
	  printf("open_dlfld_read) error reading TBP entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	  delete data;
	  return NULL;
	}
	data->ntbp++;
      } else if(section==13) { // FBP
	scancount = sscanf(fbuffer, "%s %s %s %s %s %f %f %f", buf1, buf2, buf3, buf4, buf5, &r1, &r2, &r3);
	if(scancount < 6) {
	  printf("open_dlfld_read) error reading FBP entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	  delete data;
	  return NULL;
	}
	data->ntbp++;
      } else if(section==14) { // METAL
	scancount = sscanf(fbuffer, "%s %s %s %f %f %f %f %f %f %f", buf1, buf2, buf3, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
	if(scancount < 4) {
	  printf("open_dlfld_read) error reading METAL entry in FIELD file in line:\nopen_dlfld_read) %s\n",fbuffer);
	  delete data;
	  return NULL;
	}
	data->nmetal++;
      } else if(section==15) { // TERSOFF
	printf("open_dlfld_read) keyword TERSOFF is not implemented yet!\n");
	delete data;
	return NULL;
      } else if(section==16) { // EXTERN
	scancount = sscanf(fbuffer, "%s %d", buf1, &i1);
	i2 = (int) (i1/5.0+0.9);
	for(j=0;j<i2;j++) {
	  if(NULL == fgets(fbuffer, 1024, fd)) {
	    printf("open_dlfld_read) unexpected end of FIELD file\n");
	    delete data;
	    return NULL;
	  }
	}
      }
    }
  }
  
//   printf("**** read: natoms=%d nbonds=%d nangles=%d ndiheds=%d\n",data->natoms,data->nbonds,data->nangles,data->ndiheds);
//   printf("**** done with first iteration of reading\n");
  
  *natoms = data->natoms;
  data->bdat1          = new int[data->nbonds];
  data->bdat2          = new int[data->nbonds];
  data->angleatoms     = new int[data->nangles*3];
  data->dihedralatoms  = new int[data->ndiheds*4];
  data->inversionatoms = new int[data->ninv*4];
  data->atmass   = new float[data->natoms];
  data->atcharge = new float[data->natoms];
  data->molid    = new int[data->natoms];
  data->resid    = new int[data->natoms];
  data->atname   = new char[data->natoms][LENATNAME];
  data->molname         = new char[data->nmols][LENMOLNAME];
  data->molnents        = new int[data->nmols];
  data->molnatoms       = new int[data->nmols];
  data->molnshells      = new int[data->nmols];
  data->molnconstraints = new int[data->nmols];
  data->molnbonds       = new int[data->nmols];
  data->molnangles      = new int[data->nmols];
  data->molndiheds      = new int[data->nmols];
  data->molnrigid       = new int[data->nmols];
  data->molntether      = new int[data->nmols];
  data->molninv         = new int[data->nmols];
  
  // now start over and actually read the data
  // no syntax checks beyond this point as we do the same thing again
  rewind(data->file);
  
  while(NULL != fgets(fbuffer, 1024, fd)) {
//     printf("**** read: %s",fbuffer);
    for(i = 0; fbuffer[i]; i++){
      fbuffer[i] = toupper(fbuffer[i]);
    }
    if(strstr(fbuffer,"MOLECULES")) {
      break;
    }
  }
  
  int iat, ibd=0, iang=0, idih=0, io=0, iinv=0, mo=0, m;
  
  for(t=0;t<data->nmols;t++) {
    // molecule name
    if(NULL==fgets(molname, LENMOLNAME, fd));
//     printf("**** read: %s",molname);
    strncpy(data->molname[t],molname,LENMOLNAME);
    // number of molecule entities
    if(NULL==fgets(fbuffer, 1024, fd));
//     printf("**** read: %s",fbuffer);
    scancount = sscanf(fbuffer, "%s %d", buf1, &nents);
    data->molnents[t] = nents;
    while(NULL != fgets(fbuffer, 1024, fd)) {
      for(i = 0; fbuffer[i]; i++){
	fbuffer[i] = toupper(fbuffer[i]);
      }
//       printf("**** read: %s",fbuffer);
      scancount = sscanf(fbuffer, "%s %d", buf1, &numentries);
      if(strstr(fbuffer,"FINISH")) {
	break;
      } else if(strstr(fbuffer,"ATOMS")) {
	section=1;
	data->molnatoms[t]=numentries;
      } else if(strstr(fbuffer,"SHELL")) {
	section=2;
	data->molnshells[t]=numentries;
      } else if(strstr(fbuffer,"BONDS")) {
	section=3;
	data->molnbonds[t]=numentries;
      } else if(strstr(fbuffer,"CONSTRAINTS")) {
	section=4;
	data->molnconstraints[t]=numentries;
      } else if(strstr(fbuffer,"PMF")) {
	section=5;
	numentries=2;
      } else if(strstr(fbuffer,"ANGLES")) {
	section=6;
	data->molnangles[t]=numentries;
      } else if(strstr(fbuffer,"DIHEDRALS")) {
	section=7;
	data->molndiheds[t]=numentries;
      } else if(strstr(fbuffer,"INVERSIONS")) {
	section=8;
	data->molninv[t]=numentries;
      } else if(strstr(fbuffer,"RIGID")) {
	section=9;
	data->molnrigid[t]=numentries;
      } else if(strstr(fbuffer,"TETH")) {
	section=10;
	data->molntether[t]=numentries;
      } else {
	printf("open_dlfld_read) unknown keyword in FIELD file:\nopen_dlfld_read) %s\n",fbuffer);
	delete data;
	return NULL;
      }
      for(i=0;i<numentries;i++) {
	if(NULL==fgets(fbuffer, 1024, fd));
	if(section==1) { // ATOMS
// 	  printf("**** read: %s",fbuffer);
	  scancount = sscanf(fbuffer, "%s %f %f %d %d %d", buf1, &r1, &r2, &i1, &i2, &i3);
	  for(m=0;m<nents;m++) {
	    for(i4=0;i4<i1;i4++) {
	      iat = io+m*data->molnatoms[t]+i+i4;
// 	      printf("**** before: %d\n",iat);
	      strncpy(data->atname[iat],buf1,LENATNAME);
// 	      printf("**** after: %s\n",data->atname[iat]);
	      data->atmass[iat]   = r1;
	      data->atcharge[iat] = r2;
	      data->resid[iat]    = t;
	      data->molid[iat]    = mo+m;
	    }
	  }
	  i += i1-1;
	} else if(section==2) { // SHELL
	  scancount = sscanf(fbuffer, "%d %d %f %f", &i1, &i2, &r1, &r2);
	} else if(section==3) { // BONDS
	  scancount = sscanf(fbuffer, "%s %d %d %f %f %f %f", buf1, &i1, &i2, &r1, &r2, &r3, &r4);
	  for(m=0;m<nents;m++) {
	    iat = io+m*data->molnatoms[t];
	    data->bdat1[ibd] = iat+i1;
	    data->bdat2[ibd] = iat+i2;
	    ibd++;
	  }
	} else if(section==4) { // CONSTRAINTS
	  scancount = sscanf(fbuffer, "%d %d %f", &i1, &i2, &r1);
	  for(m=0;m<nents;m++) {
	    iat = io+m*data->molnatoms[t];
	    data->bdat1[ibd] = iat+i1;
	    data->bdat2[ibd] = iat+i2;
	    ibd++;
	  }
	} else if(section==5) { // PMF
	    printf("open_dlfld_read) keyword PMF is not implemented yet!\n");
	    delete data;
	    return NULL;
	} else if(section==6) { // ANGLES
	  scancount = sscanf(fbuffer, "%s %d %d %d %f %f %f %f", buf1, &i1, &i2, &i3, &r1, &r2, &r3, &r4);
	  for(m=0;m<nents;m++) {
	    iat = io+m*data->molnatoms[t];
	    data->angleatoms[iang]   = iat+i1;
	    data->angleatoms[iang+1] = iat+i2;
	    data->angleatoms[iang+2] = iat+i3;
	    iang += 3;
	  }
	} else if(section==7) { // DIHEDRALS
	  scancount = sscanf(fbuffer, "%s %d %d %d %d %f %f %f %f %f", buf1, &i1, &i2, &i3, &i4, &r1, &r2, &r3, &r4, &r5);
	  for(m=0;m<nents;m++) {
	    iat = io+m*data->molnatoms[t];
	    data->dihedralatoms[idih]   = iat+i1;
	    data->dihedralatoms[idih+1] = iat+i2;
	    data->dihedralatoms[idih+2] = iat+i3;
	    data->dihedralatoms[idih+3] = iat+i4;
	    idih += 4;
	  }
	} else if(section==8) { // INVERSIONS
	  scancount = sscanf(fbuffer, "%s %d %d %d %d %f %f", buf1, &i1, &i2, &i3, &i4, &r1, &r2);
	  for(m=0;m<nents;m++) {
	    iat = io+m*data->molnatoms[t];
	    data->inversionatoms[iinv]   = iat+i1;
	    data->inversionatoms[iinv+1] = iat+i2;
	    data->inversionatoms[iinv+2] = iat+i3;
	    data->inversionatoms[iinv+3] = iat+i4;
	    iinv += 4;
	  }
	} else if(section==9) { // RIGID
	    printf("open_dlfld_read) keyword RIGID is not implemented yet!\n");
	    delete data;
	    return NULL;
	} else if(section==10) { // TETH
	  scancount = sscanf(fbuffer, "%s %d %f %f %f %f", buf1, &i1, &r1, &r2, &r3, &r4);
	}
      }
    }
    io += data->molnents[t] * data->molnatoms[t];
    mo += data->molnents[t];
  }
//   printf("**** reading nonbonded interactions\n");
  // read non-bonded interactions
  while(NULL != fgets(fbuffer, 1024, fd)) {
    for(i = 0; fbuffer[i]; i++){
      fbuffer[i] = toupper(fbuffer[i]);
    }
    if(strstr(fbuffer,"CLOSE")) {
      break;
    } else if(strstr(fbuffer,"VDW")) {
      section=11;
    } else if(strstr(fbuffer,"TBP")) {
      section=12;
    } else if(strstr(fbuffer,"FBP")) {
      section=13;
    } else if(strstr(fbuffer,"METAL")) {
      section=14;
    } else if(strstr(fbuffer,"TERSOFF")) {
      section=15;
    } else if(strstr(fbuffer,"EXTERN")) {
      section=16;
    } else {
      printf("open_dlfld_read) unknown keyword in FIELD file:\nopen_dlfld_read) %s\n",fbuffer);
      delete data;
      return NULL;
    }
    if(section!=16) {
      scancount = sscanf(fbuffer, "%s %d", buf1, &numentries);
    } else {
      numentries=1;
    }
    for(i=0;i<numentries;i++) {
      if(NULL==fgets(fbuffer, 1024, fd));
      if(section==11) { // VDW
	scancount = sscanf(fbuffer, "%s %s %s %f %f %f %f %f", buf1, buf2, buf3, &r1, &r2, &r3, &r4, &r5);
      } else if(section==12) { // TPB
	scancount = sscanf(fbuffer, "%s %s %s %s %f %f %f %f %f", buf1, buf2, buf3, buf4, &r1, &r2, &r3, &r4, &r5);
      } else if(section==13) { // FBP
	scancount = sscanf(fbuffer, "%s %s %s %s %s %f %f %f", buf1, buf2, buf3, buf4, buf5, &r1, &r2, &r3);
      } else if(section==14) { // METAL
	scancount = sscanf(fbuffer, "%s %s %s %f %f %f %f %f %f %f", buf1, buf2, buf3, &r1, &r2, &r3, &r4, &r5, &r6, &r7);
      } else if(section==15) { // TERSOFF
	printf("open_dlfld_read) keyword TERSOFF is not implemented yet!\n");
	delete data;
	return NULL;
      } else if(section==16) { // EXTERN
	scancount = sscanf(fbuffer, "%s %d", buf1, &i1);
	i2 = (int) (i1/5.0+0.9);
	for(j=0;j<i2;j++) {
	  if(NULL==fgets(fbuffer, 1024, fd));
	}
      }
    }
  }
//   printf("**** done with open_dlfld_read\n");
//   printf("**** last atom name (id %d): %s\n",data->natoms,data->atname[data->natoms]);
  return data;
}

static int read_dlfld_structure(void *mydata, int *optflags,
    molfile_atom_t *atoms) {
//   printf("**** read_dlfld_structure\n");
  dlflddata *data = (dlflddata *)mydata;
  int i;
  
  *optflags = MOLFILE_ATOMICNUMBER | MOLFILE_CHARGE | MOLFILE_MASS;
  for (i=0; i<data->natoms; i++) {
//     printf("**** atom: %d\n",i);
    molfile_atom_t *atom = atoms+i;
    // segmentation faults occur when atom->name/type/segid/segname/resname are overfilled!
    strncpy(atom->name,data->atname[i],LENATNAME);
    strncpy(atom->type,data->atname[i],LENATNAME);
    strncpy(atom->resname,data->molname[data->resid[i]],8);
    sprintf(atom->segid,"%d",data->molid[i]);
    atom->mass   = data->atmass[i];
    atom->charge = data->atcharge[i];
    atom->resid  = data->resid[i];
    atom->chain[0] = '\0';
    int idx = get_pte_idx_from_mass(atom->mass);
    atom->radius = get_pte_vdw_radius(idx);
    atom->atomicnumber = idx;
  }
//   printf("**** done with read_dlfld_structure\n");
  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int read_dlfld_bonds(void *mydata, int *nbonds, int **from, int **to, 
                           float **bondorder,  int **bondtype, 
                           int *nbondtypes, char ***bondtypename) {
#else
static int read_dlfld_bonds(void *mydata, int *nbonds, int **from, int **to, float **bondorder) {
#endif

//   printf("**** read_dlfld_bonds\n");
  dlflddata *data = (dlflddata *)mydata;

  *nbonds  = data->nbonds;
  *from = data->bdat1;
  *to   = data->bdat2;
  *bondorder = NULL;
#if vmdplugin_ABIVERSION > 14
  *bondtype     = NULL;
  *nbondtypes   = 0;
  *bondtypename = NULL;
#endif
  
//   printf("**** done with read_dlfld_bonds\n");
  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int read_dlfld_angles(void *mydata, int *numangles, int **angles, int **angletypes,
                      int *numangletypes, char ***angletypenames, int *numdihedrals,
                      int **dihedrals, int **dihedraltypes, int *numdihedraltypes,
                      char ***dihedraltypenames, int *numimpropers, int **impropers,        
                      int **impropertypes, int *numimpropertypes, char ***impropertypenames,
                      int *numcterms, int **cterms, int *ctermcols, int *ctermrows) {
#else
static int read_dlfld_angles(void *mydata,
                int *numangles,    int **angles,    double **angleforces,
                int *numdihedrals, int **dihedrals, double **dihedralforces,
                int *numimpropers, int **impropers, double **improperforces,
                int *numcterms,    int **cterms, 
                int *ctermcols,    int *ctermrows,  double **ctermforces);
#endif
  
//   printf("**** read_dlfld_angles\n");
  dlflddata *data = (dlflddata *)mydata;
  
  *numangles        = data->nangles;
  *numdihedrals     = data->ndiheds;
  *numimpropers     = data->ninv;
  *angles           = data->angleatoms;
  *dihedrals        = data->dihedralatoms;
  *impropers        = data->inversionatoms;
  *numcterms        = 0;
  *ctermcols        = 0;
  *ctermrows        = 0;
  *cterms            = NULL;
  
#if vmdplugin_ABIVERSION > 14
  *numangletypes    = 0;
  *numdihedraltypes = 0;
  *numimpropertypes = 0;
  *angletypenames    = NULL;
  *dihedraltypenames = NULL;
  *impropertypenames = NULL;
  *angletypes        = NULL;
  *dihedraltypes     = NULL;
  *impropertypes     = NULL;
#else
  *numangletypes    = 0;
  *numdihedraltypes = 0;
  *numimpropertypes = 0;
  *angleforces    = NULL;
  *dihedralforces = NULL;
  *improperforces = NULL;
  *ctermforces    = NULL;
#endif
  
//   printf("**** done with read_dlfld_angles\n");
  return VMDPLUGIN_SUCCESS;
}

static void close_dlfld_read(void *mydata) {
//   printf("**** close_dlfld_read\n");
  dlflddata *data = (dlflddata *)mydata;
  delete data;
//   printf("**** done with close_dlfld_read\n");
}


/* registration stuff */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "dlpolyfld";
  plugin.prettyname = "DLPOLY FIELD file";
  plugin.author = "Hanno Dietrich";
  plugin.majorv = 0;
  plugin.minorv = 5;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "dlfld,field";
  plugin.open_file_read = open_dlfld_read;
  plugin.read_structure = read_dlfld_structure;
  plugin.read_bonds = read_dlfld_bonds;
  plugin.read_angles = read_dlfld_angles;
  plugin.close_file_read = close_dlfld_read;
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

#ifndef WRITER_H
#define WRITER_H
#ifdef __cplusplus
//extern "C" {
#endif
#include <petscmat.h> 
#include "fem_inf.h"
#include "oce1.h"

struct _p_WFWriter {
  MPI_Comm comm;
  PetscBool active;
  PetscInt num;
  PetscReal xmax;
  FILE *fp;
};
typedef struct _p_WFWriter* WFWriter;

PetscErrorCode WFWriterCreate(WFWriter *writer, MPI_Comm comm);
PetscErrorCode WFWriterSet(WFWriter self, PetscInt num, PetscReal xmax);
PetscErrorCode WFWriterSetPath(WFWriter self, char path[]);
PetscErrorCode WFWriterSetFromOptions(WFWriter self);
//PetscErrorCode WFWriterCreateFromOptions(WFWriter *p_self, MPI_Comm comm);
PetscErrorCode WFWriterDestroy(WFWriter *p_self);
PetscErrorCode WFWriterView(WFWriter self);
PetscBool WFWriterIsActive(WFWriter self);

PetscErrorCode WFWriterWrite(WFWriter self, FILE *fp, FEMInf fem, Vec c);
PetscErrorCode WFWriterWriteFile(WFWriter self, char *fn, FEMInf fem, Vec c);

PetscErrorCode WFWriterWriteOCE1(WFWriter self, FILE *fp, OCE1 oce, Vec c);
PetscErrorCode WFWriterWriteFileOCE1(WFWriter self, char *fn, OCE1 oce, Vec c);

PetscErrorCode WFWriterWritePOT(WFWriter self, FILE *fp, POT pot, int J, PetscReal mu);
PetscErrorCode WFWriterWriteFilePOT(WFWriter self, char *fn, POT pot, int J, PetscReal mu);

#ifdef __cplusplus
//}
#endif
#endif

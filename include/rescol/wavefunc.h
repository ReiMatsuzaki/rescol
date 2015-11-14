#ifndef WAVEFUNC_H
#define WAVEFUNC_H
#ifdef __cplusplus
extern "C" {
#endif 
#include <petscpf.h>
#include <rescol/fem_inf.h>
typedef PF WaveFunc;
PetscErrorCode WaveFuncCreate(MPI_Comm comm, WaveFunc *p_self);
PetscErrorCode WaveFuncSetHEig(WaveFunc self, int n, int L, PetscScalar z);
PetscErrorCode WaveFuncSetFromFEM(WaveFunc self, FEMInf fem, Vec c);
PetscErrorCode WaveFuncSetFromOptions(WaveFunc self);
#ifdef __cplusplus
}
#endif 
#endif

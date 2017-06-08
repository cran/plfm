
#include <R.h>
#include <Rinternals.h>

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


extern "C" void PlFm_XZ_Y_DC(int * , int * , int * , int * , int * , int * , int * , 
						  double * ,
						  double * , double * , double *,
						  double * , double * ,
						  double * , double * , double * ,
						  double * , double * , double * , double * , double * , double * ,double *, int *, 
						  int *, double *, double * , double * , double * );
extern "C" void PlFm_X_YZ_DC(int * , int * , int * , int * , int * , int * , int * , 
						  double * ,
						  double * , double * , double *,
						  double * , double * ,
						  double * , double * , double * ,
						  double * , double * , double * , double * , double * , double * ,double *, int *,
						  int *, double *, double * , double * , double * );
extern "C" void PlFm_XZ_YZ_DC(int * , int * , int * , int * , int * , int * , int * , 
						  double * ,
						  double * , double * , double *,
						  double * , double * ,
						  double * , double * , double * ,
						  double * , double * , double * , double * , double * , double * ,double *, int *,
						  int *, double *, double * , double * , double * );
extern "C" void PlFm_XZ_Y_ADD(int * , int * , int * , int * , int * , int * , int * , 
						  double * ,
						  double * , double * , double *,
						  double * , double * ,
						  double * , double * , double * ,
						  double * , double * , double * , double * , double * , double * ,double *, int *, 
						  int *, double *, double * , double * , double * );
extern "C" void PlFm_X_YZ_ADD(int * , int * , int * , int * , int * , int * , int * , 
						  double * ,
						  double * , double * , double *,
						  double * , double * ,
						  double * , double * , double * ,
						  double * , double * , double * , double * , double * , double * ,double *, int *,
						  int *, double *, double * , double * , double * );
extern "C" void PlFm_XZ_YZ_ADD(int * , int * , int * , int * , int * , int * , int * , 
						  double * ,
						  double * , double * , double *,
						  double * , double * ,
						  double * , double * , double * ,
						  double * , double * , double * , double * , double * , double * ,double *, int *,
						  int *, double *, double * , double * , double * );


static const R_CMethodDef CEntries[] = {
    {"PlFm_X_YZ_ADD",  (DL_FUNC) &PlFm_X_YZ_ADD,  29},
    {"PlFm_X_YZ_DC",   (DL_FUNC) &PlFm_X_YZ_DC,   29},
    {"PlFm_XZ_Y_ADD",  (DL_FUNC) &PlFm_XZ_Y_ADD,  29},
    {"PlFm_XZ_Y_DC",   (DL_FUNC) &PlFm_XZ_Y_DC,   29},
    {"PlFm_XZ_YZ_ADD", (DL_FUNC) &PlFm_XZ_YZ_ADD, 29},
    {"PlFm_XZ_YZ_DC",  (DL_FUNC) &PlFm_XZ_YZ_DC,  29},
    {NULL, NULL, 0}
};

extern "C" void R_init_plfm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}





/*
 * MATLAB Compiler: 3.0
 * Date: Thu Aug 28 15:51:05 2003
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-B" "sgl" "-m" "-W"
 * "main" "-L" "C" "-t" "-T" "link:exe" "-h" "libmmfile.mlib" "-W" "mainhg"
 * "libmwsglm.mlib" "process_rowdata" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __detectstepsplot_h
#define __detectstepsplot_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_detectstepsplot(void);
extern void TerminateModule_detectstepsplot(void);
extern _mexLocalFunctionTable _local_function_table_detectstepsplot;

extern mxArray * mlfDetectstepsplot(mxArray * * figh,
                                    mxArray * stpfr,
                                    mxArray * mdf);
extern void mlxDetectstepsplot(int nlhs,
                               mxArray * plhs[],
                               int nrhs,
                               mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif

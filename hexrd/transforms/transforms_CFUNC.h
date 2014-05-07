/******************************************************************************/
/** The functions declared in this header use standard C to implement the    **/
/** various computations needed by the functions declared in the file        **/
/** "transforms_CAPI.h".                                                     **/
/**                                                                          **/
/** This separation of Python C API calls and C implementations allows the C **/
/** functions to call eachother without the unnecessary overhead of passing  **/
/** arguments and return values as Python objects.                           **/
/**                                                                          **/
/******************************************************************************/

#include <stdbool.h>

/******************************************************************************/
/* Funtions */

void gvecToDetectorXY_cfunc(long int npts, double * gVec_c_Ptr,
			    double * rMat_d_Ptr, double * rMat_s_Ptr, double * rMat_c_Ptr,
			    double * tVec_d_Ptr, double * tVec_s_Ptr, double * tVec_c_Ptr,
			    double * beamVec_Ptr,
			    double * result_Ptr);

void detectorXYToGvec_cfunc(long int npts, double * xy,
			    double * rMat_d, double * rMat_s,
			    double * tVec_d, double * tVec_s, double * tVec_c,
			    double * beamVec, double * etaVec,
			    double * tTh, double * eta, double * gVec_l);

void oscillAnglesOfHKLs_cfunc(long int npts, double * hkls, double chi,
			      double * rMat_c, double * bMat, double wavelength,
			      double * beamVec, double * etaVec,
			      double * oangs0, double * oangs1);

/******************************************************************************/
/* Utility Funtions */

void unitRowVector_cfunc(int n, double * cIn, double * cOut);

void unitRowVectors_cfunc(int m, int n, double * cIn, double * cOut);

void makeDetectorRotMat_cfunc(double * tPtr, double * rPtr);

void makeOscillRotMat_cfunc(double * oPtr, double * rPtr);

void makeRotMatOfExpMap_cfunc(double * ePtr, double * rPtr);

void makeRotMatOfQuat_cfunc(int nq, double * qPtr, double * rPtr);

void makeBinaryRotMat_cfunc(double * aPtr, double * rPtr);

void makeEtaFrameRotMat_cfunc(double * bPtr, double * ePtr, double * rPtr);

void validateAngleRanges_cfunc(int na, double * aPtr, int nr, double * minPtr, double * maxPtr, bool * rPtr, int ccw);

//void rotateVectorAboutAxis_cfunc(double tTh, double * n_g, double * dHat_l, double * gVec_l);

void rotate_vecs_about_axis_cfunc(long int na, double * angles,
				  long int nax, double * axes,
				  long int nv, double * vecs,
				  double * rVecs);

double quat_distance_cfunc(int nsym, double * q1, double * q2, double * qsym);

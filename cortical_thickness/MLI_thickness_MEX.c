/* MLI_thickness_MEX.c
 *
 * Computes the tissue thickness from a probability image, using minimum 
 * line integrals. See the help of MLI_thickness.m for more information.
 *
 * This file must be compiled using "mex MLI_thickness_MEX.c".
 * Alternatively, the compiled files can be downloaded from:
 * www.nitrc.org/projects/thickness
 *
 * Reference:
 * I. Aganj, G. Sapiro, N. Parikshak, S. K. Madsen, and P. Thompson,
 * "Measurement of cortical thickness from MRI by minimum line integrals on
 * soft-classified tissue," Human Brain Mapping, vol. 30, no. 10,
 * pp. 3188-3199, 2009. http://doi.org/10.1002/hbm.20740
 *
 * Codes by Iman Aganj.
 *
 */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int nI, nAng, nL, xM, yM, zM, x, y, z, m, n, t;
    int nR, nTr, c, nr;
    double tr, ntr, I0, voxI, J1, Ja1;
    int count = 0;
    double *I, *J, *A, *Ja, *L;
    const mwSize *sI;
    const double *lm;
    bool down;
    
    nI = mxGetNumberOfElements(prhs[0]);
    sI = mxGetDimensions(prhs[0]);
    I = (double *)mxGetPr(prhs[0]);
    nAng = mxGetNumberOfElements(prhs[1]);
    nR = (int) mxGetScalar(prhs[2]);
    tr = mxGetScalar(prhs[3]);
    nTr = (double) mxGetScalar(prhs[4]);
    if (nrhs>5)
        lm = (double *) mxGetPr(prhs[5]);
    
    plhs[0] = mxCreateNumericArray(3, sI, mxDOUBLE_CLASS, false);
    J = mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericArray(3, sI, mxDOUBLE_CLASS, false);
    A = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(3, sI, mxDOUBLE_CLASS, false);
    Ja = mxGetPr(plhs[2]);
    
    for (m = 0; m < nI; m++){
        *(A+m) = 0;
        *(Ja+m)= 0;
        *(J+m) = 1e10;
    }
    for (n = 0; n < nAng; n++){
        L = mxGetPr(mxGetCell(prhs[1],n));
        nL = mxGetM(mxGetCell(prhs[1],n));
        mexPrintf("Direction #%d out of %d\n", n+1, nAng);
        mexEvalString("drawnow");
        for (m = 0; m < nI; m++){
            if (nrhs<=5 || (m+1>=*lm && m+1<=*(lm+1))){
                xM = (m % (*sI * *(sI+1))) % *sI;
                yM = (m % (*sI * *(sI+1))) / *sI;
                zM = m / (*sI * *(sI+1));
                J1 = 0;
                for (t = 1; t >= -1; t-=2){
                    down = false;
                    I0 = 0;
                    c = 0;
                    nr = 0;
                    ntr = 0;
                    while ((c<nL) && (J1<*(J+m))){
                        x = xM + t * *(L + c);
                        y = yM + t * *(L + nL+c);
                        z = zM + t * *(L + 2*nL+c);
                        if (x>=0 && x<*sI && y>=0 && y<*(sI+1) && z>=0 && z<*(sI+2))
                            voxI = *(I + z * (*sI * *(sI+1)) + y * *sI + x);
                        else
                            voxI = 0;
                        J1 += voxI * *(L + 3*nL+c);
                        if (voxI<tr){
                            ntr += *(L + 3*nL+c);
                            if (ntr>nTr)
                                break;
                        }
                        if (nR>0){
                            if (down){
                                if (voxI>I0){
                                    nr++;
                                    if (nr==nR)
                                        break;
                                }
                            } else if (voxI<I0){
                                nr++;
                                if (nr==nR){
                                    down=true;
                                    nr=0;
                                }
                            }
                        }
                        I0 = voxI;
                        c++;
                    }
                    if (t==1)
                        Ja1 = J1;
                }
                if (J1 < *(J+m)){
                    *(J+m) = J1;
                    *(A+m) = n+1;
                    *(Ja+m) = Ja1;
                }
            } else
                *(J+m) = 0;
        }
    }
}


#include <mex.h>
#include<iostream>
#include<iomanip>
#include <omp.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Dense>
#include<vector>

using namespace std;
using namespace Eigen;

//#define PRINT_DEBUGGING_OUTPUT
#define USE_SCALAR_IMPLEMENTATION
// #define USE_SSE_IMPLEMENTATION
// #define USE_AVX_IMPLEMENTATION

#define COMPUTE_V_AS_MATRIX
//#define COMPUTE_V_AS_QUATERNION
#define COMPUTE_U_AS_MATRIX
//#define COMPUTE_U_AS_QUATERNION

#include "Singular_Value_Decomposition_Preamble.hpp"
#define e 2.71828182846

void GetDiff(double* Eh,double* Eg,double* E,double* Js,double *et,int nt)
{
    double sum=0;
    int energyType=*et;

    

#pragma omp parallel for reduction(+:sum)
    for(int ii=0;ii<nt;ii++)
    {
        Map<MatrixXd> Block_Eg(Eg + ii * 9, 9, 1);
        Map<MatrixXd> Block_Eh(Eh+ii*9*9, 9, 9);
        float A11, A21, A31, A12, A22, A32, A13, A23, A33;
        float u11, u12, u13, u21, u22, u23, u31, u32, u33, v11, v12, v13, v21, v22, v23, v31, v32, v33;
        float sigma[3];
        
        int Adress_begin = 9 * ii;
        A11 = (float)Js[Adress_begin + 0];
        A12 = (float)Js[Adress_begin + 1];
        A13 = (float)Js[Adress_begin + 2];
        A21 = (float)Js[Adress_begin + 3];
        A22 = (float)Js[Adress_begin + 4];
        A23 = (float)Js[Adress_begin + 5];
        A31 = (float)Js[Adress_begin + 6];
        A32 = (float)Js[Adress_begin + 7];
        A33 = (float)Js[Adress_begin + 8];
        
        float norm_inverse = (float)(1. / sqrt((double)A11*(double)A11 + (double)A21*(double)A21 + (double)A31*(double)A31
                + (double)A12*(double)A12 + (double)A22*(double)A22 + (double)A32*(double)A32
                + (double)A13*(double)A13 + (double)A23*(double)A23 + (double)A33*(double)A33));
        
        A11 *= norm_inverse;
        A21 *= norm_inverse;
        A31 *= norm_inverse;
        A12 *= norm_inverse;
        A22 *= norm_inverse;
        A32 *= norm_inverse;
        A13 *= norm_inverse;
        A23 *= norm_inverse;
        A33 *= norm_inverse;
        
        
#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"
        ENABLE_SCALAR_IMPLEMENTATION(Sa11.f = A11;)                                        //ENABLE_SSE_IMPLEMENTATION(Va11 = _mm_set1_ps(A11);)                                    ENABLE_AVX_IMPLEMENTATION(Va11 = _mm256_set1_ps(A11);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa21.f = A21;)                                    //    ENABLE_SSE_IMPLEMENTATION(Va21 = _mm_set1_ps(A21);)                                    ENABLE_AVX_IMPLEMENTATION(Va21 = _mm256_set1_ps(A21);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa31.f = A31;)                                   //     ENABLE_SSE_IMPLEMENTATION(Va31 = _mm_set1_ps(A31);)                                    ENABLE_AVX_IMPLEMENTATION(Va31 = _mm256_set1_ps(A31);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa12.f = A12;)                                   //     ENABLE_SSE_IMPLEMENTATION(Va12 = _mm_set1_ps(A12);)                                    ENABLE_AVX_IMPLEMENTATION(Va12 = _mm256_set1_ps(A12);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa22.f = A22;)                                   //     ENABLE_SSE_IMPLEMENTATION(Va22 = _mm_set1_ps(A22);)                                    ENABLE_AVX_IMPLEMENTATION(Va22 = _mm256_set1_ps(A22);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa32.f = A32;)                                   //     ENABLE_SSE_IMPLEMENTATION(Va32 = _mm_set1_ps(A32);)                                    ENABLE_AVX_IMPLEMENTATION(Va32 = _mm256_set1_ps(A32);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa13.f = A13;)                                    //    ENABLE_SSE_IMPLEMENTATION(Va13 = _mm_set1_ps(A13);)                                    ENABLE_AVX_IMPLEMENTATION(Va13 = _mm256_set1_ps(A13);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa23.f = A23;)                                   //     ENABLE_SSE_IMPLEMENTATION(Va23 = _mm_set1_ps(A23);)                                    ENABLE_AVX_IMPLEMENTATION(Va23 = _mm256_set1_ps(A23);)
        ENABLE_SCALAR_IMPLEMENTATION(Sa33.f = A33;)                                   //     ENABLE_SSE_IMPLEMENTATION(Va33 = _mm_set1_ps(A33);)  									ENABLE_AVX_IMPLEMENTATION(Va33 = _mm256_set1_ps(A33);)
        
        
#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"
        ENABLE_SCALAR_IMPLEMENTATION(u11 = Su11.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u21 = Su21.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u31 = Su31.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u12 = Su12.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u22 = Su22.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u32 = Su32.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u13 = Su13.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u23 = Su23.f;)
        ENABLE_SCALAR_IMPLEMENTATION(u33 = Su33.f;)
        
        ENABLE_SCALAR_IMPLEMENTATION(v11 = Sv11.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v21 = Sv21.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v31 = Sv31.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v12 = Sv12.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v22 = Sv22.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v32 = Sv32.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v13 = Sv13.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v23 = Sv23.f;)
        ENABLE_SCALAR_IMPLEMENTATION(v33 = Sv33.f;)
        
        ENABLE_SCALAR_IMPLEMENTATION(sigma[0] = Sa11.f;)
        ENABLE_SCALAR_IMPLEMENTATION(sigma[1] = Sa22.f;)
        ENABLE_SCALAR_IMPLEMENTATION(sigma[2] = Sa33.f;)
        
        Matrix3d U;
        Matrix3d V;
        U << u11, u12, u13, u21, u22, u23, u31, u32, u33;
        V << v11, v12, v13, v21, v22, v23, v31, v32, v33;
        Matrix<double, 9, 9> p;
        p << u11*V, u12*V, u13*V, u21*V, u22*V, u23*V, u31*V, u32*V, u33*V;
        Matrix<double, 9, 9> pq;
        pq << p.col(0), p.col(4), p.col(8), p.col(1), p.col(3), p.col(2), p.col(6), p.col(5), p.col(7);
        Matrix<double, 9, 9>leftMat;
        leftMat =pq;
        for (int j = 0; j<3; j++)
        {
            sigma[j] /= norm_inverse;
        }

        Matrix<double, 9, 9> K = MatrixXd::Zero(9, 9);
        Vector3d partialE;
        switch(energyType)
        {
            case 1:    //Eiso
            {
                for (int j = 0; j<3; j++)
                {
                    K(j, j) = 2 * (1 + 3 / pow(sigma[j], 4));
                    int id1, id2;
                    id1 = (j * 2) % 3; id2 = (id1 + 1) % 3;
                    double a = 1 + 1 / pow((sigma[id1] * sigma[id2]), 2);
                    double b = (sigma[id1] * sigma[id1] + sigma[id2] * sigma[id2]) / pow((sigma[id1] * sigma[id2]), 3);
                    if (a>b)
                        K.block<2, 2>(3 + 2 * j, 3 + 2 * j) << 2 * a, 2 * b, 2 * b, 2 * a;
                    else
                        K.block<2, 2>(3 + 2 * j, 3 + 2 * j) << a + b, a + b, a + b, a + b;
                }
                partialE << sigma[0] - 1 / pow(sigma[0], 3), sigma[1] - 1 / pow(sigma[1], 3), sigma[2] - 1 / pow(sigma[2], 3);
                partialE=2 * partialE;
                sum += (sigma[0] * sigma[0] + 1 / sigma[0] / sigma[0] + sigma[1] * sigma[1] + 1 / sigma[1] / sigma[1] + sigma[2] * sigma[2] + 1 / sigma[2] / sigma[2] - 6);
                break;
            }
            case 5:
            {
                for (int j = 0; j<3; j++)
                {
                    K(j, j) = 2;
                    int id1, id2;
                    id1 = (j * 2) % 3; id2 = (id1 + 1) % 3;
                    double a = 1 - 1 / (sigma[id1] + sigma[id2]);
                    double b = 1 / (sigma[id1] + sigma[id2]);
                    if (a>b)
                        K.block<2, 2>(3 + 2 * j, 3 + 2 * j) << 2 * a, 2 * b, 2 * b, 2 * a;
                    else
                        K.block<2, 2>(3 + 2 * j, 3 + 2 * j) << a + b, a + b, a + b, a + b;
                }
                partialE << 2*(sigma[0] - 1),2*(sigma[1] - 1),2*(sigma[2] - 1);
                sum +=((sigma[0] -1)* (sigma[0] -1) + (sigma[1] -1)* (sigma[1] -1)+(sigma[2] -1)* (sigma[2] -1));
                break;
            }
        }
        Block_Eg = leftMat.block<9, 3>(0, 0)*partialE;
        Block_Eh = leftMat*K*leftMat.transpose();

        
        
    }

    
    
    *E=sum;
}
// [Eh,Eg,E]=cppGetJacobDiff(Js,et,temp)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    int nv = mxGetN(prhs[0])/3;

    
    plhs[0] = mxDuplicateArray(prhs[2]);
    plhs[1]=mxCreateDoubleMatrix(9*nv,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
    double* E=mxGetPr(plhs[2]);
    *E=0;
    
//     Map<MatrixXd> Js(mxGetPr(prhs[0]), 3, 3*nt);
    
    GetDiff(mxGetPr(plhs[0]),mxGetPr(plhs[1]),E,mxGetPr(prhs[0]),mxGetPr(prhs[1]),nv);
}
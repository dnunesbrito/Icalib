/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OptFunctions.h
 * Author: darlan
 *
 * Created on 11 de Maio de 2016, 09:35
 */

#ifndef OPTFUNCTIONS_H
#define OPTFUNCTIONS_H
#include <cstdlib>
#include <sstream>
#include "Constants.h"
#include "IntegerVector.h"
#include "Interval.h"
#include <iostream>
#include <Vector.h>
#include <IntervalVector.h>
#include <Matrix.h>
#include <IntervalMatrix.h>
#include <Functions.h>
#include <Utilities.h>
#include <vector>
#include <LSS.h>
#include <GlobalOpt/Expand.h>
#include <GlobalOpt/AppList.h>
#include <GlobalOpt/VecList.h>
#include <GlobalOpt/VecUtils.h>
#include <GlobalOpt/UnconstrainedOpt.h>
#include <time.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/Jacobi>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/src/Core/VectorBlock.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>
#include <dirent.h>
#include <AutoDiff/MATRIX_AUTODIFF.hpp>
#include <BIAS/BiasType.h>
#include <math.h>
#include <ga/GASimpleGA.h>
#include <ga/GABin2DecGenome.h>
#include <ga/std_stream.h>
#include <LinOpt/LinOpt.hpp>
#include <NonLinOpt/NonLinOpt.hpp>
#include <EFunctions.h>
#include "EINTERVAL.hpp"
#include "EINTERVALVECTOR.hpp"
#include "EINTERVAL_MATRIX.hpp"

 

extern struct OPTMICALIBDATA{
    INTERVAL_MATRIX Xw;
    INTERVAL_MATRIX xc;
    INTERVAL_MATRIX K;
    INTERVAL_MATRIX RRt;
    INTERVAL_MATRIX max_min_xc;
    INTERVAL_VECTOR t;
    INTERVAL_VECTOR EulAngles;
    INTERVAL_VECTOR XVals;
    INTERVAL_VECTOR X;
    EINTERVAL_VECTOR EIVX;
    EINTERVAL_VECTOR grd_lim;
    REAL UpperBound;
}optmicalibdata;

/*Função utilizada na otimização do sistema Ap=0 
 O parâmetro de entrada é o vetor p composto da seguinte
 forma:
 p(1-3) = Angulos de Euler da matriz de rotação sendo
          (Angulo x,Angulo y, Angulo z)
 p(4-6) = Vetor de translação da câmera.
 p(7-8) = O fm_x e fm_y.
 Para se utilizar esta função a matriz A_T1C deve ser 
 inicializada em algum ponto do programa principal com 
 a matriz A do sistema.
 A matriz IK também deve ser inicializada sendo que 
 desta matriz são utilizados os valores de u0 e v0.
 Este é a matriz de parâmetros intrínsecos.*/
INTERVAL_AUTODIFF ApFuncAutoDiff (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF ApFuncAutoDiffFF (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF ApFuncAutoDiffFA (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF AxXFuncAD (CONST INTERVAL_AUTODIFF & X,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF cFunction(CONST INTERVAL_AUTODIFF& c);
INTERVAL AxXFuncADI (CONST INTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
INTERVAL AxXFunc2ADI (CONST INTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF PXFuncAD (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF KRtFuncAD (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata);
INTERVAL PXFuncADI (CONST INTERVAL_VECTOR & x,std::shared_ptr<VOID> userdata);
INTERVAL PXFuncADII (CONST INTERVAL_VECTOR & x,std::shared_ptr<VOID> userdata);
REAL PXFuncR (CONST VECTOR & x,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF kFunc(CONST INTERVAL_AUTODIFF & x);
REAL AxXFuncADR (CONST VECTOR & X,std::shared_ptr<VOID> userdata);
INTERVAL_MATRIX X_Y_triangula(CONST INTERVAL_MATRIX& K,CONST INTERVAL_MATRIX& R,CONST INTERVAL_VECTOR& t,CONST INTERVAL_VECTOR& x);
INTERVAL xPXHC(CONST INTERVAL_VECTOR& X,std::shared_ptr<VOID> userdata);
INTERVAL_MATRIX rtm2Eul2(CONST INTERVAL_MATRIX & R);
INTERVAL_VECTOR EulRotmHC(INTERVAL& psi,INTERVAL& phi,CONST INTERVAL_MATRIX &r);
EINTERVAL More_Garbow_Hillstrom_func(CONST EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF More_Garbow_Hillstrom_AD(CONST INTERVAL_AUTODIFF& X, std::shared_ptr<VOID> userdata);
EINTERVAL More_Garbow_Hillstrom_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
EINTERVAL_VECTOR More_Garbow_Hillstrom_grd(CONST EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF More_Garbow_Hillstrom_grd_AD(CONST INTERVAL_AUTODIFF & X,INT i,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF More_Garbow_Hillstrom_Bgrd_AD(CONST INTERVAL_AUTODIFF & X,INT i,std::shared_ptr<VOID> userdata);
EINTERVAL More_Garbow_Hillstrom_grd_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
EINTERVAL More_Garbow_Hillstrom_Bgrd_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
EINTERVAL More_Garbow_Hillstrom_hess_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata);
INTERVAL_AUTODIFF More_Garbow_Hillstrom_grd(CONST INTERVAL_AUTODIFF & X,INT i,std::shared_ptr<VOID> userdata);
#endif /* OPTFUNCTIONS_H */


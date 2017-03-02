/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <float.h>

#include "OptFunctions.h"
#include <EFunctions.h>

extern struct OPTMICALIBDATA optmicalibdata;
struct AbLinsys{
    INTERVAL_MATRIX A;
    INTERVAL_VECTOR b;
    INTERVAL_VECTOR X;
    INTERVAL_MATRIX IK;
};

INTERVAL_AUTODIFF ApFuncAutoDiff (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata);
    INTERVAL_MATRIX A_T1C = Albl.get()->A,IK = Albl.get()->IK;
    MATRIX_AUTODIFF_INTERVAL RMAD,RMADt,RRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian),
            KRRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL tAD(3,1,Dimension(xIAD),xIAD.ComputeHessian);
    INTERVAL_AUTODIFF Vf(Dimension(xIAD),xIAD.ComputeHessian);
    RMAD = Eul2rtmAD<INTERVAL,INTERVAL_VECTOR,INTERVAL_MATRIX>(xIAD);
    MATRIX_AUTODIFF_INTERVAL KMAD(3,3,Dimension(xIAD),xIAD.ComputeHessian);
    KMAD(1,1)=xIAD(7);
    KMAD(1,3).value = IK(1,3);
    KMAD(2,2)=xIAD(8);
    KMAD(2,3).value = IK(2,3);
    KMAD(3,3).value = 1.0;
    tAD(1,1)=xIAD(4);
    tAD(2,1)=xIAD(5);
    tAD(3,1)=xIAD(6);
    RMADt = -RMAD;
    RMADt = RMADt*tAD;
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRtMAD(i,j)=RMAD(i,j);
    for(int i = 1;i <= 3;i++)
        RRtMAD(i,4)=RMADt(i,1);
    KRRtMAD = KMAD*RRtMAD;
    MATRIX_AUTODIFF_INTERVAL p(12,1,Dimension(xIAD),xIAD.ComputeHessian);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp,1)=KRRtMAD(i,j);
            contp++;
        }
    }
    MATRIX_AUTODIFF_INTERVAL A(RowDimension(A_T1C),ColDimension(A_T1C),A_T1C,Dimension(xIAD),xIAD.ComputeHessian),VfMAD(12,1,Dimension(xIAD),xIAD.ComputeHessian),
            SumVfMAD(1,1,Dimension(xIAD),xIAD.ComputeHessian);
    VfMAD = A*p;
    for(int i = 1;i <= VfMAD.nRows;i++)
        SumVfMAD(1,1)=SumVfMAD(1,1)+Sqr(VfMAD(i,1));
    SumVfMAD(1,1) = Sqrt(SumVfMAD(1,1));
    Vf.fkt() = SumVfMAD(1,1).value;
    Vf.grd() = SumVfMAD(1,1).grd;
    if(xIAD.ComputeHessian)
        Vf.hessian() = SumVfMAD(1,1).hess;
    return Vf;
}

INTERVAL_AUTODIFF ApFuncAutoDiffFF (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata);
    INTERVAL_MATRIX A_T1C = Albl.get()->A,IK = Albl.get()->IK;
    MATRIX_AUTODIFF_INTERVAL RMAD,RMADt,RRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian),
            KRRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL tAD(3,1,Dimension(xIAD),xIAD.ComputeHessian);
    INTERVAL_AUTODIFF Vf(Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL EulAng(1,3,Dimension(xIAD),xIAD.ComputeHessian);
    EulAng(1,1) = xIAD(1);
    EulAng(1,2) = xIAD(2);
    EulAng(1,3) = xIAD(3);
    RMAD = Eul2rtmAD<INTERVAL,INTERVAL_VECTOR,INTERVAL_MATRIX>(EulAng);
    MATRIX_AUTODIFF_INTERVAL KMAD(3,3,Dimension(xIAD),xIAD.ComputeHessian);
    KMAD(1,1).value = IK(1,1);
    KMAD(1,3).value = IK(1,3);
    KMAD(2,2).value = IK(2,2);
    KMAD(2,3).value = IK(2,3);
    KMAD(3,3).value = 1.0;
    tAD(1,1)=xIAD(4);
    tAD(2,1)=xIAD(5);
    tAD(3,1)=xIAD(6);
    RMADt = -RMAD;
    RMADt = RMADt*tAD;
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRtMAD(i,j)=RMAD(i,j);
    for(int i = 1;i <= 3;i++)
        RRtMAD(i,4)=RMADt(i,1);
    KRRtMAD = KMAD*RRtMAD;
    MATRIX_AUTODIFF_INTERVAL p(12,1,Dimension(xIAD),xIAD.ComputeHessian);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp,1)=KRRtMAD(i,j);
            contp++;
        }
    }
    MATRIX_AUTODIFF_INTERVAL A(RowDimension(A_T1C),ColDimension(A_T1C),A_T1C,Dimension(xIAD),xIAD.ComputeHessian),VfMAD(12,1,Dimension(xIAD),xIAD.ComputeHessian),
            SumVfMAD(1,1,Dimension(xIAD),xIAD.ComputeHessian);
    VfMAD = A*p;
    for(int i = 1;i <= VfMAD.nRows;i++)
        SumVfMAD(1,1)=SumVfMAD(1,1)+Abs(VfMAD(i,1));
    Vf.fkt() = SumVfMAD(1,1).value;
    Vf.grd() = SumVfMAD(1,1).grd;
    if(xIAD.ComputeHessian)
        Vf.hessian() = SumVfMAD(1,1).hess;
    return Vf;
}

INTERVAL_AUTODIFF ApFuncAutoDiffFA (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata){
    std::shared_ptr<INTERVAL_MATRIX> A_T1C = std::static_pointer_cast<INTERVAL_MATRIX>(userdata);    
    MATRIX_AUTODIFF_INTERVAL RMAD,RMADt,RRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian),
            KRRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL tAD(3,1,Dimension(xIAD),xIAD.ComputeHessian);
    INTERVAL_AUTODIFF Vf(Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL EulAng(1,3,Dimension(xIAD),xIAD.ComputeHessian);
//    EulAng(1,1).value = FixedAngles(1);
//    EulAng(1,2).value = FixedAngles(2);
//    EulAng(1,3).value = FixedAngles(3);
    RMAD = Eul2rtmAD<INTERVAL,INTERVAL_VECTOR,INTERVAL_MATRIX>(EulAng);
    MATRIX_AUTODIFF_INTERVAL KMAD(3,3,Dimension(xIAD),xIAD.ComputeHessian);
    KMAD(1,1)=xIAD(4);
//    KMAD(1,3).value = IK(1,3);
    KMAD(2,2)=xIAD(5);
//    KMAD(2,3).value = IK(2,3);
    KMAD(3,3).value = 1.0;
    tAD(1,1)=xIAD(1);
    tAD(2,1)=xIAD(2);
    tAD(3,1)=xIAD(3);
    RMADt = -RMAD;
    RMADt = RMADt*tAD;
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRtMAD(i,j)=RMAD(i,j);
    for(int i = 1;i <= 3;i++)
        RRtMAD(i,4)=RMADt(i,1);
    KRRtMAD = KMAD*RRtMAD;
    MATRIX_AUTODIFF_INTERVAL p(12,1,Dimension(xIAD),xIAD.ComputeHessian);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp,1)=KRRtMAD(i,j);
            contp++;
        }
    }
    MATRIX_AUTODIFF_INTERVAL A(RowDimension((*A_T1C.get())),ColDimension((*A_T1C.get())),(*A_T1C.get()),Dimension(xIAD),xIAD.ComputeHessian),VfMAD(12,1,Dimension(xIAD),xIAD.ComputeHessian),
            SumVfMAD(1,1,Dimension(xIAD),xIAD.ComputeHessian);
    VfMAD = A*p;
    for(int i = 1;i <= VfMAD.nRows;i++)
        SumVfMAD(1,1)=SumVfMAD(1,1)+Abs(VfMAD(i,1));
    MATRIX_AUTODIFF_ELEMENT_INTERVAL KMAD1_2;
    KMAD1_2 = KMAD(1,1)-KMAD(2,2);
    SumVfMAD(1,1)=SumVfMAD(1,1)+Abs(KMAD1_2);
    Vf.fkt() = SumVfMAD(1,1).value;
    Vf.grd() = SumVfMAD(1,1).grd;
    if(xIAD.ComputeHessian)
        Vf.hessian() = SumVfMAD(1,1).hess;
    return Vf;
}

INTERVAL_AUTODIFF AxXFuncAD (CONST INTERVAL_AUTODIFF & X,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata);
    INTERVAL_MATRIX A_T1C = Albl.get()->A;
    INTERVAL_VECTOR b = Albl.get()->b;
    INTERVAL_AUTODIFF Vf;
    MATRIX_AUTODIFF_INTERVAL A(RowDimension(A_T1C),ColDimension(A_T1C),A_T1C,Albl.get()->A.cols(),X.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL VfMAD;
    MATRIX_AUTODIFF_INTERVAL SumVfMAD(1,1,Albl.get()->A.cols(),X.ComputeHessian);
    VfMAD = A*X-b;
    for(int i = 1;i <= VfMAD.nRows;i++){
        SumVfMAD(1,1)=SumVfMAD(1,1)+Sqr(VfMAD(i,1));
    }
    SumVfMAD(1,1) = Sqrt(SumVfMAD(1,1));
    Vf.fkt() = SumVfMAD(1,1).value;
    Vf.grd() = SumVfMAD(1,1).grd;
    if(X.ComputeHessian)
        Vf.hessian() = SumVfMAD(1,1).hess;
    return Vf;
}
INTERVAL_AUTODIFF cFunction(CONST INTERVAL_AUTODIFF& c,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata); 
    INTERVAL_MATRIX A_T1C = Albl.get()->A;
    INTERVAL_AUTODIFF A(1,c.ComputeHessian),B(1,c.ComputeHessian),C(1,c.ComputeHessian);
    A.fkt() = A_T1C(1,1);
    B.fkt() = A_T1C(1,2);
    C.fkt() = A_T1C(1,3);
    C = Sqr(A*c(1))+Sqr(B*c(1))+Sqr(C*c(1))-1.0;
    return Sqr(A*c(1))+Sqr(B*c(1))+Sqr(C*c(1))-1.0;
}
INTERVAL AxXFuncADI (CONST INTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata); 
    INTERVAL_MATRIX A_T1C = Albl.get()->A;
    INTERVAL_VECTOR b = Albl.get()->b;
    INTERVAL_VECTOR IVVf;
    IVVf = A_T1C*X-b;
    return Norm(IVVf);
}
INTERVAL AxXFunc2ADI (CONST INTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata);
    INTERVAL_VECTOR IVVf;
    IVVf = (*Albl.get()).A*X-(*Albl.get()).b;
    return Norm(IVVf);
}
REAL AxXFuncADR (CONST VECTOR & X,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata);
    INTERVAL_MATRIX A_T1C = Albl.get()->A;
    INTERVAL_VECTOR b = Albl.get()->b;
    VECTOR IVVf;
    IVVf = Mid(A_T1C)*X-Mid(b);
    return Norm(IVVf);
}
INTERVAL PXFuncADI (CONST INTERVAL_VECTOR & x,std::shared_ptr<VOID> userdata){
    INTERVAL_MATRIX R(3,3),RRt(3,4),P(3,4);
    INTERVAL_VECTOR t(3);
    double Vf = 0;
    INTERVAL_VECTOR EulAng(3);
    EulAng(1) = x(1);
    EulAng(2) = x(2);
    EulAng(3) = x(3);
    R = Eul2rtm(EulAng);
    cout << "x(PXFuncADI)=" << x << endl;
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRt(i,j)=R(i,j);
    INTERVAL_MATRIX K(3,3);
    Clear(K);
    K(1,1) = x(7);
    K(1,3) = optmicalibdata.K(1,3);
    K(2,2) = x(8);
    K(2,3) = optmicalibdata.K(2,3);
    K(3,3) = 1.0;
    t(1)=x(4);
    t(2)=x(5);
    t(3)=x(6);
    t = -R*t;
    for(int i = 1;i <= 3;i++)
        for(int j = 4;j <= 4;j++)
            RRt(i,j)=t(i);
    P = K*RRt;
    INTERVAL_MATRIX xccalc;
    xccalc = P*optmicalibdata.Xw;
    GetNonHomogeneous(xccalc,xccalc);
    INTERVAL ptinter;
    for(int j = 1;j <= ColDimension(xccalc);j++){
        if(AnyNAN(Col(xccalc,j))){
            Vf = -FLT_MAX;
            continue;
        }
        for(int i = 1;i <= RowDimension(xccalc);i++){
            if(Inf(xccalc(i,j)) == Machine::NegInfinity || Sup(xccalc(i,j)) == Machine::PosInfinity){
                Vf = -FLT_MAX;
                continue;
            }
            if(Intersection(ptinter,xccalc(i,j),optmicalibdata.xc(i,j)) == 1){
                Vf += 1000*ptinter.diam()-(xccalc(i,j).diam());
            }else{
                Vf -= (Diam(xccalc(i,j))+Diam(optmicalibdata.xc(i,j)));
            }
            if(Sup(Sqr(xccalc(i,j)-optmicalibdata.xc(i,j))) != Machine::NaN && Sup(Sqr(xccalc(i,j)-optmicalibdata.xc(i,j))) != Machine::PosInfinity
                    && Sup(Sqr(xccalc(i,j)-optmicalibdata.xc(i,j))) != Machine::NegInfinity)
                Vf -= Sup(Sqr(xccalc(i,j)-optmicalibdata.xc(i,j)));
        }
    }
    return Vf;
}
INTERVAL PXFuncADII (CONST INTERVAL_VECTOR & x,std::shared_ptr<VOID> userdata){
    std::shared_ptr<INTERVAL_MATRIX> A_T1C = std::static_pointer_cast<INTERVAL_MATRIX>(userdata);    
    INTERVAL_MATRIX R(3,3),RRt(3,4),P(3,4);
    INTERVAL_VECTOR t(3);
    INTERVAL_VECTOR Vf;
    INTERVAL_VECTOR EulAng(3);
    EulAng(1) = x(1);
    EulAng(2) = x(2);
    EulAng(3) = x(3);
    R = Eul2rtm(EulAng);
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRt(i,j)=R(i,j);
    INTERVAL_MATRIX K(3,3);
    Clear(K);
    K(1,1) = x(7);
    K(1,3) = optmicalibdata.K(1,3);
    K(2,2) = x(8);
    K(2,3) = optmicalibdata.K(2,3);
    K(3,3) = 1.0;
    t(1)=x(4);
    t(2)=x(5);
    t(3)=x(6);
    t = -R*t;
    for(int i = 1;i <= 3;i++)
        for(int j = 4;j <= 4;j++)
            RRt(i,j)=t(i);
    P = K*RRt;
    INTERVAL_VECTOR p(12);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp)=P(i,j);
            contp++;
        }
    }
    Vf = (*A_T1C.get())*p;
    return Norm(Vf);
}

INTERVAL_AUTODIFF PXFuncAD (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata){
    std::shared_ptr<INTERVAL_MATRIX> A_T1C = std::static_pointer_cast<INTERVAL_MATRIX>(userdata);    
    MATRIX_AUTODIFF_INTERVAL RMAD,RMADt(3,1,Dimension(xIAD),xIAD.ComputeHessian),RRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian),
            KRRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL tAD(3,1,Dimension(xIAD),xIAD.ComputeHessian);
    INTERVAL_AUTODIFF Vf(Dimension(xIAD),xIAD.ComputeHessian);
    RMAD = Eul2rtmAD<INTERVAL,INTERVAL_VECTOR,INTERVAL_MATRIX>(xIAD);
    MATRIX_AUTODIFF_INTERVAL KMAD(3,3,Dimension(xIAD),xIAD.ComputeHessian);
    KMAD(1,1)=xIAD(7);
    KMAD(1,3).value = optmicalibdata.K(1,3);
    KMAD(2,2)=xIAD(8);
    KMAD(2,3).value = optmicalibdata.K(2,3);
    KMAD(3,3).value = 1.0;
    tAD(1,1)=xIAD(4);
    tAD(2,1)=xIAD(5);
    tAD(3,1)=xIAD(6);
    RMADt = -RMAD;
    RMADt = RMADt*tAD;
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRtMAD(i,j)=RMAD(i,j);
    for(int i = 1;i <= 3;i++)
        RRtMAD(i,4)=RMADt(i,1);
    KRRtMAD = KMAD*RRtMAD;
    MATRIX_AUTODIFF_INTERVAL p(12,1,Dimension(xIAD),xIAD.ComputeHessian);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp,1)=KRRtMAD(i,j);
            contp++;
        }
    }
    MATRIX_AUTODIFF_INTERVAL A(RowDimension((*A_T1C.get())),ColDimension((*A_T1C.get())),(*A_T1C.get()),Dimension(xIAD),xIAD.ComputeHessian),
            VfMAD(12,1,Dimension(xIAD),xIAD.ComputeHessian),
            SumVfMAD(1,1,Dimension(xIAD),xIAD.ComputeHessian);
    VfMAD = A*p;
    for(int i = 1;i <= VfMAD.nRows;i++){
        SumVfMAD(1,1)=SumVfMAD(1,1)+Sqr(VfMAD(i,1));
    }
    SumVfMAD(1,1) = Sqrt(SumVfMAD(1,1));
    Vf.fkt() = SumVfMAD(1,1).value;
    Vf.grd() = SumVfMAD(1,1).grd;
    if(xIAD.ComputeHessian)
        Vf.hessian() = SumVfMAD(1,1).hess;
    return Vf;
}

INTERVAL_AUTODIFF KRtFuncAD (CONST INTERVAL_AUTODIFF & xIAD,std::shared_ptr<VOID> userdata){
    std::shared_ptr<INTERVAL_MATRIX> A_T1C = std::static_pointer_cast<INTERVAL_MATRIX>(userdata);    
    MATRIX_AUTODIFF_INTERVAL RMAD,RMADt(3,1,Dimension(xIAD),xIAD.ComputeHessian),RRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian),
            KRRtMAD(3,4,Dimension(xIAD),xIAD.ComputeHessian);
    MATRIX_AUTODIFF_INTERVAL tAD(3,1,Dimension(xIAD),xIAD.ComputeHessian);
    INTERVAL_AUTODIFF Vf(Dimension(xIAD),xIAD.ComputeHessian);
    RMAD = Eul2rtmAD<INTERVAL,INTERVAL_VECTOR,INTERVAL_MATRIX>(xIAD);
    MATRIX_AUTODIFF_INTERVAL KMAD(3,3,Dimension(xIAD),xIAD.ComputeHessian);
    KMAD(1,1)=xIAD(7);
    KMAD(1,3)=xIAD(8);
    KMAD(2,2)=xIAD(9);
    KMAD(2,3)=xIAD(10);
    KMAD(3,3).value = 1.0;
    tAD(1,1)=xIAD(4);
    tAD(2,1)=xIAD(5);
    tAD(3,1)=xIAD(6);
    RMADt = -RMAD;
    RMADt = RMADt*tAD;
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRtMAD(i,j)=RMAD(i,j);
    for(int i = 1;i <= 3;i++)
        RRtMAD(i,4)=RMADt(i,1);
    KRRtMAD = KMAD*RRtMAD;
    MATRIX_AUTODIFF_INTERVAL p(12,1,Dimension(xIAD),xIAD.ComputeHessian);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp,1)=KRRtMAD(i,j);
            contp++;
        }
    }
    MATRIX_AUTODIFF_INTERVAL A(RowDimension((*A_T1C.get())),ColDimension((*A_T1C.get())),(*A_T1C.get()),Dimension(xIAD),xIAD.ComputeHessian),
            VfMAD(12,1,Dimension(xIAD),xIAD.ComputeHessian),
            SumVfMAD(1,1,Dimension(xIAD),xIAD.ComputeHessian);
    VfMAD = A*p;
    for(int i = 1;i <= VfMAD.nRows;i++){
        SumVfMAD(1,1)=SumVfMAD(1,1)+Sqr(VfMAD(i,1));
    }
    SumVfMAD(1,1) = Sqrt(SumVfMAD(1,1));
    Vf.fkt() = SumVfMAD(1,1).value;
    Vf.grd() = SumVfMAD(1,1).grd;
    if(xIAD.ComputeHessian)
        Vf.hessian() = SumVfMAD(1,1).hess;
    return Vf;
}

REAL PXFuncR (CONST VECTOR & Vx,std::shared_ptr<VOID> userdata){
    std::shared_ptr<AbLinsys> Albl = std::static_pointer_cast<AbLinsys>(userdata);
    INTERVAL_MATRIX A_T1C = Albl.get()->A;
    INTERVAL_MATRIX R(3,3),RRt(3,4),P(3,4);
    INTERVAL_VECTOR t(3);
    INTERVAL_VECTOR Vf;
    INTERVAL_VECTOR EulAng(3);
    INTERVAL_VECTOR x(Albl.get()->X.nrows());
    x = Hull(Vx);
    EulAng(1) = x(1);
    EulAng(2) = x(2);
    EulAng(3) = x(3);
    R = Eul2rtm(EulAng);
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRt(i,j)=R(i,j);
    INTERVAL_MATRIX K(3,3);
    K(1,1) = x(7);
    K(1,3) = optmicalibdata.K(1,3);
    K(2,2) = x(8);
    K(2,3) = optmicalibdata.K(2,3);
    K(3,3) = 1.0;
    t(1)=x(4);
    t(2)=x(5);
    t(3)=x(6);
    t = -R*t;
    for(int i = 1;i <= 3;i++)
        for(int j = 4;j <= 4;j++)
            RRt(i,j)=t(i);
    P = K*RRt;
    INTERVAL_VECTOR p(12);
    int contp = 1;
    for(int i = 1;i <= 3;i++){
        for(int j = 1;j <= 4;j++){
            p(contp)=P(i,j);
            contp++;
        }
    }
    Vf = A_T1C*p;
    INTERVAL IVf;
    IVf = Norm(Vf);
    return Mid(IVf);
}

INTERVAL_AUTODIFF kFunc(CONST INTERVAL_AUTODIFF & x){
    INTERVAL_AUTODIFF out(1);
    out.fkt() = optmicalibdata.XVals(3);
    Clear(out.grd());
    out.fkt() = out.fkt()+x(1).fkt()*optmicalibdata.XVals(2);
    out.grd()(1) = optmicalibdata.XVals(2);
    out.fkt() = Sqr(x(1).fkt())*optmicalibdata.XVals(3);
    out.grd()(1) = out.grd()(1)+2*optmicalibdata.XVals(3)*x.grd()(1);
    return out;
}

INTERVAL_MATRIX X_Y_triangula(CONST INTERVAL_MATRIX& K,CONST INTERVAL_MATRIX& R,CONST INTERVAL_VECTOR& t,CONST INTERVAL_VECTOR& x){
    INTERVAL_MATRIX X(3,2);
    INTERVAL f = K(1,1);
    X(2,1) = (t(3)*x(2) + x(2)*((R(3,1)*(x(1)*(t(2)) - x(2)*(t(1))))/
            (R(1,1)*x(2) - R(1,2)*x(1))) - f*(t(2) + (R(2,1)*(x(1)*(t(2)) - x(2)*(t(1))))/
            (R(1,1)*x(2) - R(1,2)*x(1))))/(f*(R(2,2) - (R(2,1)*(R(1,2)*x(2) - R(2,2)*x(1)))/(R(1,1)*x(2) - R(1,2)*x(1))) 
            - x(2)*(R(3,2) - (R(3,1)*(R(1,2)*x(2) - R(2,2)*x(1)))/(R(1,1)*x(2) - R(1,2)*x(1))));
    X(2,2) = (t(3)*x(2) + x(2)*(2*R(3,3) + (R(3,1)*(x(1)*(t(2) + 2*R(1,3)) - x(2)*(t(1) + 2*R(1,3))))/
            (R(1,1)*x(2) - R(1,2)*x(1))) - f*(t(2) + 2*R(2,3) + (R(2,1)*(x(1)*(t(2) + 2*R(1,3)) - x(2)*(t(1) + 2*R(1,3))))/
            (R(1,1)*x(2) - R(1,2)*x(1))))/(f*(R(2,2) - (R(2,1)*(R(1,2)*x(2) - R(2,2)*x(1)))/(R(1,1)*x(2) - R(1,2)*x(1)))
            - x(2)*(R(3,2) - (R(3,1)*(R(1,2)*x(2) - R(2,2)*x(1)))/(R(1,1)*x(2) - R(1,2)*x(1))));
    X(1,1) = -(x(2)*(t(1) + X(2,1)*R(1,2)) - x(1)*(t(2) + X(2,1)*R(2,2)))/(R(1,1)*x(2) - R(1,2)*x(1));
    X(1,2) = -(x(2)*(t(1) + X(2,2)*R(1,2) + 2*R(1,3)) - x(1)*(t(2) + 2*R(1,3) + X(2,2)*R(2,2)))/(R(1,1)*x(2) - R(1,2)*x(1));
    X(3,1) = 0;
    X(3,2) = 2;
    return X;
}
NONLINSTATUS KHC(CONST EINTERVAL_VECTOR& XIV,
                    CONST EINTERVAL_VECTOR xIV,
                    CONST INTERVAL& fbar,
                    CONST EINTERVAL_VECTOR& EulAng,
                    CONST EINTERVAL_VECTOR t,
                    EINTERVAL_MATRIX& Kin){
    EINTERVAL a = EulAng(1);
    EINTERVAL b = EulAng(2);
    EINTERVAL g = EulAng(3);
    EINTERVAL A,B,C,D;
    INTERVAL f_bar = fbar;
    EINTERVAL_MATRIX K(Kin);
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)));
    B = (t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    C = xIV(2) + (XIV(3)*(K(2,3)*Cos(a)*Cos(b) - K(2,2)*Cos(b)*Sin(a)) - K(2,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)) - K(2,2)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - Cos(b)*Sin(a)*t(3)) + XIV(1)*(K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) + XIV(2)*(K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))))/(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    K(1,1) = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2));
    EINTERVAL Kinter;
    if(EIntervalIntersection(Kinter,K(1,1),Kin(1,1)) == 1){
        K(1,1) = Kinter;
        Kin(1,1) = Kinter;
    }else{
        return EMPTY;
    }
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)));
    D =  - (t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3));
    K(1,3) = (B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/D;
    if(EIntervalIntersection(Kinter,K(1,3),Kin(1,3)) == 1){
        K(1,3) = Kinter;
        Kin(1,3) = Kinter;
    }else{
        return EMPTY;
    }
    A = (XIV(3)*(K(2,3)*Cos(a)*Cos(b) - K(2,2)*Cos(b)*Sin(a)) - K(2,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)) + XIV(1)*(K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) + XIV(2)*(K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))));
    B = (t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    f_bar = 10*f_bar;
    C = xIV(1) + (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)))/(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    D =   (t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - Cos(b)*Sin(a)*t(3));//- (t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g))
    K(2,2) = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(2)-A)/D;
    if(EIntervalIntersection(Kinter,K(2,2),Kin(2,2)) == 1){
        K(2,2) = Kinter;
        Kin(2,2) = Kinter;
    }else{
        return EMPTY;
    }
    A = (XIV(3)*(K(2,3)*Cos(a)*Cos(b) - K(2,2)*Cos(b)*Sin(a))  - K(2,2)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - Cos(b)*Sin(a)*t(3)) + XIV(1)*(K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) + XIV(2)*(K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))));
    D = -(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))+ t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3));
    K(2,3) = (B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(2)-A)/D;
    if(EIntervalIntersection(Kinter,K(2,3),Kin(2,3)) == 1){
        K(2,3) = Kinter;
        Kin(2,3) = Kinter;
    }else{
        return EMPTY;
    }
    return DEFAULTOPT;
}
NONLINSTATUS tHC(CONST EINTERVAL_VECTOR& XIV,
                    CONST EINTERVAL_VECTOR xIV,
                    CONST INTERVAL& fbar,
                    CONST EINTERVAL_VECTOR& EulAng,
                    EINTERVAL_VECTOR tin,
                    CONST EINTERVAL_MATRIX& K){
    EINTERVAL a = EulAng(1);
    EINTERVAL b = EulAng(2);
    EINTERVAL g = EulAng(3);
    INTERVAL f_bar = fbar;
    EINTERVAL A,B,C,D;
    EINTERVAL_VECTOR t(3);
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)));
    B = (t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    C = xIV(2) + (XIV(3)*(K(2,3)*Cos(a)*Cos(b) - K(2,2)*Cos(b)*Sin(a)) - K(2,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)) - K(2,2)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - Cos(b)*Sin(a)*t(3)) + XIV(1)*(K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) + XIV(2)*(K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))))/(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    D = K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) ;
    t(1) = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/D;
    EINTERVAL tinter;
    if(EIntervalIntersection(tinter,t(1),tin(1)) == 1){
        tin(1) = tinter;
        t(1) = tinter;
    }else{
        return EMPTY;
    }
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + Cos(a)*Cos(b)*t(3)));
    D = K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g));
    t(2) = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/D;
    if(EIntervalIntersection(tinter,t(2),tin(2)) == 1){
        tin(2) = tinter;
        t(2) = tinter;
    }else{
        return EMPTY;
    }
    A = (XIV(3)*(K(2,3)*Cos(a)*Cos(b) - K(2,2)*Cos(b)*Sin(a)) - K(2,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))) - K(2,2)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - Cos(b)*Sin(a)*t(3)) + XIV(1)*(K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) + XIV(2)*(K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))));
    B = (t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    f_bar = 10*f_bar;
    C = xIV(1) + (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)))/(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    D = K(2,3)*Cos(a)*Cos(b);
    t(3) = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(2)-A)/D;
    if(EIntervalIntersection(tinter,t(3),tin(3)) == 1){
        tin(3) = tinter;
        t(3) = tinter;
    }else{
        return EMPTY;
    }
    return DEFAULTOPT;
}
NONLINSTATUS AngHC(CONST EINTERVAL_VECTOR& XIV,
                    CONST EINTERVAL_VECTOR xIV,
                    CONST INTERVAL& fbar,
                    EINTERVAL_VECTOR& EulAng,
                    CONST EINTERVAL_VECTOR t,
                    CONST EINTERVAL_MATRIX & K){
    EINTERVAL a = EulAng(1);
    EINTERVAL b = EulAng(2);
    EINTERVAL g = EulAng(3);
    INTERVAL f_bar = fbar;
    EINTERVAL A,B,C,D;
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)));
    B = (t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    C = xIV(2) + (XIV(3)*(K(2,3)*Cos(a)*Cos(b) - K(2,2)*Cos(b)*Sin(a)) - K(2,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)) - K(2,2)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - Cos(b)*Sin(a)*t(3)) + XIV(1)*(K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) + XIV(2)*(K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))))/(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - XIV(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - XIV(3)*Cos(a)*Cos(b) + Cos(a)*Cos(b)*t(3));
    D = K(1,1)*t(3);
    EINTERVAL T1;
    T1 = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/D;
    EINTERVAL Anginter;
    if(T1.maxabs() <= 1){
        b = ArcSin(T1);
        if(EIntervalIntersection(Anginter,b,EulAng(1)) == 1){
            b = Anginter;
            EulAng(1) = Anginter;
        }else{
            return EMPTY;
        }
    }
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) ) + XIV(3)*(K(1,1)*Sin(b) + K(1,3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)));
    D = -K(1,1)*Cos(b)*t(2);
    T1 = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/D;
    if(T1.maxabs() <= 1){
        g = ArcSin(T1);
        if(EIntervalIntersection(Anginter,g,EulAng(2)) == 1){
            g = Anginter;
            EulAng(2) = Anginter;
        }else{
            return EMPTY;
        }
    }
    A = (XIV(1)*(K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g)) + XIV(2)*(K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g)) - K(1,1)*(Sin(b)*t(3) + Cos(b)*Cos(g)*t(1) - Cos(b)*Sin(g)*t(2)) + XIV(3)*(K(1,1)*Sin(b)) - K(1,3)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + Cos(a)*Cos(b)*t(3)));
    D = XIV(3)*K(1,3)*Cos(b);
    T1 = -(B*Sqrt(Sqr(f_bar)-Sqr(C))-B*xIV(1)-A)/D;
    if(Abs(T1(1)) <= 1){
        a = ArcCos(T1);
        if(EIntervalIntersection(Anginter,a,EulAng(3)) == 1){
            a = Anginter;
            EulAng(3) = Anginter;
        }else{
           return EMPTY;
        }
    }
    return REDUCED;
}
INTERVAL xPXHC(CONST EINTERVAL_VECTOR& X,std::shared_ptr<VOID> userdata){
    assert(X.rows() == 10);
    assert(!optmicalibdata.Xw.empty());
    assert(!optmicalibdata.xc.empty());
    EINTERVAL_VECTOR EulAng(3),t(3);
    EINTERVAL_MATRIX K(3,3);
    EulAng(1) = X(1);
    EulAng(2) = X(2);
    EulAng(3) = X(3);
    t(1) = X(4);
    t(2) = X(5);
    t(3) = X(6);
    K(1,1) = X(7);
    K(2,2) = X(8);
    K(1,3) = X(9);
    K(2,3) = X(10);
    for(int i = 1;i <= optmicalibdata.Xw.cols();i++){
        if(KHC(optmicalibdata.Xw.col(i),optmicalibdata.xc.col(i),optmicalibdata.UpperBound,EulAng,t,K) == EMPTY){
            return INTERVAL(Machine::NaN);
        }
        if(tHC(optmicalibdata.Xw.col(i),optmicalibdata.xc.col(i),optmicalibdata.UpperBound,EulAng,t,K) == EMPTY){
            return INTERVAL(Machine::NaN);
        }
        if(AngHC(optmicalibdata.Xw.col(i),optmicalibdata.xc.col(i),optmicalibdata.UpperBound,EulAng,t,K) == EMPTY){
            return INTERVAL(Machine::NaN);
        }
    }
    Resize(optmicalibdata.X,10);
//    nonlinopt.X(1) = EulAng(1);
//    optmicalibdata.X(2) = EulAng(2);
//    optmicalibdata.X(3) = EulAng(3);
//    optmicalibdata.X(4) = t(1);
//    optmicalibdata.X(5) = t(2);
//    optmicalibdata.X(6) = t(3);
//    optmicalibdata.X(7) = K(1,1);
//    optmicalibdata.X(8) = K(1,3);
//    optmicalibdata.X(9) = K(2,2);
//    optmicalibdata.X(10) = K(2,3);
//    optmicalibdata.K = K;
    return PXFuncADII(optmicalibdata.X,userdata);
}

INTERVAL_MATRIX rtm2Eul2(CONST INTERVAL_MATRIX & R){
    INTERVAL_VECTOR Psy(2),Theta(2),Phi(2),tmpAng(3);
    INTERVAL_MATRIX EulAngles(3,2);
    INTERVAL_MATRIX r = R;
//    for(int i = 1;i <= 3;i++){
//        INTERVAL_VECTOR Lnorm;
//        Lnorm = R.row(i)/Mid(Norm(R.row(i)));
//        SetRow(r,i,Lnorm);
//    }
    Theta(1) = ArcSin(r(1,3));
    if(!(1 <= r(1,3)) && !(-1 <= r(1,3))){
        Theta(2) = Constant::Pi - Theta(1);
        if(0.0 <= Cos(Theta(1))){
            Psy(1) = INTERVAL(-Constant::Pi/2,Constant::Pi/2);
            Phi(1) = INTERVAL(-Constant::Pi/2,Constant::Pi/2);
        }else{
            Psy(1) = ArcTan2(-r(1,2)/Cos(Theta(1)),r(1,1)/Cos(Theta(1)));
            Phi(1) = ArcTan2(-r(2,3)/Cos(Theta(1)),r(3,3)/Cos(Theta(1)));
        }
        if(0.0 <= Cos(Theta(2))){
            Psy(2) = INTERVAL(-Constant::Pi/2,Constant::Pi/2);
            Phi(2) = INTERVAL(-Constant::Pi/2,Constant::Pi/2);
        }else{
            Psy(2) = ArcTan2(-r(1,2)/Cos(Theta(2)),r(1,1)/Cos(Theta(2)));
            Phi(2) = ArcTan2(-r(2,3)/Cos(Theta(2)),r(3,3)/Cos(Theta(2)));
        }
    }else{
        if(Diam(r(1,3)) == 0){
            Phi(1) = INTERVAL(0.0);
            if(-1 <= r(1,3)){
                Psy(1)=ArcTan2(-r(2,1),-r(3,1));
            }else{
                Psy(1)=ArcTan2(r(2,1),r(3,1));            
            }
        }else if(Diam(r(1,2)) == 0 && Inf(r(1,2)) == 0){
            Psy(1)=ArcTan2(r(2,1),r(3,1));
            cout << "Psy=" << Psy << endl;
            Phi(1) = ArcTan2(-Sup(r(2,3))/Cos(Sup(Theta(1))),Sup(r(3,3))/Cos(Sup(Theta(1))));
            cout << "Phi(1)=" << Phi(1) << endl;
        }else{
            Psy(1) = ArcTan2(-Inf(r(1,2))/Cos(Sup(Theta(1))),Sup(r(1,1))/Cos(Sup(Theta(1))));
            Psy(1) = Hull(Psy(1),Constant::Pi);
            cout << "Psy=" << Psy << endl;
            Phi(1) = ArcTan2(-Inf(r(2,3))/Cos(Sup(Theta(1))),Sup(r(3,3))/Cos(Sup(Theta(1))));
            Phi(1) = Hull(Phi(1),Constant::Pi);
            cout << "Phi(1)=" << Phi(1) << endl;
            EulRotmHC(Psy(1),Phi(1),r);
        }
    }
    tmpAng(3)=Psy(1);
    tmpAng(2)=Theta(1);
    tmpAng(1)=Phi(1);
    SetCol(EulAngles,1,tmpAng);
    tmpAng(3)=Psy(2);
    tmpAng(2)=Theta(2);
    tmpAng(1)=Phi(2);
    SetCol(EulAngles,2,tmpAng);
    cout << "EulAngles: " << EulAngles << endl;
    return EulAngles;
}

INTERVAL_VECTOR EulRotmHC(INTERVAL& psi,INTERVAL& phi,CONST INTERVAL_MATRIX &r){
    INTERVAL_VECTOR Psi(4),Phi(4);
    if(!(0.0 <= Sin(phi)) || !(0.0 <= Cos(phi))){
        Psi(1) = ArcSin((Cos(psi)*Sin(phi)+r(2,1))/(-r(1,3)*Cos(phi)));
        Psi(2) = ArcSin((Cos(psi)*Cos(phi)+r(2,2))/(-r(1,3)*Sin(phi)));
        Psi(3) = ArcSin((Sin(psi)*Sin(phi)+r(3,1))/(-r(1,3)*Cos(phi)));
        Psi(4) = ArcSin((Sin(psi)*Cos(phi)+r(3,2))/(-r(1,3)*Sin(phi)));
        cout << "Psi=" << Psi << endl;
        INTERVAL tmppsi;
        for(int i = 1;i <= Psi.nrows();i++){
            if(Intersection(tmppsi,psi,Psi(i)) == 1){
                psi = tmppsi;
            }
        }
    }else if(!(0.0 <= Sin(phi))){
        Psi(1) = ArcCos((r(1,3)*Sin(psi)*Cos(phi)+r(2,1))/Sin(phi));
        Psi(2) = ArcSin((Cos(psi)*Cos(phi)+r(2,2))/(-r(1,3)*Sin(phi)));
        Psi(3) = ArcSin((r(1,3)*Cos(psi)*Cos(phi)+r(3,1))/Sin(phi));
        Psi(4) = ArcSin((Sin(psi)*Cos(phi)+r(3,2))/(-r(1,3)*Sin(phi)));
    }else if(!(0.0 <= Cos(phi))){
        Psi(1) = ArcSin((Cos(psi)*Sin(phi)+r(2,1))/(-r(1,3)*Cos(phi)));
        Psi(2) = ArcCos((r(1,3)*Sin(psi)*Sin(phi)+r(2,2))/Cos(phi));
        Psi(3) = ArcSin((Sin(psi)*Sin(phi)+r(3,1))/(-r(1,3)*Cos(phi)));
        Psi(4) = ArcCos((r(1,3)*Cos(psi)*Sin(phi)+r(3,2))/Cos(phi));
    }
    if(!(0.0 <= Sin(psi)) || !(0.0 <= Cos(psi))){
        Phi(1) = ArcSin((Cos(psi)*Sin(phi)+r(2,1))/(-r(1,3)*Sin(psi)));
        Phi(2) = ArcSin((Cos(psi)*Cos(phi)+r(2,2))/(-r(1,3)*Sin(psi)));
        Phi(3) = ArcSin((Sin(psi)*Sin(phi)+r(3,1))/(-r(1,3)*Cos(psi)));
        Phi(4) = ArcSin((Sin(psi)*Cos(phi)+r(3,2))/(-r(1,3)*Cos(psi)));
        cout << "Phi=" << Phi << endl;
        INTERVAL tmpphi;
        for(int i = 1;i <= Phi.nrows();i++){
            if(Intersection(tmpphi,phi,Phi(i)) == 1){
                phi = tmpphi;
            }
        }
    }
    INTERVAL_VECTOR Angout(2);
    Angout(1) = psi;
    Angout(2) = phi;
}
INTERVAL grdHC(CONST EINTERVAL_VECTOR& p,std::shared_ptr<VOID> userdata){
    std::shared_ptr<INTERVAL_MATRIX> A_T1C = std::static_pointer_cast<INTERVAL_MATRIX>(userdata);    
    EINTERVAL_MATRIX K(3,3);
    EINTERVAL a,g,b;
    EINTERVAL_VECTOR t(3);
    EINTERVAL_MATRIX dp_dk1_1(3,4,{Cos(b)*Cos(g), -Cos(b)*Sin(g), Sin(b), t(2)*Cos(b)*Sin(g) - t(1)*Cos(b)*Cos(g) - t(3)*Sin(b),
                                0,              0,      0,                                               0,
                                0,              0,      0,                                               0});
    EINTERVAL_MATRIX dp_dk1_2(3,4,{Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b), Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g), -Cos(b)*Sin(a), t(3)*Cos(b)*Sin(a) - t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)),
                                    0,                                    0,              0,                                                                                                        0,
                                    0,                                    0,              0,                                                                                                        0});
    EINTERVAL_MATRIX dp_dk1_3(3,4,{ Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b), Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g), Cos(a)*Cos(b), - t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - t(3)*Cos(a)*Cos(b),
                                    0,                                    0,             0,                                                                                                          0,
                                    0,                                    0,             0,                                                                                                          0});
    EINTERVAL_MATRIX dp_dk2_2(3,4,{                                    0,                                    0,              0,                                                                                                        0,
                                    Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b), Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g), -Cos(b)*Sin(a), t(3)*Cos(b)*Sin(a) - t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)),
                                                                    0,                                    0,              0,                                                                                                        0});
    EINTERVAL_MATRIX dp_dk2_3(3,4,{                                  0,                                    0,             0,                                                                                                          0,
                                    Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b), Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g), Cos(a)*Cos(b), - t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - t(3)*Cos(a)*Cos(b),
                                                                    0,                                    0,             0,                                                                                                          0});
    EINTERVAL_MATRIX dp_da(3,4,{K(1,3)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) - K(1,2)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)), K(1,3)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - K(1,2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)), - K(1,2)*Cos(a)*Cos(b) - K(1,3)*Cos(b)*Sin(a), K(1,2)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + t(3)*Cos(a)*Cos(b)) - K(1,3)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(3)*Cos(b)*Sin(a)),
                               K(2,3)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) - K(2,2)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)), K(2,3)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - K(2,2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)), - K(2,2)*Cos(a)*Cos(b) - K(2,3)*Cos(b)*Sin(a), K(2,2)*(t(1)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + t(2)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) + t(3)*Cos(a)*Cos(b)) - K(2,3)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(3)*Cos(b)*Sin(a)),
                                               K(3,3)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)),                                               K(3,3)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)),                       -K(3,3)*Cos(b)*Sin(a),                                                                                                                  -K(3,3)*(t(1)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + t(2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(3)*Cos(b)*Sin(a))}); 
    EINTERVAL_MATRIX dp_db(3,4,{K(1,2)*Cos(b)*Cos(g)*Sin(a) - K(1,1)*Cos(g)*Sin(b) - K(1,3)*Cos(a)*Cos(b)*Cos(g), K(1,1)*Sin(b)*Sin(g) + K(1,3)*Cos(a)*Cos(b)*Sin(g) - K(1,2)*Cos(b)*Sin(a)*Sin(g), K(1,1)*Cos(b) - K(1,3)*Cos(a)*Sin(b) + K(1,2)*Sin(a)*Sin(b), K(1,3)*(t(3)*Cos(a)*Sin(b) + t(1)*Cos(a)*Cos(b)*Cos(g) - t(2)*Cos(a)*Cos(b)*Sin(g)) - K(1,1)*(t(3)*Cos(b) - t(1)*Cos(g)*Sin(b) + t(2)*Sin(b)*Sin(g)) - K(1,2)*(t(3)*Sin(a)*Sin(b) + t(1)*Cos(b)*Cos(g)*Sin(a) - t(2)*Cos(b)*Sin(a)*Sin(g)),
                                K(2,2)*Cos(b)*Cos(g)*Sin(a) - K(2,3)*Cos(a)*Cos(b)*Cos(g),                      K(2,3)*Cos(a)*Cos(b)*Sin(g) - K(2,2)*Cos(b)*Sin(a)*Sin(g),               K(2,2)*Sin(a)*Sin(b) - K(2,3)*Cos(a)*Sin(b),                                                          K(2,3)*(t(3)*Cos(a)*Sin(b) + t(1)*Cos(a)*Cos(b)*Cos(g) - t(2)*Cos(a)*Cos(b)*Sin(g)) - K(2,2)*(t(3)*Sin(a)*Sin(b) + t(1)*Cos(b)*Cos(g)*Sin(a) - t(2)*Cos(b)*Sin(a)*Sin(g)),
                                                 -K(3,3)*Cos(a)*Cos(b)*Cos(g),                                                  K(3,3)*Cos(a)*Cos(b)*Sin(g),                                   -K(3,3)*Cos(a)*Sin(b),                                                                                                                                        K(3,3)*(t(3)*Cos(a)*Sin(b) + t(1)*Cos(a)*Cos(b)*Cos(g) - t(2)*Cos(a)*Cos(b)*Sin(g))});
    EINTERVAL_MATRIX dp_dg(3,4,{K(1,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,1)*Cos(b)*Sin(g), - K(1,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) - K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - K(1,1)*Cos(b)*Cos(g), 0, K(1,1)*(t(2)*Cos(b)*Cos(g) + t(1)*Cos(b)*Sin(g)) - K(1,3)*(t(1)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - t(2)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))) - K(1,2)*(t(1)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b))),
                      K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) + K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)),                      - K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) - K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)), 0,                                            - K(2,2)*(t(1)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - t(2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b))) - K(2,3)*(t(1)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - t(2)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))),
                                                                    K(3,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)),                                                                     -K(3,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)), 0,                                                                                                                                            -K(3,3)*(t(1)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - t(2)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)))});
    EINTERVAL_MATRIX dp_dt1(3,4,{0, 0, 0, - K(1,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) - K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) - K(1,1)*Cos(b)*Cos(g),
                                0, 0, 0,                      - K(2,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) - K(2,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)),
                                0, 0, 0,                                                                     -K(3,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b))});
    EINTERVAL_MATRIX dp_dt2(3,4,{0, 0, 0, K(1,1)*Cos(b)*Sin(g) - K(1,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)) - K(1,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)),
                                0, 0, 0,                    - K(2,2)*(Cos(a)*Cos(g) - Sin(a)*Sin(b)*Sin(g)) - K(2,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g)),
                                0, 0, 0,                                                                   -K(3,3)*(Cos(g)*Sin(a) + Cos(a)*Sin(b)*Sin(g))});
    EINTERVAL_MATRIX dp_dt3(3,4,{0, 0, 0, K(1,2)*Cos(b)*Sin(a) - K(1,3)*Cos(a)*Cos(b) - K(1,1)*Sin(b),
                                0, 0, 0,               K(2,2)*Cos(b)*Sin(a) - K(2,3)*Cos(a)*Cos(b),
                                0, 0, 0,                                   -K(3,3)*Cos(a)*Cos(b)});
    EINTERVAL_MATRIX A((*A_T1C.get()));
    EINTERVAL_VECTOR S(3);
    S.reset();
    EINTERVAL M1(0);
    EINTERVAL_VECTOR R;
    R = A*p;
    S = A*dp_dk1_1.torowvector();
    for(int o = 1;o <= A.cols();o++){
        switch(o){
            case 1:
                for(int j = 1;j <= A.rows();j++){
                    for(int k = 1;k <= A.rows();k++){
                        if(k == j)continue;
                        M1 += R(k)*S(k);
                    }
                    EINTERVAL_VECTOR L;
                    L = A.row(j);
                    L(o) = 0;
                    M1 -= L*p*S(j);
                    M1 = M1/(A(j,o)*S(j));
                    EINTERVAL IK;
                    IK = M1 - (K(1,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)))/(Cos(b)*Cos(g));
                    if(EIntervalIntersection(IK,IK,K(1,1)) == 1)
                        K(1,1) = IK;
                    else
                        return INTERVAL(Machine::NaN);
                    IK = M1 - (K(1,3)*(Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)) + K(1,1)*Cos(b)*Cos(g))/((Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)));
                    if(EIntervalIntersection(IK,IK,K(1,1)) == 1)
                        K(1,2) = IK;
                    else
                        return INTERVAL(Machine::NaN);
                    IK = M1 - (K(1,2)*(Cos(a)*Sin(g) + Cos(g)*Sin(a)*Sin(b)) + K(1,1)*Cos(b)*Cos(g))/((Sin(a)*Sin(g) - Cos(a)*Cos(g)*Sin(b)));
                    if(EIntervalIntersection(IK,IK,K(1,1)) == 1)
                        K(1,3) = IK;
                    else
                        return INTERVAL(Machine::NaN);
                }
                break;
        }
    }
}
EINTERVAL More_Garbow_Hillstrom_func(CONST EINTERVAL_VECTOR & X, std::shared_ptr<VOID> userdata){
    return Sqr(1.5 - X(1)*(1 - X(2))) + Sqr(2.25 - X(1)*(1 - Sqr(X(2)))) + Sqr(2.625 - X(1)*(1 - Power(X(2),3)));
}
INTERVAL_AUTODIFF More_Garbow_Hillstrom_AD(CONST INTERVAL_AUTODIFF& X, std::shared_ptr<VOID> userdata){
    return Sqr(1.5 - X(1)*(1 - X(2))) + Sqr(2.25 - X(1)*(1 - Sqr(X(2)))) + Sqr(2.625 - X(1)*(1 - X(2)*Sqr(X(2))));
}
EINTERVAL More_Garbow_Hillstrom_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    EINTERVAL tR;
    EINTERVAL A,B,x,y;
    x = X(1);y = X(2);
    A =  Sqr(2.25 - x*(1 - Sqr(y))) + Sqr(2.625 - x*(1 - Power(y,3)));
    B = Sqr(1.5) + Sqr(x) + Sqr(x*y) + 3*x*y - 2*Sqr(x)*y+A;
    tR = -(f_bar - B)/3;
    if(EIntervalIntersection(tR,tR,x) == 1){
        x = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    B = Sqr(1.5)/(3*x) + x/3 + x*Sqr(y)/3 - 1 - 2*x*y/3 + A/3*x;
    tR = f_bar/(3*x)-B;
    if(EIntervalIntersection(tR,tR,y) == 1){
        y = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    X(1) = x;
    X(2) = y;
    return More_Garbow_Hillstrom_func(X,NULL);
}
EINTERVAL_VECTOR More_Garbow_Hillstrom_grd(CONST EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    EINTERVAL_VECTOR grd(2);
    grd(1) = 2*(Sqr(X(2)) - 1)*(X(1)*(Sqr(X(2)) - 1) + 2.25) + 2*(Power(X(2),3) - 1)*
            (X(1)*(Power(X(2),3) - 1) + 2.625) + 2*(X(1)*(X(2) - 1) + 1.5)*(X(2) - 1);
    grd(2) = 2*X(1)*(X(1)*(X(2) - 1) + 1.5) + 6*X(1)*Sqr(X(2))*(X(1)*(Power(X(2),3) - 1) + 2.625) + 
            4*X(1)*X(2)*(X(1)*(Sqr(X(2)) - 1) + 2.25);
    return grd;
}
//d2f_dx2 = 2*(y^2 - 1)*(x*(y^2 - 1) + 9/4) + 2*(y^3 - 1)*(x*(y^3 - 1) + 21/8) + 2*(x*(y - 1) + 3/2)*(y - 1)
//d2f_dy2 = 2*x*(x*(y - 1) + 3/2) + 6*x*y^2*(x*(y^3 - 1) + 21/8) + 4*x*y*(x*(y^2 - 1) + 9/4)
INTERVAL_AUTODIFF More_Garbow_Hillstrom_grd_AD(CONST INTERVAL_AUTODIFF & X,INT i,std::shared_ptr<VOID> userdata){
    EINTERVAL_VECTOR grd(2);
    if(i == 1)
        return 2*(Sqr(X(2)) - 1)*(X(1)*(Sqr(X(2)) - 1) + 2.25) + 2*(X(2)*Sqr(X(2)) - 1)*
            (X(1)*(X(2)*Sqr(X(2)) - 1) + 2.625) + 2*(X(1)*(X(2) - 1) + 1.5)*(X(2) - 1);
    if(i == 2)
        return 2*X(1)*(X(1)*(X(2) - 1) + 1.5) + 6*X(1)*Sqr(X(2))*(X(1)*(X(2)*Sqr(X(2)) - 1) + 2.625) + 
            4*X(1)*X(2)*(X(1)*(Sqr(X(2)) - 1) + 2.25);
}
INTERVAL_AUTODIFF More_Garbow_Hillstrom_Bgrd_AD(CONST INTERVAL_AUTODIFF & X,INT i,std::shared_ptr<VOID> userdata){
    EINTERVAL_VECTOR grd(2);
    std::shared_ptr<MATRIX> B = std::static_pointer_cast<MATRIX>(userdata);
    if(i == 1){
        INTERVAL_AUTODIFF B11(2,X.ComputeHessian),B12(2,X.ComputeHessian);
        B11.fkt() = B.get()[0](1,1);
        Initialize(B11.grd(),0.0);
        B12.fkt() = B.get()[0](1,2);
        Initialize(B12.grd(),0.0);
        return B11*(2*(Sqr(X(2)) - 1)*(X(1)*(Sqr(X(2)) - 1) + 2.25) + 2*(Power(X(2),3) - 1)*(X(1)*(Power(X(2),3) - 1) + 2.625) + 2*(X(1)*(X(2) - 1) + 1.5)*(X(2) - 1)) + 
                B12*(2*X(1)*(X(1)*(X(2) - 1) + 1.5) + 6*X(1)*Sqr(X(2))*(X(1)*(Power(X(2),3) - 1) + 2.625) + 4*X(1)*X(2)*(X(1)*(Sqr(X(2)) - 1) + 2.25));
    }
    if(i == 2){
        INTERVAL_AUTODIFF B21(2,X.ComputeHessian),B22(2,X.ComputeHessian);
        B21.fkt() = B.get()[0](2,1);
        Initialize(B21.grd(),0);
        B22.fkt() = B.get()[0](2,2);
        Initialize(B22.grd(),0);
        return B21(1)*(2*(Sqr(X(2)) - 1)*(X(1)*(Sqr(X(2)) - 1) + 2.25) + 2*(Power(X(2),3) - 1)*(X(1)*(Power(X(2),3) - 1) + 2.625) + 2*(X(1)*(X(2) - 1) + 1.5)*(X(2) - 1)) + 
                B22(1)*(2*X(1)*(X(1)*(X(2) - 1) + 1.5) + 6*X(1)*Sqr(X(2))*(X(1)*(Power(X(2),3) - 1) + 2.625) + 4*X(1)*X(2)*(X(1)*(Sqr(X(2)) - 1) + 2.25));
    }
}

EINTERVAL More_Garbow_Hillstrom_grd_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    if(userdata)
        return More_Garbow_Hillstrom_Bgrd_HC(f_bar,X,userdata);
    EINTERVAL A;
    A = 6*X(1)*Sqr(X(2))*(X(1)*(Power(X(2),3) - 1) + 21.0/8.0) + 
            4*X(1)*X(2)*(X(1)*(Sqr(X(2)) - 1) + 9.0/4.0);
    EINTERVAL tR;
    tR = -(2.0/3.0*Sqr(X(1))*(X(2)-1)+A/3);
    if(EIntervalIntersection(tR,tR,X(1)) == 1){
        X(1) = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    tR = -(-1+1.5/X(1)+A/(2*X(1)));
    if(EIntervalIntersection(tR,tR,X(2)) == 1){
        X(2) = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    A = 2*(Power(X(2),3) - 1)* (X(1)*(Power(X(2),3) - 1) + 21.0/8.0) + 2*(X(1)*(X(2) - 1) + 3.0/2.0)*(X(2) - 1);
    tR = -(Sqr(X(2))*(X(1)*(Sqr(X(2))-2)+9.0/4.0)-9.0/4.0+A/2);
    if(EIntervalIntersection(tR,tR,X(1)) == 1){
        X(1) = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    tR = Sqrt(-(X(1)*(4.0/9.0*Sqr(X(2))*(Sqr(X(2))-2)+4.0/9.0)-1+4.0/18.0*A));
    if(EIntervalIntersection(tR,tR,X(2)) == 1){
        X(2) = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    return More_Garbow_Hillstrom_func(X,NULL);
}
//B1_1*(2*(y^2 - 1)*(x*(y^2 - 1) + 9/4) + 2*(y^3 - 1)*(x*(y^3 - 1) + 21/8) + 2*(x*(y - 1) + 3/2)*(y - 1)) + B1_2*(2*x*(x*(y - 1) + 3/2) + 6*x*y^2*(x*(y^3 - 1) + 21/8) + 4*x*y*(x*(y^2 - 1) + 9/4))
//B2_1*(2*(y^2 - 1)*(x*(y^2 - 1) + 9/4) + 2*(y^3 - 1)*(x*(y^3 - 1) + 21/8) + 2*(x*(y - 1) + 3/2)*(y - 1)) + B2_2*(2*x*(x*(y - 1) + 3/2) + 6*x*y^2*(x*(y^3 - 1) + 21/8) + 4*x*y*(x*(y^2 - 1) + 9/4))
EINTERVAL More_Garbow_Hillstrom_Bgrd_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    EINTERVAL T1,x,y;
    x = X(1);
    y = X(2);
    std::shared_ptr<INTERVAL_MATRIX> B = std::static_pointer_cast<INTERVAL_MATRIX>(userdata);
    T1 = B.get()[0](1,2)*(x*(x*(y - 1) + 3/2) + 3*x*Sqr(y)*(x*(Power(y,3) - 1) + 2.625) + x*y*(x*(Sqr(y) - 1) + 2.25));
    EINTERVAL T2;
    T2 = (Sqr(y) - 1)*(x*(Sqr(y) - 1) + 2.25) + (Power(y,3) - 1)*(x*(Power(y,3) - 1) + 2.625);
    EINTERVAL tR;
    tR = -(x*y*(y-2)+1.5*(y-1)+(T1/B.get()[0](1,1))+T2);
    if(EIntervalIntersection(x,tR,x) == 0){
        return EINTERVAL(Machine::NaN);
    }
    tR = -(x*(y*(y-2)+1)-1.5+(T1/B.get()[0](1,1))+T2);
    if(EIntervalIntersection(y,tR,y) == 1){
        return EINTERVAL(Machine::NaN);
    }
    T1 = B.get()[0](2,2)*(3*x*Sqr(y)*(x*(Power(y,3) - 1) + 2.625) + 2*x*y*(x*(Sqr(y) - 1) + 2.25));
    tR = -(x*y*(y-2)+1.5*(y-1)+(T1/B.get()[0](2,1))+T2);
    if(EIntervalIntersection(x,tR,x) == 1){
        X(1) = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    tR = -(x*(y*(y-2)+1)-1.5+(T1/B.get()[0](2,1))+T2);
    if(EIntervalIntersection(y,tR,y) == 1){
        return EINTERVAL(Machine::NaN);
    }
    X(1) = x;
    X(2) = y;
    return More_Garbow_Hillstrom_func(X,NULL);
}
EINTERVAL More_Garbow_Hillstrom_hess_HC(INTERVAL f_bar,EINTERVAL_VECTOR & X,std::shared_ptr<VOID> userdata){
    EINTERVAL A;
    EINTERVAL x,y;
    x = X(1);y = X(2);
    A = Sqr(Sqr(y)-1)+Sqr(Power(y,3)-1);
    EINTERVAL  tR;
    tR = -(EINTERVAL(0,Machine::PosInfinity)-(A+Sqr(y)+1));
    if(EIntervalIntersection(tR,tR,y) == 1){
        y = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    A = 2*Sqr(x*y) + 4.5*Sqr(x)*Power(y,4) + 0.5*Sqr(x) + 3*x*y*(x*(Power(y,3) - 1) + 2.625);
    tR = -(A/(x*(Sqr(y) - 1) + 2.25))+EINTERVAL(0,Machine::PosInfinity)/1.5;
    if(EIntervalIntersection(tR,tR,x) == 1){
        x = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    A = 2*Sqr(x*y) + 4.5*Sqr(x)*Power(y,4) + x*(x*(Sqr(y) - 1) + 2.25) + 0.5*Sqr(x);
    tR = -(A/(3*x*(x*(Power(y,3) - 1) + 2.625)))+EINTERVAL(0,Machine::PosInfinity);
    if(EIntervalIntersection(tR,tR,y) == 1){
        y = tR;
    }else{
        return EINTERVAL(Machine::NaN);
    }
    X(1) = x;
    X(2) = y;
    return More_Garbow_Hillstrom_func(X,NULL);
}
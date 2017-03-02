 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Icalib.cpp
 * Author: darlan
 * 
 * Created on 24 de Maio de 2016, 09:18
 */

#include <Icalib.h>
#include <FileOp.h>
#include <float.h>
#include <limits.h>

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
struct OPTMICALIBDATA  optmicalibdata;

Icalib::Icalib() {

}

Icalib::Icalib(const Icalib& orig) {
}
std::unique_ptr<FUNCTION2> IGAObjFunc(new FUNCTION2);

float
GADotObjFunc(GAGenome & c){
    GABin2DecGenome & genome = (GABin2DecGenome &)c;
    AbLinsys *Albl = (AbLinsys *)(c.userData());
    INTERVAL_VECTOR xin = Albl->X;
    float Vf = 0;
    INTERVAL_VECTOR x(xin.nrows());
    INTERVAL IVf=0;
    for(int i = 1;i <= x.nrows();i++){
        x(i) = xin(i)+SymHull(genome.phenotype(i-1));
        Vf -= Diam(x(i));
    }
    IVf = x.Box(1,3)*x.Box(4,6);
    if(std::isnan(Inf(IVf)) || std::isnan(Sup(IVf)))Vf = -FLT_MAX;
    Vf -= (Diam(IVf)+Abs(IVf)+Inf(IVf));
    return Vf;
}
float
GAxPXPObjFunc(GAGenome & c)
{
    GABin2DecGenome & genome = (GABin2DecGenome &)c;
    AbLinsys *Albl = (AbLinsys *)(c.userData());
    INTERVAL IVf;
    INTERVAL_VECTOR x(Albl->X.nrows());
    for(int i = 1;i <= x.nrows();i++){
        x(i) = Albl->X(i)+SymHull(genome.phenotype(i-1));
    }
    MATRIX Delta,Amais,Amenos,Ac;
    Ac = Mid(Albl->A);
    Delta = Sup(Albl->A)-Ac;
    for(int iz = 1;iz <= Albl->z.cols();iz++){
        MATRIX Tz;
        Tz = Diag(Albl->z.col(iz));
        Amais = Ac-Delta*Tz;
        Amenos = -(Ac+Delta*Tz);
        IVf += Norm(Amais*x-Sup(Albl->b));
        IVf += Norm(Amenos*x-Inf(Albl->b));
    }
    if(IVf.inf() < 1e-30)
        IVf.ival.inf = 0;
    if(IVf.inf() > 1e30)
        IVf.ival.inf = 1e30;
    if(IVf.sup() < 1e-30)
        IVf.ival.inf = 0;
    if(IVf.sup() > 1e30)
        IVf.ival.inf = 1e30;
    return -(abs(IVf.inf())+1e3*abs(IVf.sup()));
}
float
GAftzexpObjFunc(GAGenome & c)
{
    GABin2DecGenome & genome = (GABin2DecGenome &)c;
    AbLinsys *Albl = (AbLinsys *)(c.userData());
    INTERVAL IVf;
    INTERVAL_VECTOR x(Albl->X.nrows());
    for(int i = 1;i <= x.nrows();i++){
        x(i) = Albl->X(i)+SymHull(genome.phenotype(i-1));
    }
    MATRIX Delta,Amais,Amenos,Ac;
    Ac = Mid(Albl->A);
    Delta = Sup(Albl->A)-Ac;
    for(int iz = 1;iz <= Albl->z.cols();iz++){
        MATRIX Tz;
        Tz = Diag(Albl->z.col(iz));
        Amais = Ac-Delta*Tz;
        Amenos = -(Ac+Delta*Tz);
        IVf += Norm(Amais*x-Sup(Albl->b));
        IVf += Norm(Amenos*x-Inf(Albl->b));
        if(std::isnan(IVf.inf())){
            IVf.ival.inf = 0;
        }
    }
    return -(abs(IVf.inf())+1e3*abs(IVf.sup()));
}
float
GAIniTsaiObjFunc(GAGenome & c)
{
    GABin2DecGenome & genome = (GABin2DecGenome &)c;
    AbLinsys *Albl = (AbLinsys *)(c.userData());
    INTERVAL IVf(0);
    INTERVAL_VECTOR x(Albl->X.nrows());
    for(int i = 1;i <= x.nrows();i++){
        x(i) = Mid(Albl->X(i))+SymHull(genome.phenotype(i-1));
    }
    MATRIX Delta,Amais,Amenos,Ac;
    Ac = Mid(Albl->A);
    Delta = Sup(Albl->A)-Ac;
    for(int iz = 1;iz <= Albl->z.cols();iz++){
        MATRIX Tz;
        Tz = Diag(Albl->z.col(iz));
        Amais = Ac-Delta*Tz;
        Amenos = -(Ac+Delta*Tz);
        IVf += Norm(Amais*x-Sup(Albl->b));
        IVf += Norm(Amenos*x-Inf(Albl->b));
    }
    return -(abs(IVf.inf())+1e3*abs(IVf.sup()));
}
float
GAftzObjFunc(GAGenome & c)
{
    GABin2DecGenome & genome = (GABin2DecGenome &)c;
    AbLinsys *Albl = (AbLinsys *)(c.userData());
    INTERVAL IVf;
    INTERVAL_VECTOR x(Albl->X.nrows());
    for(int i = 1;i <= x.nrows();i++){
        x(i) = Mid(Albl->X(i))+SymHull(genome.phenotype(i-1));
    }
    MATRIX Delta,Amais,Amenos,Ac;
    Ac = Mid(Albl->A);
    Delta = Sup(Albl->A)-Ac;
    for(int iz = 1;iz <= Albl->z.cols();iz++){
        MATRIX Tz;
        Tz = Diag(Albl->z.col(iz));
        Amais = Ac-Delta*Tz;
        Amenos = -(Ac+Delta*Tz);
        IVf += Norm(Amais*x-Sup(Albl->b));
        IVf += Norm(Amenos*x-Inf(Albl->b));
    }
    return -(abs(IVf.inf())+1e3*abs(IVf.sup()));
}

Icalib::~Icalib() {
    calccameradata.clear();
    truecameradata.clear();
    targets.clear();
}
VOID Icalib::erase(){
    calccameradata.clear();
    truecameradata.clear();
    targets.clear();    
}
VOID Icalib::cleanpointdata(){
    for(int contCam = 0;contCam < Ncamera;contCam++){
        calccameradata[contCam].t.clear();
        calccameradata[contCam].VEMxc.clear();
        calccameradata[contCam].VEMXw.clear();
        calccameradata[contCam].PtsIdx.clear();
        calccameradata[contCam].vIXw.clear();
        calccameradata[contCam].vxci.clear();
    }
}
VOID Icalib::LoadCameraMatrix(const std::vector<CAMERADATA>& camdata){
    Eigen::Matrix3d K;
    CAMERA tmpCamera;
    Ncamera = camdata.size();
//    ofstream csvcamdata("camdata.csv");
    for(int contCam = 0;contCam < camdata.size();contCam++){
        //csvcamdata << contCam << endl;
        tmpCamera.EMR = Eigen::AngleAxisd(-(camdata[contCam].rotX(3,0)/180)*M_PI, Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-(camdata[contCam].rotY(3,0)/180)*M_PI,  Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd(-(camdata[contCam].rotZ(3,0)/180)*M_PI, Eigen::Vector3d::UnitZ());
        tmpCamera.TrueEulAng << -(camdata[contCam].rotX(3,0)/180)*M_PI,-(camdata[contCam].rotY(3,0)/180)*M_PI,
                -(camdata[contCam].rotZ(3,0)/180)*M_PI;
//        csvcamdata << -(camdata[contCam].rotX(3,0)/180)*M_PI << "," << -(camdata[contCam].rotY(3,0)/180)*M_PI << ","
//                   << -(camdata[contCam].rotZ(3,0)/180)*M_PI;
        tmpCamera.EVt = camdata[contCam].local;
//        csvcamdata << "," << camdata[contCam].local(0) << "," << camdata[contCam].local(1) << "," << camdata[contCam].local(2);
        float fm_x,fm_y;
        fm_x = camdata[contCam].f*camdata[contCam].m_x;
        fm_y = fm_x;
        K  << fm_x, 0 , camdata[contCam].ImageSize(0)/2,
              0, fm_y, camdata[contCam].ImageSize(1)/2,
              0, 0, 1;
//        csvcamdata << "," << K(1,1) << "," << camdata[contCam].ImageSize(0)/2 << endl;
        tmpCamera.EMK = K;
        Resize(tmpCamera.K,3,3);
        tmpCamera.K = MatrixXd2MATRIX(K);
        tmpCamera.ImageSize = camdata[contCam].ImageSize;
        tmpCamera.ID = camdata[contCam].ID;
        tmpCamera.calibrated = false;
        Eigen::Matrix3Xd RRt;
        RRt.resize(Eigen::NoChange,4);
        RRt.block(0,0,3,3) = tmpCamera.EMR;
        RRt.col(3) = -tmpCamera.EMR*tmpCamera.EVt;
        tmpCamera.RRt = MatrixXd2MATRIX(RRt);
        truecameradata.push_back(tmpCamera);

        tmpCamera.EulAngl = INTERVAL_VECTOR(3,{INTERVAL(-(camdata[contCam].rotX(3,0)/180)*M_PI),
                                               INTERVAL(-(camdata[contCam].rotY(3,0)/180)*M_PI),
                                               INTERVAL(-(camdata[contCam].rotZ(3,0)/180)*M_PI)});
        tmpCamera.R = Eul2rtm(tmpCamera.EulAngl);
        tmpCamera.If = INTERVAL(-0.3,-0.016);
        tmpCamera.Im_x = tmpCamera.ImageSize(0)/INTERVAL(0.04,0.015);
        tmpCamera.Im_y = tmpCamera.ImageSize(0)/INTERVAL(0.04,0.015);
        INTERVAL Ifm_x = tmpCamera.If*tmpCamera.Im_x;
        INTERVAL Ifm_y = tmpCamera.If*tmpCamera.Im_y;
        INTERVAL x0,y0;
        x0 = tmpCamera.ImageSize(0)/2;
        y0 = tmpCamera.ImageSize(1)/2;
        tmpCamera.K = INTERVAL_MATRIX(3,3,{Ifm_x,0,   x0,
                                            0   ,Ifm_y,y0,
                                            0   ,0   ,1});
        Resize(tmpCamera.RRt,3,4);
        Initialize(tmpCamera.RRt,INTERVAL(-40,40));
        tmpCamera.RRt(2,4)=INTERVAL(1,1);
        calccameradata.push_back(tmpCamera);
    }
//    csvcamdata.close();
}

VOID Icalib::LoadTargetMatrix(CONST std::vector<FORMDATA>& formdata){
    TARGET tmpTarget;
    Ntarget = formdata.size();
    for(int countTarg = 0;countTarg < formdata.size();countTarg++){
        tmpTarget.RelativePos.resize(Eigen::NoChange,formdata[countTarg].vertex.cols());
        Eigen::MatrixXd AbsolutePosH,AbsolutePosNH,AbsoluteNormH,AbsoluteNormNH;
        AbsolutePosH = formdata[countTarg].TransMatrix*formdata[countTarg].vertex.colwise().homogeneous();
        AbsolutePosNH = AbsolutePosH.colwise().hnormalized();
        AbsoluteNormH = formdata[countTarg].TransMatrix*formdata[countTarg].normals.colwise().homogeneous();
        AbsoluteNormNH = AbsoluteNormH.colwise().hnormalized();
        tmpTarget.Npoints = AbsolutePosNH.cols();
        tmpTarget.RelativePos = AbsolutePosNH;
        tmpTarget.ID = formdata[countTarg].ID;
        tmpTarget.AbsolutePos = formdata[countTarg].vertex;
        tmpTarget.local = formdata[countTarg].local;
        targets.push_back(tmpTarget);
    }
}

VOID Icalib::PointProjection(REAL imgdiam,REAL imgruido,bool savedata){
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,imgruido);
    YAML::Node elementnode,groupnode;
    WriteNode(elementnode,imgdiam);
    groupnode["Diameter"].push_back(elementnode);
    elementnode.reset();
    WriteNode(elementnode,imgruido);
    groupnode["Noise"].push_back(elementnode);
    elementnode.reset();
    Eigen::Matrix3Xd tmpx;
    Eigen::Matrix2Xd x,tmpruidox,tmpdiam;
    Eigen::Matrix3d K,R;
    Eigen::Vector3d t;
    ofstream csvpmatrix(filenames.prjmatrx.c_str());
    ofstream csvpixeldata(filenames.csvpxlocation);
    for(int contCam = 0;contCam < truecameradata.size();contCam++){
        calccameradata[contCam].VEMxc.clear();
        truecameradata[contCam].VEMxc.clear();
        truecameradata[contCam].vxci.clear();
        calccameradata[contCam].vxci.clear();
        truecameradata[contCam].VEMXw.clear();
        calccameradata[contCam].VEMXw.clear();
        K = truecameradata[contCam].EMK;
        R = truecameradata[contCam].EMR;
        t = truecameradata[contCam].EVt;
        Eigen::Matrix<double,3,4> Rt;
        Rt.block(0,0,3,3) = R;
        Rt.col(3) = -R*t;
        Resize(truecameradata[contCam].RRt,3,4);
        truecameradata[contCam].RRt = MatrixXd2MATRIX(Rt);
        Initialize(calccameradata[contCam].RRt,INTERVAL(-40,40));
        calccameradata[contCam].RRt(2,4) = 1;
        truecameradata[contCam].TargIdx.resize(targets.size());
        calccameradata[contCam].TargIdx.resize(targets.size());
        truecameradata[contCam].EMP = K*Rt;
        Resize(truecameradata[contCam].P,3,4);
        truecameradata[contCam].P = MatrixXd2MATRIX(truecameradata[contCam].EMP);
        int contTrgPrj = 0;
        INTERVAL_MATRIX IMx;
        for(int contTrg = 0;contTrg < targets.size();contTrg++){
            tmpx = truecameradata[contCam].EMP*targets[contTrg].RelativePos.colwise().homogeneous();
            x = tmpx.colwise().hnormalized();
            if(x.row(0).minCoeff() < 0 || x.row(0).maxCoeff() > truecameradata[contCam].ImageSize(0))continue;
            if(x.row(1).minCoeff() < 0 || x.row(1).maxCoeff() > truecameradata[contCam].ImageSize(1))continue;
            tmpruidox.resize(Eigen::NoChange,tmpx.cols());
            tmpdiam.resize(Eigen::NoChange,tmpx.cols());
            generator.seed(time(NULL));
            for(int i = 0;i < tmpruidox.rows();i++)
                for(int j = 0;j < tmpruidox.cols();j++){
                    tmpruidox(i,j)=distribution(generator);
                }
            for(int i = 0;i < x.cols();i++)
                csvpixeldata << contTrg << "," << i << "," << x(0,i)+tmpruidox(0,i) << "," << 
                        x(1,i)+tmpruidox(1,i) << "," << imgdiam << imgruido << endl;
            tmpdiam.setConstant(imgdiam);
            calccameradata[contCam].VEMxc.push_back(x);
            truecameradata[contCam].VEMxc.push_back(x);
            truecameradata[contCam].xc = Hull(MatrixXd2MATRIX(x));
            Eigen::Matrix2Xd xil;
            xil.resize(Eigen::NoChange,x.cols());
            xil.setOnes();
            xil.row(0) *= truecameradata[contCam].ImageSize(0)/2;
            xil.row(1) *= truecameradata[contCam].ImageSize(1)/2;
            xil = x - xil;
            IMx = MatrixXd2MATRIX(xil);
            truecameradata[contCam].vxci.push_back(IMx);
            IMx = IMx + SymHull(MatrixXd2MATRIX(tmpdiam)) + MatrixXd2MATRIX(tmpruidox);
            WriteNode(elementnode,IMx);
            groupnode[calccameradata[contCam].ID][targets[contTrg].ID].push_back(elementnode);
            elementnode.reset();
            calccameradata[contCam].vxci.push_back(IMx);
            truecameradata[contCam].VEMXw.push_back(targets[contTrg].RelativePos);
            calccameradata[contCam].VEMXw.push_back(targets[contTrg].RelativePos);
            truecameradata[contCam].TargIdx(contTrgPrj) = contTrg;
            calccameradata[contCam].TargIdx(contTrgPrj) = contTrg;
            contTrgPrj++;
        }
        truecameradata[contCam].TargIdx.conservativeResize(Eigen::NoChange,contTrgPrj);
        calccameradata[contCam].TargIdx.conservativeResize(Eigen::NoChange,contTrgPrj);
        VECTOR TargIdx(calccameradata[contCam].TargIdx.rows());
        for(int i = 1;i <= Dimension(TargIdx);i++)
            TargIdx(i) = calccameradata[contCam].TargIdx(i-1);
        WriteNode(elementnode,TargIdx);
        groupnode[calccameradata[contCam].ID]["PrjtedTarg"].push_back(elementnode);
        elementnode.reset();
    }
    if(savedata){
        FileOp CameraData(filenames.pxlocation,ios::in);
        CameraData.node.push_back(groupnode);
        FileOp SaveData(filenames.pxlocation,ios::out);
        SaveData.node = CameraData.node;
    }
    csvpixeldata.close();
    csvpmatrix.close();
}

VOID Icalib::XYZINTAbsolLoad(const string& filename){
    FileOp XYZdataLoad(filename,ios::in);
    INTERVAL_MATRIX IXw;
    if(XYZdataLoad.IsOpen()){
        if(!XYZdataLoad.node.IsNull()){
            for(int contTar = 0;contTar < targets.size();contTar++){
                ReadNode(XYZdataLoad.node[targets[contTar].ID][0],0,IXw);
                targets[contTar].IMAbsolutePos = IXw;
            }
            return;
        }
        XYZdataLoad.SaveFile.close();
        XYZdataLoad.node.reset();
    }
    FileOp XYZdataSave(filename,ios::out);
    INTERVAL_VECTOR Xw(4);
    Resize(IXw,4,targets[0].Npoints);
    for(int contPts = 0;contPts < targets[0].Npoints;contPts++){
        Xw(4) = 1;
        if(targets[0].AbsolutePos(2,contPts) < 2){
            Xw(3) = 0;
        }else{
            Xw(3) = INTERVAL(1,2.1);
        }
        if(targets[0].AbsolutePos(0,contPts) < 2){
            Xw(1) = 0;
        }else{
            Xw(1) = INTERVAL(0.45,0.9);
        }
        if(targets[0].AbsolutePos(1,contPts) < 2){
            Xw(2) = 0;
        }else{
            Xw(2) = INTERVAL(0.3,0.45);
        }
        SetCol(IXw,contPts+1,Xw);
    }
    for(int contTar = 0;contTar < targets.size();contTar++){
        YAML::Node elementnode,tmpnode;
        WriteNode(elementnode,IXw);
        tmpnode[targets[contTar].ID].push_back(elementnode);
        XYZdataSave.node["Absolute"].push_back(tmpnode);
        targets[contTar].IMAbsolutePos = IXw;
    }
}
VOID Icalib::XYZINTRelatLoad(CONST double diameter,CONST double noise,bool savedata){
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,noise);
    YAML::Node elementnode,groupnode;
    WriteNode(elementnode,diameter);
    groupnode["Diameter"].push_back(elementnode);
    elementnode.reset();
    WriteNode(elementnode,noise);
    groupnode["Noise"].push_back(elementnode);
    elementnode.reset();
    ofstream csvfile(filenames.csvptlocation);
    for(int contTar = 0;contTar < targets.size();contTar++){
        Eigen::MatrixXd EMnoise;
        EMnoise.resize(targets[contTar].RelativePos.rows(),targets[contTar].RelativePos.cols());
        generator.seed(time(NULL));
        for(int i = 0;i < EMnoise.rows();i++){
            for(int j = 0;j < EMnoise.cols();j++){
                EMnoise(i,j)=distribution(generator);
            }
        }
        INTERVAL_MATRIX MDiam(targets[contTar].RelativePos.rows(),targets[contTar].RelativePos.cols());
        Initialize(MDiam,SymHull(diameter));
        MATRIX Datapoint;
        Resize(targets[contTar].IMRelativePos,targets[contTar].RelativePos.rows(),targets[contTar].RelativePos.cols());
        targets[contTar].IMRelativePos = MatrixXd2MATRIX(targets[contTar].RelativePos);
        targets[contTar].IMRelativePos = targets[contTar].IMRelativePos + MDiam + MatrixXd2MATRIX(EMnoise);
        Datapoint = Mid(targets[contTar].IMRelativePos) + MatrixXd2MATRIX(EMnoise);
        for(int i = 1;i <= targets[contTar].IMRelativePos.cols();i++)
            csvfile << contTar << "," << i << "," << noise << "," << Datapoint(1,i) << ","
                                                    << Datapoint(2,i) << ","
                                                    << Datapoint(3,i) << endl;
        if(contTar == 0){
            for(int contPts = 1;contPts <= targets[contTar].Npoints;contPts++){
                if(targets[contTar].AbsolutePos(0,contPts-1) < 2 &&
                        targets[contTar].AbsolutePos(1,contPts-1) < 2 &&
                        targets[contTar].AbsolutePos(2,contPts-1) < 2){
                    targets[contTar].IMRelativePos(1,contPts)=0;
                    targets[contTar].IMRelativePos(2,contPts)=0;
                    targets[contTar].IMRelativePos(3,contPts)=0;
                    break;
                }
            }
        }
        WriteNode(elementnode,targets[contTar].IMRelativePos);
        groupnode[targets[contTar].ID].push_back(elementnode);
        elementnode.reset();
    }
    csvfile.close();
    if(savedata){
        FileOp TranslData(filenames.ptlocation,ios::in);
        TranslData.node["Relative"].push_back(groupnode);
        FileOp SaveData(filenames.ptlocation,ios::out);
        SaveData.node = TranslData.node;
    }
}
VOID Icalib::IniTsai(INT cameraidx){
    if(Ntarget == 0){
        cout << "Não dados de alvos carregados" << endl;
        EXIT_FAILURE;
    }
    if(Ncamera == 0){
        cout << "Não há dados de cameras carregados" << endl;
        EXIT_FAILURE;
    }
    cout << endl << "Calculating (r1,r2,tx,ty) to camera: " << calccameradata[cameraidx].ID << endl;
    cout.flush();
    INTERVAL_VECTOR Xw,xc;
    INTERVAL_MATRIX RRt(3,4),IMxcnormal,T;
    YAML::Node tmpnode;
    if(file_exists(filenames.rrtcalc)){
        FileOp RRtfilenode(filenames.rrtcalc,ios::in);
        if(RRtfilenode.IsOpen() && calccameradata[cameraidx].RRtFromFile){
                tmpnode = RRtfilenode.node[calccameradata[cameraidx].ID][0]["RRtFinal"];
                if(!ReadNode(tmpnode,0,calccameradata[cameraidx].RRt)){
                    RRtfilenode.SaveFile.close();
                    RRtfilenode.node.reset();
                    tmpnode.reset();
                    FileOp RRtfile(filenames.rrtcalc);
                    WriteNode(tmpnode,calccameradata[cameraidx].RRt);
                    RRtfile.node[calccameradata[cameraidx].ID][0]["RRtFinal"].push_back(tmpnode);
                    tmpnode.reset();
                }
        }    
    }
//    INT Nrows = 0;
//    for(int i = 0;i < calccameradata[cameraidx].vxci.size();i++)
//        Nrows += ColDimension(calccameradata[cameraidx].vxci[i]);
//    Resize(A_T1C,Nrows,Dimension_X);
//    Initialize(A_T1C,0);
//    Resize(b,Nrows);
//    if(usetruedata && !Fixedtx_ty){
//        optmiset.TheDomain = optmiset.TheDomain/truecameradata[cameraidx].RRt(2,4);
//    }
    FileOp PixelData(filenames.pxlocation,ios::in),TranslationData(filenames.ptlocation,ios::in);
    int contColxc = 1;
    int nDiam = PixelData.node.size();
    INTERVAL_MATRIX vxcin;
    for(int contDiam = 0;contDiam < nDiam;contDiam++){
        INTERVAL_MATRIX IMxc,IMXw;
        for(int contTar = 0;contTar < calccameradata[cameraidx].vxci.size();contTar++){
            if(!ReadNode(
             PixelData.node[contDiam][calccameradata[cameraidx].ID][targets[calccameradata[cameraidx].TargIdx(contTar)].ID][0],
                    0,
                    vxcin))
                continue;
            IMxc.hcat(vxcin);
        }
        Normalize2Dpoints(IMxcnormal,IMxc,T);
        int nTransl = TranslationData.node["Relative"].size();
        double CamDiameter,CamNoise;
        if(!ReadNode(PixelData.node[contDiam]["Diameter"][0],0,CamDiameter)){
            cout << "Can't read cam diameter data" << endl;
            return;
        }
        if(!ReadNode(PixelData.node[contDiam]["Noise"][0],0,CamNoise)){
            cout << "Can't read cam noise data" << endl;
            return;
        }
        INT Dimension_X = 7;
        for(int contTrans = 0;contTrans < nTransl;){
            Resize(optmiset.TheDomain,Dimension_X);
            Initialize(RRt,INTERVAL(-40,40));
            int contDomain = 1;
//            for(int i = 1;i <= 2;i++){
//                for(int j = 1;j <= 4;j++){
//                    if(i == 2 && j == 4)continue;
//                    if(j == 4){
//                        if(!Fixedtx_ty){
//                            optmiset.TheDomain(contDomain) = calccameradata[cameraidx].RRt(i,j);
//                            contDomain++;
//                        }else{
//                            continue;
//                        }
//                    }else{
//                        optmiset.TheDomain(contDomain) = calccameradata[cameraidx].RRt(i,j);
//                        contDomain++;
//                    }
//                }
//            }
            double TarDiameter,TarNoise;
            if(!ReadNode(TranslationData.node["Relative"][contTrans]["Diameter"][0],0,TarDiameter)){
                cout << "Can't read target diameter data" << endl;
                return;
            }
            if(!ReadNode(TranslationData.node["Relative"][contTrans]["Noise"][0],0,TarNoise)){
                cout << "Can't read target noise data" << endl;
                return;
            }
            int contRow = 1;
            INTERVAL_MATRIX IMRelativePos;
            IMxc.Delete();
            AbLinsys albl;
            for(int contTar = 0;contTar < calccameradata[cameraidx].vxci.size();contTar++){
                if(!ReadNode(TranslationData.node["Relative"][contTrans][targets[calccameradata[cameraidx].TargIdx(contTar)].ID][0],
                        0,
                        IMRelativePos))
                    continue;
                if(!ReadNode(PixelData.node[contDiam][calccameradata[cameraidx].ID][targets[calccameradata[cameraidx].TargIdx(contTar)].ID][0],0,vxcin))
                    continue;
                GetHomogeneous(vxcin,vxcin);
                vxcin = T*vxcin;
                for(int contPts = 1;contPts <= ColDimension(calccameradata[cameraidx].vxci[contTar]);contPts++){
                    Xw = Col(IMRelativePos,contPts);
                    xc = Col(vxcin,contPts);
                    contColxc++;
                    INTERVAL_VECTOR LA(7),Lb(1);
                    LA(1) = Xw(1)*xc(2);LA(2) = Xw(2)*xc(2);LA(3) = Xw(3)*xc(2);
                    LA(4) = xc(2);
                    LA(5) = -Xw(1)*xc(1);LA(6) = -Xw(2)*xc(1);LA(7) = -Xw(3)*xc(1);
                    Lb(1) = xc(1);
                    contRow++;
                    albl.A.vcat(LA);
                    albl.b.vcat(Lb);
                }
            }
            INTERVAL_VECTOR Xlin;
            Xlin = optmiset.TheDomain;
            LinOpt linopt(albl.A,albl.b,Xlin);
            linopt.optmiset = optmiset;
            linopt.RohnOptimalLinSys();
            MATRIX Ac,Delta;
            VECTOR bc,delta;
            Ac = Mid(albl.A);
            bc = Mid(albl.b);
            Delta = Sup(albl.A)-Mid(albl.A);
            delta = Sup(albl.b)-Mid(albl.b);
            REAL alpha = 1e-5;
            bool passou = false;
            if(linopt.X.max() > 1e10){
                Initialize(linopt.X,0);
                while(linopt.X.max() < 1e10 && alpha <= 1){
                    linopt.A = Ac+alpha*SymHull(Delta);
                    linopt.B = bc+alpha*SymHull(delta);
                    linopt.RohnOptimalLinSys();
                    alpha += 1e-2;
                    passou = true;
                }
            }
            if(passou){
                alpha -= 1e-2;
                linopt.A = Ac+alpha*SymHull(Delta);
                linopt.B = bc+alpha*SymHull(delta);
                linopt.RohnOptimalLinSys();
            }
//            cout << "X=" << linopt.X ;
            std::shared_ptr<VOID> userdata;// = make_shared<AbLinsys>(albl);
            INTERVAL MinValue;// =  AxXFunc2ADI(linopt.X,userdata);
//            cout << " = (" << MinValue << ")" << endl;
            std::unique_ptr<FUNCTION2> OptFunction(new FUNCTION2(Dimension_X,
                            NULL,
                            AxXFuncAD,
                            AxXFuncADR,
                            AxXFuncADI,
                            NULL,
                            NULL,
                            NULL,
                            NULL));
            GAGenome::Evaluator GAFunc;//(GAIniTsaiObjFunc);
            AbLinsys Albl;
            Albl.X = linopt.X;
            Albl.A = linopt.A;
            Albl.b = linopt.B;
            optmiset.TheDomain = Albl.X;
            if(Inf(MinValue) < 1e-1){
                GAFunc = GAGenome::Evaluator(GAIniTsaiObjFunc);
                optmiset.expansion = false;
            }else{
                GAFunc = GAGenome::Evaluator(GAxPXPObjFunc);
            }
            Albl.z = ynset(Albl.A.cols());
            linopt.X = HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset);
            optmiset.expansion = true;
            cout << "X=" << linopt.X ;
            std::shared_ptr<VOID> userdata2 = make_shared<AbLinsys>(albl);
            MinValue =  AxXFunc2ADI(linopt.X,userdata2);
            cout << " = (" << MinValue << ")" << endl;
            INT contHeuBox = 1;
            if(AnyNAN(linopt.X)){
                Initialize(linopt.X,SymHull(Machine::PosInfinity));
            }
            for(int i = 1;i <= 2;i++){
                for(int j = 1;j <= 4;j++){
                    if(i == 2 && j == 4)continue;
                    if(j == 4){
                        RRt(i,j)=linopt.X(contHeuBox);
                        contHeuBox++;
                    }else{
                        RRt(i,j)=linopt.X(contHeuBox);
                        contHeuBox++;
                    }
                }
            }
            INTERVAL c,d,s;
            d = Sqrt(Sqr(RRt(2,1))+
                Sqr(RRt(2,2))+
                Sqr(RRt(2,3)));
            EINTERVAL tyFinal(Machine::NegInfinity,Machine::PosInfinity);
            if(d.inf() != Machine::NegInfinity && d.sup() != Machine::PosInfinity){
                if(!(0.0 <= d)){
                    for(int i = 1;i <= albl.A.rows();i++){
                        EINTERVAL ty(0,0);
                        INT contX = 1;
                        for(int j = 1;j <= albl.A.cols();j++){
                            ty += linopt.X(contX)*albl.A(i,j);
                            contX++;
                        }
                        ty = ty/albl.b(i);
                        EIntervalIntersection(tyFinal,ty,tyFinal);
                    }
                    if(!tyFinal.empty())
                        c = tyFinal(1)/d;
                    else
                        c = 1/d;
                }else{
                    for(int i = 1;i <= albl.A.rows();i++){
                        EINTERVAL ty(0,0);
                        INT contX = 1;
                        for(int j = 1;j <= albl.A.cols();j++){
                            ty += linopt.X(contX)*albl.A(i,j);
                            contX++;
                        }
                        ty = ty/albl.b(i);
                        EIntervalIntersection(tyFinal,ty,tyFinal);
                    }
                    if(!tyFinal.empty())
                        c = tyFinal(1)/Mid(d);
                    else
                        c = 1/d;
                }
            }else{
                c = SymHull(Machine::PosInfinity);
            }
            YAML::Node elementnode,tmpnode;
            INTERVAL_VECTOR tmpR;
//            WriteNode(elementnode,calccameradata[cameraidx].RRt);
//            tmpnode["RRtBcBO"].push_back(elementnode);
//            elementnode.reset();
//            WriteNode(elementnode,c);
//            tmpnode["c"].push_back(elementnode);
//            elementnode.reset();
//            WriteNode(elementnode,s);
//            tmpnode["s"].push_back(elementnode);
//            elementnode.reset();
            WriteNode(elementnode,CamDiameter);
            tmpnode["imgdiam"].push_back(elementnode);
            elementnode.reset();
            WriteNode(elementnode,TarDiameter);
            tmpnode["tardiam"].push_back(elementnode);
            elementnode.reset();
            WriteNode(elementnode,CamNoise);
            tmpnode["imgnoise"].push_back(elementnode);
            elementnode.reset();
            WriteNode(elementnode,TarNoise);
            tmpnode["tarnoise"].push_back(elementnode);
            elementnode.reset();
//            WriteNode(elementnode,ApproximationList.front().MinPoint);
//            tmpnode["ApproxMinPoint"].push_back(elementnode);
//            elementnode.reset();
            tmpR = c*Row(RRt,2);
            SetRow(RRt,2,tmpR);
            tmpR = c*Row(RRt,1);
            SetRow(RRt,1,tmpR);
            RRt(2,4) = c;
//            s = c*Sqrt(Sqr(RRt(1,1))+
//                Sqr(RRt(1,2))+
//                Sqr(RRt(1,3)));
//            cout << "s=" << s << endl;
            OrthonormalizeR1X_R2XOptm(cameraidx,RRt);
//            cout << "RRt=" << RRt << endl;
//            cout << "TrueRRt=" << truecameradata[cameraidx].RRt << endl;
            if(RRt.hasNan()){
                Initialize(RRt,SymHull(Machine::PosInfinity));
            }
            WriteNode(elementnode,RRt);
            tmpnode["RRtAcAO"].push_back(elementnode);
            elementnode.reset();
            FileOp SaveData(filenames.rrtcalc);
            SaveData.node[calccameradata[cameraidx].ID].push_back(tmpnode);
            contTrans++;
        }
    }
}

VOID Icalib::OrthonormalizeR1X_R2XOptm(INT cameraidx,INTERVAL_MATRIX& RRt){
    INTERVAL k;
    INTERVAL_MATRIX R;
    INTERVAL_VECTOR R1X,R2X,R3X;
    R = RRt.Box(1,3,1,3);
    R1X = Row(R,1);
    R2X = Row(R,2);
    INTERVAL dotproduct;
    dotproduct = R1X*R2X;
    REAL maxval = 1;
    if(!(0.0 <= dotproduct)) {
        INTERVAL_VECTOR coeffs(3);
        INT sgn = 1;
        coeffs(1) = R1X*R2X;
        coeffs(2) = R1X*R1X+R2X*R2X;
        coeffs(3) = R1X*R2X;
        if(Sup(coeffs(1)) < 0)sgn = -1;
        coeffs = sgn*coeffs;
        INTERVAL_VECTOR roots;
        roots = Roots(coeffs);
        int i;
        for(i = 1;i <= Dimension(roots) && !std::isnan(Inf(roots(i)));i++);
        if(i == 1)
            k = -(1.0/2.0)*R1X*R2X;
        else if(i == 2)
            k = sgn*roots(1);
        else if(i > 2)
            k = sgn*roots(2);
        if(std::isinf(Inf(k)) || std::isinf(Sup(k)))
            k = -(1.0/2.0)*R1X*R2X;
        R1X = R1X+k*R2X;
        R2X = R2X+k*R1X;
        dotproduct = R1X*R2X;
        if(!(0.0 <= dotproduct)){
            INTERVAL_VECTOR tmpR1_2X;
            tmpR1_2X = R1X;
            tmpR1_2X.vcat(R2X);
            GAGenome::Evaluator GAFunc(GADotObjFunc);
            Resize(optmiset.TheDomain,6);
            Initialize(optmiset.TheDomain,SymHull(1));
            std::unique_ptr<FUNCTION2> NULLFUNC;
            AbLinsys Albl;
            Albl.X = tmpR1_2X;
            tmpR1_2X = HeuristicExpansion(NULLFUNC,Albl,GAFunc,optmiset);
            R1X = tmpR1_2X.Box(1,3);
            R2X = tmpR1_2X.Box(4,6);
        }
    }
    R3X = Cross(R1X,R2X);
//    if(R1X.max() > maxval)
//        R1X = R1X/R1X.max();
//    if(R2X.max() > maxval)
//        R2X = R2X/R2X.max();
//    if(R3X.max() > maxval)
//        R3X = R3X/R3X.max();
    SetRow(R,1,R1X);
    SetRow(R,2,R2X);
    SetRow(R,3,R3X);
//    R = R/maxval;
//    INTERVAL_MATRIX Rdiam(3,3);
//    MATRIX MRdiam(3,3);
//    Initialize(MRdiam,1e-3);
//    Rdiam = R;
//    if(!isRegular(R)){
//        cout << "R=" << R << endl;
//        R = Mid(R);
//        while(1){
//            R = R + Rdiam;
//            if(!isRegular(R)){
//                R = R - Rdiam;
//                break;
//            }
//        }
//        cout << "R=" << R << endl;
//        return;
//    }
    for(int i = 1;i <= 3;i++)
        for(int j = 1;j <= 3;j++)
            RRt(i,j)=R(i,j);
}
VOID Icalib::f_t_zest(INT cameraidx){
    cout << endl << "Calculating (f,tz) to camera: " << calccameradata[cameraidx].ID << endl;
    cout.flush();
    INTERVAL_VECTOR Xw,xc;
    INTERVAL_MATRIX RRt,K(3,3);
    Resize(optmiset.TheDomain,2);
    INT contDomain = 1;
    contDomain++;
    FileOp TranslationData(filenames.ptlocation,ios::in);
    if(!TranslationData.IsOpen()){
        cout << "No translation data found" << endl;
        return;
    }
    FileOp PixelData(filenames.pxlocation,ios::in);
    if(!PixelData.IsOpen()){
        cout << "No camera pixel data found" << endl;
        return;
    }
    FileOp RRtCalcData(filenames.rrtcalc,ios::in);
    if(!RRtCalcData.IsOpen()){
        cout << "No camera pixel data found" << endl;
        return;
    }
    int contRRtDiam;
    Node elementnode;
//    INTERVAL s;
    for(int contPxDiam = 0;contPxDiam < PixelData.node.size();contPxDiam++){
        for(int contTrans = 0;contTrans < TranslationData.node["Relative"].size();contTrans++){
            INTERVAL_MATRIX vxcin;
            INTERVAL_MATRIX IMxc,T,IMxcnormal,IMXw;
            double PxDiam,PxNoise;
            if(!usetruedata){
                elementnode = PixelData.node[contPxDiam]["Diameter"][0];
                if(!ReadNode(elementnode,0,PxDiam))continue;
                elementnode.reset();
                elementnode = PixelData.node[contPxDiam]["Noise"][0];
                if(!ReadNode(elementnode,0,PxNoise))continue;
                elementnode.reset();
                for(int contTar = 0;contTar < targets.size();contTar++){
                    elementnode = PixelData.node[contPxDiam][calccameradata[cameraidx].ID][targets[contTar].ID][0];
                    if(!ReadNode(elementnode,0,vxcin)){
                        elementnode.reset();
                        continue;
                    }
                    elementnode.reset();
                    IMxc.hcat(vxcin);
                }
                Normalize2Dpoints(IMxcnormal,IMxc,T);
                double RRtDiam;
                bool foundRRt = false;
                double TarDiameter,TarNoise;
                if(!ReadNode(TranslationData.node["Relative"][contTrans]["Diameter"][0],0,TarDiameter)){
                    cout << "Can't read target diameter data" << endl;
                    return;
                }
                if(!ReadNode(TranslationData.node["Relative"][contTrans]["Noise"][0],0,TarNoise)){
                    cout << "Can't read target noise data" << endl;
                    return;
                }
                for(contRRtDiam = 0;contRRtDiam < RRtCalcData.node[calccameradata[cameraidx].ID].size();contRRtDiam++){
                    elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["imgdiam"][0];
                    if(!ReadNode(elementnode,0,RRtDiam)){
                        elementnode.reset();
                        continue;
                    }
                    elementnode.reset();
                    if((abs(PxDiam-RRtDiam)) > 1e-5){
                        continue;
                    }
                    elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["tardiam"][0];
                    if(!ReadNode(elementnode,0,RRtDiam)){
                        elementnode.reset();
                        continue;                        
                    }
                    elementnode.reset();
                    if((abs(TarDiameter-RRtDiam)) > 1e-5){
                        continue;
                    }
                    foundRRt = true;
                    break;
                }
                if(!foundRRt)continue;

                elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["RRtAcAO"][0];
                if(!ReadNode(elementnode,0,RRt)){
                    elementnode.reset();
                    continue;
                }
                optmiset.TheDomain(contDomain) = RRt(3,4);
            }else{
                RRt = truecameradata[cameraidx].RRt;
                optmiset.TheDomain(contDomain) = truecameradata[cameraidx].RRt(3,4);
            }
            INTERVAL_MATRIX IMRelativePos;
            INTERVAL_MATRIX R,EulAng;
            INTERVAL_VECTOR b;
            INTERVAL_MATRIX A_T1C;
            for(int targetidx = 0;targetidx < calccameradata[cameraidx].vxci.size();targetidx++){
                if(!ReadNode(TranslationData.node["Relative"][contTrans][targets[calccameradata[cameraidx].TargIdx(targetidx)].ID][0],
                        0,
                        IMRelativePos))
                    continue;
                elementnode = PixelData.node[contPxDiam][calccameradata[cameraidx].ID][targets[targetidx].ID][0];
                if(!ReadNode(elementnode,0,vxcin)){
                    elementnode.reset();
                    continue;
                }
                elementnode.reset();
                GetHomogeneous(vxcin,vxcin);
                vxcin = T*vxcin;
                for(int contPts = 1;contPts <= ColDimension(calccameradata[cameraidx].vxci[targetidx]);contPts++){
                    Xw = Col(IMRelativePos,contPts);
                    IMXw.hcat(Xw);
                    xc = Col(vxcin,contPts);
                    INTERVAL_VECTOR LA(2),Lb(1);
                    LA(1) = RRt(1,1)*Xw(1)+RRt(1,2)*Xw(2)+RRt(1,3)*Xw(3)+RRt(1,4);
                    LA(2)  = -xc(1);
                    A_T1C.vcat(LA);
                    Lb(1) = (RRt(3,1)*Xw(1)+RRt(3,2)*Xw(2)+RRt(3,3)*Xw(3))*xc(1)+xc(1);
                    b.vcat(Lb);
                    LA(1) = RRt(2,1)*Xw(1)+RRt(2,2)*Xw(2)+RRt(2,3)*Xw(3)+RRt(2,4);
                    LA(2) = -xc(2);
                    A_T1C.vcat(LA);
                    Lb(1) = (RRt(3,1)*Xw(1)+RRt(3,2)*Xw(2)+RRt(3,3)*Xw(3))*xc(2)+xc(2);
                    b.vcat(Lb);
                }
            }
            INTERVAL_VECTOR Xlin;
            Xlin = optmiset.TheDomain;
            Initialize(Xlin,SymHull(Machine::PosInfinity));
            LinOpt linopt(A_T1C,b,Xlin);
            linopt.optmiset = optmiset;
            linopt.RohnOptimalLinSys();
            MATRIX Ac,Delta;
            VECTOR bc,delta;
            Ac = Mid(A_T1C);
            bc = Mid(b);
            Delta = Sup(A_T1C)-Mid(A_T1C);
            delta = Sup(b)-Mid(b);
            REAL alpha = 1e-5;
            bool passou = false;
            if(linopt.X.max() > 1e10){
                Initialize(linopt.X,0);
                while(linopt.X.max() < 1e10 && alpha <= 1){
                    linopt.A = Ac+alpha*SymHull(Delta);
                    linopt.B = bc+alpha*SymHull(delta);
                    linopt.RohnOptimalLinSys();
                    alpha += 1e-2;
                    passou = true;
                }
            }
            if(passou){
                alpha -= 1e-2;
                linopt.A = Ac+alpha*SymHull(Delta);
                linopt.B = bc+alpha*SymHull(delta);
                linopt.RohnOptimalLinSys();
            }
            AbLinsys albl;
            albl.A = linopt.A;
            albl.b = linopt.B;
            albl.X = linopt.X;
            cout << "X=" << linopt.X ;
            std::shared_ptr<VOID> userdata = make_shared<AbLinsys>(albl);
            INTERVAL MinValue =  AxXFunc2ADI(linopt.X,userdata);
            cout << " = (" << MinValue << ")" << endl;
//            std::unique_ptr<FUNCTION2> OptFunction(new FUNCTION2(2,
//                            NULL,
//                            AxXFuncAD,
//                            AxXFuncADR,
//                            AxXFuncADI,
//                            NULL,
//                            NULL,
//                            NULL,
//                            NULL));
//            GAGenome::Evaluator GAFunc;
//            AbLinsys Albl;
//            Albl.X = linopt.X;
//            Albl.A = linopt.A;
//            Albl.b = linopt.B;
//            optmiset.TheDomain = Albl.X;
//            if(Inf(MinValue) < 1){
//                GAFunc = GAGenome::Evaluator(GAftzObjFunc);
//                optmiset.expansion = false;
//            }else{
//                GAFunc = GAGenome::Evaluator(GAftzexpObjFunc);
//            }
//            Albl.z = ynset(Albl.A.cols());
//            linopt.X = HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset);
//            optmiset.expansion = true;
//            cout << "X=" << linopt.X ;
//            std::shared_ptr<VOID> userdata2 = make_shared<AbLinsys>(albl);
//            MinValue =  AxXFunc2ADI(linopt.X,userdata2);
//            cout << " = (" << MinValue << ")" << endl;
            RRt(3,4)=linopt.X(2);
            Initialize(K,0);
            K(1,1) = linopt.X(1);
            K(2,2) = linopt.X(1);
            K(1,3)=0;
            K(2,3)=0;
            K(3,3)=1;
            FileOp LoadData(filenames.rrtcalc,ios::in);
            YAML::Node elementnode,tmpnode,loadnode;
            WriteNode(elementnode,RRt);
            loadnode = LoadData.node;
            tmpnode = loadnode[calccameradata[cameraidx].ID][contRRtDiam];
            LoadData.SaveFile.close();
            LoadData.node.reset();
            tmpnode["RRtFinal"].push_back(elementnode);
            elementnode.reset();
            WriteNode(elementnode,K);
            tmpnode["K"].push_back(elementnode);
            elementnode.reset();
            FileOp SaveData(filenames.rrtcalc,ios::out);
            SaveData.node=loadnode;
        }
    }
//    char chcomputername[512];
//    gethostname(&chcomputername[0],512);
//    string computername;
//    computername = "-";
//    computername += chcomputername;
//    string date;
//    time_t tme = time(0);   // get time now
//    struct tm * now = localtime( & tme );
//    date = "-" + to_string(now->tm_mday) + "-" +
//    to_string(now->tm_mon + 1) + "-" +
//    to_string(now->tm_year + 1900);
//    string filename;
//    filename = filenames.rrtcalc;
//    int dotidx = filename.find_last_of(".");
//    filename.insert(dotidx,computername+date);
//    std::ifstream  src(filenames.rrtcalc, std::ios::binary);
//    std::ofstream  dst(filename,   std::ios::binary);
//    dst << src.rdbuf();
//    src.close();
//    dst.close();
}

bool readStdIn(std::string& line){
    struct pollfd pfd = { STDIN_FILENO, POLLIN, 0 };

    int ret = poll(&pfd, 1, 10);  // timeout of 1000ms
    if(ret == 1) // there is something to read
    {
        std::getline(std::cin, line);
        return true;
    }
    else if(ret == -1){
            std::cout << "Error: " << strerror(errno) << std::endl;
            return false;
    }
    return false;
}
INTERVAL_MATRIX Icalib::InitializeA_T1C(const INTERVAL_MATRIX& Xw, const INTERVAL_MATRIX& xc){
    INTERVAL_MATRIX A_T1C;
    if(ColDimension(Xw) != ColDimension(xc)){
        cerr << "Xw and xc dimension differ" << endl;
        return A_T1C;
    }
    Resize(A_T1C,ColDimension(xc)*2,12);
    Initialize(A_T1C,INTERVAL(0));
    int contColXw = 1;
    int contLA = 1;
    for(int contPts = 1;contPts <= ColDimension(xc);contPts++){
        for(int contColA = 5;contColA <= 12;contColA++){
            if(contColA < 9)
                A_T1C(contLA,contColA)=-xc(3,contColXw)*Xw(contColA-4,contColXw);
            else
                A_T1C(contLA,contColA)=xc(2,contColXw)*Xw(contColA-8,contColXw);
        }
        contLA++;
        for(int contColA = 1;contColA <= 12;contColA++){
            if(contColA < 5)
                A_T1C(contLA,contColA)=xc(3,contColXw)*Xw(contColA,contColXw);
            if(contColA >= 5 && contColA < 9)
                A_T1C(contLA,contColA)=0.0;
            if(contColA > 8)
                A_T1C(contLA,contColA)=-xc(1,contColXw)*Xw(contColA-8,contColXw);
        }
        contLA++;
        contColXw++;
    }
    return A_T1C;
}
VOID Icalib::Triangulation( ){
    cout << "Starting triangulation ..." << endl;
    if(!file_exists(filenames.config)){
        cerr << "Config file doesn't exist" << endl;
        return;
    }
    if(!file_exists(filenames.rrtcalc)){
        cerr << "Data file doesn't exist" << endl;
        return;
    }
    std::vector<std::vector<INTERVAL_VECTOR> > p;
    std::vector<INTERVAL_VECTOR> tmpp;
    FileOp PixelData(filenames.pxlocation,ios::in),TranslationData(filenames.ptlocation,ios::in);
    int nDiam = PixelData.node.size();
    INTERVAL_MATRIX vxcin,IMxcnormal,T;
    vector<INTERVAL_MATRIX> vT;
    Node elementnode,cameranode;
    int nTransl = TranslationData.node["Relative"].size();
    fstream DataCSV(filenames.csvtriangresult,ios::in|ios::out|ios::app);
    
    for(int contDiam = 0;contDiam < nDiam;contDiam++){
        for(int contTrans = 0;contTrans < nTransl;contTrans++){
            INTERVAL_MATRIX IMXw,IMXwcalc;
            vector<INTERVAL_MATRIX> vIMxctar;
            vector<INTERVAL_MATRIX> vIMXw;
            double PxDiam,PxNoise;
            if(!ReadNode(PixelData.node[contDiam]["Diameter"][0],0,PxDiam)){
                cout << "Can't read cam diameter data" << endl;
                return;
            }
            if(!ReadNode(PixelData.node[contDiam]["Noise"][0],0,PxNoise)){
                cout << "Can't read cam noise data" << endl;
                return;
            }
            INTERVAL_MATRIX IMRelativePos;
            for(int contTar = 0;contTar < targets.size();contTar++){
                if(!ReadNode(TranslationData.node["Relative"][contTrans][targets[contTar].ID][0],
                        0,
                        IMRelativePos))
                    continue;
                vIMXw.push_back(IMRelativePos);
            }
            FileOp RRtCalcData(filenames.rrtcalc,ios::in);
            double TarDiameter,TarNoise;
            for(int cameraidx = 0;cameraidx < calccameradata.size();cameraidx++){
                INTERVAL_MATRIX P(3,4);
                INTERVAL_MATRIX IMxc,IMxcnormal;
                vector<INTERVAL_MATRIX> vIMxc;
                for(int contTar = 0;contTar < targets.size();contTar++){
                    if(!ReadNode(
                        PixelData.node[contDiam][calccameradata[cameraidx].ID][targets[contTar].ID][0],
                        0,
                        vxcin))
                        continue;
                    IMxc.hcat(vxcin);
                    vIMxc.push_back(vxcin);
                }
                Normalize2Dpoints(IMxcnormal,IMxc,T);
                vT.push_back(T);
                vIMxctar.push_back(IMxcnormal);
                vIMxc.clear();
                if(!ReadNode(TranslationData.node["Relative"][contTrans]["Diameter"][0],0,TarDiameter)){
                    cout << "Can't read target diameter data" << endl;
                    return;
                }
                if(!ReadNode(TranslationData.node["Relative"][contTrans]["Noise"][0],0,TarNoise)){
                    cout << "Can't read target noise data" << endl;
                    return;
                }
                int contRRtDiam;
                bool foundRRt = false;
                bool foundP = false;
                REAL imgdiam,tardiam;
                for(contRRtDiam = 0;contRRtDiam < RRtCalcData.node[calccameradata[cameraidx].ID].size();contRRtDiam++){
                    elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["imgdiam"][0];
                    if(!ReadNode(elementnode,0,imgdiam)){
                        elementnode.reset();
                        continue;
                    }
                    elementnode.reset();
                    if((abs(PxDiam-imgdiam)) > 1e-5){
                        continue;
                    }
                    elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["tardiam"][0];
                    if(!ReadNode(elementnode,0,tardiam)){
                        elementnode.reset();
                        continue;                        
                    }
                    elementnode.reset();
                    if((abs(TarDiameter-tardiam)) > 1e-5){
                        continue;
                    }
                    foundRRt = true;
                    break;
                }
                if(!foundRRt)continue;
                INTERVAL_MATRIX RRt,K;
                INTERVAL_VECTOR EulAngl,t;
                cameranode = RRtCalcData.node[calccameradata[cameraidx].ID];
                if(!usetruedata){
                    elementnode = cameranode[contRRtDiam]["EulAnglAAdj"];
                    if(ReadNode(elementnode[0],0,EulAngl)){
                        elementnode.reset();
                        foundRRt = true;
                        elementnode = cameranode[contRRtDiam]["KAAdj"];
                        if(ReadNode(elementnode[0],0,K)){
                            elementnode.reset();
                            foundRRt = true;
                            elementnode = cameranode[contRRtDiam]["tAAdj"];
                            if(ReadNode(elementnode[0],0,t)){
                                elementnode.reset();
                                foundRRt = true;
                            }else{
                                foundRRt = false;
                            }
                        }else{
                            foundRRt = false;
                        }
                    }else{
                        foundRRt = false;
                    }
                    if(foundRRt){
//                        cout << "K=" << K << endl;
//                        cout << "EulAngl=" << EulAngl << endl;
//                        cout << "t=" << t << endl;
                        elementnode.reset();
                        RRt = Eul2rtm(EulAngl);
//                        cout << "RRt=" << RRt << endl;
                        t = -RRt*t;
//                        cout << "t=" << t << endl;
                        RRt.hcat(t);
                    }else{
                        foundP = false;
                        elementnode = cameranode[contRRtDiam]["RRtFinal"];
                        if(ReadNode(elementnode[0],0,P)){
                            elementnode.reset();
                            foundRRt = true;
                            foundP = true;
                        }else{
                            foundRRt = false;
                        }
                        elementnode.reset();                       
                    }
                }else{
                    imgnoise = 0;
                    tarnoise = 0;
                }
                if(foundRRt && !foundP){
                    P = K*RRt;
                }else if(!foundRRt && !foundP){
                    break;
                }                
//                cout << "P=" << P << endl;
                INTERVAL_VECTOR prow(4);
                for(int i = 1;i <= RowDimension(P);i++){
                    prow = Row(P,i);
                    tmpp.push_back(prow);
                }
                p.push_back(tmpp);
                tmpp.clear();
            }
            cout << "imgnoise: " << PxNoise << " - ";
            cout << "tarnoise: " << TarNoise << " - ";
            cout << "imgdiam: " << PxDiam << " - ";
            cout << "tardiam: " << TarDiameter << endl;
            Node tmpnode;
            WriteNode(elementnode,tarnoise);
            tmpnode["Diameter"].push_back(elementnode);
            elementnode.reset();
            INT nRows = pow(2,calccameradata.size());
            INT Dimension_X = 3;
            INTERVAL_MATRIX A_T1C;
            INTERVAL_VECTOR b;
            for(int contTar = 0;contTar < targets.size();contTar++){
                    for(int contPts = 1;contPts <= ColDimension(vIMXw[contTar]);contPts++){
                        INTERVAL_VECTOR xc;
                        INT contRowA = 1;
                        INT IniX = 1;
                        Resize(A_T1C,nRows,Dimension_X);
                        Initialize(A_T1C,0);
                        Resize(b,nRows);
                        Initialize(b,0);
                        xc = Col(vIMxctar[0],(contTar*8)+contPts);
                        GetHomogeneous(xc,xc);
                        //xc = vT[0]*xc;
                        for(int contCam = 1;contCam < calccameradata.size();contCam++){
                            INTERVAL_VECTOR xcl,L;
                            xcl = Col(vIMxctar[contCam],(contTar*8)+contPts);
                            GetHomogeneous(xcl,xcl);
                            //xcl = vT[contCam]*xcl;
                            L = (xc(1)*p[0][2]-p[0][0]);
                            SetRow(A_T1C,contRowA,L.Box(IniX,Dimension_X));
                            b(contRowA) = -L(4);
                            contRowA++;
                            L = (xc(2)*p[0][2]-p[0][1]);
                            SetRow(A_T1C,contRowA,L.Box(IniX,Dimension_X));
                            b(contRowA) = -L(4);
                            contRowA ++;
                            L = (xcl(1)*p[contCam][2]-p[contCam][0]);
                            SetRow(A_T1C,contRowA,L.Box(IniX,Dimension_X));
                            b(contRowA) = -L(4);
                            contRowA ++;
                            L = (xcl(2)*p[contCam][2]-p[contCam][1]);
                            SetRow(A_T1C,contRowA,L.Box(IniX,Dimension_X));
                            b(contRowA) = -L(4);
                            contRowA ++;
                        }
                        std::unique_ptr<FUNCTION2> OptFunction(new FUNCTION2(Dimension_X,
                                        NULL,
                                        AxXFuncAD,
                                        AxXFuncADR,
                                        AxXFuncADI,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL));
//                        cout << "A=" << A_T1C << endl;
//                        cout << "B=" << b << endl;
                        LinOpt linopt(A_T1C,b);
                        linopt.optmiset.MaxIter = optmiset.MaxIter;
                        linopt.RohnOptimalLinSys();
                        MATRIX Ac,Delta;
                        VECTOR bc,delta;
                        Ac = Mid(A_T1C);
                        bc = Mid(b);
                        Delta = Sup(A_T1C)-Mid(A_T1C);
                        delta = Sup(b)-Mid(b);
                        REAL alpha = 1e-5;
                        bool passou = false;
                        if(linopt.X.max() > 1e10){
                            Initialize(linopt.X,0);
                            while(linopt.X.max() < 1e10 && alpha <= 1){
                                linopt.A = Ac+alpha*SymHull(Delta);
                                linopt.B = bc+alpha*SymHull(delta);
                                linopt.RohnOptimalLinSys();
                                alpha += 1e-2;
                                passou = true;
                            }
                        }
                        if(passou){
                            alpha -= 1e-2;
                            linopt.A = Ac+alpha*SymHull(Delta);
                            linopt.B = bc+alpha*SymHull(delta);
                            linopt.RohnOptimalLinSys();
                        }
                        AbLinsys albl;
                        albl.A = linopt.A;
                        albl.b = linopt.B;
                        albl.X = linopt.X;
                        std::shared_ptr<VOID> userdata = make_shared<AbLinsys>(albl);
                        INTERVAL MinValue =  AxXFunc2ADI(linopt.X,userdata);
                        GAGenome::Evaluator GAFunc;
                        AbLinsys Albl;
                        Albl.X = linopt.X;
                        Albl.A = linopt.A;
                        Albl.b = linopt.B;
                        optmiset.TheDomain.clear();
                        optmiset.TheDomain = Albl.X;
                        if(Inf(MinValue) < 1e-1){
                            GAFunc = GAGenome::Evaluator(GAIniTsaiObjFunc);
                            optmiset.expansion = false;
                        }else{
                            GAFunc = GAGenome::Evaluator(GAxPXPObjFunc);
                        }
                        //linopt.X = HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset);
                        optmiset.expansion = true;
                        IMXwcalc.hcat(linopt.X);
                    }
                    INTERVAL_MATRIX X,Xl;
                    XYZINTRelatLoad(0,0);
                    PointProjection(0,0);
                    Xl = targets[contTar].IMRelativePos;
                    INTERVAL_MATRIX Xlcalc;
                    Xlcalc = IMXwcalc;
                    if(Xlcalc.rows() > 3)
                        GetNonHomogeneous(Xlcalc,Xlcalc);
                    if(Xl.rows() > 3)
                        GetNonHomogeneous(Xl,Xl);
                    INTERVAL diffnorm(0);
                    INTERVAL_MATRIX diffXl;
                    diffXl = Xlcalc-Xl;
                    int contbelongs = 0;
                    for(int i = 1;i <= ColDimension(Xl);i++){
                        diffnorm += Norm(diffXl.col(i));
                        if(Col(Xl,i) <= Col(Xlcalc,i))
                            contbelongs++;
                    }
                    double belongspercent = 1.0*contbelongs/ColDimension(Xl);
                    cout << "Error=" << diffnorm << endl;
                    DataCSV << PxNoise << "," << TarNoise << "," << PxDiam << "," << TarDiameter << ","
                            << Inf(diffnorm) << "," << Sup(diffnorm) << "," << belongspercent << endl;
//                    for(int i = 1;i <= RowDimension(Xlcalc);i++)
//                        for(int j = 1;j <= ColDimension(Xlcalc);j++)
//                            if(Intersection(Xlcalc(i,j),Xlcalc(i,j),vIMXw[contTar](i,j)) != 1)
//                                Xlcalc(i,j) = vIMXw[contTar](i,j);
                    FileOp TranslData(filenames.calcptlocation);
                    WriteNode(elementnode,Xlcalc);
                    IMXwcalc.Delete();
                    tmpnode[targets[contTar].ID].push_back(elementnode);
                    elementnode.reset();
                    TranslData.node["Relative"].push_back(tmpnode);
                    tmpnode.reset();
            }
            p.clear();
        }
    }
    DataCSV.close();
}
VOID Icalib::ExpansionP(INT cameraidx){
    if(!file_exists(filenames.config))return;
    if(!file_exists(filenames.rrtcalc))return;
    cout << endl << "Expanding P to camera: " << calccameradata[cameraidx].ID << endl;
    YAML::Node elementnode;
    double imgnoise,tarnoise,imgdiam,tardiam;
    FileOp RRtCalcData(filenames.rrtcalc,ios::in);
    YAML::Node cameranode;

    cameranode = RRtCalcData.node[calccameradata[cameraidx].ID];
    Emitter emitter;
    emitter << cameranode;
    if(emitter.size() == 0)return;
    INT Dimension_X = 12;
    Resize(optmiset.TheDomain,Dimension_X);
    if(!RRtCalcData.IsOpen() && !calccameradata[cameraidx].RRtFromFile)return;
    INTERVAL_MATRIX R(3,3),RRt(3,4),K(3,3);
    INTERVAL_MATRIX BothEulAngles;
    FileOp PixelData(filenames.pxlocation,ios::in),TranslationData(filenames.ptlocation,ios::in);
    int nDiam = PixelData.node.size();
    INTERVAL_MATRIX vxcin,IMxcnormal,T;
    for(int contDiam = 0;contDiam < nDiam;contDiam++){
        INTERVAL_MATRIX IMxc,IMxcnormal,TrueIMxcnormal;
        XYZINTRelatLoad(0,0);
        PointProjection(0,0);
        for(int contTar = 0;contTar < targets.size();contTar++){
            if(!ReadNode(
             PixelData.node[contDiam][calccameradata[cameraidx].ID][targets[contTar].ID][0],
                    0,
                    vxcin))
                continue;
            IMxc.hcat(vxcin);
            vxcin = truecameradata[cameraidx].vxci[contTar];
            TrueIMxcnormal.hcat(vxcin);
        }
        Normalize2Dpoints(IMxcnormal,IMxc,T);
        if(RowDimension(TrueIMxcnormal) < 3)
            GetHomogeneous(TrueIMxcnormal,TrueIMxcnormal);
        TrueIMxcnormal = T*TrueIMxcnormal;
        GetNonHomogeneous(TrueIMxcnormal,TrueIMxcnormal);
        int nTransl = TranslationData.node["Relative"].size();
        double PxDiam,PxNoise;
        if(!ReadNode(PixelData.node[contDiam]["Diameter"][0],0,PxDiam)){
            cout << "Can't read cam diameter data" << endl;
            return;
        }
        if(!ReadNode(PixelData.node[contDiam]["Noise"][0],0,PxNoise)){
            cout << "Can't read cam noise data" << endl;
            return;
        }
        for(int contTrans = 0;contTrans < nTransl;contTrans++){
            INTERVAL_MATRIX Xw,TrueXw;
            double TarDiameter,TarNoise;
            if(!ReadNode(TranslationData.node["Relative"][contTrans]["Diameter"][0],0,TarDiameter)){
                cout << "Can't read target diameter data" << endl;
                return;
            }
            if(!ReadNode(TranslationData.node["Relative"][contTrans]["Noise"][0],0,TarNoise)){
                cout << "Can't read target noise data" << endl;
                return;
            }
            INTERVAL_MATRIX IMRelativePos;
            int contRRtDiam;
            bool foundRRt = false;
            for(contRRtDiam = 0;contRRtDiam < RRtCalcData.node[calccameradata[cameraidx].ID].size();contRRtDiam++){
                elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["imgdiam"][0];
                if(!ReadNode(elementnode,0,imgdiam)){
                    elementnode.reset();
                    continue;
                }
                elementnode.reset();
                if((abs(PxDiam-imgdiam)) > 1e-5){
                    continue;
                }
                elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["tardiam"][0];
                if(!ReadNode(elementnode,0,tardiam)){
                    elementnode.reset();
                    continue;                        
                }
                elementnode.reset();
                if((abs(TarDiameter-tardiam)) > 1e-5){
                    continue;
                }
                foundRRt = true;
                break;
            }
            if(!foundRRt)continue;
                elementnode = cameranode[contRRtDiam]["RRtFinal"];
                if(!ReadNode(elementnode[0],0,RRt)){
                    elementnode.reset();
                    continue;
                }
                elementnode.reset();
                elementnode = cameranode[contRRtDiam]["K"];
                if(!ReadNode(elementnode[0],0,K))continue;
                elementnode.reset();
            cout << "imgnoise: " << PxNoise << " - ";
            cout << "tarnoise: " << TarNoise << " - ";
            cout << "imgdiam: " << PxDiam << " - ";
            cout << "tardiam: " << TarDiameter << " - " << endl;
            INTERVAL_MATRIX P(3,4);
            K(1,3)=0;
            K(2,3)=0;
            if(usetruedata){
                IMxc.Delete();
                for(int contTar = 0;contTar < calccameradata[cameraidx].vxci.size();contTar++){
                    IMxc.hcat(truecameradata[cameraidx].vxci[contTar]);
                    Xw.hcat(targets[truecameradata[cameraidx].TargIdx(contTar)].IMRelativePos);
                }
                Normalize2Dpoints(IMxcnormal,IMxc,T);
            }else{
                for(int contTar = 0;contTar < targets.size();contTar++){
                    if(!ReadNode(TranslationData.node["Relative"][contTrans][targets[contTar].ID][0],
                            0,
                            IMRelativePos))
                        continue;
                    Xw.hcat(IMRelativePos);
                    IMRelativePos = MatrixXd2MATRIX(targets[contTar].RelativePos);
                    TrueXw.hcat(IMRelativePos);
                }
                Normalize2Dpoints(IMxcnormal,IMxc,T);
                P = K*RRt;
            }
            INTERVAL_MATRIX b;
            Initialize(b,0);
            optmicalibdata.K = K;
            optmicalibdata.RRt = RRt;
            if(RowDimension(Xw) < 4)
                GetHomogeneous(Xw,Xw);
            optmicalibdata.Xw = Xw;
            Resize(optmicalibdata.max_min_xc,3,2);
            INTERVAL_VECTOR max_xc(3),min_xc(3),t;
            min_xc = INTERVAL_VECTOR(3,{INTERVAL(-calccameradata[cameraidx].ImageSize(0)/2),
                                        INTERVAL(-calccameradata[cameraidx].ImageSize(1)/2),
                                        1});
            max_xc = INTERVAL_VECTOR(3,{INTERVAL(calccameradata[cameraidx].ImageSize(0)/2),
                                        INTERVAL(calccameradata[cameraidx].ImageSize(1)/2),
                                        1});
            min_xc = T*min_xc;
            max_xc = T*max_xc;
            SetCol(optmicalibdata.max_min_xc,1,Col(max_xc,1));
            SetCol(optmicalibdata.max_min_xc,2,Col(min_xc,1));
            Resize(optmicalibdata.t,3);
            optmicalibdata.t = INTERVAL_VECTOR(3,{Hull(100),INTERVAL(-100,0),1});
            R = RRt.Box(1,3,1,3);
            INTERVAL_MATRIX Xwfiltered,xcfiltered;
            filter_points(Xw,IMxcnormal,RRt,K,Xwfiltered,xcfiltered);
            bool appfirstline = true;
            if(file_exists(filenames.csvpxresult))
                appfirstline = false;
            fstream csvfile(filenames.csvpxresult,ios::in|ios::out|ios::app);
            if(appfirstline)
                csvfile << "CAMID,PxDiam,TarDiam,PxNoise,TarNoise,"
                        << "Inf(Xw-xc),Sup(Xw-xc),NBlngs,Mid(A)-Mid(B),"
                        << "Inf(TrueXw-xc),Sup(TrueXw-xc),NBlngs,Mid(A)-Mid(B),"
                        << "Inf(Xw-Truexc),Sup(Xw-Truexc),NBlngs,Mid(A)-Mid(B),"
                        << "Inf(TrueXw-Truexc),Sup(TrueXw-Truexc),NBlngs,Mid(A)-Mid(B)," 
                        << "Inf(TrueXwUnnorm-TruexcUnnorm),Sup(TrueXwUnnorm-TruexcUnnorm),NBlngs,Mid(A)-Mid(B),Regular" << ","
                        << "InfAng1" << "," << "SupAng1" << "," << "InfAng2" << "," << "SupAng2" << "," << "InfAng3" << "," 
                        << "SupAng3" << "," << "Inft1" << "," << "Supt1" << "," << "Inft2" << "," << "Supt2" << ","
                        << "Inft3" << "," << "Supt3" << "," << "InfK1" << "," << "SupK1" << "," << "InfK2" << "," << "SupK2" << endl;
            INTERVAL_VECTOR MinValue(3);
            if(!isRegular(R)){
                cout << "Not regular" << endl;
                csvfile << cameraidx << "," << PxDiam << "," << TarDiameter << "," << PxNoise << "," << TarNoise << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";        
                MinValue = ShowProjectionUnnorm(RRt,K,Xw,IMxc,T,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << "," 
                        << 0 << ",";
                INTERVAL_MATRIX SUPP,INFP,PDiff;
                if(Inf(P(3,4)) != 0 && Sup(P(3,4)) != 0){
                    SUPP = P/Inf(P(3,4));
                    INFP = P/Sup(P(3,4));
                    P = Hull(INFP,SUPP);
                }else if(Inf(P(3,4)) != 0){
                    P = P/Inf(P(3,4));
                }else if(Sup(P(3,4)) != 0){
                    P = P/Sup(P(3,4));
                }
                truecameradata[cameraidx].P = T*truecameradata[cameraidx].P;
                PDiff = P - truecameradata[cameraidx].P;
                cout << "PDiff=" << PDiff << endl;
                INTERVAL FroNorm(0);
                for(int i = 1;i <= PDiff.rows();i++){
                    FroNorm += Norm(PDiff.row(i));
                }
                csvfile << FroNorm.inf()/3 << "," << FroNorm.sup()/3 << endl;
                csvfile.close();
                continue;
            }
//            ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
            BothEulAngles = rtm2Eul(R);
            INTERVAL_VECTOR EulAngles(3);
            EulAngles = Col(BothEulAngles,1);
//            filter_points(Xw,IMxcnormal,RRt,K,Xwfiltered,xcfiltered);
            if(ColDimension(Xwfiltered) == 0)continue;
            t = LocalTriangulation(P,xcfiltered);
            if(AnyNAN(t)){
                cout << "Can't do the local triangulation" << endl;
                csvfile << cameraidx << "," << PxDiam << "," << TarDiameter << "," << PxNoise << "," << TarNoise << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
                MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";        
                MinValue = ShowProjectionUnnorm(RRt,K,Xw,IMxc,T,false);
                csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << "," 
                        << 0 << ",";
                INTERVAL_MATRIX SUPP,INFP,PDiff;
                if(Inf(P(3,4)) != 0 && Sup(P(3,4)) != 0){
                    SUPP = P/Inf(P(3,4));
                    INFP = P/Sup(P(3,4));
                    P = Hull(INFP,SUPP);
                }else if(Inf(P(3,4)) != 0){
                    P = P/Inf(P(3,4));
                }else if(Sup(P(3,4)) != 0){
                    P = P/Sup(P(3,4));
                }
                truecameradata[cameraidx].P = T*truecameradata[cameraidx].P;
                PDiff = P - truecameradata[cameraidx].P;
                INTERVAL FroNorm(0);
                for(int i = 1;i <= PDiff.rows();i++){
                    FroNorm += Norm(PDiff.row(i));
                }
                csvfile << FroNorm.inf() << "," << FroNorm.sup() << endl;
                csvfile.close();
                continue;
            }
            R = Eul2rtm(EulAngles);
            for(int i = 1;i <= 3;i++)
                for(int j = 1;j <= 3;j++)
                    RRt(i,j)=R(i,j);
            INTERVAL_VECTOR Rt;
            Rt = -R*t;
            for(int i = 1;i <= 3;i++)
                for(int j = 4;j <= 4;j++)
                    RRt(i,j)=Rt(i);
            INTERVAL_VECTOR xIAD(8);
            xIAD(1)=EulAngles(1);
            xIAD(2)=EulAngles(2);
            xIAD(3)=EulAngles(3);
            xIAD(4)=t(1);
            xIAD(5)=t(2);
            xIAD(6)=t(3);
            xIAD(7)=K(1,1);
            xIAD(8)=K(2,2);
            INTERVAL_VECTOR tmp(3);
            tmp = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
//            if(tmp(1).inf()/Xwfiltered.cols() > 1){
//                Dimension_X = 8;
//                cout << "Starting Heuristic expasion for P                             " << '\r';
//                cout.flush();
//                GAGenome::Evaluator GAFunc(GAxPXPObjFunc);
//                std::unique_ptr<FUNCTION2> OptFunction(new FUNCTION2(Dimension_X,
//                                NULL,
//                                PXFuncAD,
//                                PXFuncR,
//                                PXFuncADI,
//                                NULL,
//                                NULL,
//                                NULL,
//                                NULL));
//                Resize(optmiset.TheDomain,Dimension_X);
//                Initialize(optmiset.TheDomain,INTERVAL(-15.0*M_PI/180,15.0*M_PI/180));
//                for(int i = 4;i <= 6;i++)
//                    optmiset.TheDomain(i) = SymHull(15.0);
//                optmiset.TheDomain(7) = SymHull(5);
//                optmiset.TheDomain(8) = SymHull(5);
//                optmicalibdata.K = K;
//                AbLinsys Albl;
//                if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                    Albl.A = InitializeA_T1C(Xwfiltered,xcfiltered);
//                    optmicalibdata.xc = xcfiltered;
//                    optmicalibdata.Xw = Xwfiltered;
//                    Albl.X = xIAD;
////                    xIAD = HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset);
////                    MATRIX A = Xwfiltered.diam();
////                    VECTOR B = A.max();
////                    REAL C = B.max();
//                    tmp = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
////                    if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered) && C > 1e-15 && tmp(1).inf()/Xwfiltered.cols() > 1){
////                        optmicalibdata.xc = xcfiltered;
////                        optmicalibdata.Xw = Xwfiltered;
////                        Albl.A = InitializeA_T1C(Hull(Inf(Xwfiltered)),Hull(Mid(xcfiltered)));
////                        Albl.X = xIAD;
////                        Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
////                        tmp = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
////                        if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered) && tmp(1).inf()/Xwfiltered.cols() > 1){
////                            optmicalibdata.xc = xcfiltered;
////                            optmicalibdata.Xw = Xwfiltered;
////                            Albl.A = InitializeA_T1C(Hull(Sup(Xwfiltered)),Hull(Mid(xcfiltered)));
////                            Albl.X = xIAD;
////                            Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
////                            tmp = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
////                            if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered) && tmp(1).inf()/Xwfiltered.cols() > 1){
////                                optmicalibdata.xc = xcfiltered;
////                                optmicalibdata.Xw = Xwfiltered;
////                                Albl.A = InitializeA_T1C(Hull(Mid(Xwfiltered)),Hull(Inf(xcfiltered)));
////                                Albl.X = xIAD;
////                                Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
////                                tmp = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
////                                if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered) && tmp(1).inf()/Xwfiltered.cols() > 1){
////                                    optmicalibdata.xc = xcfiltered;
////                                    optmicalibdata.Xw = Xwfiltered;
////                                    Albl.A = InitializeA_T1C(Hull(Mid(Xwfiltered)),Hull(Sup(xcfiltered)));
////                                    Albl.X = xIAD;
////                                    Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
////                                    tmp = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
////                                    if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered) && tmp(1).inf()/Xwfiltered.cols() > 1){
////                                        optmicalibdata.xc = xcfiltered;
////                                        optmicalibdata.Xw = Xwfiltered;
////                                        Albl.A = InitializeA_T1C(Hull(Inf(Xwfiltered)),Hull(Inf(xcfiltered)));
////                                        Albl.X = xIAD;
////                                        Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
////                                        tmp = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
////                                        if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered) && tmp(1).inf()/Xwfiltered.cols() > 1){
////                                            optmicalibdata.xc = xcfiltered;
////                                            optmicalibdata.Xw = Xwfiltered;
////                                            Albl.A = InitializeA_T1C(Hull(Sup(Xwfiltered)),Hull(Sup(xcfiltered)));
////                                            Albl.X = xIAD;
////                                            Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
////                                            filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered);
////                                        }
////                                    }
////                                }
////                            }
////                        }
////                    }
//                }
//            }else{
//                
//            }
            EulAngles(1) = xIAD(1);
            EulAngles(2) = xIAD(2);
            EulAngles(3) = xIAD(3);
            cout << "EulAngles=" << EulAngles << endl;
            EulAngles(3) = Abs(EulAngles(3)) - Constant::Pi;
            INTERVAL_VECTOR TrueEulAng = Hull(VectorXd2VECTOR(truecameradata[cameraidx].TrueEulAng));
            cout << "EulAngles=" << EulAngles << " <-> ";
            cout << "TrueAngles=" << TrueEulAng << endl;
            INTERVAL_VECTOR EulAngDiff(3);
            for(int g = 1;g <= 3;g++){
                if(!(0.0 <= TrueEulAng(g))){
                    if(!(0.0 <= EulAngles(g))){
                        EulAngDiff(g).ival.inf = abs(EulAngles(g).ival.inf)-abs(TrueEulAng(g).ival.inf);
                        EulAngDiff(g).ival.sup = abs(EulAngles(g).ival.sup)-abs(TrueEulAng(g).ival.sup);
                        EulAngDiff(g) = EulAngDiff(g)/TrueEulAng(g);
                    }
                }else
                    EulAngDiff(g) = EulAngles(g);
            }
            cout << "EulAngDiff=" << EulAngDiff << endl;
            Resize(t,3);
            t(1) = xIAD(4);
            t(2) = xIAD(5);
            t(3) = xIAD(6);
            cout << "t=" << t << " <-> ";
            INTERVAL_VECTOR tDiff,Truet = Hull(VectorXd2VECTOR(truecameradata[cameraidx].EVt));
            cout << "Truet=" << Truet << endl;
            PointwiseDiv(tDiff,t-Truet,Truet);
            cout << "tDiff=" << tDiff << endl;
            K(1,1) = xIAD(7);
            K(1,1) = xIAD(8);
            K(1,3) = 0;
            K(2,3) = 0;
            cout << "K=" << K << " <-> ";
            INTERVAL_MATRIX TrueK;
            TrueK = T*MatrixXd2MATRIX(truecameradata[cameraidx].EMK);
            cout << "TrueK=" << TrueK << endl;
            INTERVAL_VECTOR KDiff(2);
            KDiff(1) = (K(1,1)+TrueK(1,1))/TrueK(1,1);
            KDiff(2) = (K(2,2)+TrueK(2,2))/TrueK(2,2);
            cout << "KDiff=" << KDiff << endl;
            YAML::Node elementnode,loadnode;
            FileOp LoadData(filenames.rrtcalc,ios::in);
            WriteNode(elementnode,EulAngles);
            loadnode = LoadData.node;
            LoadData.SaveFile.close();
            loadnode[calccameradata[cameraidx].ID][contRRtDiam]["EulAnglAAdj"].push_back(elementnode);
            elementnode.reset();
            WriteNode(elementnode,t);
            loadnode[calccameradata[cameraidx].ID][contRRtDiam]["tAAdj"].push_back(elementnode);
            elementnode.reset();
            WriteNode(elementnode,K);
            loadnode[calccameradata[cameraidx].ID][contRRtDiam]["KAAdj"].push_back(elementnode);
            elementnode.reset();
            FileOp SaveData(filenames.rrtcalc,ios::out);
            SaveData.node.reset();
            SaveData.node = loadnode;
            csvfile << cameraidx << "," << PxDiam << "," << TarDiameter << "," << PxNoise << "," << TarNoise << ",";
            MinValue = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
            filter_points(TrueXw,IMxcnormal,xIAD,Xwfiltered,xcfiltered);
            MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
            filter_points(Xw,TrueIMxcnormal,xIAD,Xwfiltered,xcfiltered);
            MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
            filter_points(TrueXw,TrueIMxcnormal,xIAD,Xwfiltered,xcfiltered);
            MinValue = ShowProjection(RRt,K,Xwfiltered,xcfiltered,false);
            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";        
            MinValue = ShowProjectionUnnorm(xIAD,Xw,IMxc,T,false);
            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << "," <<
                    1 << ",";
            csvfile << EulAngDiff(1).inf() << "," << EulAngDiff(1).sup() << "," 
                    << EulAngDiff(2).inf() << "," << EulAngDiff(2).sup() << "," 
                    << EulAngDiff(3).inf() << "," << EulAngDiff(3).sup() << ","; 
            csvfile << tDiff(1).inf() << "," << tDiff(1).sup() << ","
                    << tDiff(2).inf() << "," << tDiff(2).sup() << ","
                    << tDiff(3).inf() << "," << tDiff(3).sup() << ",";
            csvfile << KDiff(1).inf() << "," << KDiff(1).sup() << "," << KDiff(2).inf() << "," << KDiff(2).sup() << endl;
            csvfile.close();
        }
    }
//    char chcomputername[512];
//    gethostname(&chcomputername[0],512);
//    string computername;
//    computername = "-";
//    computername += chcomputername;
//    string date;
//    time_t tme = time(0);   // get time now
//    struct tm * now = localtime( & tme );
//    date = "-" + to_string(now->tm_mday) + "-" +
//    to_string(now->tm_mon + 1) + "-" +
//    to_string(now->tm_year + 1900);
//    string filename;
//    filename = RRtFileName;
//    int dotidx = filename.find_last_of(".");
//    filename.insert(dotidx,computername+date);
//    std::ifstream  src(RRtFileName, std::ios::binary);
//    std::ofstream  dst(filename,   std::ios::binary);
//    dst << src.rdbuf();
//    src.close();
//    dst.close();
}

//VOID Icalib::BundleAdjustment(INT cameraidx,string RRtFileName,string ConfigFileName){
//    if(!file_exists(ConfigFileName))return;
//    if(!file_exists(RRtFileName))return;
//    cout << endl << "Expanding P to camera: " << calccameradata[cameraidx].ID << endl;
//    YAML::Node elementnode;
//    double imgnoise,tarnoise;
//    FileOp RRtCalcData(RRtFileName,ios::in);
//    YAML::Node cameranode;
//    if(usetruedata){
//            calccameradata[cameraidx].RRt = truecameradata[cameraidx].RRt;
//    }else if(file_exists(RRtFileName)){
//        cameranode = RRtCalcData.node[calccameradata[cameraidx].ID];
//        Emitter emitter;
//        emitter << cameranode;
//        if(emitter.size() == 0)return;
//    }else{
//        cerr << "RRt file not exist" << endl;
//        return;
//    }
//    INT Dimension_X = 12;
//    Resize(optmiset.TheDomain,Dimension_X);
//    if(!RRtCalcData.IsOpen() && !calccameradata[cameraidx].RRtFromFile)return;
//    INTERVAL_MATRIX R(3,3),RRt(3,4),K(3,3);
//    INTERVAL_MATRIX BothEulAngles;
//    FileOp PixelData("PixelData.yml",ios::in),TranslationData("TranslationData.yml",ios::in);
//    int nDiam = PixelData.node.size();
//    INTERVAL_MATRIX vxcin,IMxcnormal,T;
//    bool appfirstline = true;
//    if(file_exists("resuldata.csv"))
//        appfirstline = false;
//    fstream csvfile("resuldata.csv",ios::in|ios::out|ios::app);
//    for(int contDiam = 0;contDiam < nDiam;contDiam++){
//        INTERVAL_MATRIX IMxc,IMxcnormal,TrueIMxcnormal;
//        XYZINTRelatLoad(0,0);
//        PointProjection(0,0);
//        for(int contTar = 0;contTar < targets.size();contTar++){
//            if(!ReadNode(
//             PixelData.node[contDiam][calccameradata[cameraidx].ID][targets[contTar].ID][0],
//                    0,
//                    vxcin))
//                continue;
//            IMxc.hcat(vxcin);
//            vxcin = truecameradata[cameraidx].vxci[contTar];
//            TrueIMxcnormal.hcat(vxcin);
//        }
//        Normalize2Dpoints(IMxcnormal,IMxc,T);
//        if(RowDimension(TrueIMxcnormal) < 3)
//            GetHomogeneous(TrueIMxcnormal,TrueIMxcnormal);
//        TrueIMxcnormal = T*TrueIMxcnormal;
//        GetNonHomogeneous(TrueIMxcnormal,TrueIMxcnormal);
//        int nTransl = TranslationData.node["Relative"].size();
//        double PxDiam;
//        if(!ReadNode(PixelData.node[contDiam]["Diameter"][0],0,PxDiam)){
//            cout << "Can't read cam diameter data" << endl;
//            return;
//        }
//        for(int contTrans = 0;contTrans < nTransl;contTrans++){
//            INTERVAL_MATRIX Xw,TrueXw;
//            double TarDiameter;
//            if(!ReadNode(TranslationData.node["Relative"][contTrans]["Diameter"][0],0,TarDiameter)){
//                cout << "Can't read target diameter data" << endl;
//                return;
//            }
//            INTERVAL_MATRIX IMRelativePos;
//            int contRRtDiam;
//            bool foundRRt = false;
//            for(contRRtDiam = 0;contRRtDiam < RRtCalcData.node[calccameradata[cameraidx].ID].size();contRRtDiam++){
//                elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["imgnoise"][0];
//                if(!ReadNode(elementnode,0,imgnoise)){
//                    elementnode.reset();
//                    continue;
//                }
//                elementnode.reset();
//                if((abs(PxDiam-imgnoise)) > 1e-5){
//                    continue;
//                }
//                elementnode = RRtCalcData.node[calccameradata[cameraidx].ID][contRRtDiam]["tarnoise"][0];
//                if(!ReadNode(elementnode,0,tarnoise)){
//                    elementnode.reset();
//                    continue;                        
//                }
//                elementnode.reset();
//                if((abs(TarDiameter-tarnoise)) > 1e-5){
//                    continue;
//                }
//                foundRRt = true;
//                break;
//            }
//            if(!foundRRt)continue;
//            if(!usetruedata){
//                elementnode = cameranode[contRRtDiam]["RRtFinal"];
//                if(!ReadNode(elementnode[0],0,RRt)){
//                    elementnode.reset();
//                    continue;
//                }
//                elementnode.reset();
//                elementnode = cameranode[contRRtDiam]["K"];
//                if(!ReadNode(elementnode[0],0,K))continue;
//                elementnode.reset();
//            }else{
//                imgnoise = 0;
//                tarnoise = 0;
//            }
//            cout << "imgnoise: " << imgnoise << " - ";
//            cout << "tarnoise: " << tarnoise << endl;
//            INTERVAL_MATRIX P(3,4);
//            calccameradata[cameraidx].K(1,3)=0;
//            calccameradata[cameraidx].K(2,3)=0;
//            if(usetruedata){
//                IMxc.Delete();
//                for(int contTar = 0;contTar < calccameradata[cameraidx].vxci.size();contTar++){
//                    IMxc.hcat(truecameradata[cameraidx].vxci[contTar]);
//                    Xw.hcat(targets[truecameradata[cameraidx].TargIdx(contTar)].IMRelativePos);
//                }
//                Normalize2Dpoints(IMxcnormal,IMxc,T);
//                truecameradata[cameraidx].K(1,3)=0;
//                truecameradata[cameraidx].K(2,3)=0;
//                RRt = truecameradata[cameraidx].RRt;
//                P = truecameradata[cameraidx].K*RRt;
//            }else{
//                calccameradata[cameraidx].K = K;
//                calccameradata[cameraidx].RRt = RRt;
//                for(int contTar = 0;contTar < targets.size();contTar++){
//                    if(!ReadNode(TranslationData.node["Relative"][contTrans][targets[contTar].ID][0],
//                            0,
//                            IMRelativePos))
//                        continue;
//                    Xw.hcat(IMRelativePos);
//                    IMRelativePos = MatrixXd2MATRIX(targets[contTar].RelativePos);
//                    TrueXw.hcat(IMRelativePos);
//                }
//                Normalize2Dpoints(IMxcnormal,IMxc,T);
//                P = K*RRt;
//            }
//            optmicalibdata.K = K;
//            optmicalibdata.RRt = RRt;
//            if(RowDimension(Xw) < 4)
//                GetHomogeneous(Xw,Xw);
//            optmicalibdata.Xw = Xw;
//            Resize(optmicalibdata.max_min_xc,3,2);
//            INTERVAL_VECTOR max_xc(3),min_xc(3),t;
//            min_xc = INTERVAL_VECTOR(3,{INTERVAL(-calccameradata[cameraidx].ImageSize(0)/2),
//                                        INTERVAL(-calccameradata[cameraidx].ImageSize(1)/2),
//                                        1});
//            max_xc = INTERVAL_VECTOR(3,{INTERVAL(calccameradata[cameraidx].ImageSize(0)/2),
//                                        INTERVAL(calccameradata[cameraidx].ImageSize(1)/2),
//                                        1});
//            min_xc = T*min_xc;
//            max_xc = T*max_xc;
//            SetCol(optmicalibdata.max_min_xc,1,Col(max_xc,1));
//            SetCol(optmicalibdata.max_min_xc,2,Col(min_xc,1));
//            Resize(optmicalibdata.t,3);
//            optmicalibdata.t = INTERVAL_VECTOR(3,{Hull(100),INTERVAL(-100,0),1});
//            Resize(optmicalibdata.EulAngles,3);
//            optmicalibdata.EulAngles = calccameradata[cameraidx].EulAngl;
//            R = RRt.Box(1,3,1,3);
//            BothEulAngles = rtm2Eul(R);
//            INTERVAL_VECTOR EulAngles(3);
//            EulAngles = Col(BothEulAngles,1);
//            INTERVAL_MATRIX Xwfiltered,xcfiltered;
//            filter_points(Xw,IMxcnormal,RRt,K,Xwfiltered,xcfiltered);
//            if(ColDimension(Xwfiltered) == 0)continue;
//            t = LocalTriangulation(P,xcfiltered);
//            if(AnyNAN(t)){
//                cout << "Can't do the local triangulation" << endl;
//                continue;
//            }
//            INTERVAL_VECTOR xIAD(8);
//            xIAD(1)=EulAngles(1);
//            xIAD(2)=EulAngles(2);
//            xIAD(3)=EulAngles(3);
//            xIAD(4)=t(1);
//            xIAD(5)=t(2);
//            xIAD(6)=t(3);
//            xIAD(7)=K(1,1);
//            xIAD(8)=K(2,2);
//            Dimension_X = 8;
//            cout << "Starting Heuristic expasion for P                             " << '\r';
//            cout.flush();
//            GAGenome::Evaluator GAFunc(GAxPXPObjFunc);
//            std::unique_ptr<FUNCTION2> OptFunction(new FUNCTION2(Dimension_X,
//                            NULL,
//                            PXFuncAD,
//                            PXFuncR,
//                            PXFuncADI,
//                            NULL,
//                            NULL,
//                            NULL,
//                            NULL));
//            Resize(optmiset.TheDomain,Dimension_X);
//            Initialize(optmiset.TheDomain,INTERVAL(-15.0*M_PI/180,15.0*M_PI/180));
//            for(int i = 4;i <= 6;i++)
//                optmiset.TheDomain(i) = SymHull(15.0);
//            optmiset.TheDomain(7) = SymHull(5);
//            optmiset.TheDomain(8) = SymHull(5);
//            optmicalibdata.K = K;
//            AbLinsys Albl;
//            if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                Albl.A = InitializeA_T1C(Hull(Mid(Xwfiltered)),Hull(Mid(xcfiltered)));
//                optmicalibdata.xc = xcfiltered;
//                optmicalibdata.Xw = Xwfiltered;
//                Albl.X = xIAD;
//                xIAD = HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset);
//                if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                    optmicalibdata.xc = xcfiltered;
//                    optmicalibdata.Xw = Xwfiltered;
//                    Albl.A = InitializeA_T1C(Hull(Inf(Xwfiltered)),Hull(Mid(xcfiltered)));
//                    Albl.X = xIAD;
//                    Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
//                    if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                        optmicalibdata.xc = xcfiltered;
//                        optmicalibdata.Xw = Xwfiltered;
//                        Albl.A = InitializeA_T1C(Hull(Sup(Xwfiltered)),Hull(Mid(xcfiltered)));
//                        Albl.X = xIAD;
//                        Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
//                        if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                            optmicalibdata.xc = xcfiltered;
//                            optmicalibdata.Xw = Xwfiltered;
//                            Albl.A = InitializeA_T1C(Hull(Mid(Xwfiltered)),Hull(Inf(xcfiltered)));
//                            Albl.X = xIAD;
//                            Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
//                            if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                                optmicalibdata.xc = xcfiltered;
//                                optmicalibdata.Xw = Xwfiltered;
//                                Albl.A = InitializeA_T1C(Hull(Mid(Xwfiltered)),Hull(Sup(xcfiltered)));
//                                Albl.X = xIAD;
//                                Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
//                                if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                                    optmicalibdata.xc = xcfiltered;
//                                    optmicalibdata.Xw = Xwfiltered;
//                                    Albl.A = InitializeA_T1C(Hull(Inf(Xwfiltered)),Hull(Inf(xcfiltered)));
//                                    Albl.X = xIAD;
//                                    Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
//                                    if(filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered)){
//                                        optmicalibdata.xc = xcfiltered;
//                                        optmicalibdata.Xw = Xwfiltered;
//                                        Albl.A = InitializeA_T1C(Hull(Sup(Xwfiltered)),Hull(Sup(xcfiltered)));
//                                        Albl.X = xIAD;
//                                        Hull(xIAD,HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset));
//                                        filter_points(Xwfiltered,xcfiltered,xIAD,Xwfiltered,xcfiltered);
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//            cout << "xIAD=" << xIAD << endl;
//            calccameradata[cameraidx].EulAngl(1) = xIAD(1);
//            calccameradata[cameraidx].EulAngl(2) = xIAD(2);
//            calccameradata[cameraidx].EulAngl(3) = xIAD(3);
//            Resize(t,3);
//            t(1) = xIAD(4);
//            t(2) = xIAD(5);
//            t(3) = xIAD(6);
//            calccameradata[cameraidx].K(1,1) = xIAD(7);
//            calccameradata[cameraidx].K(1,1) = xIAD(8);
//            calccameradata[cameraidx].K(1,3) = 0;
//            calccameradata[cameraidx].K(2,3) = 0;
//            R = Eul2rtm(calccameradata[cameraidx].EulAngl);
//            for(int i = 1;i <= 3;i++)
//                for(int j = 1;j <= 3;j++)
//                    RRt(i,j)=R(i,j);
//            INTERVAL_VECTOR Rt;
//            Rt = -R*t;
//            for(int i = 1;i <= 3;i++)
//                for(int j = 4;j <= 4;j++)
//                    RRt(i,j)=Rt(i);
//            YAML::Node elementnode,loadnode;
//            FileOp LoadData(RRtFileName,ios::in);
//            WriteNode(elementnode,calccameradata[cameraidx].EulAngl);
//            loadnode = LoadData.node;
//            LoadData.SaveFile.close();
//            loadnode[calccameradata[cameraidx].ID][contRRtDiam]["EulAnglAAdj"].push_back(elementnode);
//            elementnode.reset();
//            WriteNode(elementnode,t);
//            loadnode[calccameradata[cameraidx].ID][contRRtDiam]["tAAdj"].push_back(elementnode);
//            elementnode.reset();
//            WriteNode(elementnode,calccameradata[cameraidx].K);
//            loadnode[calccameradata[cameraidx].ID][contRRtDiam]["KAAdj"].push_back(elementnode);
//            elementnode.reset();
//            FileOp SaveData(RRtFileName,ios::out);
//            SaveData.node.reset();
//            SaveData.node = loadnode;
//            INTERVAL_VECTOR MinValue(3);
//            if(appfirstline)
//                csvfile << "PxDiam,TarDiam,"
//                        << "Inf(Xw-xc),Sup(Xw-xc),NBlngs,Mid(A)-Mid(B),"
//                        << "Inf(TrueXw-xc),Sup(TrueXw-xc),NBlngs,Mid(A)-Mid(B),"
//                        << "Inf(Xw-Truexc),Sup(Xw-Truexc),NBlngs,Mid(A)-Mid(B),"
//                        << "Inf(TrueXw-Truexc),Sup(TrueXw-Truexc),NBlngs,Mid(A)-Mid(B)," 
//                        << "Inf(TrueXwUnnorm-TruexcUnnorm),Sup(TrueXwUnnorm-TruexcUnnorm),NBlngs,Mid(A)-Mid(B)," << endl;
//            csvfile << PxDiam << "," << TarDiameter << ",";
//            MinValue = ShowProjection(xIAD,Xwfiltered,xcfiltered,false);
//            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
//            filter_points(TrueXw,IMxcnormal,xIAD,Xwfiltered,xcfiltered);
//            MinValue = ShowProjection(RRt,calccameradata[cameraidx].K,Xwfiltered,xcfiltered,false);
//            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
//            filter_points(Xw,TrueIMxcnormal,xIAD,Xwfiltered,xcfiltered);
//            MinValue = ShowProjection(RRt,calccameradata[cameraidx].K,Xwfiltered,xcfiltered,false);
//            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";
//            filter_points(TrueXw,TrueIMxcnormal,xIAD,Xwfiltered,xcfiltered);
//            MinValue = ShowProjection(RRt,calccameradata[cameraidx].K,Xwfiltered,xcfiltered,false);
//            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << ",";        
//            MinValue = ShowProjectionUnnorm(xIAD,Xw,IMxc,T,false);
//            csvfile << Inf(MinValue(1)) << "," << Sup(MinValue(1)) << "," << Inf(MinValue(2))/ColDimension(IMxcnormal) << "," << Sup(MinValue(3)) << endl;        
//        }
//        csvfile.close();
//    }
//    char chcomputername[512];
//    gethostname(&chcomputername[0],512);
//    string computername;
//    computername = "-";
//    computername += chcomputername;
//    string date;
//    time_t tme = time(0);   // get time now
//    struct tm * now = localtime( & tme );
//    date = "-" + to_string(now->tm_mday) + "-" +
//    to_string(now->tm_mon + 1) + "-" +
//    to_string(now->tm_year + 1900);
//    string filename;
//    filename = RRtFileName;
//    int dotidx = filename.find_last_of(".");
//    filename.insert(dotidx,computername+date);
//    std::ifstream  src(RRtFileName, std::ios::binary);
//    std::ofstream  dst(filename,   std::ios::binary);
//    dst << src.rdbuf();
//}


INTERVAL_VECTOR Icalib::LocalTriangulation(const INTERVAL_MATRIX & P,const INTERVAL_MATRIX & xc){
    vector<INTERVAL_VECTOR> p;
    for(int i = 1;i <= RowDimension(P);i++){
        p.push_back(Row(P,i));
    }
    INTERVAL_VECTOR x,L;
    INTERVAL_MATRIX A_T1C;
    INTERVAL_VECTOR b;
    Resize(A_T1C,4,3);
    Initialize(A_T1C,0);
    Resize(b,4);
    Initialize(b,0);
    INT contRowA = 1;
    INT Dimension_X = 3;
    INTERVAL_VECTOR t(3);
    bool firsttime = true;
    INT contTimes = 1;
    for(int contPts1 = 1;contPts1 <= ColDimension(xc);contPts1++){
        x = Col(xc,contPts1);
        for(int contPts2 = contPts1+1;contPts2 <= ColDimension(xc);contPts2++){
            if(contTimes > 5)break;
            contRowA = 1;
            INTERVAL_VECTOR xl;
            xl = Col(xc,contPts2);
            if(Abs(Norm(x-xl)) < 
                    0.4*Abs(Norm(Col(optmicalibdata.max_min_xc,1)-Col(optmicalibdata.max_min_xc,2))))continue;
            L = (x(1)*p[2]-p[0]);
            SetRow(A_T1C,contRowA,L.Box(1,3));
            b(contRowA) = -L(4);
            contRowA++;
            L = (x(2)*p[2]-p[1]);
            SetRow(A_T1C,contRowA,L.Box(1,3));
            b(contRowA) = -L(4);
            contRowA ++;
            L = (xl(1)*p[2]-p[0]);
            SetRow(A_T1C,contRowA,L.Box(1,3));
            b(contRowA) = -L(4);
            contRowA ++;
            L = (xl(2)*p[2]-p[1]);
            SetRow(A_T1C,contRowA,L.Box(1,3));
            b(contRowA) = -L(4);
            contRowA ++;
            LinOpt linopt(A_T1C,b);
            linopt.optmiset.MaxIter = optmiset.MaxIter;
            linopt.RohnOptimalLinSys();
            std::unique_ptr<FUNCTION2> OptFunction(new FUNCTION2(Dimension_X,
                            NULL,
                            AxXFuncAD,
                            AxXFuncADR,
                            AxXFuncADI,
                            NULL,
                            NULL,
                            NULL,
                            NULL));
            AbLinsys Albl;
            Albl.A = A_T1C;
            Albl.b = b;
            Albl.X = linopt.X;
            std::shared_ptr<VOID> userdata;
            userdata = make_shared<AbLinsys>(Albl);
            INTERVAL MinValue =  AxXFunc2ADI(linopt.X,userdata);
            GAGenome::Evaluator GAFunc(GAIniTsaiObjFunc);
            optmiset.TheDomain.clear();
            optmiset.TheDomain = Albl.X;
            if(Inf(MinValue) < 1e-1){
                GAFunc = GAGenome::Evaluator(GAIniTsaiObjFunc);
                optmiset.expansion = false;
            }else{
                GAFunc = GAGenome::Evaluator(GAxPXPObjFunc);
            }
            linopt.X = HeuristicExpansion(OptFunction,Albl,GAFunc,optmiset);
            optmiset.expansion = true;
            if(firsttime){
                t = linopt.X;
                firsttime = false;
                contTimes++;
            }else{
                t = Hull(t,linopt.X);
                contTimes++;
            }
        }
    }
    if(firsttime){
        Initialize(t,Machine::NaN);
    }
    return t;
}
INTERVAL_VECTOR Icalib::ShowProjection(const INTERVAL_MATRIX& RRt, 
        const INTERVAL_MATRIX& K, 
        const INTERVAL_MATRIX& Xw, 
        const INTERVAL_MATRIX& xcreal,
        const bool showpoints){
    INTERVAL_VECTOR out(3);
    if(Xw.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    if(xcreal.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    INTERVAL_MATRIX xccalc(RowDimension(xcreal),ColDimension(xcreal)),
            xclocal(2,ColDimension(xcreal)),
            diffxc(2,ColDimension(xcreal));
    MATRIX diffxc2(2,ColDimension(xcreal));
    INTERVAL_MATRIX localXw(4,ColDimension(Xw));
    INTERVAL_MATRIX P(3,4);
    if(RowDimension(Xw) < 4)
        GetHomogeneous(localXw,Xw);
    else
        localXw = Xw;
    P = K*RRt;
    xccalc = P*localXw;
    GetNonHomogeneous(xccalc,xccalc);
    if(RowDimension(xcreal) > 2)
        GetNonHomogeneous(xclocal,xcreal);
    else
        xclocal = xcreal;
    int totalBelongs = 0;
    if(showpoints){
        cout << "Calculated xc               Real xc            Real belongs to calc?" << endl;
        for(int j = 1;j <= ColDimension(xccalc);j++){
            cout << Col(xccalc,j) << " <-> ";
            cout << Col(xclocal,j) << " <-> ";
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                cout << (xclocal(i,j) <= xccalc(i,j)) << " - ";
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            cout << endl;
            if(contBelongs >= 2)totalBelongs++;
        }
    }else{
        for(int j = 1;j <= ColDimension(xccalc);j++){
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            if(contBelongs >= 2)totalBelongs++;
        }        
    }
    diffxc = xclocal-xccalc;
    diffxc2 = Mid(xclocal)-Mid(xccalc);
    INTERVAL normdiffxc(0,0);
    REAL normdiffxc2 = 0;
    for(int j = 1;j <= ColDimension(diffxc);j++){
        normdiffxc = normdiffxc+Norm(Col(diffxc,j));
        normdiffxc2 = normdiffxc2+Norm(Col(diffxc2,j));
    }
    cout << "Norm(diffxc)=" << Sqrt(normdiffxc) << endl;
    out(1) = Sqrt(normdiffxc);
    out(2) = totalBelongs;
    out(3) = Sqrt(normdiffxc2);
    return out;
}
INTERVAL_VECTOR Icalib::ShowProjection(INTERVAL_VECTOR xIAD, 
        const INTERVAL_MATRIX& Xw, 
        const INTERVAL_MATRIX& xcreal,
        const bool showpoints){
    INTERVAL_VECTOR out(3),EulAngl;
    if(Xw.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    if(xcreal.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    INTERVAL_MATRIX xccalc(RowDimension(xcreal),ColDimension(xcreal)),
            xclocal(2,ColDimension(xcreal)),
            diffxc(2,ColDimension(xcreal));
    MATRIX diffxc2(2,ColDimension(xcreal));
    INTERVAL_MATRIX localXw(4,ColDimension(Xw));
    INTERVAL_MATRIX P(3,4),K(3,3),RRt,t,R,Rt;
    if(RowDimension(Xw) == 0 || ColDimension(Xw) == 0){
        Initialize(out,Machine::NaN);
    }
    if(RowDimension(Xw) < 4)
        GetHomogeneous(localXw,Xw);
    else
        localXw = Xw;
    EulAngl = xIAD.Box(1,3);
    t = xIAD.Box(4,6);
    Clear(K);
    K(1,1)=xIAD(7);
    K(2,2)=xIAD(8);
    K(3,3)=1;
    R = Eul2rtm(EulAngl);
    Rt = -R*t;
    RRt = R;
    RRt.hcat(Rt);    
    P = K*RRt;
    xccalc = P*localXw;
    GetNonHomogeneous(xccalc,xccalc);
    if(RowDimension(xcreal) > 2)
        GetNonHomogeneous(xclocal,xcreal);
    else
        xclocal = xcreal;
    int totalBelongs = 0;
    if(showpoints){
        cout << "Calculated xc               Real xc            Real belongs to calc?" << endl;
        for(int j = 1;j <= ColDimension(xccalc);j++){
            cout << Col(xccalc,j) << " <-> ";
            cout << Col(xclocal,j) << " <-> ";
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                cout << (xclocal(i,j) <= xccalc(i,j)) << " - ";
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            cout << endl;
            if(contBelongs >= 2)totalBelongs++;
        }
    }else{
        for(int j = 1;j <= ColDimension(xccalc);j++){
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            if(contBelongs >= 2)totalBelongs++;
        }        
    }
    diffxc = xclocal-xccalc;
    diffxc2 = Mid(xclocal)-Mid(xccalc);
    INTERVAL normdiffxc(0,0);
    REAL normdiffxc2 = 0;
    for(int j = 1;j <= ColDimension(diffxc);j++){
        normdiffxc = normdiffxc+Norm(Col(diffxc,j));
        normdiffxc2 = normdiffxc2+Norm(Col(diffxc2,j));
    }
    cout << "Norm(diffxc)=" << Sqrt(normdiffxc) << endl;
    out(1) = Sqrt(normdiffxc);
    out(2) = totalBelongs;
    out(3) = Sqrt(normdiffxc2);
    return out;
}
INTERVAL_VECTOR Icalib::ShowProjectionUnnorm(INTERVAL_VECTOR xIAD, 
        const INTERVAL_MATRIX& Xw, 
        const INTERVAL_MATRIX& xcreal,
        INTERVAL_MATRIX& T,
        const bool showpoints){
    INTERVAL_VECTOR out(3),EulAngl;
    if(Xw.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    if(xcreal.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    INTERVAL_MATRIX xccalc(RowDimension(xcreal),ColDimension(xcreal)),
            xclocal(2,ColDimension(xcreal)),
            diffxc(2,ColDimension(xcreal));
    MATRIX diffxc2(2,ColDimension(xcreal));
    INTERVAL_MATRIX localXw(4,ColDimension(Xw));
    INTERVAL_MATRIX P(3,4),K(3,3),RRt,t,R,Rt;
    if(RowDimension(Xw) == 0 || ColDimension(Xw) == 0){
        Initialize(out,Machine::NaN);
    }
    if(RowDimension(Xw) < 4)
        GetHomogeneous(localXw,Xw);
    else
        localXw = Xw;
    EulAngl = xIAD.Box(1,3);
    t = xIAD.Box(4,6);
    Clear(K);
    K(1,1)=xIAD(7);
    K(2,2)=xIAD(8);
    K(3,3)=1;
    R = Eul2rtm(EulAngl);
    Rt = -R*t;
    RRt = R;
    RRt.hcat(Rt);
    INTERVAL_MATRIX invT;
    Inverse_Interval_Matrix(3,0,T,invT);
    P = invT*K*RRt;
    xccalc = P*localXw;
    GetNonHomogeneous(xccalc,xccalc);
    if(RowDimension(xcreal) > 2)
        GetNonHomogeneous(xclocal,xcreal);
    else
        xclocal = xcreal;
    int totalBelongs = 0;
    if(showpoints){
        cout << "Calculated xc               Real xc            Real belongs to calc?" << endl;
        for(int j = 1;j <= ColDimension(xccalc);j++){
            cout << Col(xccalc,j) << " <-> ";
            cout << Col(xclocal,j) << " <-> ";
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                cout << (xclocal(i,j) <= xccalc(i,j)) << " - ";
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            cout << endl;
            if(contBelongs >= 2)totalBelongs++;
        }
    }else{
        for(int j = 1;j <= ColDimension(xccalc);j++){
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            if(contBelongs >= 2)totalBelongs++;
        }        
    }
    diffxc = xclocal-xccalc;
    diffxc2 = Mid(xclocal)-Mid(xccalc);
    INTERVAL normdiffxc(0,0);
    REAL normdiffxc2 = 0;
    for(int j = 1;j <= ColDimension(diffxc);j++){
        normdiffxc = normdiffxc+Norm(Col(diffxc,j));
        normdiffxc2 = normdiffxc2+Norm(Col(diffxc2,j));
    }
    cout << "Norm(diffxc)=" << Sqrt(normdiffxc) << endl;
    out(1) = Sqrt(normdiffxc);
    out(2) = totalBelongs;
    out(3) = Sqrt(normdiffxc2);
    return out;
}
INTERVAL_VECTOR Icalib::ShowProjectionUnnorm(const INTERVAL_MATRIX& RRt, 
        const INTERVAL_MATRIX& K,  
        const INTERVAL_MATRIX& Xw, 
        const INTERVAL_MATRIX& xcreal,
        INTERVAL_MATRIX& T,
        const bool showpoints){
    INTERVAL_VECTOR out(3);
    if(Xw.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    if(xcreal.empty()){
        Initialize(out,SymHull(Machine::PosInfinity));
        return out;
    }
    INTERVAL_MATRIX xccalc(RowDimension(xcreal),ColDimension(xcreal)),
            xclocal(2,ColDimension(xcreal)),
            diffxc(2,ColDimension(xcreal));
    MATRIX diffxc2(2,ColDimension(xcreal));
    INTERVAL_MATRIX localXw(4,ColDimension(Xw));
    INTERVAL_MATRIX P(3,4);
    if(RowDimension(Xw) == 0 || ColDimension(Xw) == 0){
        Initialize(out,Machine::NaN);
    }
    if(RowDimension(Xw) < 4)
        GetHomogeneous(localXw,Xw);
    else
        localXw = Xw;
    INTERVAL_MATRIX invT;
    //Inverse_Interval_Matrix(3,0,T,invT);
    invT = Inverse(Mid(T));
    P = invT*K*RRt;
    xccalc = P*localXw;
    GetNonHomogeneous(xccalc,xccalc);
    if(RowDimension(xcreal) > 2)
        GetNonHomogeneous(xclocal,xcreal);
    else
        xclocal = xcreal;
    int totalBelongs = 0;
    if(showpoints){
        cout << "Calculated xc               Real xc            Real belongs to calc?" << endl;
        for(int j = 1;j <= ColDimension(xccalc);j++){
            cout << Col(xccalc,j) << " <-> ";
            cout << Col(xclocal,j) << " <-> ";
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                cout << (xclocal(i,j) <= xccalc(i,j)) << " - ";
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            cout << endl;
            if(contBelongs >= 2)totalBelongs++;
        }
    }else{
        for(int j = 1;j <= ColDimension(xccalc);j++){
            int contBelongs = 0;
            for(int i = 1;i <= 2;i++){
                if(xclocal(i,j) <= xccalc(i,j))contBelongs++;
            }
            if(contBelongs >= 2)totalBelongs++;
        }        
    }
    diffxc = xclocal-xccalc;
    diffxc2 = Mid(xclocal)-Mid(xccalc);
    INTERVAL normdiffxc(0,0);
    REAL normdiffxc2 = 0;
    for(int j = 1;j <= ColDimension(diffxc);j++){
        normdiffxc = normdiffxc+Norm(Col(diffxc,j));
        normdiffxc2 = normdiffxc2+Norm(Col(diffxc2,j));
    }
    cout << "Norm(diffxc)=" << Sqrt(normdiffxc) << endl;
    out(1) = Sqrt(normdiffxc);
    out(2) = totalBelongs;
    out(3) = Sqrt(normdiffxc2);
    return out;
}


INTERVAL_VECTOR Icalib::HeuristicExpansion(std::unique_ptr<FUNCTION2> & f,
                                            CONST AbLinsys & Albl,
                                            CONST GAGenome::Evaluator& GAFunc,
                                            CONST OPTMISET& optmiset){
    if(f)
        IGAObjFunc = std::move(f);
    CONST AbLinsys *userdata = &Albl;
// See if we've been given a seed to use (for testing purposes).  When you
// specify a random seed, the evolution will be exactly the same each time
// you use that seed number.

    unsigned int seed = 0;

// Create a phenotype for two variables.  The number of bits you can use to
// represent any number is limited by the type of computer you are using.  In
// this case, we use 16 bits to represent a floating point number whose value
// can range from -5 to 5, inclusive.  The bounds on x1 and x2 can be applied
// here and/or in the objective function.
    INTERVAL_VECTOR xout(Albl.X.nrows());
    xout = Albl.X;
//    for(int j = 1;j <= 2;j++){
        GABin2DecPhenotype map;
        for(int i = 1;i <= Albl.X.nrows();i++){
                map.add(24, 0,Sup(optmiset.TheDomain(i)));
        }
  // Create the template genome using the phenotype map we just made.
        GABin2DecGenome genome(map, GAFunc, (void *)userdata);

      // Now create the GA using the genome and run it.  We'll use sigma truncation
      // scaling so that we can handle negative objective scores.

        GASimpleGA ga(genome);
        GASigmaTruncationScaling scaling;
        ga.populationSize(optmiset.popsize);
        ga.nGenerations(optmiset.ngen);
        ga.pMutation(optmiset.pmut);
        ga.pCrossover(optmiset.pcross);
        ga.scaling(scaling);
        ga.scoreFilename("bog.dat");
        ga.scoreFrequency(10);
        ga.flushFrequency(50);
        ga.evolve(seed);
    // Dump the results of the GA to the screen.
        genome = ga.statistics().bestIndividual();
        if(optmiset.expansion){
            for(int j = 1;j <= Albl.X.nrows();j++){
                xout(j) = xout(j)+SymHull(genome.phenotype(j-1));
            }
        }else{
            for(int j = 1;j <= Albl.X.nrows();j++){
                xout(j) = Mid(xout(j))+SymHull(genome.phenotype(j-1));
            }
        }
//    }
    return xout;
}

float
GAExpObjFunc(GAGenome & c)
{
    GABin2DecGenome & genome = (GABin2DecGenome &)c;
    AbLinsys *Albl = (AbLinsys *)(c.userData());   
    INTERVAL_VECTOR xin;
    xin = Albl->X;
    float Vf = 0.0;
    INTERVAL_VECTOR x(xin.nrows());
//    x = Mid((*xin));
    INTERVAL IVf=0;
    for(int i = 1;i <= x.nrows();i++){
        x(i) = xin(i)+SymHull(genome.phenotype(i-1));
    }
    std::shared_ptr<VOID> userdata;
    userdata = make_shared<AbLinsys>((*Albl));
    IVf = Function(*IGAObjFunc,x,userdata);
    if(std::isnan(Inf(IVf)) || std::isnan(Sup(IVf))){
        return -FLT_MAX;
    }
    if(Diam(IVf) != 0){
        Vf += 1/(5*Diam(IVf));
        Vf += 1/(10*Sup(IVf));
    }else{
        if(Inf(IVf) != 0)
            Vf += 1/Inf(IVf);
        else
            Vf += 1/10*Sup(IVf);
    }
    if(Inf(IVf) == 0)Vf *= 1e1;
    return Vf;
}
void
RinvRExpFunc(INTERVAL_MATRIX& R)
{
    float Vf = 0.0;
    INTERVAL_VECTOR x(9);
    INTERVAL IVf=0;
}
bool Icalib::OptimidataConsistence(CONST INTERVAL_VECTOR& X,bool test_b,bool test_K){
}
bool filter_points(CONST INTERVAL_MATRIX& Xw,
                                CONST INTERVAL_MATRIX& xc,
                                CONST INTERVAL_MATRIX& RRt,
                                CONST INTERVAL_MATRIX& K,
                                INTERVAL_MATRIX& Xwout,
                                INTERVAL_MATRIX& xcout){
    INTERVAL_MATRIX Xwlocal,xclocal,Xwin,xcin;
    INTERVAL_MATRIX P(3,4);
    Xwin = Xw;
    xcin = xc;
    Xwout.Delete();
    xcout.Delete();
    if(RowDimension(Xw) < 4)
        GetHomogeneous(Xwlocal,Xw);
    else
        Xwlocal = Xw;
    P = K*RRt;
    xclocal = P*Xwlocal;
    GetNonHomogeneous(xclocal,xclocal);
    for(int j = 1;j <= ColDimension(xclocal);j++){
        if(AnyNAN(Col(xclocal,j)))continue;
        bool hasinf = false;
        for(int i = 1;i <= RowDimension(xclocal);i++){
            if(Inf(xclocal(i,j)) == Machine::NegInfinity || Sup(xclocal(i,j)) == Machine::PosInfinity){
                hasinf = true;
                break;
            }
        }
        if(hasinf)continue;
        Xwout.hcat(Col(Xwin,j));
        xcout.hcat(Col(xcin,j));
    }
}
bool filter_points(CONST INTERVAL_MATRIX& Xw,
                                CONST INTERVAL_MATRIX& xc,
                                INTERVAL_VECTOR& xIAD,  
                                INTERVAL_MATRIX& Xwout,
                                INTERVAL_MATRIX& xcout){
    INTERVAL_MATRIX Xwlocal,xclocal,Xwin,xcin;
    INTERVAL_MATRIX P(3,4),R,RRt,K(3,3);
    INTERVAL_VECTOR Rt,EulAngl,t;
    if(RowDimension(Xw) < 4)
        GetHomogeneous(Xwlocal,Xw);
    else
        Xwlocal = Xw;
    Xwin = Xw;
    xcin = xc;
    EulAngl = xIAD.Box(1,3);
    t = xIAD.Box(4,6);
    Clear(K);
    K(1,1)=xIAD(7);
    K(2,2)=xIAD(8);
    K(3,3)=1;
    R = Eul2rtm(EulAngl);
    Rt = -R*t;
    RRt = R;
    RRt.hcat(Rt);
    P = K*RRt;
    xclocal = P*Xwlocal;
    GetNonHomogeneous(xclocal,xclocal);
    xcout.Delete();
    Xwout.Delete();
    for(int j = 1;j <= ColDimension(xclocal);j++){
        if(AnyNAN(Col(xclocal,j)))continue;
        bool hasinf = false;
        for(int i = 1;i <= RowDimension(xclocal);i++){
            if(Inf(xclocal(i,j)) == Machine::NegInfinity || Sup(xclocal(i,j)) == Machine::PosInfinity){
                hasinf = true;
                break;
            }
        }
        if(hasinf)continue;
        Xwout.hcat(Col(Xwin,j));
        xcout.hcat(Col(xcin,j));
    }
}
void Icalib::SaveData(REAL imgnoise,REAL tarnoise){
    std::string dirname;
    string date;
    time_t tme = time(0);   // get time now
    struct tm * now = localtime( & tme );
    date = to_string(now->tm_mday) + "-" +
    to_string(now->tm_mon + 1) + "-" +
    to_string(now->tm_year + 1900);
    int convimgnoise = floor(imgnoise*10);
    int convtarnoise = floor(tarnoise*10);
    char *cdirname;
    cdirname = new char[100];
    snprintf(cdirname,100,"./Data/%s-imgnoise-%02i-tarnoise%02i",date.c_str(),convimgnoise,convtarnoise);
    mkdir(cdirname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    dirname = cdirname ;
    dirname += "/";
    string filename;
    
    filename = dirname+filenames.rrtcalc;
    std::ifstream  src(filenames.rrtcalc, std::ios::binary);
    if(src.is_open()){
        std::ofstream  dst(filename,   std::ios::binary);
        dst << src.rdbuf();
        dst.close();
        src.close();
        remove(filenames.rrtcalc.c_str());
    }
    filename = dirname+filenames.pxlocation;
    std::ifstream  src1(filenames.pxlocation, std::ios::binary);
    if(src1.is_open()){
        std::ofstream  dst1(filename,   std::ios::binary);
        dst1 << src1.rdbuf();
        dst1.close();
        src1.close();
        remove(filenames.pxlocation.c_str());
    }
    filename = dirname+filenames.csvpxlocation;
    std::ifstream  src2(filenames.csvpxlocation, std::ios::binary);
    if(src2.is_open()){
        std::ofstream  dst2(filename,   std::ios::binary);
        dst2 << src2.rdbuf();
        dst2.close();
        src2.close();
        remove(filenames.csvpxlocation.c_str());
    }
    filename = dirname+filenames.csvptlocation;
    std::ifstream  src3(filenames.csvptlocation, std::ios::binary);
    if(src3.is_open()){
        std::ofstream  dst3(filename,   std::ios::binary);
        dst3 << src3.rdbuf();
        dst3.close();
        src3.close();
        remove(filenames.csvptlocation.c_str());
    }
    
    filename = dirname+filenames.prjmatrx;
    std::ifstream  src4(filenames.prjmatrx, std::ios::binary);
    if(src4.is_open()){
        std::ofstream  dst4(filename,   std::ios::binary);
        dst4 << src4.rdbuf();
        dst4.close();
        src4.close();
        remove(filenames.prjmatrx.c_str());
    }
    filename = dirname+filenames.csvpxresult;
    std::ifstream  src5(filenames.csvpxresult, std::ios::binary);
    if(src5.is_open()){
        std::ofstream  dst5(filename,   std::ios::binary);
        dst5 << src5.rdbuf();
        dst5.close();
        src5.close();
        remove(filenames.csvpxresult.c_str());
    }
    filename = dirname+filenames.csvtriangresult;
    std::ifstream  src6(filenames.csvtriangresult, std::ios::binary);
    if(src6.is_open()){
        std::ofstream  dst6(filename,   std::ios::binary);
        dst6 << src6.rdbuf();
        dst6.close();
        src6.close();
        remove(filenames.csvtriangresult.c_str());
    }
    filename = dirname+filenames.ptlocation;
    std::ifstream  src7(filenames.ptlocation, std::ios::binary);
    if(src7.is_open()){
        std::ofstream  dst7(filename,   std::ios::binary);
        dst7 << src7.rdbuf();
        dst7.close();
        src7.close();
        remove(filenames.ptlocation.c_str());
    }
    filename = dirname+filenames.refptlocation;
    std::ifstream  src8(filenames.refptlocation, std::ios::binary);
    if(src8.is_open()){
        std::ofstream  dst8(filename,   std::ios::binary);
        dst8 << src8.rdbuf();
        dst8.close();
        src8.close();
        remove(filenames.refptlocation.c_str());
    }
}
VOID Icalib::LoadFileNames(const string& namesfile){
    FileOp fop(namesfile,ios::in);
    filenames.collada = fop.node["collada"].as<string>();
    filenames.config = fop.node["config"].as<string>();
    filenames.calcptlocation = fop.node["calcptlocation"].as<string>();
    filenames.csvptlocation = fop.node["csvptlocation"].as<string>();
    filenames.csvpxlocation = fop.node["csvpxlocation"].as<string>();
    filenames.csvpxresult = fop.node["csvpxresult"].as<string>();
    filenames.csvtriangresult = fop.node["csvtriangresult"].as<string>();
    filenames.prjmatrx = fop.node["prjmatrx"].as<string>();
    filenames.ptlocation = fop.node["ptlocation"].as<string>();
    filenames.pxlocation = fop.node["pxlocation"].as<string>();
    filenames.refptlocation = fop.node["refptlocation"].as<string>();
    filenames.rrtcalc = fop.node["rrtcalc"].as<string>();
}
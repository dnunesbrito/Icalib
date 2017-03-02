/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Icalib.h
 * Author: darlan
 *
 * Created on 24 de Maio de 2016, 09:18
 */

#ifndef ICALIB_H
#define ICALIB_H
#include <FileOp.h>
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
#include "OptFunctions.h"
#include <DOMPrint.hpp>
#include <unistd.h>
#include <poll.h>
#include <signal.h>
#include <errno.h>
#include <NLF/NLF.hpp>
#include <NLF/NLF2.hpp>
#include <math.h>
#include <ga/GASimpleGA.h>
#include <ga/GABin2DecGenome.h>
#include <ga/std_stream.h>
#include <chrono>
#include <time.h>
#include <PolinomialRoots/PolinomialRoots.hpp>
#include <opencv2/opencv.hpp>
#include <random>
#include <LinOpt/LinOpt.hpp>
#include <NonLinOpt/NonLinOpt.hpp>

struct TARGET{
    /*Todo objeto é representado por uma caixa onde a 
     * ordem das coordenadas dos lados é
    /*     7    8
    /*   5/_____/
         |     |6
      y  |     |
         | 3   | 4
         |/____|/
         1     2
      (0,0) x
     */
    INT Npoints;
    /*Pontos do alvo em coordenadas nãohomogêneas, em metros*/
    Eigen::Matrix3Xd RelativePos;
    /*Intervalo dos pontos do alvo em coordenadas nãohomogêneas, em metros*/
    INTERVAL_MATRIX IMRelativePos;
    /*Intervalo dos pontos do alvo em coordenadas nãohomogêneas, em metros
     absolutas*/
    INTERVAL_MATRIX IMAbsolutePos;
    /*Pontos do alvo em coordenadas nãohomogêneas 
     * , em metros, normalizadas [-1,1]*/
    Eigen::Matrix3Xd AbsolutePos;
    /*Identificação do alvo no blender;*/
    Eigen::Vector3d local;
    std::string ID;
};
struct CAMERA{
    bool RRtFromFile;
    bool UseSolListFile;
    bool DelListFolder;
    bool Initfixtxty;
    bool repeatsamedata;
    double imgnoise;
    double tarnoise;
    INTERVAL s;
    /*Matriz de rotação Eigen*/
    Eigen::Matrix3d EMR;
    /*Matriz intervalar de rotação*/
    INTERVAL_MATRIX R;
    /*Matriz de parâmetros intrínsecos*/
    Eigen::Matrix3d EMK;
    /*Matriz intervalar de parâmetros intrínsecos 
     * estimada*/
    INTERVAL_MATRIX K;
    /*Vetor de translação da câmera*/
    Eigen::Vector3d EVt;
    /*Vetor intervalar de translação da câmera*/
    std::vector<INTERVAL_VECTOR> t;
    /*Tamanho da imagem*/
    Eigen::Vector2i ImageSize;
    /*Vetor de matrizes 2xN contendo os pontos 
     * simulados projetados para cada alvo*/
    std::vector<Eigen::Matrix2Xd> VEMxc;
    /*Matriz intervalar 2xN contendo pontos do SCCG 
     * projetados no plano de   *
     * imagem estimado para cada alvo*/
    INTERVAL_MATRIX xc;
    /*Vetor de matrizes 3xN contendo os pontos 
     * simulados no SCCG para cada alvo*/
    std::vector<Eigen::Matrix3Xd> VEMXw;
    /*Matriz intervalar 4xN contendo pontos no SCCG 
     * estimado para cada câmera e cada ponto do
     alvo*/
    INTERVAL_MATRIX Xw;
    /*Matriz contendo os índices dos alvos para 
     * cada conjunto de pontos xccalc e Xwcalc*/
    Eigen::VectorXi TargIdx;
    /*Índices dos pontos relativos a cada alvo projeto.             
     *Cada linha equivale a um alvo e cada coluna ao *
     * índice do ponto*/
    std::vector<Eigen::VectorXi> PtsIdx;
    /*Vetor de matriz intervalar 4xN contendo pontos no SCCG 
     * estimado para cada câmera e cada ponto de cada
     alvo*/
    std::vector<INTERVAL_MATRIX> vIXw;
    /*Matriz intervalar 2xN contendo pontos do SCCG 
     * projetados no plano de   *
     * imagem estimado deslocado do centro da imagem
     para cada alvo.*/
    vector<INTERVAL_MATRIX> vxci;
    /*Identificação da câmera no blender;*/
    std::string ID;
    /*Matriz de projeção da câmera 3x4*/
    Eigen::Matrix<double,3,4> EMP;
    /*Intervalo para o valor da distância focal*/
    INTERVAL If;
    /*Intervalo para relação px/m em x*/
    INTERVAL Im_x;
    /*Intervalo para relação px/m em y*/
    INTERVAL Im_y;
    /*Matriz intervalar com os intervalos para os 
     *valores da matriz projeção*/
    INTERVAL_MATRIX P;
    /*Indica se a câmera já foi calibrada*/
    bool calibrated;
    /*Intervalo contendo o ângulo mínimo e 
     * máximo de rotação da camera para * 
     * os eixos*/
    INTERVAL_VECTOR EulAngl;
    /*Vetor contendo os ângulos da
     camera corretos*/
    Eigen::Vector3d TrueEulAng;
    /*Matriz intervalar [R -Rt] real*/
    INTERVAL_MATRIX RRt;
    /*Matriz intervalar [R -Rt] real*/
    INTERVAL_MATRIX RRtAc;
    /*Approximation list*/
    list<APPROX_ELEMENT> Approx_cam;
};
struct SETINGS{
    /*Nome do arquivo gerado com os dados da cena
     gerado pelo blender*/
    string blender_filename;
    string tar_loc_INTERVAL_filename;
    string opt_setting_filename;
};
struct AbLinsys{
    INTERVAL_MATRIX A;
    INTERVAL_VECTOR b;
    INTERVAL_VECTOR X;
    INTERVAL_MATRIX IK;
    MATRIX z;
};
struct FILENAMES{
    string config;
    string collada;
    string rrtcalc;
    string prjmatrx;
    string pxlocation;
    string ptlocation;
    string csvpxlocation;
    string csvptlocation;
    string csvpxresult;
    string csvtriangresult;
    string refptlocation;
    string calcptlocation;
    FILENAMES(){
        config = "ConfigFile.yml";
    }
};
class Icalib {
public:
    bool usetruedata;
    bool usealltar;
    bool lsq;
    double imgnoise_min,imgnoise_max,imgnoisestep;
    double tarnoise_min,tarnoise_max,tarnoisestep;
    INT Cameraidx,Targetidx;
    double imgnoise,tarnoise;
    bool findtbytriangulation;
    Icalib();
    Icalib(const Icalib& orig);
    virtual ~Icalib();
    std::vector<CAMERA> truecameradata, calccameradata;
    std::vector<TARGET> targets;
    FILENAMES filenames;
    /*Função utilizada para projetar cada ponto 
     * de cada alvo no plano de*
     *imagem. Os pontos devem ser inicializados 
     * antes da função ser     *
     *chamada.                                                          *
     *Como parâmetro de entrada a variável LoadPar 
     * indica ser serão     *
     *utilizados parâmetros padrão (0) ou parâmetros 
     * definidos pelo     *
     *usuário (1).*/
    VOID PointProjection(REAL imgdiam,REAL imgruido,bool savedata = false);
    /* Função utilizada para carregar os dados da câmera lidos do arquivos xml lido 
     do blender para o formato local*/
    VOID LoadCameraMatrix(CONST std::vector<CAMERADATA>& camdata);
    /* Função utilizada para carregar os dados dos alvos lidos do arquivos xml lido 
     do blender para o formato local*/
    VOID LoadTargetMatrix(CONST std::vector<FORMDATA>& formdata);
    VOID LoadFileNames(const string& namesfile);
    INT Ncamera;
    INT Ntarget;
    /*Função utilizada para calcular o intervalo
     no SCCG que contém o X e o Y para cada ponto
     projetado */
    VOID XYZINTAbsolLoad(const string& filename);
    /*Função utilizada para calcular o intervalo
     no SCCG que contém o X e o Y para cada ponto
     projetado */
    VOID XYZINTRelatLoad(CONST double diameter,CONST double noise,bool savedata = false);
    /*Função utilizada para realizar a calibração utilizando o método
     de Tsai.*/
    VOID CalibTsai();
    /*Função utilizada para carregar o arquivo de configurações
     o arquivo de configurações deve se chamar configcalib.yml
     e possuir os seguintes campos
     entrada_manual: [true ou false]
     arquivo_blender: <caminho completo para arquivo>*/
    VOID ReadSettings( );
    /*Função utilizada para inicializar os parâmetros
     da otimização global intervalar*/
    VOID ReadIniOptPar(CONST string& camID = string(),CONST string& tarID = string(),CONST string& fase = string());
    /*Retorna o número de cameras*/
    INT camsize(){
        return Ncamera;
    }
    /*Retorna o número de alvos*/
    INT tarsize(){
        return Ntarget;
    }
    VOID IniTsai(INT cameraidx,INT targetidx,bool Fixedtx_ty);
    void SaveData(REAL imgnoise,REAL tarnoise);
    VOID OrthonormalizeR1X_R2X(INT cameraidx);
    /*Faz a ortonormalização da matriz de rotação*/
    VOID OrthonormalizeR1X_R2XOptm(INT cameraidx,INTERVAL_MATRIX& RRt);
    VOID f_t_zest(INT cameraidx,INT targetidx,bool Fixed);
    VOID f_t_zest(INT cameraidx);
    VOID f_est(INT cameraidx,bool Fixed,double imgnoise = 0,double tarnoise = 0);
    VOID ShowCompareRRt(INT Divided,INT cameraidx,INT targetidx,INTERVAL_MATRIX RRt);
    VOID IniTsai(INT cameraidx);
//    VOID BundleAdjustment(INT cameraidx,string RRtFileName,string ConfigFileName);
    INTERVAL_MATRIX InitializeA_T1C(const INTERVAL_MATRIX& Xw, const INTERVAL_MATRIX& xc);
    VOID LinTsai(INT cameraidx,bool Fixedtx_ty,string RRtFileName,double imgnoise,double tarnoise);
    OPTMISET optmiset;
    VOID Triangulation( );
    INTERVAL_VECTOR LocalTriangulation(const INTERVAL_MATRIX & P,const INTERVAL_MATRIX & xc);
    VOID ExpansionP(INT cameraidx);
    INTERVAL_VECTOR ShowProjection(const INTERVAL_MATRIX& RRt, 
                                            const INTERVAL_MATRIX& K, 
                                            const INTERVAL_MATRIX& Xw, 
                                            const INTERVAL_MATRIX& xcreal,
                                            const bool showpoints = true);
    INTERVAL_VECTOR ShowProjection(INTERVAL_VECTOR xIAD, 
                                            const INTERVAL_MATRIX& Xw, 
                                            const INTERVAL_MATRIX& xcreal,
                                            const bool showpoints);    
    INTERVAL_VECTOR ShowProjectionUnnorm(INTERVAL_VECTOR xIAD, 
                                            const INTERVAL_MATRIX& Xw, 
                                            const INTERVAL_MATRIX& xcreal,
                                            INTERVAL_MATRIX& T,
                                            const bool showpoints);   
    INTERVAL_VECTOR ShowProjectionUnnorm(const INTERVAL_MATRIX& RRt, 
        const INTERVAL_MATRIX& K,  
        const INTERVAL_MATRIX& Xw, 
        const INTERVAL_MATRIX& xcreal,
        INTERVAL_MATRIX& T,
        const bool showpoints);    
    INTERVAL_VECTOR HeuristicExpansion(std::unique_ptr<FUNCTION2> & f,
                                            CONST AbLinsys & Albl,
                                            CONST GAGenome::Evaluator& GAFunc,
                                            CONST OPTMISET& optmiset);
    bool OptimidataConsistence(CONST INTERVAL_VECTOR& X,bool test_b = true,bool test_K = true);
    void erase();
    void cleanpointdata();
    friend bool filter_points(CONST INTERVAL_MATRIX& Xw,
                                CONST INTERVAL_MATRIX& xc,
                                CONST INTERVAL_MATRIX& RRt,
                                CONST INTERVAL_MATRIX& K,
                                INTERVAL_MATRIX& Xwout,
                                INTERVAL_MATRIX& xcout);
    friend bool filter_points(CONST INTERVAL_MATRIX& Xw,
                                CONST INTERVAL_MATRIX& xc,
                                INTERVAL_VECTOR& xIAD,  
                                INTERVAL_MATRIX& Xwout,
                                INTERVAL_MATRIX& xcout);
private:
};
bool readStdIn(std::string& line);

float
GAExpObjFunc(GAGenome & c);

void RinvRExpFunc(INTERVAL_MATRIX& R);
bool filter_points(CONST INTERVAL_MATRIX& Xw,
                                CONST INTERVAL_MATRIX& xc,
                                CONST INTERVAL_MATRIX& RRt,
                                CONST INTERVAL_MATRIX& K,
                                INTERVAL_MATRIX& Xwout,
                                INTERVAL_MATRIX& xcout);
bool filter_points(CONST INTERVAL_MATRIX& Xw,
                                CONST INTERVAL_MATRIX& xc,
                                INTERVAL_VECTOR& xIAD,  
                                INTERVAL_MATRIX& Xwout,
                                INTERVAL_MATRIX& xcout);
#endif /* ICALIB_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: darlan
 *
 * Created on 1 de Fevereiro de 2017, 19:46
 */

#include <cstdlib>
#include <Icalib.h>

using namespace std;

/*
 * 
 */
std::vector<double> imgnoise,tarnoise;
std::vector<double> imgdiam,tardiam;
void LoadNoiseDiam(const std::string& filename,std::vector<float>& noisediamvals){
    FileOp noisediamfile(filename,ios::in);
    noisediamvals.clear();
    Emitter emitter;
    emitter << noisediamfile.node;
    if(emitter.size() == 0)exit(EXIT_FAILURE);
    noisediamvals.push_back(noisediamfile.node["imgminnoise"].as<float>());
    noisediamvals.push_back(noisediamfile.node["imgmaxnoise"].as<float>());
    noisediamvals.push_back(noisediamfile.node["imgstepnoise"].as<float>());

    noisediamvals.push_back(noisediamfile.node["tarminnoise"].as<float>());
    noisediamvals.push_back(noisediamfile.node["tarmaxnoise"].as<float>());
    noisediamvals.push_back(noisediamfile.node["tarstepnoise"].as<float>());

    noisediamvals.push_back(noisediamfile.node["imgmindiam"].as<float>());
    noisediamvals.push_back(noisediamfile.node["imgmaxdiam"].as<float>());
    noisediamvals.push_back(noisediamfile.node["imgstepdiam"].as<float>());

    noisediamvals.push_back(noisediamfile.node["tarmindiam"].as<float>());
    noisediamvals.push_back(noisediamfile.node["tarmaxdiam"].as<float>());
    noisediamvals.push_back(noisediamfile.node["tarstepdiam"].as<float>());
}
void LoadData(const std::string& noisediamfile){
    double imgminnoise,imgmaxnoise,imgnoisestep;
    double tarminnoise,tarmaxnoise,tarnoisestep;
    double imgmindiam,imgmaxdiam,imgdiamstep;
    double tarmindiam,tarmaxdiam,tardiamstep;
    std::vector<float> noisediamvals;
    LoadNoiseDiam(noisediamfile,noisediamvals);
    imgminnoise = noisediamvals[0];
    imgmaxnoise = noisediamvals[1];
    imgnoisestep = noisediamvals[2];
    if(imgnoisestep == 0 && imgmaxnoise != imgminnoise){
        cerr << "O passo do ruído na imagem errado igual a zero." << endl;
    }
    tarminnoise = noisediamvals[3];
    tarmaxnoise = noisediamvals[4];
    tarnoisestep = noisediamvals[5];
    if(tarnoisestep == 0 && tarmaxnoise != tarminnoise){
        cerr << "O passo do ruído na imagem errado igual a zero." << endl;
    }
    imgmindiam = noisediamvals[6];
    imgmaxdiam = noisediamvals[7];
    imgdiamstep = noisediamvals[8];
    if(imgdiamstep == 0 && imgmaxdiam != imgmindiam){
        cerr << "O passo do ruído na imagem errado igual a zero." << endl;
    }

    tarmindiam = noisediamvals[9];
    tarmaxdiam = noisediamvals[10];
    tardiamstep = noisediamvals[11];
    if(tardiamstep == 0 && tarmaxdiam != tarmindiam){
        cerr << "O passo do ruído na imagem errado igual a zero." << endl;
    }

    if(imgminnoise > imgmaxnoise || imgmaxnoise < (imgminnoise+imgnoisestep)){
        cerr << "Error on noises limits to image data" << endl;
        return;
    }
    if(tarminnoise > tarmaxnoise || tarmaxnoise < (tarminnoise+tarnoisestep)){
        cerr << "Error on tarnoise limits" << endl;
        return;
    }
    if(imgminnoise == imgmaxnoise)
        imgnoise.push_back(imgminnoise);
    else
        for(double noise = imgminnoise;noise < imgmaxnoise+imgnoisestep;noise+=imgnoisestep)imgnoise.push_back(noise);
    if(tarminnoise == tarmaxnoise)tarnoise.push_back(tarminnoise);
    else
        for(double noise = tarminnoise;noise < tarmaxnoise+imgnoisestep;noise+=tarnoisestep)tarnoise.push_back(noise);

    if(imgmindiam == imgmaxdiam)
        imgdiam.push_back(imgmindiam);
    else
        for(double diam = imgmindiam;diam < imgmaxdiam+imgdiamstep;diam+=imgdiamstep)imgdiam.push_back(diam);
    if(tarmindiam == tarmaxdiam)tardiam.push_back(tarmindiam);
    else
        for(double diam = tarmindiam;diam < tarmaxdiam+tardiamstep;diam+=tardiamstep)tardiam.push_back(diam);
}
int main(int argc, char** argv) {
    string filename,configfilename,noisediamfile,filenames;
    if(argc > 1){
        for(int i = 0;argv[i];i++){
            if(strcmp(argv[i],"-colladafilename") == 0){
                i++;
                filename = argv[i];
            }
            if(strcmp(argv[i],"-configfilename") == 0){
                i++;
                configfilename = argv[i];
            }
            if(strcmp(argv[i],"-noisediamfile") == 0){
                i++;
                noisediamfile = argv[i];
                cout << noisediamfile << endl;
            }
            if(strcmp(argv[i],"-filenames") == 0){
                i++;
                filenames = argv[i];
            }
        }
    }else{
        cout << "argv deu errado" << endl;
        return -1;
    }
    int MaxTar;
    LoadData(noisediamfile);
    for(int contimgnoise = 0;contimgnoise < imgnoise.size();contimgnoise++){
        for(int conttarnoise = 0;conttarnoise < tarnoise.size();conttarnoise++){
            Icalib calibdata;
            calibdata.LoadFileNames(filenames);
            FileOp confnode(calibdata.filenames.config,ios::in);
            if(!confnode.IsOpen()) return EXIT_FAILURE;
            std::vector<CAMERADATA> camdata;
            std::vector<FORMDATA> formdata;
            XMLAcces(calibdata.filenames.collada,camdata,formdata);
            calibdata.LoadCameraMatrix(camdata);
            calibdata.LoadTargetMatrix(formdata);
            calibdata.usetruedata = false;
            if(file_exists(calibdata.filenames.ptlocation))
                remove(calibdata.filenames.ptlocation.c_str());
            if(file_exists(calibdata.filenames.pxlocation))
                remove(calibdata.filenames.pxlocation.c_str());
            for(int contimgdiam = 0;contimgdiam < imgdiam.size();contimgdiam++){
                calibdata.PointProjection(imgdiam[contimgdiam],imgnoise[contimgnoise],true);
            }
            for(int conttardiam = 0;conttardiam < tardiam.size();conttardiam++){
                calibdata.XYZINTRelatLoad(tardiam[conttardiam],tarnoise[conttarnoise],true);
            }
            for(int etapa = 0;etapa < 2;etapa++){
                for(int contCam = 0;contCam < calibdata.camsize();contCam++){
                    MaxTar = 1;
                    for(int contTar = 0;contTar < MaxTar;contTar++){
                        YAML::Node tmpnode = confnode.node[etapa][calibdata.calccameradata[contCam].ID]
                                                        [calibdata.targets[calibdata.calccameradata[contCam].TargIdx(contTar)].ID];
                        YAML::Emitter emitter;
                        emitter << tmpnode;
                        if(emitter.size() == 0)continue;
                        calibdata.optmiset.MaxIter = tmpnode["MaxIter"].as<int>();
                        calibdata.optmiset.Iterations = tmpnode["Iterations"].as<int>();
                        calibdata.optmiset.Eps = tmpnode["Eps"].as<double>();
                        calibdata.optmiset.ngen = tmpnode["ngen"].as<int>();
                        calibdata.optmiset.popsize = tmpnode["popsize"].as<int>();
                        if(etapa == 0){
                            calibdata.IniTsai(contCam);
                        }else if(etapa == 1){
                            calibdata.f_t_zest(contCam);                            
                            calibdata.ExpansionP(contCam);
                        }
                        continue;
                    }
                }
            }
            calibdata.SaveData(imgnoise[contimgnoise],tarnoise[conttarnoise]);
        }
    }
    cout << "Calibration finish" << endl;
    return 0;
}

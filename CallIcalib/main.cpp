#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for fork */
#include <sys/types.h> /* for pid_t */
#include <sys/wait.h> /* for wait */
#include <Icalib.h>
#include <sys/wait.h>
#include <sys/syscall.h>
#include <math.h>

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
        for(double noise = tarminnoise;noise < tarmaxnoise+tarnoisestep;noise+=tarnoisestep)tarnoise.push_back(noise);

    if(imgmindiam == imgmaxdiam)
        imgdiam.push_back(imgmindiam);
    else
        for(double diam = imgmindiam;diam < imgmaxdiam+imgdiamstep;diam+=imgdiamstep)imgdiam.push_back(diam);
    if(tarmindiam == tarmaxdiam)tardiam.push_back(tarmindiam);
    else
        for(double diam = tarmindiam;diam < tarmaxdiam+tardiamstep;diam+=tardiamstep)tardiam.push_back(diam);
}
void modififilenames(const string& namefiles,float imgnoise,float tarnoise){
    FileOp fop(namefiles,ios::in);
    string outfilename = namefiles;
    int Idx = outfilename.find_last_of(".");
    int iimgnoise = ceil(imgnoise);
    int itarnoise = ceil(tarnoise*10);
    string strnoises = to_string(iimgnoise) + "_" + to_string(itarnoise);
    outfilename.insert(Idx,strnoises);
    outfilename = outfilename;
    FileOp fopout(outfilename,ios::out);
    vector<string> filenames(12);
    filenames[0] = fop.node["collada"].as<string>();
    filenames[1] = fop.node["calcptlocation"].as<string>();
    filenames[2] = fop.node["config"].as<string>();
    filenames[3] = fop.node["csvptlocation"].as<string>();
    filenames[4] = fop.node["csvpxlocation"].as<string>();
    filenames[5] = fop.node["csvpxresult"].as<string>();
    filenames[6] = fop.node["csvtriangresult"].as<string>();
    filenames[7] = fop.node["prjmatrx"].as<string>();
    filenames[8] = fop.node["ptlocation"].as<string>();
    filenames[9] = fop.node["pxlocation"].as<string>();
    filenames[10] = fop.node["refptlocation"].as<string>();
    filenames[11] = fop.node["rrtcalc"].as<string>();
    for(int i = 1;i <= 11;i++){
        if(i == 2)continue;
        string strnoises;
        Idx = filenames[i].find_last_of(".");
        strnoises = to_string(iimgnoise) + "_" + to_string(itarnoise);
        filenames[i].insert(Idx,strnoises);
    }
    fopout.node["collada"] = filenames[0];
    fopout.node["calcptlocation"] = filenames[1];
    fopout.node["config"] = filenames[2];
    fopout.node["csvptlocation"] = filenames[3];
    fopout.node["csvpxlocation"] = filenames[4];
    fopout.node["csvpxresult"] = filenames[5];
    fopout.node["csvtriangresult"] = filenames[6];
    fopout.node["prjmatrx"] = filenames[7];
    fopout.node["ptlocation"] = filenames[8];
    fopout.node["pxlocation"] = filenames[9];
    fopout.node["refptlocation"] = filenames[10];
    fopout.node["rrtcalc"] = filenames[11];
}
int main(int argc, char** argv)
{
    string filename,configfilename,noisediamfile,filenames;
    int nprocess = 3;
    if(argc > 1){
        for(int i = 0;argv[i];i++){
            if(strcmp(argv[i],"-colladafilename") == 0){
                i++;
                filename = argv[i];
                continue;
            }
            if(strcmp(argv[i],"-configfilename") == 0){
                i++;
                configfilename = argv[i];
                continue;
            }
            if(strcmp(argv[i],"-noisediamfile") == 0){
                i++;
                noisediamfile = argv[i];
                continue;
            }
            if(strcmp(argv[i],"-filenames") == 0){
                i++;
                filenames = argv[i];
                continue;
            }
            if(strcmp(argv[i],"-nprocess") == 0){
                i++;
                nprocess = atoi(argv[i]);
            }
        }
    }else{
        cout << "argv deu errado" << endl;
        return -1;
    }
    LoadData(noisediamfile);
    pid_t pid[5];
    int continiprocess = 0;
    string noisefile;
    int pidx = 0;
    int parent = 0;
    string outfilename;
    int status;
    for(int contimgnoise = 0;contimgnoise < imgnoise.size();contimgnoise++){
        for(int conttarnoise = 0;conttarnoise < tarnoise.size();conttarnoise++){
            pid[pidx] = fork();
            if(pid[pidx] == 0){
                if(!noisefile.empty()){
                    char noisefilename[200];
                    sprintf(noisefilename,"%s",noisefile.c_str());
                    char chfilenames[200];
                    sprintf(chfilenames,"%s",outfilename.c_str());
                    static char *argvi[]={"-filenames", chfilenames,"-noisediamfile",noisefilename,NULL};
                    execv("icalibprg",argvi);
                    exit(EXIT_FAILURE);
                }
                return EXIT_SUCCESS;
            }
            modififilenames(filenames,imgnoise[contimgnoise],tarnoise[conttarnoise]);
            outfilename = filenames;
            int Idx = outfilename.find_last_of(".");
            int iimgnoise = floor(imgnoise[contimgnoise]);
            int itarnoise = floor(tarnoise[conttarnoise]*10);
            string strnoises = to_string(iimgnoise) + "_" + to_string(itarnoise);
            outfilename.insert(Idx,strnoises);
            noisefile = noisediamfile;
            Idx = noisefile.find_last_of(".");
            noisefile.insert(Idx,strnoises);
            string tmpsave = noisefile;
            FileOp noisediamsave(tmpsave,ios::out);
            FileOp noisediam(noisediamfile,ios::in);
            noisediamsave.node = noisediam.node;
                    noisediamsave.node["imgminnoise"] = imgnoise[contimgnoise];
                    noisediamsave.node["imgmaxnoise"] = imgnoise[contimgnoise];
            noisediamsave.node["imgstepnoise"] = 0;

                    noisediamsave.node["tarminnoise"] = tarnoise[conttarnoise];
                    noisediamsave.node["tarmaxnoise"] = tarnoise[conttarnoise];
            noisediamsave.node["tarstepnoise"] = 0;
            noisediamsave.close();
            noisediam.close();
            bool stop = false;
            if(continiprocess >= nprocess){
                while(!stop){
                    for(int idxp = 0;idxp < nprocess;idxp++){
                        pid_t result = waitpid(pid[idxp],&status,WNOHANG);
                        if(result == -1){
                            pidx = idxp;
                            stop = true;
                            break;
                        }else if(result > 0){
                            pidx = idxp;
                            stop = true;
                            break;
                        }
                    }
                }
            }else{
                continiprocess++;
                if(continiprocess > 1)
                    pidx++;
            }
        }
    }
    wait(NULL);
    pid[0] = fork();
    if(pid[0] == 0){
        return EXIT_SUCCESS;
    }else{
        wait(NULL);
    }
    return 0;
}



/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GTKFileOp.hpp"

void LoadGtkListStore(const std::string& filename,
                        GtkListStore *liststore,
                        GtkListStore *colunmnames,
                        int CameIDColumn,
                        int TargetIDColumn,
                        int EtapaColumn){
    GtkTreeIter iter,colnameiter;
    FileOp filelist(filename,ios::in);
    if(!filelist.IsOpen())return;
    YAML::Emitter emit;
    emit << filelist.node;
    if(emit.size() == 0)return;
    gint ncolumns = gtk_tree_model_get_n_columns(GTK_TREE_MODEL(liststore));
    gint nrows = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(liststore),NULL);
    string strcam,strtarg,strcolname;
    string strrow,strcol;
    gint etapa;
    for(int i = 0;i < nrows;){
        GValue value = G_VALUE_INIT;
        GValue gvcam = G_VALUE_INIT;
        GValue gvtar = G_VALUE_INIT;
        GValue gvetapa = G_VALUE_INIT;
        strrow=to_string(i);
        gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(liststore),&iter,strrow.c_str());
        gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,CameIDColumn,&gvcam);
        strcam = g_value_get_string(&gvcam);
        g_value_reset(&gvcam);
        g_value_unset(&gvcam);
        gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,TargetIDColumn,&gvtar);
        strtarg = g_value_get_string(&gvtar);
        g_value_reset(&gvtar);
        g_value_unset(&gvtar);
        gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,EtapaColumn,&gvetapa);
        etapa = g_value_get_int(&gvetapa);
        g_value_reset(&gvetapa);
        g_value_unset(&gvetapa);
        bool hasline = false;
        for(int j = 2;j < ncolumns;j++){
            bool hascolumn = false;
            GValue gcolname = G_VALUE_INIT;
            strcol = to_string(j);
            if(gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(colunmnames),&colnameiter,strcol.c_str())){
                gtk_tree_model_get_value(GTK_TREE_MODEL(colunmnames),&colnameiter,0,&gcolname);
                strcolname = g_value_get_string(&gcolname);
            }
            YAML::Emitter emitt;
            emitt << filelist.node[etapa][strcam][strtarg][strcolname];
            if(emitt.size() == 0)continue;
            string s = filelist.node[etapa][strcam][strtarg][strcolname].as<string>();
            bool isint = (s.find_first_not_of( "0123456789" ) == string::npos);
            bool isdouble = (s.find_first_not_of( "0123456789.Ee-" ) == string::npos);
            if(strcmp(strcolname.c_str(),"Etapa") == 0)continue;
            if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),j) == G_TYPE_STRING && !isdouble && !isint){
                YAML::Emitter emitter;
                emitter << filelist.node[etapa][strcam][strtarg][strcolname];
                if(emitter.size() != 0){
                    g_value_init(&value,G_TYPE_STRING);
                    g_value_set_string(&value,filelist.node[etapa][strcam][strtarg][strcolname].as<string>().c_str());
                    hasline = true;
                    hascolumn = true;
                }else{
                    continue;
                }
            }
            if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),j) == G_TYPE_DOUBLE && isdouble){
                YAML::Emitter emitter;
                emitter << filelist.node[etapa][strcam][strtarg][strcolname];
                if(emitter.size() != 0){
                    g_value_init(&value,G_TYPE_DOUBLE);
                    g_value_set_double(&value,filelist.node[etapa][strcam][strtarg][strcolname].as<double>());
                    hasline = true;
                    hascolumn = true;
                }else{
                    continue;
                }

            }
            if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),j) == G_TYPE_INT && isint){
                YAML::Emitter emitter;
                emitter << filelist.node[etapa][strcam][strtarg][strcolname];
                if(emitter.size() != 0){
                    g_value_init(&value,G_TYPE_INT);
                    g_value_set_int(&value,filelist.node[etapa][strcam][strtarg][strcolname].as<int>());
                    hasline = true;
                    hascolumn = true;
                }else{
                    continue;
                }

            }
            if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),j) == G_TYPE_BOOLEAN && isint){
                YAML::Emitter emitter;
                emitter << filelist.node[etapa][strcam][strtarg][strcolname];
                if(emitter.size() != 0){
                    g_value_init(&value,G_TYPE_BOOLEAN);
                    g_value_set_boolean(&value,filelist.node[etapa][strcam][strtarg][strcolname].as<int>());
                    hasline = true;
                    hascolumn = true;
                }else{
                    continue;
                }

            }
            if(hascolumn){
                gtk_list_store_set_value(liststore,&iter,j,&value);
                g_value_unset(&value);
            }
        }
        if(!hasline){
            gtk_list_store_remove(liststore,&iter);
            nrows--;
            continue;
        }
        i++;
    }
}
void LoadGtkListStoreIDs(const std::string& filename,
                            GtkListStore *liststore,
                            int CameIDColumn,
                            int TargetIDColumn,
                            int EtapaColumn){
    GtkTreeIter iter;
    gtk_list_store_clear(liststore);
    std::vector<CAMERADATA> camdata;
    std::vector<FORMDATA> targets;
    XMLAcces(filename,camdata,targets);
    gint ncam = camdata.size();
    gint ntarg = targets.size();
    GValue value = G_VALUE_INIT;
    GValue gvetapa = G_VALUE_INIT;
    for(int etapa = 0;etapa < 2;etapa++){
        for(int contCam = 0;contCam < ncam;contCam++){
            for(int contTarg = 0;contTarg < ntarg;contTarg++){
                gtk_list_store_append (liststore, &iter);
                g_value_init(&value,G_TYPE_STRING);
                g_value_set_string(&value,camdata[contCam].ID.c_str());
                gtk_list_store_set_value(liststore,&iter,CameIDColumn,&value);
                g_value_set_string(&value,targets[contTarg].ID.c_str());
                gtk_list_store_set_value(liststore,&iter,TargetIDColumn,&value);
                g_value_unset(&value);
                g_value_init(&gvetapa,G_TYPE_INT);
                g_value_set_int(&gvetapa,etapa);
                gtk_list_store_set_value(liststore,&iter,EtapaColumn,&gvetapa);
                g_value_unset(&gvetapa);
            }
        }
    }
}

void SaveGtkListStore(CONST std::string& filename,GtkListStore *liststore,GtkListStore *columnnames){
    gint ncolumns = gtk_tree_model_get_n_columns(GTK_TREE_MODEL(liststore));
    GtkTreeIter iter,columnnameiter;
    FileOp filelist(filename,ios::out);
    if(!filelist.IsOpen())return;
    string strcam,strtarg,strcolname,strvalue;
    double doubcell;
    int intcell,etapa;
    if(!gtk_tree_model_get_iter_first(GTK_TREE_MODEL(liststore),&iter))
        return;
    do{
        GValue gvcam = G_VALUE_INIT;
        GValue gvtar = G_VALUE_INIT;
        GValue gvetapa = G_VALUE_INIT;
        if(!gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(columnnames),&columnnameiter,"2"))
            return;
        gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,0,&gvcam);
        strcam = g_value_get_string(&gvcam);
        gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,1,&gvtar);
        strtarg = g_value_get_string(&gvtar);
        gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,14,&gvetapa);
        etapa = g_value_get_int(&gvetapa);
        for(int i = 2;i < ncolumns;i++){
            GValue value = G_VALUE_INIT;
            GValue gcolname = G_VALUE_INIT;
            gtk_tree_model_get_value(GTK_TREE_MODEL(columnnames),&columnnameiter,0,&gcolname);
            strcolname = g_value_get_string(&gcolname);
            if(strcmp(strcolname.c_str(),"Etapa") == 0)continue;
            gtk_tree_model_get_value(GTK_TREE_MODEL(liststore),&iter,i,&value);
            if(G_VALUE_HOLDS_STRING(&value)){
                strvalue = g_value_get_string(&value);
                filelist.node[etapa][strcam][strtarg][strcolname] = strvalue;
            }
            if(G_VALUE_HOLDS_DOUBLE(&value)){
                doubcell = g_value_get_double(&value);
                filelist.node[etapa][strcam][strtarg][strcolname] = doubcell;
            }
            if(G_VALUE_HOLDS_INT(&value)){
                intcell = g_value_get_int(&value);
                filelist.node[etapa][strcam][strtarg][strcolname] = intcell;
            }
            if(G_VALUE_HOLDS_BOOLEAN(&value)){
                intcell = g_value_get_boolean(&value);
                filelist.node[etapa][strcam][strtarg][strcolname] = intcell;
            }
            if(!gtk_tree_model_iter_next(GTK_TREE_MODEL(columnnames),&columnnameiter))
                break;
            g_value_reset(&value);
            g_value_reset(&gcolname);
            g_value_unset(&value);
            g_value_unset(&gcolname);
        }
        g_value_reset(&gvcam);
        g_value_reset(&gvtar);
        g_value_reset(&gvetapa);
        g_value_unset(&gvcam);
        g_value_unset(&gvtar);
        g_value_unset(&gvetapa);
    }while(gtk_tree_model_iter_next(GTK_TREE_MODEL(liststore),&iter));
}
void SaveNoiseDiam(const std::string& filename,std::vector<float>& noisediamvals){
    FileOp noisediamfile(filename,ios::out);
    noisediamfile.node["imgminnoise"] = noisediamvals[0];
    noisediamfile.node["imgmaxnoise"] = noisediamvals[1];
    noisediamfile.node["imgstepnoise"] = noisediamvals[2];

    noisediamfile.node["tarminnoise"] = noisediamvals[3];
    noisediamfile.node["tarmaxnoise"] = noisediamvals[4];
    noisediamfile.node["tarstepnoise"] = noisediamvals[5];

    noisediamfile.node["imgmindiam"] = noisediamvals[6];
    noisediamfile.node["imgmaxdiam"] = noisediamvals[7];
    noisediamfile.node["imgstepdiam"] = noisediamvals[8];

    noisediamfile.node["tarmindiam"] = noisediamvals[9];
    noisediamfile.node["tarmaxdiam"] = noisediamvals[10];
    noisediamfile.node["tarstepdiam"] = noisediamvals[11];
}
void LoadNoiseDiam(const std::string& filename,std::vector<float>& noisediamvals){
    FileOp noisediamfile(filename,ios::in);
    noisediamvals.clear();
    
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
void SaveData(REAL imgnoise,REAL tarnoise){
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
    
    filename = dirname+"RRtCalc.yml";
    std::ifstream  src("RRtCalc.yml", std::ios::binary);
    if(src.is_open()){
        std::ofstream  dst(filename,   std::ios::binary);
        dst << src.rdbuf();
        dst.close();
        src.close();
        remove("RRtCalc.yml");
    }
    filename = dirname+"PixelData.yml";
    std::ifstream  src1("PixelData.yml", std::ios::binary);
    if(src1.is_open()){
        std::ofstream  dst1(filename,   std::ios::binary);
        dst1 << src1.rdbuf();
        dst1.close();
        src1.close();
        remove("PixelData.yml");
    }
    filename = dirname+"pixellocation.csv";
    std::ifstream  src2("pixellocation.csv", std::ios::binary);
    if(src2.is_open()){
        std::ofstream  dst2(filename,   std::ios::binary);
        dst2 << src2.rdbuf();
        dst2.close();
        src2.close();
        remove("pixellocation.csv");
    }
    filename = dirname+"pointlocation.csv";
    std::ifstream  src3("pointlocation.csv", std::ios::binary);
    if(src3.is_open()){
        std::ofstream  dst3(filename,   std::ios::binary);
        dst3 << src3.rdbuf();
        dst3.close();
        src3.close();
        remove("pointlocation.csv");
    }
    
    filename = dirname+"projectionmatrix.csv";
    std::ifstream  src4("projectionmatrix.csv", std::ios::binary);
    if(src4.is_open()){
        std::ofstream  dst4(filename,   std::ios::binary);
        dst4 << src4.rdbuf();
        dst4.close();
        src4.close();
        remove("projectionmatrix.csv");
    }
    filename = dirname+"resuldata.csv";
    std::ifstream  src5("resuldata.csv", std::ios::binary);
    if(src5.is_open()){
        std::ofstream  dst5(filename,   std::ios::binary);
        dst5 << src5.rdbuf();
        dst5.close();
        src5.close();
        remove("resuldata.csv");
    }
    filename = dirname+"TriangulationResult.csv";
    std::ifstream  src6("TriangulationResult.csv", std::ios::binary);
    if(src6.is_open()){
        std::ofstream  dst6(filename,   std::ios::binary);
        dst6 << src6.rdbuf();
        dst6.close();
        src6.close();
        remove("TriangulationResult.csv");
    }
    filename = dirname+"TranslatioData.yml";
    std::ifstream  src7("TranslatioData.yml", std::ios::binary);
    if(src7.is_open()){
        std::ofstream  dst7(filename,   std::ios::binary);
        dst7 << src7.rdbuf();
        dst7.close();
        src7.close();
        remove("TranslatioData.yml");
    }
    filename = dirname+"TranslationData.yml";
    std::ifstream  src8("TranslationData.yml", std::ios::binary);
    if(src8.is_open()){
        std::ofstream  dst8(filename,   std::ios::binary);
        dst8 << src8.rdbuf();
        dst8.close();
        src8.close();
        remove("TranslationData.yml");
    }
    filename = dirname+"camdata.csv";
    std::ifstream  src9("camdata.csv", std::ios::binary);
    if(src9.is_open()){
        std::ofstream  dst9(filename,   std::ios::binary);
        dst9 << src9.rdbuf();
        dst9.close();
        src9.close();
        remove("camdata.csv");
    }
}
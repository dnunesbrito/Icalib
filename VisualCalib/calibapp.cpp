/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//glib-compile-schemas ./
//glib-compile-resources calibapp.gresource.xml --target=resources.c --generate-source
#include <gtk/gtk.h>
#include <glib-2.0/gobject/gtype.h>

#include "calibapp.hpp"
#include "calibappwin.hpp"
#include "OptSettings.hpp"

struct _CalibApp
{
  GtkApplication parent;
};

string filename;

struct _CalibAppClass
{
  GtkApplicationClass parent_class;
};

G_DEFINE_TYPE(CalibApp, calib_app, GTK_TYPE_APPLICATION);
static std::vector<string> colnames;
static string RRtFileName;

static void
calib_app_init (CalibApp *app)
{
}


static void
calib_app_startup (GApplication *app)
{

}
static void
calib_app_activate (GApplication *app)
{
  GList *windows;
  CalibAppWindow *win;
  windows = gtk_application_get_windows (GTK_APPLICATION (app));
  if (windows)
    win = CALIB_APP_WINDOW (windows->data);
  else
    win = calib_app_window_new (CALIB_APP (app));
  colnames = calib_app_window_get_columnnames(win);
  gtk_window_present (GTK_WINDOW (win));
}

static void
calib_app_open (GApplication  *app,
                  GFile        **files,
                  gint           n_files,
                  const gchar   *hint)
{
  GList *windows;
  CalibAppWindow *win;
  int i;
  windows = gtk_application_get_windows (GTK_APPLICATION (app));
  if (windows)
    win = CALIB_APP_WINDOW (windows->data);
  else
    win = calib_app_window_new (CALIB_APP (app));
for (i = 0; i < n_files; i++)
    calib_app_window_open (win, files[i]);
}

static void
calib_app_class_init (CalibAppClass *Class)
{
  G_APPLICATION_CLASS (Class)->activate = calib_app_activate;
  G_APPLICATION_CLASS (Class)->open = calib_app_open;
}

CalibApp *
calib_app_new (void)
{
  return (CalibApp *)g_object_new (CALIB_APP_TYPE,
                       "application-id", "org.gtk.calibapp",
                       "flags", G_APPLICATION_HANDLES_OPEN,
                       NULL);
}

extern "C"{
    void
    on_bttnLoadCalibData_clicked(GtkButton *button,
               GtkWidget   *entrfile){
//        if(strcmp(gtk_widget_get_name(entrfile),"chcktruedata") == 0)
//            calibdata.usetruedata = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(entrfile));
//        if(strcmp(gtk_widget_get_name(entrfile),"chckbttnalltar") == 0)
//            calibdata.usealltar = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(entrfile));
//        if(strcmp(gtk_widget_get_name(entrfile),"chckbttnlinlsq") == 0){
//            calibdata.lsq = !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(entrfile));
//        }
        if(strcmp(gtk_widget_get_name(entrfile),"camfile") == 0){
            filename = gtk_entry_get_text(GTK_ENTRY(entrfile));
        }
        if(strcmp(gtk_widget_get_name(entrfile),"translation") == 0){
//            filename = gtk_entry_get_text(GTK_ENTRY(entrfile));
//            calibdata.XYZINTAbsolLoad(filename);
        }
        if(strcmp(gtk_widget_get_name(entrfile),"rrtfilename") == 0){
            RRtFileName = gtk_entry_get_text(GTK_ENTRY(entrfile));
        }
        if(strcmp(gtk_widget_get_name(entrfile),"BttnCalib") == 0){
            gtk_widget_set_sensitive(entrfile,true);
        }
        if(strcmp(gtk_widget_get_name(entrfile),"bttbundeladj") == 0){
            gtk_widget_set_sensitive(entrfile,true);
        }
        if(strcmp(gtk_widget_get_name(entrfile),"bttntriangulation") == 0){
            gtk_widget_set_sensitive(entrfile,true);
            gtk_button_set_label(button,"Dados de calibração carregados");
        }
    }
    void
    Message_end_calib(GtkButton *button,
               CalibAppWindow   *win){
        GtkWidget *dialog;
        GtkDialogFlags flags = GTK_DIALOG_MODAL;
        dialog = gtk_message_dialog_new (GTK_WINDOW(win),
                                 flags,
                                 GTK_MESSAGE_OTHER,
                                 GTK_BUTTONS_OK,
                                 "%s",
                                 "Finish Calibration");
        gtk_dialog_run (GTK_DIALOG(dialog));
        gtk_widget_destroy (dialog);
    }
}
void  CalibCameras(std::vector<double>& imgnoise,std::vector<double>& tarnoise,
                    std::vector<double>& imgdiam,std::vector<double>& tardiam,
                    const string& ConfigFileName,bool usealltar){
    FileOp confnode(ConfigFileName,ios::in);
    if(!confnode.IsOpen())return;
    int MaxTar;
    for(int contimgnoise = 0;contimgnoise < imgnoise.size();contimgnoise++){
        for(int conttarnoise = 0;conttarnoise < tarnoise.size();conttarnoise++){
            Icalib calibdata;
            std::vector<CAMERADATA> camdata;
            std::vector<FORMDATA> formdata;
            XMLAcces(filename,camdata,formdata);
            calibdata.LoadCameraMatrix(camdata);
            calibdata.LoadTargetMatrix(formdata);
            calibdata.usetruedata = false;
                if(file_exists("TranslationData.yml"))
                    remove("TranslationData.yml");
                if(file_exists("PixelData.yml"))
                    remove("PixelData.yml");
                for(int contimgdiam = 0;contimgdiam < imgdiam.size();contimgdiam++){
                    calibdata.PointProjection(imgdiam[contimgdiam],imgnoise[contimgnoise],true);
                }
                for(int conttardiam = 0;conttardiam < tardiam.size();conttardiam++){
                    calibdata.XYZINTRelatLoad(tardiam[conttardiam],tarnoise[conttarnoise],"TranslationData.yml",true);
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
                        for(int c = 0;c < colnames.size();c++){
//                            if(strcmp(colnames[c].c_str(),"Iterations") == 0)
//                                calibdata.optmiset.Iterations = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"MaxBranches") == 0)
//                                calibdata.optmiset.MaxBranches = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"Eps") == 0)
//                                calibdata.optmiset.Eps = tmpnode[colnames[c]].as<double>();
//                            if(strcmp(colnames[c].c_str(),"RaiseFactor") == 0)
//                                calibdata.optmiset.RaiseFactor = tmpnode[colnames[c]].as<double>();
//                            if(strcmp(colnames[c].c_str(),"HeuristicExpansioFactor") == 0)
//                                calibdata.optmiset.HeuristicExpandFactor = tmpnode[colnames[c]].as<double>();
                            if(strcmp(colnames[c].c_str(),"MaxIter") == 0)
                                calibdata.optmiset.MaxIter = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"RRtFromFile") == 0)
//                                calibdata.calccameradata[contCam].RRtFromFile = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"UseSolListFile") == 0)
//                                calibdata.calccameradata[contCam].UseSolListFile = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"Dellistfolder") == 0)
//                                calibdata.calccameradata[contCam].DelListFolder = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"Initfixtxty") == 0)
//                                calibdata.calccameradata[contCam].Initfixtxty = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"repeatsamedata") == 0)
//                                calibdata.calccameradata[contCam].repeatsamedata = tmpnode[colnames[c]].as<int>();
//                            if(strcmp(colnames[c].c_str(),"BranchLevels") == 0)
//                                calibdata.optmiset.BranchLevels = tmpnode[colnames[c]].as<int>();
                        }
                        if(etapa == 0){
                            if(calibdata.lsq)
                                calibdata.IniTsai(contCam,
                                    calibdata.calccameradata[contCam].Initfixtxty,
                                    RRtFileName);
                        }else if(etapa == 1){
                            calibdata.f_t_zest(contCam);                            
                            calibdata.ExpansionP(contCam,
                                                calibdata.calccameradata[contCam].Initfixtxty,
                                                "RRtCalc.yml",
                                                ConfigFileName);
                        }
                        continue;
                    }
                }
            }
//            calibdata.Triangulation(ConfigFileName,"RRtCalc.yml","TranslatioData.yml");
            SaveData(imgnoise[contimgnoise],tarnoise[conttarnoise]);
        }
    }
    cout << "Calibration finish" << endl;
}
void BundleAdjustment(const string& ConfigFileName,const string& DataFileName){
    Icalib calibdata;
    std::vector<CAMERADATA> camdata;
    std::vector<FORMDATA> formdata;
    XMLAcces(filename,camdata,formdata);
    calibdata.LoadCameraMatrix(camdata);
    calibdata.LoadTargetMatrix(formdata);
    for(int contCam = 0;contCam < calibdata.calccameradata.size();contCam++){
        calibdata.ExpansionP(contCam,
                                    calibdata.calccameradata[contCam].Initfixtxty,
                                    DataFileName,
                                    ConfigFileName);
    }
}
void Triangulation(const string ConfigFileName,const string DataFileName,const string TranslFileName){
    Icalib calibdata;
    std::vector<CAMERADATA> camdata;
    std::vector<FORMDATA> formdata;
    XMLAcces(filename,camdata,formdata);
    calibdata.LoadCameraMatrix(camdata);
    calibdata.LoadTargetMatrix(formdata);
    calibdata.Triangulation(ConfigFileName,DataFileName,TranslFileName);
}
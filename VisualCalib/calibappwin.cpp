/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <gtk/gtk.h>

#include "calibapp.hpp"
#include "calibappwin.hpp"
#include "OptSettings.hpp"
#include "NoiseDiamwin.hpp"

using namespace std;

struct _CalibAppWindow
{
  GtkApplicationWindow parent;
  CalibAppWindowPrivate *priv;
};

struct _CalibAppWindowClass
{
  GtkApplicationWindowClass parent_class;
};

typedef struct _CalibAppWindowPrivate CalibAppWindowPrivate;

struct _CalibAppWindowPrivate
{
  GSettings *settings;
  GtkWidget *entrcamfile;
  GtkWidget *entrconfigfilename;
  GtkWidget *entrtranslation;
  GtkWidget *Images;
  GtkWidget *Text;
  GtkWidget *YAML;
  GtkWidget *Blender;
  GtkWidget *entrRRtFileName;
  GtkWidget *chckbttntrue;
  GtkWidget *chckbttnalltar;
  GtkListStore *liststorecolname;
  GtkWidget *chckbttnlinlsq;
  GtkListStore *liststore;
};

struct _GtkCellRendererTextPrivate{
    int columid;
};

G_DEFINE_TYPE_WITH_PRIVATE(CalibAppWindow, calib_app_window, GTK_TYPE_APPLICATION_WINDOW);

static void
calib_app_window_init (CalibAppWindow *win)
{
  win->priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (win);
  gtk_widget_init_template (GTK_WIDGET (win));
  win->priv->settings = g_settings_new ("org.gtk.calibapp");

  g_settings_bind (win->priv->settings, "configfilename",
                   win->priv->entrconfigfilename, "text",
                   G_SETTINGS_BIND_DEFAULT);

  g_settings_bind (win->priv->settings, "camfile",
                   win->priv->entrcamfile, "text",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (win->priv->settings, "translation",
                   win->priv->entrtranslation, "text",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (win->priv->settings, "rrtfilename",
                   win->priv->entrRRtFileName, "text",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (win->priv->settings, "chcktruedata",
                   win->priv->chckbttntrue, "active",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (win->priv->settings, "chckbttnalltar",
                   win->priv->chckbttnalltar, "active",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (win->priv->settings, "chckbttnlinlsq",
                   win->priv->chckbttnlinlsq, "active",
                   G_SETTINGS_BIND_DEFAULT);
}

static void
calib_app_window_dispose (GObject *object)
{
  CalibAppWindow *win;

  win = CALIB_APP_WINDOW (object);
  win->priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (win);
  g_clear_object (&win->priv->settings);

  G_OBJECT_CLASS (calib_app_window_parent_class)->dispose (object);
}



static void
calib_app_window_class_init (CalibAppWindowClass *Class)
{
    G_OBJECT_CLASS (Class)->dispose = calib_app_window_dispose;
  gtk_widget_class_set_template_from_resource (GTK_WIDGET_CLASS (Class),
                                               "/org/gtk/calibapp/window.ui");
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, entrcamfile);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, entrconfigfilename);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, entrtranslation);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, Images);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, Text);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, YAML);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, Blender);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, entrRRtFileName);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, chckbttntrue);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, liststorecolname);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, chckbttnalltar);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), CalibAppWindow, chckbttnlinlsq);
}

CalibAppWindow *
calib_app_window_new (CalibApp *app)
{
  return (CalibAppWindow *)g_object_new (CALIB_APP_WINDOW_TYPE, "application", app, NULL);
}
INT ClickedGAColumn;
extern "C"{
    void
    find_files (GtkWidget *widget,
                 CalibAppWindow   *app)
    {
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (app);
      GtkWidget *dialog;
      GtkFileChooserAction act = GTK_FILE_CHOOSER_ACTION_OPEN;
      gint res;
      dialog = gtk_file_chooser_dialog_new ("Open File",
                                          GTK_WINDOW(app),
                                          act,
                                          "_Cancel",
                                          GTK_RESPONSE_CANCEL,
                                          "_Open",
                                          GTK_RESPONSE_ACCEPT,
                                          NULL);
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      gtk_file_chooser_add_filter(chooser,GTK_FILE_FILTER(priv->Blender));
      gtk_file_chooser_add_filter(chooser,GTK_FILE_FILTER(priv->Images));
      gtk_file_chooser_add_filter(chooser,GTK_FILE_FILTER(priv->YAML));
      gtk_file_chooser_add_filter(chooser,GTK_FILE_FILTER(priv->Text));
      res = gtk_dialog_run (GTK_DIALOG (dialog));
      gchar *filename;
      if (res == GTK_RESPONSE_ACCEPT)
      {
        filename = gtk_file_chooser_get_filename (chooser);
        if(strcmp(gtk_widget_get_name(widget),gtk_widget_get_name(priv->entrcamfile)) == 0)
            gtk_entry_set_text(GTK_ENTRY(priv->entrcamfile),filename);
        if(strcmp(gtk_widget_get_name(widget),gtk_widget_get_name(priv->entrconfigfilename)) == 0)
            gtk_entry_set_text(GTK_ENTRY(priv->entrconfigfilename),filename);
        if(strcmp(gtk_widget_get_name(widget),gtk_widget_get_name(priv->entrtranslation)) == 0)
            gtk_entry_set_text(GTK_ENTRY(priv->entrtranslation),filename);
        if(strcmp(gtk_widget_get_name(widget),gtk_widget_get_name(priv->entrRRtFileName)) == 0)
            gtk_entry_set_text(GTK_ENTRY(priv->entrRRtFileName),filename);
        g_free (filename);
      }
      gtk_widget_destroy (dialog);  
    }
    void
    on_bttncalib_clicked(GtkWidget *widget,
                CalibAppWindow *win){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (win);
        double imgminnoise,imgmaxnoise,imgnoisestep;
        double tarminnoise,tarmaxnoise,tarnoisestep;
        double imgmindiam,imgmaxdiam,imgdiamstep;
        double tarmindiam,tarmaxdiam,tardiamstep;
        std::vector<float> noisediamvals;
        LoadNoiseDiam("noisediam.yml",noisediamvals);
        imgminnoise = noisediamvals[0];
        imgmaxnoise = noisediamvals[1];
        imgnoisestep = noisediamvals[2];
        if(imgnoisestep == 0 && imgmaxnoise != imgminnoise){
            GtkWidget *dialog;
            GtkDialogFlags flags = GTK_DIALOG_MODAL;
            dialog = gtk_message_dialog_new (GTK_WINDOW(win),
                                 flags,
                                 GTK_MESSAGE_OTHER,
                                 GTK_BUTTONS_OK,
                                 "%s",
                                 "O passo do ruído na imagem errado igual a zero.");
            gtk_dialog_run (GTK_DIALOG(dialog));
            gtk_widget_destroy (dialog);
        }
        tarminnoise = noisediamvals[3];
        tarmaxnoise = noisediamvals[4];
        tarnoisestep = noisediamvals[5];
        if(tarnoisestep == 0 && tarmaxnoise != tarminnoise){
            GtkWidget *dialog;
            GtkDialogFlags flags = GTK_DIALOG_MODAL;
            dialog = gtk_message_dialog_new (GTK_WINDOW(win),
                                 flags,
                                 GTK_MESSAGE_OTHER,
                                 GTK_BUTTONS_OK,
                                 "%s",
                                 "O passo do ruído na imagem errado igual a zero.");
            gtk_dialog_run (GTK_DIALOG(dialog));
            gtk_widget_destroy (dialog);
        }

        imgmindiam = noisediamvals[6];
        imgmaxdiam = noisediamvals[7];
        imgdiamstep = noisediamvals[8];
        if(imgdiamstep == 0 && imgmaxdiam != imgmindiam){
            GtkWidget *dialog;
            GtkDialogFlags flags = GTK_DIALOG_MODAL;
            dialog = gtk_message_dialog_new (GTK_WINDOW(win),
                                 flags,
                                 GTK_MESSAGE_OTHER,
                                 GTK_BUTTONS_OK,
                                 "%s",
                                 "O passo do ruído na imagem errado igual a zero.");
            gtk_dialog_run (GTK_DIALOG(dialog));
            gtk_widget_destroy (dialog);
        }
        
        tarmindiam = noisediamvals[9];
        tarmaxdiam = noisediamvals[10];
        tardiamstep = noisediamvals[11];
        if(tardiamstep == 0 && tarmaxdiam != tarmindiam){
            GtkWidget *dialog;
            GtkDialogFlags flags = GTK_DIALOG_MODAL;
            dialog = gtk_message_dialog_new (GTK_WINDOW(win),
                                 flags,
                                 GTK_MESSAGE_OTHER,
                                 GTK_BUTTONS_OK,
                                 "%s",
                                 "O passo do ruído na imagem errado igual a zero.");
            gtk_dialog_run (GTK_DIALOG(dialog));
            gtk_widget_destroy (dialog);
        }

        std::vector<double> imgnoise,tarnoise;
        std::vector<double> imgdiam,tardiam;
        if(imgminnoise > imgmaxnoise || imgmaxnoise < (imgminnoise+imgnoisestep)){
            g_print("%s","Error on noises limits to image data");
            return;
        }
        if(tarminnoise > tarmaxnoise || tarmaxnoise < (tarminnoise+tarnoisestep)){
            GtkWidget *dialog;
            GtkDialogFlags flags = GTK_DIALOG_MODAL;
            dialog = gtk_message_dialog_new (GTK_WINDOW(win),
                                     flags,
                                     GTK_MESSAGE_ERROR,
                                     GTK_BUTTONS_OK,
                                     "%s",
                                     "Error on tarnoise limits");
            gtk_dialog_run (GTK_DIALOG(dialog));
            gtk_widget_destroy (dialog);
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

        CalibCameras(imgnoise,tarnoise,imgdiam,tardiam,gtk_entry_get_text(GTK_ENTRY(priv->entrconfigfilename)));
    }
    void
    on_bttntriangulation_clicked(GtkWidget *widget,
                CalibAppWindow *app){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (app);
        Triangulation(gtk_entry_get_text(GTK_ENTRY(priv->entrconfigfilename)),
                        gtk_entry_get_text(GTK_ENTRY(priv->entrRRtFileName)),
                        gtk_entry_get_text(GTK_ENTRY(priv->entrtranslation)));
    }
    void
    on_bttbundeladj_clicked(GtkWidget *widget,
                CalibAppWindow *app){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (app);
        BundleAdjustment(gtk_entry_get_text(GTK_ENTRY(priv->entrconfigfilename)),
                         gtk_entry_get_text(GTK_ENTRY(priv->entrRRtFileName)));
    }

    void
    on_chckbttntrue_toggled(GtkWidget *widget,
                 GtkEntry   *entrtranslation){
        gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget));
        string strfilename = gtk_entry_get_text(entrtranslation);
        if(active){
            if(!strfilename.empty()){
                if(strfilename.find("True") != string::npos)return;
                int idx = strfilename.find_last_of(".");
                strfilename.insert(idx,"True");
            }else{
                strfilename = "TranslationTrue.yml";
            }
            gtk_entry_set_text(entrtranslation,strfilename.c_str());
        }else{
            if(!strfilename.empty()){
                int idx = strfilename.find_last_of("True");
                if(idx != string::npos){
                    strfilename.erase(idx-3,4);
                }
            }else{
                strfilename = "Translation.yml";
            }
            gtk_entry_set_text(entrtranslation,strfilename.c_str());
        }
    }
    /*Carrega dados dos arquivos para as variáveis*/
    void
    on_bttnloaddata_clicked(GtkWidget *widget,
                CalibAppWindow *app){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (app);
    }
    /*Utiliza um arquivo de configurações para carregar todos os dados*/
    void
    on_bttnloaddefault_clicked(GtkWidget *widget,
                CalibAppWindow *app){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (app);
    }
    /*Salva todas as configurações das variáveis do programa.*/
    void
    on_bttnSaveConfig_clicked(GtkWidget *widget,
                CalibAppWindow *app){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (app);
    }
    void
    on_chckbttnlin_lsq_toggled(GtkWidget *widget,
                CalibAppWindow *app){
        gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget));
        GValue value = G_VALUE_INIT;
        g_value_init(&value,G_TYPE_STRING);
        if(active){
            g_value_set_string(&value,"f-tz");
            g_object_set_property(G_OBJECT(widget),"label",&value);
        }else{
            g_value_set_string(&value,"IniTsai+f-tz");
            g_object_set_property(G_OBJECT(widget),"label",&value);
        }
    }
    void
    on_GAoptmiset_clicked(GtkWidget *widget,
                          CalibAppWindow *appwin){
        CalibAppWindowPrivate *priv;
        priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (CALIB_APP_WINDOW(appwin));
        LoadGtkListStoreIDs(gtk_entry_get_text(GTK_ENTRY(priv->entrcamfile)),priv->liststore,0,1,2);
        LoadGtkListStore(gtk_entry_get_text(GTK_ENTRY(priv->entrconfigfilename)),priv->liststore,
                        priv->liststorecolname,0,1,2);
    }
    void
    on_GAoptmiset_open_window(GtkWidget *widget,
                          GtkWindow *appwin){
        gtk_window_present(appwin);
    }
    void
    on_bttnCancelGA_clicked(GtkWidget *widget,
                          GtkWindow *GAsetwindow){
        gtk_window_close(GAsetwindow);
    }
    void
    on_cellrdrspinSel(GtkCellRendererText *renderer,
               gchar               *path,
               gchar               *new_text,
               GtkTreeViewColumn    *column){
        ClickedGAColumn = gtk_tree_view_column_get_sort_column_id(column);
        cout << ClickedGAColumn << endl;
    }
    void
    on_cellrdrspinpopsize_edited(GtkCellRendererText *renderer,
               gchar               *path,
               gchar               *new_text,
               GtkListStore        *liststore){
        GtkTreeIter iter;
        gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(liststore),&iter,path);
        INT val;
        val = std::stoi(new_text);
        gtk_list_store_set(liststore,&iter,ClickedGAColumn,val,-1);
    }

    void
    on_cellrdrspin_edited(GtkCellRendererText *renderer,
               gchar               *path,
               gchar               *new_text,
               GtkListStore        *liststore){
        GtkTreeIter iter;
        gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(liststore),&iter,path);
        INT val;
        val = std::stoi(new_text);
        gtk_list_store_set(liststore,&iter,ClickedGAColumn,val,-1);
    }
    void
    on_mnitmconfgnoisediam_activate(GtkMenuItem *menuitem,
                                    CalibAppWindow  *window){
        NoiseDiamSettings *noisediamwin;
        noisediamwin = noisediamsettings_new (window);
  
        gtk_window_present (GTK_WINDOW (noisediamwin));
    }
    void
    on_mnitmconfgopt_activate(GtkMenuItem *menuitem,
                              CalibAppWindow  *window){
        OptSettings *optwin;
        optwin = optsettings_new (window);
        gtk_window_present (GTK_WINDOW (optwin));
    }
    void
    on_bttncleanfiles_clicked(GtkWidget *widget,
                          CalibAppWindow *appwin){
        remove("RRtCalc.yml");
        remove("TranslationData.yml");
        remove("projectionmatrix.csv");
        remove("pointlocation.csv");
        remove("pixellocation.csv");
        remove("PixelData.yml");
        remove("camdata.csv");
    }
}
void
calib_app_window_open (CalibAppWindow *win,
                         GFile            *file)
{
  CalibAppWindowPrivate *priv;

  priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (win);
}
const gchar *calib_app_window_get_config_file_name(CalibAppWindow *win){
  return gtk_entry_get_text(GTK_ENTRY(win->priv->entrconfigfilename));
}
const gchar *calib_app_window_get_blender_file_name(CalibAppWindow *win){
  return gtk_entry_get_text(GTK_ENTRY(win->priv->entrcamfile));
}
const gchar *calib_app_window_get_cam_id(CalibAppWindow *win){
  return gtk_entry_get_text(GTK_ENTRY(win->priv->entrcamfile));
}
const gchar *calib_app_window_get_RRtFileName(CalibAppWindow *win){
  return gtk_entry_get_text(GTK_ENTRY(win->priv->entrRRtFileName));
}
std::vector<string> calib_app_window_get_columnnames(CalibAppWindow *win){
    CalibAppWindowPrivate *priv;

    priv = (CalibAppWindowPrivate *)calib_app_window_get_instance_private (win);
    std::vector<string> out;
    GtkTreeIter iter;
    GValue value = G_VALUE_INIT;
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(priv->liststorecolname),&iter);
    do{
        gtk_tree_model_get_value(GTK_TREE_MODEL(priv->liststorecolname),&iter,0,&value);
        if(G_VALUE_HOLDS_STRING(&value))
            out.push_back(g_value_get_string(&value));
        g_value_reset(&value);
        g_value_unset(&value);
    }while(gtk_tree_model_iter_next(GTK_TREE_MODEL(priv->liststorecolname),&iter));
    return out;
}
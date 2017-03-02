/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <gtk/gtk.h>
#include <gtkmm-3.0/gtkmm/liststore.h>
#include <glib-2.0/gobject/gvalue.h>
#include <glib-2.0/gobject/gtype.h>
#include <gtk-3.0/gtk/gtktreemodel.h>
#include <glib-2.0/gobject/gvaluetypes.h>
#include <glib-2.0/glib/glist.h>

#include "calibapp.hpp"
#include "calibappwin.hpp"
#include "OptSettings.hpp"

struct _OptSettings
{
  GtkDialog parent;
};

struct _OptSettingsClass
{
  GtkDialogClass parent_class;
};

typedef struct _OptSettingsPrivate OptSettingsPrivate;

struct _OptSettingsPrivate
{
  GSettings *settings;
  GtkWidget *font;
  GtkWidget *transition;
  GtkTreeView *treeview;
  GtkListStore *liststore;
  GtkWidget *spnbttn_rowconftab;
  GtkListStore *liststorecolname;
};

OptSettingsPrivate *priv;
static   gint ClickedColumn;
const gchar *configfilename;
const gchar *blenderfilename;

G_DEFINE_TYPE_WITH_PRIVATE(OptSettings, optsettings, GTK_TYPE_DIALOG)

static void
optsettings_init (OptSettings *prefs)
{

  priv = (OptSettingsPrivate *)optsettings_get_instance_private (prefs);
  gtk_widget_init_template (GTK_WIDGET (prefs));
  LoadGtkListStoreIDs(blenderfilename,priv->liststore);
  LoadGtkListStore(configfilename,priv->liststore,priv->liststorecolname);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(priv->spnbttn_rowconftab),gtk_tree_model_iter_n_children(GTK_TREE_MODEL(priv->liststore),NULL));
}

static void
optsettings_dispose (GObject *object)
{

  priv = (OptSettingsPrivate *)optsettings_get_instance_private (OPT_SETTINGS (object));

  G_OBJECT_CLASS (optsettings_parent_class)->dispose (object);
}

static void
optsettings_class_init (OptSettingsClass *Class)
{
  G_OBJECT_CLASS (Class)->dispose = optsettings_dispose;

  gtk_widget_class_set_template_from_resource (GTK_WIDGET_CLASS (Class),
                                               "/org/gtk/calibapp/OptSettings.ui");
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),OptSettings,treeview);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),OptSettings,liststore);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),OptSettings,liststorecolname);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),OptSettings,spnbttn_rowconftab);
}

OptSettings *
optsettings_new (CalibAppWindow *win)
{
    ClickedColumn = 0;
    configfilename = calib_app_window_get_config_file_name(win);
    blenderfilename = calib_app_window_get_blender_file_name(win);
  return (OptSettings *)g_object_new (OPT_SETTINGS_TYPE, "transient-for", win, NULL);
}
extern "C"{
    void
    on_renderer_edited(GtkCellRendererText *renderer,
               gchar               *path,
               gchar               *new_text,
               GtkListStore        *liststore){
        GtkTreeIter iter;
        gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(liststore),&iter,path);
        if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),ClickedColumn) == G_TYPE_STRING)
            gtk_list_store_set(liststore,&iter,ClickedColumn,new_text,-1);
        if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),ClickedColumn) == G_TYPE_DOUBLE){
            double value;
            stringstream(new_text) >> value;
            gtk_list_store_set(liststore,&iter,ClickedColumn,value,-1);
        }
        if(gtk_tree_model_get_column_type(GTK_TREE_MODEL(liststore),ClickedColumn) == G_TYPE_INT){
            int value;
            stringstream(new_text) >> value;
            gtk_list_store_set(liststore,&iter,ClickedColumn,value,-1);
        }
        ClickedColumn = 0;
    }
    void
    on_renderer_editedSel(GtkCellRendererText *renderer,
               gchar               *path,
               gchar               *new_text,
               GtkTreeViewColumn    *column){
        GList *tree_view_columns = gtk_tree_view_get_columns(priv->treeview);
        int i = 0;
        for(;i < gtk_tree_view_get_n_columns(priv->treeview);tree_view_columns = g_list_next(tree_view_columns),i++){
            if(strcmp(gtk_tree_view_column_get_title(GTK_TREE_VIEW_COLUMN(tree_view_columns->data)),
                    gtk_tree_view_column_get_title(column)) == 0)
                break;
        }
        ClickedColumn = i;
        g_list_free(tree_view_columns);
    }
    void
    on_renderer_bool_edited(GtkCellRendererToggle *renderer,
               gchar               *path,
               GtkListStore        *liststore){
        GtkTreeIter iter;
        gtk_tree_model_get_iter_from_string(GTK_TREE_MODEL(liststore),&iter,path);
        gboolean active = gtk_cell_renderer_toggle_get_active(renderer);
        gtk_list_store_set(liststore,&iter,ClickedColumn,!active,-1);
        ClickedColumn = 0;
    }
    void
    on_renderer_bool_editedSel(GtkCellRendererToggle *renderer,
               gchar               *path,
               GtkTreeViewColumn    *column){
        GList *tree_view_columns = gtk_tree_view_get_columns(priv->treeview);
        int i = 0;
        for(;i < gtk_tree_view_get_n_columns(priv->treeview);tree_view_columns = g_list_next(tree_view_columns),i++){
            if(strcmp(gtk_tree_view_column_get_title(GTK_TREE_VIEW_COLUMN(tree_view_columns->data)),
                    gtk_tree_view_column_get_title(column)) == 0)
                break;
        }
        ClickedColumn = i;
        g_list_free(tree_view_columns);
    }
    void
    on_bttn_cancel_clicked(GtkButton *button,
               GtkDialog   *windows){
        gtk_window_close(&windows->window);
    }
    void
    on_bttn_ok_clicked(GtkButton *button,
               OptSettings   *win){
        priv = (OptSettingsPrivate *)optsettings_get_instance_private (win);
        SaveGtkListStore(configfilename,priv->liststore,priv->liststorecolname);
    }
    void
    on_spnbttn_rowconftab_output(GtkSpinButton *spin_button,
               GtkListStore       *liststore){
        GtkTreeIter iter;
        gint nrows = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(priv->liststore),NULL);
        if(nrows > 0){
            gtk_tree_model_iter_nth_child(GTK_TREE_MODEL(liststore),&iter,NULL,
                    nrows-1);
            if(gtk_spin_button_get_value(spin_button) > nrows)
                gtk_list_store_append(liststore,&iter);
            else if(gtk_spin_button_get_value(spin_button) < nrows)
                gtk_list_store_remove(liststore,&iter);
        }else{
            gtk_list_store_append(liststore,&iter);           
        }
    }
}
std::vector<string> optsettings_getcolumnnames(OptSettings *opt){
    std::vector<string> out;
    priv = (OptSettingsPrivate *)optsettings_get_instance_private (opt);
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
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "NoiseDiamwin.hpp"

struct _NoiseDiamSettings
{
  GtkWindow parent;
};

struct _NoiseDiamSettingsClass
{
  GtkWindowClass parent_class;
};

typedef struct _NoiseDiamSettingsPrivate NoiseDiamSettingsPrivate;

struct _NoiseDiamSettingsPrivate
{
  GSettings *settings;

  GtkWidget *spnbttnimgminnoise;
  GtkWidget *spnbttnimgmaxnoise;
  GtkWidget *spnbttnimgnoisestep;

  GtkWidget *spnbttntarminnoise;
  GtkWidget *spnbttntarmaxnoise;
  GtkWidget *spnbttntarnoisestep;

  GtkWidget *spnbttnimgmindia;
  GtkWidget *spnbttnimgmaxdiam;
  GtkWidget *spnbttnimgdiamstep;
  
  
  GtkWidget *spnbttntarmindiam;
  GtkWidget *spnbttntarmaxdiam;
  GtkWidget *spnbttntardiamstep;
};

NoiseDiamSettingsPrivate *noisediampriv;

G_DEFINE_TYPE_WITH_PRIVATE(NoiseDiamSettings, noisediamsettings, GTK_TYPE_WINDOW)

static void
noisediamsettings_init (NoiseDiamSettings *noisediamwin)
{
  noisediampriv = (NoiseDiamSettingsPrivate *)noisediamsettings_get_instance_private (noisediamwin);
  gtk_widget_init_template (GTK_WIDGET (noisediamwin));

  noisediampriv->settings = g_settings_new ("org.gtk.calibapp");

  g_settings_bind (noisediampriv->settings, "imgminnoise",
                   noisediampriv->spnbttnimgminnoise, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "imgmaxnoise",
                   noisediampriv->spnbttnimgmaxnoise, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "imgnoisestep",
                   noisediampriv->spnbttnimgnoisestep, "value",
                   G_SETTINGS_BIND_DEFAULT);
  
  
  g_settings_bind (noisediampriv->settings, "tarminnoise",
                   noisediampriv->spnbttntarminnoise, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "tarmaxnoise",
                   noisediampriv->spnbttntarmaxnoise, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "tarnoisestep",
                   noisediampriv->spnbttntarnoisestep, "value",
                   G_SETTINGS_BIND_DEFAULT);

  g_settings_bind (noisediampriv->settings, "imgmindiam",
                   noisediampriv->spnbttnimgmindia, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "imgmaxdiam",
                   noisediampriv->spnbttnimgmaxdiam, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "imgdiamstep",
                   noisediampriv->spnbttnimgdiamstep, "value",
                   G_SETTINGS_BIND_DEFAULT);
  
  
  g_settings_bind (noisediampriv->settings, "tarmindiam",
                   noisediampriv->spnbttntarmindiam, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "tarmaxdiam",
                   noisediampriv->spnbttntarmaxdiam, "value",
                   G_SETTINGS_BIND_DEFAULT);
  g_settings_bind (noisediampriv->settings, "tardiamstep",
                   noisediampriv->spnbttntardiamstep, "value",
                   G_SETTINGS_BIND_DEFAULT);

}

static void
noisediamsettings_dispose (GObject *object)
{

  noisediampriv = (NoiseDiamSettingsPrivate *)noisediamsettings_get_instance_private (NOISEDIAM_SETTINGS (object));

  G_OBJECT_CLASS (noisediamsettings_parent_class)->dispose (object);
}

static void
noisediamsettings_class_init (NoiseDiamSettingsClass *Class)
{
  G_OBJECT_CLASS (Class)->dispose = noisediamsettings_dispose;

  gtk_widget_class_set_template_from_resource (GTK_WIDGET_CLASS (Class),
                                               "/org/gtk/calibapp/NoiseDiamSettings.ui");
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttnimgminnoise);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttnimgmaxnoise);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttnimgnoisestep);

  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttntarminnoise);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttntarmaxnoise);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttntarnoisestep);
    
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttnimgmindia);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttnimgmaxdiam);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttnimgdiamstep);
    
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttntarmindiam);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttntarmaxdiam);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class),NoiseDiamSettings,spnbttntardiamstep);
}

NoiseDiamSettings *
noisediamsettings_new (CalibAppWindow *win)
{
  return (NoiseDiamSettings *)g_object_new (NOISEDIAM_SETTINGS_TYPE, "transient-for", win, NULL);
}

extern "C"{
    void
    on_bttnexitnoisediam_clicked(GtkButton *button,
                                 GtkWindow *win){
        std::vector<float> noisediamvals;
        float val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttnimgminnoise));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttnimgmaxnoise));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttnimgnoisestep));
        noisediamvals.push_back(val);

        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttntarminnoise));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttntarmaxnoise));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttntarnoisestep));
        noisediamvals.push_back(val);

        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttnimgmindia));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttnimgmaxdiam));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttnimgdiamstep));
        noisediamvals.push_back(val);

        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttntarmindiam));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttntarmaxdiam));
        noisediamvals.push_back(val);
        val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(noisediampriv->spnbttntardiamstep));
        noisediamvals.push_back(val);
        
        SaveNoiseDiam("noisediam.yml",noisediamvals);
        gtk_window_close(win);
    }
}
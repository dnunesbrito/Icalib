/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <gtk/gtk.h>

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
};

G_DEFINE_TYPE_WITH_PRIVATE(OptSettings, optsettings, GTK_TYPE_DIALOG)

static void
optsettings_init (OptSettings *prefs)
{
  OptSettingsPrivate *priv;

  priv = (OptSettingsPrivate *)optsettings_get_instance_private (prefs);
  gtk_widget_init_template (GTK_WIDGET (prefs));
}

static void
optsettings_dispose (GObject *object)
{
  OptSettingsPrivate *priv;

  priv = (OptSettingsPrivate *)optsettings_get_instance_private (OPT_SETTINGS (object));

  G_OBJECT_CLASS (optsettings_parent_class)->dispose (object);
}

static void
optsettings_class_init (OptSettingsClass *Class)
{
  G_OBJECT_CLASS (Class)->dispose = optsettings_dispose;

  gtk_widget_class_set_template_from_resource (GTK_WIDGET_CLASS (Class),
                                               "/org/gtk/calibapp/OptSettings.ui");
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), OptSettings, font);
  gtk_widget_class_bind_template_child_private (GTK_WIDGET_CLASS (Class), OptSettings, transition);
}

OptSettings *
optsettings_new (CalibAppWindow *win)
{
  return (OptSettings *)g_object_new (OPT_SETTINGS_TYPE, "transient-for", win, "use-header-bar", TRUE, NULL);
}



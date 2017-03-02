/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   NoiseDiamwin.hpp
 * Author: darlan
 *
 * Created on 22 de Setembro de 2016, 11:16
 */

#ifndef NOISEDIAMWIN_HPP
#define NOISEDIAMWIN_HPP


#include <gtk/gtk.h>
#include <gtkmm-3.0/gtkmm/liststore.h>
#include <glib-2.0/gobject/gvalue.h>
#include <glib-2.0/gobject/gtype.h>
#include <gtk-3.0/gtk/gtktreemodel.h>
#include <glib-2.0/gobject/gvaluetypes.h>
#include <glib-2.0/glib/glist.h>

#include <FileOp.h>
#include "GTKFileOp.hpp"
#include <gtk/gtk.h>
#include "calibappwin.hpp"
#include <sstream>

#define NOISEDIAM_SETTINGS_TYPE (noisediamsettings_get_type ())
#define NOISEDIAM_SETTINGS(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), NOISEDIAM_SETTINGS_TYPE, NoiseDiamSettings))


typedef struct _NoiseDiamSettings          NoiseDiamSettings;
typedef struct _NoiseDiamSettingsClass     NoiseDiamSettingsClass;


GType                   noisediamsettings_get_type     (void);
NoiseDiamSettings        *noisediamsettings_new        (CalibAppWindow *win);


#endif /* NOISEDIAMWIN_HPP */


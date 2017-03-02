/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   OptSettings.hpp
 * Author: darlan
 *
 * Created on 31 de Maio de 2016, 14:08
 */

#ifndef OPTSETTINGS_HPP
#define OPTSETTINGS_HPP

#include <gtk/gtk.h>
#include "calibappwin.hpp"


#define OPT_SETTINGS_TYPE (optsettings_get_type ())
#define OPT_SETTINGS(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), OPT_SETTINGS_TYPE, OptSettings))


typedef struct _OptSettings          OptSettings;
typedef struct _OptSettingsClass     OptSettingsClass;


GType                   optsettings_get_type     (void);
OptSettings        *optsettings_new          (CalibAppWindow *win);


#endif /* OPTSETTINGS_HPP */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   calibappwin.hpp
 * Author: darlan
 *
 * Created on 31 de Maio de 2016, 08:45
 */

#ifndef CALIBAPPWIN_HPP
#define CALIBAPPWIN_HPP

#include <gtk/gtk.h>
#include "calibapp.hpp"
#include <DOMPrint.hpp>



#define CALIB_APP_WINDOW_TYPE (calib_app_window_get_type ())
#define CALIB_APP_WINDOW(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), CALIB_APP_WINDOW_TYPE, CalibAppWindow))


typedef struct _CalibAppWindow         CalibAppWindow;
typedef struct _CalibAppWindowClass    CalibAppWindowClass;
typedef struct _CalibAppWindowPrivate CalibAppWindowPrivate;


GType                   calib_app_window_get_type     (void);
CalibAppWindow       *calib_app_window_new          (CalibApp *app);
void                    calib_app_window_open         (CalibAppWindow *win,
                                                         GFile *file);
const gchar *calib_app_window_get_config_file_name(CalibAppWindow *win);
const gchar *calib_app_window_get_blender_file_name(CalibAppWindow *win);
const gchar *calib_app_window_get_cam_id(CalibAppWindow *win);
std::vector<string> calib_app_window_get_columnnames(CalibAppWindow *win);

#endif /* CALIBAPPWIN_HPP */


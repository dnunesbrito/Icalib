/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   calibapp.hpp
 * Author: darlan
 *
 * Created on 31 de Maio de 2016, 08:36
 */

#ifndef CALIBAPP_HPP
#define CALIBAPP_HPP

#include <gtk/gtk.h>
#include <Icalib.h>
#include <IcalibThread.hpp>
#include "GTKFileOp.hpp"

#define CALIB_APP_TYPE (calib_app_get_type ())
#define CALIB_APP(obj) (G_TYPE_CHECK_INSTANCE_CAST ((obj), CALIB_APP_TYPE, CalibApp))


typedef struct _CalibApp       CalibApp;
typedef struct _CalibAppClass  CalibAppClass;


GType           calib_app_get_type    (void);
CalibApp     *calib_app_new         (void);
void        CalibCamera(string& CamId,string& TargID);
void        CalibCameras(std::vector<double>& imgnoise,
                                std::vector<double>& tarnoise,
                                std::vector<double>& imgdiam,
                                std::vector<double>& tardiam,
                                const string& ConfigFileName,bool usealltar = true);
void        BundleAdjustment(const string& ConfigFileName,const string& DataFileName);
void        Triangulation(const string ConfigFileName,const string DataFileName,const string TranslFileName);

#endif /* CALIBAPP_HPP */


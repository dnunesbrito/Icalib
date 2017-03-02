/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GTKFileOp.hpp
 * Author: darlan
 *
 * Created on 24 de Junho de 2016, 09:29
 */

#ifndef GTKFILEOP_HPP
#define GTKFILEOP_HPP

#include <gtk-3.0/gtk/gtk.h>
#include <FileOp.h>
#include <DOMPrint.hpp>
#include <sys/stat.h>
#include <sys/types.h>

void LoadGtkListStore(const std::string& filename,
                        GtkListStore *liststore,
                        GtkListStore *colunmnames,
                            int CameIDColumn = 0,
                            int TargetIDColumn = 1,
                            int EtapaColumn = 14);
void LoadGtkListStoreIDs(const std::string& filename,
                            GtkListStore *liststore,
                            int CameIDColumn = 0,
                            int TargetIDColumn = 1,
                            int EtapaColumn = 14);
void SaveGtkListStore(const std::string& filename,GtkListStore *liststore,GtkListStore *columnnames);
void SaveNoiseDiam(const std::string& filename,std::vector<float>& noisediamvals);
void LoadNoiseDiam(const std::string& filename,std::vector<float>& noisediamvals);
void SaveData(REAL imgnoise,REAL tarnoise);

#endif /* GTKFILEOP_HPP */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: darlan
 *
 * Created on 31 de Maio de 2016, 08:30
 */

#include <cstdlib>
#include "calibapp.hpp"
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
  /* Since this example is running uninstalled,
   * we have to help it find its schema. This
   * is *not* necessary in properly installed
   * application.
   */
  g_setenv ("GSETTINGS_SCHEMA_DIR", ".", FALSE);  

  return g_application_run (G_APPLICATION (calib_app_new ()), argc, argv);
}


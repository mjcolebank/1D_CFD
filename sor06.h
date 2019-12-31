/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 1.0                                                           */
/*  Date: 30 Dec. 2019                                                     */
/*                                                                         */
/*  Primary Authors: M.S. Olufsen                                          */
/*  Key Contributers: M.U. Qureshi & M.J. Colebank                         */
/*  Institute: North Carolina State University, Raleigh, NC, USA           */
/*  Funded by: NSF-DMS # 1615820 & AHA #19PRE34380459                      */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/


#ifndef _SOR06_H
#define _SOR06_H

#include <cmath>
#include "tools.h"

int    nbrves;                            // Number of vessels in the tree.

int    tmstps = 8192,                     // The number of timesteps per period.
       plts   = 1024;                     // Number of plots per period.

const char *CO_filename = "Qin.dat";      // Input flow file at the heart for the control mouse.


         
double conv   = 1332.220,               // Conversion from mmHg to SI-units.
       rho    = 1.055,                  // Density of blood [g/cm^3].
       mu     = 0.049,                  // Viscosity of blood [g/cm/s].
       nu     = mu/rho,                 // Dynamic viscosity of blood [cm^2/s].
       Lr     = 1.0,                    // Characteristic radius of the
                                        // vessels in the tree [cm].
       Lr2    = sq(Lr),                 // The squared radius [cm2].
       Lr3    = cu(Lr),                 // The radius to the third power [cm^3].
       g      = 981.0,                  // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,               // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),      // The squared Froudes number.
       Re     = q*rho/mu/Lr,            // Reynolds number.
       p0     = 2.0/rho/g/Lr*conv,      // Ensures a certain diastolic pressure.
       *fjac[18],                       // Work space used by bound_bif.
       xr, f, df;                       // Work space used by bound_right.

int max_cycles = 40,
    cycles     = 1;
double Tper, Fcst, Period, k, Deltat;

#endif

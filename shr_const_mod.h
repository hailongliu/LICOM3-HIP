#ifndef INCLUDED_SHR_CONST
#define INCLUDED_SHR_CONST //to avoid duplicated include

//
// SVN $Id: shr_const_mod.F90 6749 2007-10-04 20:58:20Z jwolfe $
// SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/branch_tags/cesm1_0_5_rel_tags/cesm1_0_5_n01_share3_110527/shr/shr_const_mod.F90 $
//

// MODULE shr_const_mod
//    use shr_kind_mod
//    integer(SHR_KIND_IN,parameter,private :: R8  SHR_KIND// rename for local readability only

   //----------------------------------------------------------------------------
   // physical constants (all data public
   //----------------------------------------------------------------------------

#define SHR_CONST_PI       (double)3.14159265358979323846 // pi
#define SHR_CONST_CDAY     (double)86400.0     // sec in calendar day ~ sec
#define SHR_CONST_SDAY     (double)86164.0     // sec in siderial day ~ sec
#define SHR_CONST_OMEGA    ((double)2.0*SHR_CONST_PI/SHR_CONST_SDAY) // earth rot ~ rad/sec
#define SHR_CONST_REARTH   (double)6.37122e6   // radius of earth ~ m
#define SHR_CONST_G        (double)9.80616     // acceleration of gravity ~ m/s^2

#define SHR_CONST_STEBOL   (double)5.67e-8     // Stefan-Boltzmann constant ~ W/m^2/K^4
#define SHR_CONST_BOLTZ    (double)1.38065e-23 // Boltzmann's constant ~ J/K/molecule
#define SHR_CONST_AVOGAD   (double)6.02214e26  // Avogadro's number ~ molecules/kmole
#define SHR_CONST_RGAS     (SHR_CONST_AVOGAD*SHR_CONST_BOLTZ)       // Universal gas constant ~ J/K/kmole
#define SHR_CONST_MWDAIR   (double)28.966      // molecular weight dry air ~ kg/kmole
#define SHR_CONST_MWWV     (double)18.016      // molecular weight water vapor
#define SHR_CONST_RDAIR    (SHR_CONST_RGAS/SHR_CONST_MWDAIR)        // Dry air gas constant     ~ J/K/kg
#define SHR_CONST_RWV      (SHR_CONST_RGAS/SHR_CONST_MWWV)          // Water vapor gas constant ~ J/K/kg
#define SHR_CONST_ZVIR     ((SHR_CONST_RWV/SHR_CONST_RDAIR)-(double)1.0)// RWV/RDAIR - 1.0
#define SHR_CONST_KARMAN   (double)0.4         // Von Karman constant
#define SHR_CONST_PSTD     (double)101325.0    // standard pressure ~ pascals
#define SHR_CONST_PDB      (double)0.0112372   // ratio of 13C/12C in Pee Dee Belemnite (C isotope standard
 
#define SHR_CONST_TKTRIP   (double)273.16      // triple point of fresh water        ~ K
#define SHR_CONST_TKFRZ    (double)273.15      // freezing T of fresh water          ~ K 
#define SHR_CONST_TKFRZSW  (SHR_CONST_TKFRZ - (double)1.8)// freezing T of salt water  ~ K

#define SHR_CONST_RHODAIR  (SHR_CONST_PSTD/(SHR_CONST_RDAIR*SHR_CONST_TKFRZ))             // density of dry air at STP  ~ kg/m^3
#define SHR_CONST_RHOFW    (double)1.000e3     // density of fresh water     ~ kg/m^3
#define SHR_CONST_RHOSW    (double)1.026e3     // density of sea water       ~ kg/m^3
#define SHR_CONST_RHOICE   (double)0.917e3     // density of ice             ~ kg/m^3
#define SHR_CONST_CPDAIR   (double)1.00464e3   // specific heat of dry air   ~ J/kg/K
#define SHR_CONST_CPWV     (double)1.810e3     // specific heat of water vap ~ J/kg/K
#define SHR_CONST_CPVIR    ((SHR_CONST_CPWV/SHR_CONST_CPDAIR)-(double)1.0)// CPWV/CPDAIR - 1.0
#define SHR_CONST_CPFW     (double)4.188e3     // specific heat of fresh h2o ~ J/kg/K
#define SHR_CONST_CPSW     (double)3.996e3     // specific heat of sea h2o   ~ J/kg/K
#define SHR_CONST_CPICE    (double)2.11727e3   // specific heat of fresh ice ~ J/kg/K
#define SHR_CONST_LATICE   (double)3.337e5     // latent heat of fusion      ~ J/kg
#define SHR_CONST_LATVAP   (double)2.501e6     // latent heat of evaporation ~ J/kg
#define SHR_CONST_LATSUB   (SHR_CONST_LATICE + SHR_CONST_LATVAP)              // latent heat of sublimation ~ J/kg
                         
#define SHR_CONST_OCN_REF_SAL  (double)34.7    // ocn ref salinity (psu
#define SHR_CONST_ICE_REF_SAL   (double)4.0    // ice ref salinity (psu

#define SHR_CONST_SPVAL    (double)1.0e30      // special missing value

#endif //INCLUDED_SHR_CONST

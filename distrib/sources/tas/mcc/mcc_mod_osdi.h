/******************************************************************************/
/*                                                                            */
/*                      Chaine de CAO & VLSI   AVERTEC                        */
/*                                                                            */
/*    Fichier : mcc_mod_osdi.h                                                 */
/*                                                                            */
/*                                                                            */
/*    (c) copyright 2001 AVERTEC                                              */
/*    Tous droits reserves                                                    */
/*                                                                            */
/*    Auteur(s) : Gregoire AVOT                                               */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#define MCC_OSDI_TUNE_NO_JUNCTION  0x01
#define MCC_OSDI_TUNE_NO_EXTRINSIC 0x02

/* USER field of mcc_modellist* */

typedef struct {
  float qg, qb, qd, qs ;
  float qgsov, qgdov, qgb ;
  float qjbd, qjbs ;
} osdicharge ;

typedef struct scacheosdicharge {
  struct scacheosdicharge *next ;
  osdicharge               charge ;
  char                   *key ;
} cacheosdicharge ;

typedef struct {
  float charge[4] ;
} osdibrcharge ;

typedef struct {
  double ab, ls, lg ;
} osdijuncapconfig ;

void mcc_initparam_osdi( mcc_modellist *ptmodel );
void mcc_clean_osdi( mcc_modellist *ptmodel );

void mcc_calcQint_osdi( mcc_modellist *ptmodel, 
                       double L, 
                       double W,
                       double temp, 
                       double vgs, 
                       double vbs, 
                       double vds,
                       double *ptQg, 
                       double *ptQs, 
                       double *ptQd,
                       double *ptQb,
                       elp_lotrs_param *lotrsparam
                     );
double mcc_calcCGP_osdi( mcc_modellist *ptmodel,
                        elp_lotrs_param *lotrsparam, 
                        double vgx, 
                        double L, 
                        double W,
                        double temp,
                        double *ptQov 
                      );
double mcc_calcCGD_osdi( mcc_modellist *ptmodel, 
                        double L, 
                        double W, 
                        double temp, 
                        double vgs0, 
                        double vgs1, 
                        double vbs, 
                        double vds,
                        elp_lotrs_param *lotrsparam
                      );
double mcc_calcCGSI_osdi( mcc_modellist *ptmodel, 
                         double L, 
                         double W, 
                         double temp, 
                         double vgsf, 
                         double vgsi, 
                         double vbs, 
                         double vdsf,
                         double vdsi,
                         elp_lotrs_param *lotrsparam
                       );
double mcc_calcCGG_osdi( mcc_modellist *ptmodel, 
                         double L, 
                         double W, 
                         double temp, 
                         double vgsf, 
                         double vgsi, 
                         double vbs, 
                         double vdsf,
                         double vdsi,
                         elp_lotrs_param *lotrsparam
                       );
double mcc_calcVTH_osdi( mcc_modellist *ptmodel, 
                        double L, 
                        double W, 
                        double temp, 
                        double vbstrue, 
                        double vds, 
                        elp_lotrs_param *lotrsparam 
                      );
double mcc_calcIDS_osdi( mcc_modellist *ptmodel, 
                        double Vbstrue,
                        double Vgs,
                        double Vds, 
                        double W,
                        double L, 
                        double Temp,
                        elp_lotrs_param *lotrsparam
                      );
double mcc_calcDWCJ_osdi( mcc_modellist *mccmodel, 
                         elp_lotrs_param *lotrsparam, 
                         double temp,
                         double L, 
                         double W
                       );
double mcc_calcCDS_osdi( mcc_modellist  *ptmodel, 
                        elp_lotrs_param *lotrsparam,
                        double          temp, 
                        double          vbx1, 
                        double          vbx2,
                        double          L,
                        double          W
                      );
double mcc_calcCDP_osdi( mcc_modellist *ptmodel, 
                        elp_lotrs_param *lotrsparam,
                        double temp, 
                        double vbx1, 
                        double vbx2,
                        double L,
                        double W
                      );
double mcc_calcCDW_osdi( mcc_modellist *ptmodel, 
                        elp_lotrs_param *lotrsparam,
                        double temp, 
                        double vbx1, 
                        double vbx2, 
                        double L, 
                        double W 
                      );

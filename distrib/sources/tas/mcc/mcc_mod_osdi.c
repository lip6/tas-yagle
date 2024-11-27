/******************************************************************************/
/*                                                                            */
/*                      Chaine de CAO & VLSI   AVERTEC                        */
/*                                                                            */
/*    Fichier : mcc_mod_osdi.c                                                 */
/*                                                                            */
/*                                                                            */
/*    (c) copyright 2003 AVERTEC                                              */
/*    Tous droits reserves                                                    */
/*                                                                            */
/*    Auteur(s) : Gregoire AVOT                                               */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#include <alloca.h>
#include "osdi.h"

#include MCC_H
#include "mcc_util.h"
#include "mcc_mod_util.h"
#include "mcc_mod_spice.h"
#include "mcc_mod_osdi.h"
#include "mcc_mod_osdi_interface.h"

#define EPSILON 0.001

void mcc_initparam_osdi( mcc_modellist *ptmodel )
{
  static int  f = 0 ;
  ptype_list *ptl ;

  if( !f )
    if( !osdi_loaddynamiclibrary() )
      exit(1);

  ptl = getptype( ptmodel->USER, OSDICACHECHARGE );
  if( ptl ) {
    printf( "warning : previous cache charge not empty !!!\n " );
    mcc_clean_osdi( ptmodel );
  }
  
  mcc_clean_osdi_interface( ptmodel, 1 );

  f = 1 ;
  ptmodel = NULL ;
}

void mcc_clean_osdi( mcc_modellist *mccmodel )
{
  ptype_list     *ptl ;
  cacheosdicharge *cache, *next ;

  ptl = getptype( mccmodel->USER, OSDICACHECHARGE );
  if ( ptl ) {
    for( cache = (cacheosdicharge*)ptl->DATA ; cache ; cache = next ) {
      next = cache->next ;
      free( cache );
    }

    mccmodel->USER = delptype( mccmodel->USER, OSDICACHECHARGE );
  }

  mcc_clean_osdi_interface( mccmodel, 0 );
}

void osdi_mcc_addtuneeffect( mcc_modellist *mccmodel, 
                        osditunedparam *tparam, 
                        int flag, 
                        osdijuncapconfig *juncapconfig
                      )
{
  if( ( flag & MCC_OSDI_TUNE_NO_JUNCTION ) == MCC_OSDI_TUNE_NO_JUNCTION ) {
    /* desactivate junction capacitance calculation 
       see doc psp 102.0 issued 06/2006 : p9 #6            */
    tparam->param[ tparam->n   ] = namealloc("SWJUNCAP") ; 
    tparam->value[ tparam->n++ ] = 0.0 ;
  }
  else {
      int  swjuncap=3;
      char *name = namealloc("swjuncap");
      for (mcc_paramlist *p=mccmodel->PARAM; p; p=p->NEXT) {
	if(p->NAME == name) swjuncap = p->VALUE;
      }
      tparam->param[ tparam->n   ] = namealloc("SWJUNCAP") ; 
      tparam->value[ tparam->n++ ] = swjuncap ;
      if( juncapconfig ) {
        if(swjuncap == 1) {
          tparam->param[ tparam->n   ] = namealloc("ABDRAIN") ; 
          tparam->value[ tparam->n++ ] = juncapconfig->ab ;
          tparam->param[ tparam->n   ] = namealloc("LSDRAIN") ; 
          tparam->value[ tparam->n++ ] = juncapconfig->ls ;
          tparam->param[ tparam->n   ] = namealloc("LGDRAIN") ; 
          tparam->value[ tparam->n++ ] = juncapconfig->lg ;

          tparam->param[ tparam->n   ] = namealloc("ABSOURCE") ; 
          tparam->value[ tparam->n++ ] = juncapconfig->ab ;
          tparam->param[ tparam->n   ] = namealloc("LSSOURCE") ; 
          tparam->value[ tparam->n++ ] = juncapconfig->ls ;
          tparam->param[ tparam->n   ] = namealloc("LGSOURCE") ; 
          tparam->value[ tparam->n++ ] = juncapconfig->lg ;
        }
    }
  }
  
  if( ( flag & MCC_OSDI_TUNE_NO_EXTRINSIC ) == MCC_OSDI_TUNE_NO_EXTRINSIC ) {
      int  swgeo;
      char *name = namealloc("swgeo");
      for (mcc_paramlist *p=mccmodel->PARAM; p; p=p->NEXT) {
	if(p->NAME == name) swgeo = p->VALUE;
      }
      switch (swgeo) {
      case 2: 
     /* level==1038 Binning models*/
      tparam->param[ tparam->n   ] = namealloc("POCGOV") ;    /* p64 #3.253 */
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLCGOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PWCGOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLWCGOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
    
      tparam->param[ tparam->n   ] = namealloc("POCGOVD") ;    /* p64 #3.254 */
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLCGOVD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PWCGOVD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLWCGOVD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
    
      tparam->param[ tparam->n   ] = namealloc("POCFR") ;    /* p65 #3.265 */
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLCFR") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PWCFR") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLWCFR") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
     
      tparam->param[ tparam->n   ] = namealloc("POCFRD") ;    /* p65 #3.266 */
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLCFRD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PWCFRD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLWCFRD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
     
      tparam->param[ tparam->n   ] = namealloc("POCGBOV") ;    /* p56 #3.3 */
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLCGBOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PWCGBOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("PLWCGBOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      break;
    case 1: // Geometrical scaling
      tparam->param[ tparam->n   ] = namealloc("LOV") ;    /* p52 #3.108 */
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("LOVD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CGBOVL") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CFRW") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CFRDW") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      break;
    default: // Local
      tparam->param[ tparam->n   ] = namealloc("CGOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CGOVD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CGBOV") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CFR") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      tparam->param[ tparam->n   ] = namealloc("CFRD") ;
      tparam->value[ tparam->n++ ] = 0.0 ;
      break;
   }
  }
}


void osdi_copycharge( osdicharge *chargeorg, osdicharge *chargedest )
{
    chargedest->qg    = chargeorg->qg ;
    chargedest->qb    = chargeorg->qb ;
    chargedest->qd    = chargeorg->qd ;
    chargedest->qs    = chargeorg->qs ;
    chargedest->qgsov = chargeorg->qgsov ;
    chargedest->qgdov = chargeorg->qgdov ;
    chargedest->qgb   = chargeorg->qgb ;
    chargedest->qjbd  = chargeorg->qjbd ;
    chargedest->qjbs  = chargeorg->qjbs ;
}

char* osdi_cachechargegenkey( double           vgs,
                         double           vds,
                         double           vbs,
                         osdijuncapconfig *juncapconfig 
                       )
{
  char buf[1024];
  char *pt ;

  if( juncapconfig )
    sprintf( buf, "%.2f %.2f %.2f %.2e %.2e %.2e", vgs, 
                                                   vds,
                                                   vbs, 
                                                   juncapconfig->ab,
                                                   juncapconfig->ls,
                                                   juncapconfig->lg
           );
  else
    sprintf( buf, "%.2f %.2f %.2f", vgs,
                                    vds,
                                    vbs
           );

  pt = namealloc( buf );

  return pt ;
}

osdicharge* osdi_findchargeincache( mcc_modellist   *mccmodel,
                              double           vgs,
                              double           vds,
                              double           vbs,
                              osdijuncapconfig *juncapconfig
                            )
{
  char            *key ;
  ptype_list      *ptl ;
  cacheosdicharge *cache ;
  
  ptl = getptype( mccmodel->USER, OSDICACHECHARGE );
  if( !ptl )
    return NULL ;
    
  key = osdi_cachechargegenkey( vgs, vds, vbs, juncapconfig );
  for( cache = (cacheosdicharge*)ptl->DATA ; cache ; cache = cache->next ) {
    if( cache->key == key )
      break ;
  }

  if( cache )
    return &(cache->charge) ;

  return NULL ;
}

void osdi_addchargeincache( mcc_modellist   *mccmodel,
                       osdicharge       *charge,
                       double           vgs,
                       double           vds,
                       double           vbs,
                       osdijuncapconfig *juncapconfig
                     )
{
  ptype_list      *ptl ;
  cacheosdicharge *cache ;
  
  ptl = getptype( mccmodel->USER, OSDICACHECHARGE );
  if( !ptl ) {
    mccmodel->USER = addptype( mccmodel->USER, OSDICACHECHARGE, NULL );
    ptl = mccmodel->USER ;
  }

  cache = mbkalloc( sizeof( cacheosdicharge ) );
  cache->next = ptl->DATA ;
  ptl->DATA   = cache ;
  cache->key  = osdi_cachechargegenkey( vgs, vds, vbs, juncapconfig );
  
  osdi_copycharge( charge, &(cache->charge) );
    
}


void osdi_mcc_getcharge( mcc_modellist   *mccmodel,
                    double           L,
                    double           W,
                    double           temp,
                    double           vgs,
                    double           vds,
                    double           vbs,
                    elp_lotrs_param *lotrsparam,
                    osdijuncapconfig *juncapconfig,
                    osdicharge       *charge 
                  )
{
  osdi_trs            model ;
  osditunedparam  tparam ;
  int            base ;
  double         *i_charge;
  osdibrcharge   qq[3] ;
  int            calc[3] ;
  int            n ;
  int            p ;
  int            nbbr ;
  double         qint[4] ;
  osdibrcharge   qout[3] ;
  osdicharge     *cache ;

  if( V_BOOL_TAB[__AVT_USE_CACHE_PSP].VALUE )
    cache = osdi_findchargeincache( mccmodel, vgs, vds, vbs, juncapconfig );
  else
    cache = NULL ;

  if( cache )
    osdi_copycharge( cache, charge );
  else {
    /* we suppose the simkit interface use the following node order :
       0 : D
       1 : G
       2 : S
       3 : B
    */

    /* compute the 3 charge components :
       1 : intrinsic
       2 : extrinsic (gs and gd overlap)
       3 : junction
       OSDI will not provide each branch charges. Then we use the operational parameters from the instance of the model. Here we focus on getting intrinsic charges.
       Other values will be used in capacitor forms, and we will get the capacitor vlaues from operational paraeters, later.
    */
    
    tparam.param = (char**)alloca( sizeof( char* ) * 16 );
    tparam.value = (double*)alloca( sizeof( double ) * 16 );

    
    calc[0] = MCC_OSDI_TUNE_NO_JUNCTION | MCC_OSDI_TUNE_NO_EXTRINSIC ;
    calc[1] = MCC_OSDI_TUNE_NO_JUNCTION ;
    calc[2] = MCC_OSDI_TUNE_NO_EXTRINSIC ;

/* osdi psp, we get only intrinsic charge */
    for( n=0 ; n<3 ; n++) {

      tparam.n = 0 ;
      
      memset(&model, 0, sizeof(osdi_trs));
      osdi_mcc_addtuneeffect( mccmodel, &tparam, calc[n], juncapconfig );
      osdi_initialize( &model, mccmodel, lotrsparam, L, W, temp, &tparam );
      if( model.model->num_terminals != 4 ) {
        printf( "number of external nodes differs than 4. can't extract charge values\n" );
        return ;
      }

      osdi_set_polarization( &model, vgs, vds, vbs );
      i_charge = (double*)calloc(model.model->num_nodes, sizeof(double));
      model.model->load_residual_react(model.idata, model.mdata, i_charge);
      for(int i=0;i<4;i++) qq[n].charge[i] = i_charge[i] + i_charge[i+4]; 
      // collapsible node charge added.
      osdi_terminate( &model );
    }


    /* fill intrinsic component */
    charge->qd = qq[0].charge[0];
    charge->qg = qq[0].charge[1];
    charge->qs = qq[0].charge[2];
    charge->qb = qq[0].charge[3];
    /* fill overlap component */
    charge->qgdov = qq[1].charge[0] - qq[0].charge[0];
    charge->qgb   = qq[1].charge[3] - qq[0].charge[3];
    charge->qgsov = qq[1].charge[2] - qq[0].charge[2];
    /* fill junction component */
    charge->qjbd  = qq[2].charge[0] - qq[0].charge[0];
    charge->qjbs  = qq[2].charge[2] - qq[0].charge[2];

    free(i_charge);

    if( V_BOOL_TAB[__AVT_USE_CACHE_PSP].VALUE )
      osdi_addchargeincache( mccmodel, charge, vgs, vds, vbs, juncapconfig );
  }
}

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
                     )
{
  osdicharge charge ;
  double    s ;

  if( ptQg ) *ptQg = 0.0 ;
  if( ptQs ) *ptQs = 0.0 ;
  if( ptQd ) *ptQd = 0.0 ;
  if( ptQb ) *ptQb = 0.0 ;

  osdi_mcc_getcharge( ptmodel, L, W, temp, vgs, vds, vbs, lotrsparam, NULL, &charge );

  s = W*L ;

  if( ptQd ) *ptQd = charge.qd / s ;;
  if( ptQg ) *ptQg = charge.qg / s ;;
  if( ptQs ) *ptQs = charge.qs / s ;;
  if( ptQb ) *ptQb = charge.qb / s ;;
}

double mcc_calcCGP_osdi( mcc_modellist   *ptmodel,
                        elp_lotrs_param *lotrsparam, 
                        double           vgx, 
                        double           L, 
                        double           W, 
                        double           temp,
                        double          *ptQov 
                      )
{
  osdicharge charge ;
  double    cgp ;
#if 0
  osdi_trs      model ;
  uint32_t id_cgp,accflag ;
  char  *name;

  name = namealloc( "cgdol" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cgp = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cgp, &accflag, 0);
  osdi_set_polarization( &model, vgx, 0.0, 0.0 );
  cgp = *ptr ;

  osdi_terminate( &model );
  if( vgx > EPSILON || vgx < -EPSILON )
    cgp = fabs(cgp) ;
  else
    cgp = 0.0 ;
  
  if( ptQov )
    *ptQov = fabs(cgp*vgx/W);

#else  
  osdi_mcc_getcharge( ptmodel, L, W, temp, vgx, 0.0, 0.0, lotrsparam, NULL, &charge);

  if( vgx > EPSILON || vgx < -EPSILON )
    cgp = fabs(charge.qgdov/vgx) ;
  else
    cgp = 0.0 ;
  
  if( ptQov )
    *ptQov = fabs(charge.qgdov/W);

#endif  
  return cgp / W;
}

double mcc_calcCGD_osdi( mcc_modellist *ptmodel, 
                        double L, 
                        double W, 
                        double temp, 
                        double vgs0, 
                        double vgs1, 
                        double vbs, 
                        double vds0,
                        double vds1,
                        elp_lotrs_param *lotrsparam
                      )
{
  osdicharge charge0 ;
  osdicharge charge1 ;
  double    cgd, cgd0, cgd1 ;
  double    s ;
  s = W*L ;

#if 0
  osdi_trs      model ;
  uint32_t id_cgd,accflag ;
  char  *name;

  name = namealloc( "cdg" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cgd = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cgd, &accflag, 0);

  osdi_set_polarization( &model, vgs0, vds0, vbs );
  cgd0 = *ptr ;
  osdi_set_polarization( &model, vgs1, vds1, vbs );
  cgd1 = *ptr ;

  osdi_terminate( &model );
  cgd = fabs( (cgd1*vds1 - cgd0*vds0)/(vgs1 - vgs0) )/s ;
#else


  osdi_mcc_getcharge( ptmodel, L, W, temp, vgs0, vds0, vbs, lotrsparam, NULL, &charge0 );
  osdi_mcc_getcharge( ptmodel, L, W, temp, vgs1, vds1, vbs, lotrsparam, NULL, &charge1 );

  s = W*L ;
  cgd = fabs( ( charge1.qd - charge0.qd ) / ( vgs1 - vgs0 ) )/s ;
#endif
  
  return cgd ;
}

// CGG return the charge to calculate input capa
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
                       )
{
  osdicharge charge0 ;
  osdicharge charge1 ;
  double cgg, cgg0, cgg1 ;
  double s ;
  osdi_trs      model ;
  uint32_t id_cgg,accflag ;
  char  *name;

  name = namealloc( "cgg" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cgg = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cgg, &accflag, 0);
  osdi_set_polarization( &model, vgsi, vdsi, vbs );
  cgg0 = *ptr ;
  osdi_set_polarization( &model, vgsf, vdsf, vbs );
  cgg1 = *ptr ;
  cgg = fabs( (cgg0*vgsi - cgg1*vgsf) ) ;

  osdi_terminate( &model );

  return cgg ;
}

double mcc_calcCGSI_osdi( mcc_modellist *ptmodel, 
                         double L, 
                         double W, 
                         double temp, 
                         double vgs,
                         double vbs, 
                         double vds,
                         elp_lotrs_param *lotrsparam
                       )
{
  osdicharge charge0 ;
  osdicharge charge1 ;
  double cgs ,cgs0, cgs1;
  double s ;
#if 0
  osdi_trs      model ;
  uint32_t id_cgs,accflag ;
  char  *name;

  name = namealloc( "csg" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cgs = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cgs, &accflag, 0);
  osdi_set_polarization( &model, 0.0, vds, vbs );
  cgs0 = *ptr ;
  osdi_set_polarization( &model, vgs, vds, vbs );
  cgs1 = *ptr ;
  s = W*L ;
  cgs = fabs( cgs1 )/s ;

  osdi_terminate( &model );
#endif

  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vds, vbs, lotrsparam, NULL, &charge0 );
  osdi_mcc_getcharge( ptmodel, L, W, temp, vgs, vds, vbs, lotrsparam, NULL, &charge1 );

  s = W*L ;
  cgs = fabs( ( charge1.qs - charge0.qs ) / vgs )/s ;

  return cgs ;
}

double mcc_calcVTH_osdi( mcc_modellist   *mccmodel,
                        double           L,
                        double           W,
                        double           temp,
                        double           vbs,
                        double           vdd,
                        elp_lotrs_param *lotrsparam 
                      )
{
  osdi_trs      model ;
  double   vth ;
  uint32_t id_vth,accflag ;
  char  *name;

  name = namealloc( osdi_op_param_label[ OSDI_OP_PARAM_VTH ]);
 
  osdi_initialize( &model, mccmodel, lotrsparam, L, W, temp, NULL );
  id_vth = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_vth, &accflag, 0);
  osdi_set_polarization( &model, vdd, vdd, vbs );
  vth = *ptr ;

  osdi_terminate( &model );

  if( mccmodel->TYPE == MCC_TRANS_P )
    vth = -vth ;
    
  return vth ;
}

double mcc_calcIDS_osdi( mcc_modellist *mccmodel, 
                        double vbs,
                        double vgs,
                        double vds, 
                        double W,
                        double L, 
                        double temp,
                        elp_lotrs_param *lotrsparam
                      )
{
  osdi_trs      model ;
  double   ids ;
  uint32_t id_ids, accflag ;
  char  *name;

  name = namealloc( osdi_op_param_label[ OSDI_OP_PARAM_IDS ]);
 
  osdi_initialize( &model, mccmodel, lotrsparam, L, W, temp, NULL );
  id_ids = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_ids,&accflag, 0);
  osdi_set_polarization( &model, vgs, vds, vbs );
  
  ids = *ptr ;

  osdi_terminate( &model );
  return fabs(ids) ;
}

double mcc_calcDWCJ_osdi( mcc_modellist *mccmodel, 
                         elp_lotrs_param *lotrsparam, 
                         double temp,
                         double L, 
                         double W
                       )
{
  osdi_trs           model ;
  double        dwcj ;
  uint32_t           id_weff, accflag ;
  double        weff ;
  char  *name;

  name = namealloc( osdi_op_param_label[ OSDI_OP_PARAM_WEFF ]);
  osdi_initialize( &model, mccmodel, lotrsparam, L, W, temp, NULL );
  id_weff = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  osdi_set_polarization( &model, 0.0, 0.0, 0.0 );
  
  double *ptr = osdi_access_ptr(&model,id_weff,&accflag, 0);
  weff = *ptr; ;
  dwcj = weff - W ;

  osdi_terminate( &model );
  return dwcj ;
}

double mcc_calcCDS_osdi( mcc_modellist   *ptmodel, 
                        elp_lotrs_param *lotrsparam,
                        double           temp, 
                        double           vbx0, 
                        double           vbx1,
                        double           L,
                        double           W
                      )
{
  osdijuncapconfig  juncapconfig ;
  double            cds, cds0, cds1 ;
  osdicharge        charge0 ;
  osdicharge        charge1 ;
#if 1
  osdi_trs      model ;
  uint32_t id_cds,accflag ;
  char  *name;

  name = namealloc( "cjd" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cds = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cds, &accflag, 0);
  osdi_set_polarization( &model, 0.0, vbx0, 0.0 );
  cds0 = *ptr;
  osdi_set_polarization( &model, 0.0, vbx1, 0.0 );
  cds1 = *ptr;

  osdi_terminate( &model );

  juncapconfig.ab = W*W ;
  juncapconfig.ls = 0.0 ;
  juncapconfig.lg = 0.0 ;

  cds = fabs((cds1*vbx1-cds0*vbx0)/(vbx1-vbx0))/juncapconfig.ab ;
#else

  juncapconfig.ab = W*W ;
  juncapconfig.ls = 0.0 ;
  juncapconfig.lg = 0.0 ;

  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vbx0, 0.0, lotrsparam, &juncapconfig, &charge0 );
  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vbx1, 0.0, lotrsparam, &juncapconfig, &charge1 );
  
  cds = fabs( (charge1.qjbd-charge0.qjbd)/(vbx1-vbx0) )/juncapconfig.ab ;
#endif


  return cds ;
}

double mcc_calcCDP_osdi( mcc_modellist *ptmodel, 
                        elp_lotrs_param *lotrsparam,
                        double temp, 
                        double vbx0, 
                        double vbx1,
                        double L,
                        double W
                      )
{
  osdijuncapconfig  juncapconfig ;
  double           cdp, cdp0, cdp1 ;
  osdicharge        charge0 ;
  osdicharge        charge1 ;
#if 1
  osdi_trs      model ;
  uint32_t id_cdp,accflag ;
  char  *name;

  name = namealloc( "cjdgat" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cdp = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cdp, &accflag, 0);
  osdi_set_polarization( &model, 0.0, vbx0, 0.0 );
  cdp0 = *ptr;
  osdi_set_polarization( &model, 0.0, vbx1, 0.0 );
  cdp1 = *ptr;
  osdi_terminate( &model );

  juncapconfig.ab = 0.0 ;
  juncapconfig.ls = W ;
  juncapconfig.lg = 0.0 ;

  cdp = fabs((cdp1*vbx1-cdp0*vbx0)/(vbx1-vbx0))/juncapconfig.ls ;
#else
  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vbx0, 0.0, lotrsparam, &juncapconfig, &charge0 );
  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vbx1, 0.0, lotrsparam, &juncapconfig, &charge1 );
  
  cdp = fabs( (charge1.qjbd-charge0.qjbd)/(vbx1-vbx0) )/juncapconfig.ls ;
#endif


  return cdp ;
}

double mcc_calcCDW_osdi( mcc_modellist *ptmodel, 
                        elp_lotrs_param *lotrsparam,
                        double temp, 
                        double vbx0, 
                        double vbx1, 
                        double L, 
                        double W 
                      )
{
  osdicharge        charge0 ;
  osdicharge        charge1 ;
  double           cdw, cdw0, cdw1 ;
  osdijuncapconfig  juncapconfig ;

  juncapconfig.ab = 0.0 ;
  juncapconfig.ls = 0.0 ;
  juncapconfig.lg = W ;

#if 1
  osdi_trs      model ;
  uint32_t id_cdw,accflag ;
  char  *name;

  name = namealloc( "cjsgat" );
 
  osdi_initialize( &model, ptmodel, lotrsparam, L, W, temp, NULL );
  id_cdw = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_cdw, &accflag, 0);
  osdi_set_polarization( &model, 0.0, vbx0, 0.0 );
  cdw0 = *ptr ;
  osdi_set_polarization( &model, 0.0, vbx1, 0.0 );
  cdw0 = *ptr ;

  osdi_terminate( &model );

  
  cdw = fabs((cdw1*vbx0-cdw0*vbx1)/(vbx1-vbx0))/juncapconfig.lg ;
#else

  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vbx0, 0.0, lotrsparam, &juncapconfig, &charge0 );
  osdi_mcc_getcharge( ptmodel, L, W, temp, 0.0, vbx1, 0.0, lotrsparam, &juncapconfig, &charge1 );
  cdw = fabs( (charge1.qjbd-charge0.qjbd)/(vbx1-vbx0) )/juncapconfig.lg ;
#endif

  return cdw ;
}


double osdi_print_jacob ( mcc_modellist   *mccmodel,
                        double           L,
                        double           W,
                        double           temp,
                        double           vbs,
                        double           vdd,
                        elp_lotrs_param *lotrsparam 
                      )
{
  osdi_trs      model ;
  double   vth ;
  uint32_t id_vth,accflag ;
  char  *name;
  double   *jacob_resist;
  double   *jacob_react ;
  double   *residual_resist;
  double   *residual_react;
  double   *spice_rhs_tran;

  name = namealloc( osdi_op_param_label[ OSDI_OP_PARAM_VTH ]);
  osdi_initialize( &model, mccmodel, lotrsparam, L, W, temp, NULL );

  double **jacobian_ptr_resist = (double**)((char*)model.idata+model.model->jacobian_ptr_resist_offset);
  jacob_resist = (double*)calloc(model.model->num_jacobian_entries, sizeof(double));
  jacob_react  = (double*)calloc(model.model->num_jacobian_entries, sizeof(double));
  residual_resist  = (double*)calloc(model.model->num_nodes, sizeof(double));
  residual_react   = (double*)calloc(model.model->num_nodes, sizeof(double));
  spice_rhs_tran   = (double*)calloc(model.model->num_nodes, sizeof(double));
  for(int i=0; i<model.model->num_jacobian_entries; i++ ) {
     jacobian_ptr_resist[i] = &jacob_resist[i];
     if(model.model->jacobian_entries[i].react_ptr_off != UINT32_MAX) {
        double **jptr = (double**)((char*)model.idata+model.model->jacobian_entries[i].react_ptr_off);
        *jptr = &jacob_react[i];
     }
  }

  id_vth = osdi_getindexparam( &model, name, OSDI_FIND_OPARAM );
  double *ptr = (double*)osdi_access_ptr(&model,id_vth, &accflag, 0);
  osdi_set_polarization( &model, vdd, vdd, vbs );
  vth = *ptr ;

  printf("vbs=%e, vdd=%e, vth=%e\n", vbs,vdd,vth);
  model.model->load_jacobian_resist(model.idata, model.mdata);
  model.model->load_jacobian_react(model.idata, model.mdata,1.0);
  model.model->load_residual_resist(model.idata, model.mdata, residual_resist);
  model.model->load_residual_react (model.idata, model.mdata, residual_react);
  model.model->load_spice_rhs_tran (model.idata, model.mdata, spice_rhs_tran, model.simdata.prev_solve, 1.0);
  for(int i=0; i<model.model->num_nodes; i++ ) {
   printf("residual %s - %e : %e ", model.model->nodes[i].name, residual_resist[i], residual_react[i]);
   printf("spice_rhs_tran -  %e \n",   spice_rhs_tran[i]);
  }
  for(int i=0; i<model.model->num_jacobian_entries; i++ ) {
    printf("%s - %s  :",
           model.model->nodes[model.model->jacobian_entries[i].nodes.node_1].name,
           model.model->nodes[model.model->jacobian_entries[i].nodes.node_2].name);
    int flag = model.model->jacobian_entries[i].flags ;
    if(flag & JACOBIAN_ENTRY_RESIST_CONST || flag & JACOBIAN_ENTRY_RESIST) printf("%e :", jacob_resist[i]);
    else printf("    :  ");
    if(flag & JACOBIAN_ENTRY_REACT_CONST || flag & JACOBIAN_ENTRY_REACT) printf("%e \n", jacob_react[i]);
    else printf("    \n");
  }
 
  osdi_terminate( &model );

  if( mccmodel->TYPE == MCC_TRANS_P )
    vth = -vth ;
    
  return vth ;
}


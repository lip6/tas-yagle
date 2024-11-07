
#define _GNU_SOURCE 1
#include <dlfcn.h>

#include <alloca.h>
#include <stdint.h>
#include "osdi.h"

#include MCC_H
#include "mcc_mod_osdi_interface.h"

char *osdi_op_param_label[ OSDI_OP_PARAM_NB ] = { "ids", "vth", "weff", "cgd", "cgdol",
  "cjd", "cjdsti", "cjdgat", "swgeo" } ;
chain_list  *OSDI_HEAD;

OsdiDescriptor 
     *dyn_osdi_psp103va = NULL,
     *dyn_osdi_psp103tva = NULL,
     *dyn_osdi_pspnqs103va = NULL;

int addosdimodelparam( osdimodelparam *param, uint32_t index, double value ) {

  if( param->n >= param->max )
    return 0;

  param->index[ param->n ] = index ;
  param->value[ param->n ] = value ;
  (param->n)++ ;

  return 1 ;
}

uint32_t osdi_getindexparam( osdi_trs *ptr, char *name, int paramtype )
{
  uint32_t      i ;
  uint32_t      base ;
  uint32_t      idx ;
  ptype_list   *ptl ;
  ht           *hash ;
  long          n ;
 
  idx = OSDI_UNDEF ;

  /* to be optimized with an hash table : name -> index */

  if( ( paramtype & OSDI_FIND_MPARAM ) == OSDI_FIND_MPARAM ) {
  
    ptl = getptype( ptr->mccmodel->USER, OSDIHASHMODELPARAM );
    if( ptl ) {
      hash = (ht*)ptl->DATA ;
    }
    else {
    
      hash = addht( 1024 );
      base = ptr->model->num_instance_params;
      
      for( i=base ; i < ptr->model->num_params; i++ ) 
        addhtitem( hash, 
                   namealloc( *(ptr->model->param_opvar[i].name) ), 
                   i 
                 );

      ptr->mccmodel->USER = addptype( ptr->mccmodel->USER, 
                                      OSDIHASHMODELPARAM, 
                                      hash 
                                    );
    }
 
    n = gethtitem( hash, name );
    if( n == EMPTYHT || n == DELETEHT )
      idx = OSDI_UNDEF ;
    else
      idx = (uint32_t)n ;
  }
  
  if( ( paramtype & OSDI_FIND_IPARAM ) == OSDI_FIND_IPARAM ) {

    ptl = getptype( ptr->mccmodel->USER, OSDIHASHINSTANCEPARAM );
    if( ptl ) {
      hash = (ht*)ptl->DATA ;
    }
    else {
    
      hash = addht( 64 );
      
      base = 0 ;
      for( i=0 ; i < ptr->model->num_instance_params ; i++ ) 
        addhtitem( hash, 
                   namealloc( (char*) *ptr->model->param_opvar[i+base].name ),
                   i+base 
                 );

      ptr->mccmodel->USER = addptype( ptr->mccmodel->USER, 
                                      OSDIHASHINSTANCEPARAM, 
                                      hash 
                                    );
    }
 
    n = gethtitem( hash, name );
    if( n == EMPTYHT || n == DELETEHT )
      idx = OSDI_UNDEF ;
    else
      idx = (uint32_t)n ;
  }
   
  if( ( paramtype & OSDI_FIND_OPARAM ) == OSDI_FIND_OPARAM ) {

    ptl = getptype( ptr->mccmodel->USER, OSDIHASHOPVARPARAM );
    if( ptl ) {
      hash = (ht*)ptl->DATA ;
    }
    else {
    
      hash = addht( 256 );
      
      base = ptr->model->num_params ;
      for( i=0 ; i < ptr->model->num_opvars ; i++ ) 
        addhtitem( hash, 
                   namealloc( (char*) *ptr->model->param_opvar[i+base].name ),
                   i+base 
                 );

      ptr->mccmodel->USER = addptype( ptr->mccmodel->USER, 
                                      OSDIHASHOPVARPARAM, 
                                      hash 
                                    );
    }
 
    n = gethtitem( hash, name );
    if( n == EMPTYHT || n == DELETEHT )
      idx = OSDI_UNDEF ;
    else
      idx = (uint32_t)n ;
  }
   
  return idx ;
}

void osdi_loadmodelparameter( osdi_trs *ptr, osdimodelparam *param )
{
  mcc_paramlist *mccparam ;
  uint32_t       i ;
  int            r ;
 
  for( mccparam = ptr->mccmodel->PARAM ; mccparam ; mccparam = mccparam->NEXT ) {
 
    i = osdi_getindexparam( ptr, mccparam->NAME, OSDI_FIND_MPARAM );
    
    if( i != OSDI_UNDEF ) {
      r = addosdimodelparam( param, i, mccparam->VALUE );
      if( !r )
        printf( "too many model parameters ! can't set param %s\n", mccparam->NAME );
    }
    else {
      i = osdi_getindexparam( ptr, mccparam->NAME, OSDI_FIND_IPARAM );
      if( i != OSDI_UNDEF ) 
          printf( "instance parameter %s is found in osdi model parameters list\n", mccparam->NAME );
      else
          printf( "model parameter %s isn't supported in osdi model\n", mccparam->NAME );
    }
  }
}

void osdi_loadinstanceparameter( osdi_trs             *ptr,
                            double           L,
                            double           W,
                            elp_lotrs_param *lotrsparam,
                            osdimodelparam   *param 
                          )
{
  char     *label[20] ;
  double    value[20] ;
  int       n ;
  int       j ;
  uint32_t  i ;
  int       r ;

  n=0 ;
  
  label[n  ] = "L" ;
  value[n++] = L ;
  label[n  ] = "W" ;
  value[n++] = W ;

  if( lotrsparam ) {
    label[n  ] = "SA" ;
    value[n++] = lotrsparam->PARAM[elpSA] ;
    label[n  ] = "SB" ;
    value[n++] = lotrsparam->PARAM[elpSB] ;
    label[n  ] = "MULT" ;
    value[n++] = lotrsparam->PARAM[elpM] ;
  }

  for( j=0 ; j<n ; j++ ) {

    i = osdi_getindexparam( ptr, namealloc(label[j]), OSDI_FIND_IPARAM );

    if( i != OSDI_UNDEF ) {
      r = addosdimodelparam( param, i, value[j] );
      if( !r ) 
        printf( "too many instance parameters ! can't set %s\n", label[j] );
    }
    else {
      printf( "instance parameter %s isn't supported in osdi model\n", label[j] );
    }
  }
}

void osdi_tuneparam( osdi_trs *trs, osditunedparam *tuned, osdimodelparam *mparam, osdimodelparam *iparam )
{
  int      i ;
  uint32_t j ;
  uint32_t idx ;
  char     *name ;

  if( tuned ) { /* to be optimized */
  
    for( i=0 ; i<tuned->n ; i++ ) {

      name = namealloc( tuned->param[i] ) ;
    
      idx = osdi_getindexparam( trs, name, OSDI_FIND_MPARAM );
      
      if( idx != OSDI_UNDEF ) { /* model parameter */
      
        for( j = 0 ; j<mparam->n ; j++ ) {
        
          if( mparam->index[j] == idx ) {
            mparam->value[j] = tuned->value[i] ;
            break ;
          }
        }
        
        if( j >= mparam->n ) 
          addosdimodelparam( mparam, idx, tuned->value[i] );
      }
      else { /* instance parameter */
    
        idx = osdi_getindexparam( trs, name, OSDI_FIND_IPARAM );
        
        if( idx != OSDI_UNDEF ) {

          for( j = 0 ; j<iparam->n ; j++ ) {
          
            if( iparam->index[j] == idx ) {
              iparam->value[j] = tuned->value[i] ;
              break ;
            }
          }
          
          if( j >= iparam->n ) 
            addosdimodelparam( iparam, idx, tuned->value[i] );
        }
        else {
          printf( "unknown tuned parameter %s\n", name );
        }
      }
    }
  }
}

void freeosdiparam( osdimodelparam *param )
{
  mbkfree( param->index );
  mbkfree( param->value );
  mbkfree( param );
}

osdimodelparam* allocosdiparam( int size )
{
  osdimodelparam *mparam ;

  mparam = (osdimodelparam*)mbkalloc( sizeof( osdimodelparam ) );
  mparam->max   = size ;
  mparam->n     = 0 ;
  mparam->index = (uint32_t*)mbkalloc( sizeof( uint32_t ) * mparam->max );
  mparam->value = (double*) mbkalloc( sizeof( double )  * mparam->max );

  return mparam ;
}

osdimodelparam* duposdimodelparam( osdimodelparam *p )
{
  osdimodelparam *n ;
  
  n = allocosdiparam( p->max );
  n->n = p->n ;
  memcpy( n->index, p->index, sizeof( uint32_t ) * p->max );
  memcpy( n->value, p->value, sizeof( double  ) * p->max );

  return n ;
}

void *osdi_access_ptr(osdi_trs *ptr, int index, int *aflag, int write) {
    int flag = ptr->model->param_opvar[index].flags;
    int access_flg;
    if(write) access_flg = ACCESS_FLAG_SET;
    *aflag = flag; //return parameter flag
    if (index < ptr->model->num_instance_params)
      access_flg = access_flg | ACCESS_FLAG_INSTANCE;

    return ptr->model->access(ptr->idata, ptr->mdata, index, access_flg );
}

void osdi_loadmodel( osdi_trs             *ptr,
                double           L,
                double           W,
                osditunedparam   *tuned,
                elp_lotrs_param *lotrsparam
              )
{
  osdimodelparam  *mparam ;
  osdimodelparam  *iparam ;
  osdimodelparam  *tmparam ;
  osdimodelparam  *tiparam ;
  OsdiInitInfo  *info_inst, *info_model;
  ptype_list     *ptl ;
  int  flag;
  void *aptr;
  double temp = 25.0;

  ptl = getptype( ptr->mccmodel->USER, OSDICACHEMODEL );
  if( !ptl ) {
    mparam = allocosdiparam( ptr->model->num_params - ptr->model->num_instance_params );
    osdi_loadmodelparameter( ptr, mparam );
    if( V_BOOL_TAB[ __AVT_USE_CACHE_OSDI ].VALUE )
      ptr->mccmodel->USER = addptype( ptr->mccmodel->USER, 
                                      OSDICACHEMODEL, 
                                      mparam 
                                    );
  }
  else
    mparam = (osdimodelparam*)ptl->DATA ;

  ptl = getptype( ptr->mccmodel->USER, OSDICACHEINSTANCE );
  if( !ptl ) {
    iparam = allocosdiparam( ptr->model->num_instance_params );
    osdi_loadinstanceparameter( ptr, L, W, lotrsparam, iparam );
    if( V_BOOL_TAB[ __AVT_USE_CACHE_OSDI ].VALUE )
      ptr->mccmodel->USER = addptype( ptr->mccmodel->USER, 
                                      OSDICACHEINSTANCE, 
                                      iparam 
                                    );
  }
  else
    iparam = (osdimodelparam*)ptl->DATA ;

  tmparam = duposdimodelparam( mparam );
  tiparam = duposdimodelparam( iparam );

  if(!ptr->mdata) ptr->mdata = calloc(1, ptr->model->model_size);
  if(!ptr->idata) ptr->idata = calloc(1, ptr->model->instance_size);

#if 1
  int typeidx = osdi_getindexparam( ptr, namealloc("TYPE"), OSDI_FIND_MPARAM );
  aptr = osdi_access_ptr(ptr, typeidx, &flag, 1);
  if( ptr->mccmodel->TYPE == MCC_TRANS_P )
    *(int*) aptr = -1;
  else
    *(int*) aptr = +1;
#endif

  osdi_tuneparam( ptr, tuned, tmparam, tiparam );
  for(int i=0; i<tmparam->n; i++) {
    aptr = osdi_access_ptr(ptr, tmparam->index[i], &flag, 1);
    switch( flag & PARA_TY_MASK ) {
      case PARA_TY_REAL:
        *(double*) aptr = tmparam->value[i];
        break;
      case PARA_TY_INT:
        *(int*) aptr = (int)tmparam->value[i];
        break;
     default:
        printf("Not supported type detected\n");
        exit(1);
    }
  }
  info_model = (OsdiInitInfo*)calloc(1,sizeof(OsdiInitInfo));
  ptr->model->setup_model( NULL, ptr->mdata, NULL, info_model );
    if (info_model->num_errors)
    for(int i=0; i<info_model->num_errors; i++)
    switch( info_model->errors[i].code ) {
    case INIT_ERR_OUT_OF_BOUNDS :
      printf( "  ->Model Init out of bounds error %s:%d\n\n" , ptr->mccmodel->NAME, info_model->errors[i].payload.parameter_id);
      return 0 ;
      break ;
    default :
      printf( "  ->unknown error\n\n" );
      return 0 ;
      break ;
  }

  for(int i=0; i<tiparam->n; i++) {
    aptr = osdi_access_ptr(ptr, tiparam->index[i], &flag, 1);
    switch( flag & PARA_TY_MASK ) {
      case PARA_TY_REAL:
        *((double*) aptr) = tiparam->value[i];
        break;
      case PARA_TY_INT:
        *((int*) aptr) = (int)tiparam->value[i];
        break;
     default:
        printf("Not supported type detected\n");
        exit(1);
    }
  }
  info_inst = (OsdiInitInfo*)calloc(1,sizeof(OsdiInitInfo));
  ptr->model->setup_instance( NULL,
                                    ptr->idata, 
                                    ptr->mdata, 
                                    temp,
                                    ptr->model->num_nodes,
                                    NULL,
                                    info_inst
                                  );
  if (info_inst->num_errors)
  for(int i=0; i<info_inst->num_errors; i++)
  switch( info_inst->errors[i].code ) {
    case INIT_ERR_OUT_OF_BOUNDS :
      printf( "  ->Inst Init out of bounds error %s:%d\n\n" , ptr->mccmodel->NAME, info_inst->errors[i].payload.parameter_id);
      return 0 ;
      break ;
    default :
      printf( "  ->unknown error\n\n" );
      return 0 ;
      break ;
  }
  uint32_t *node_mapping = 
       (uint32_t *)((char*)ptr->idata + ptr->model->node_mapping_offset);
  for(int i=0;i<ptr->model->num_nodes;i++)
       node_mapping[i] = i;

  for(int i=0; i < ptr->model->num_collapsible; i++) {
      int to,from,tmp;
      if (0&& ! (bool*)((char*) ptr->idata + ptr->model->collapsed_offset)[i] )
        continue;
      from = ptr->model->collapsible[i].node_1;
      to   = ptr->model->collapsible[i].node_2;

      if (to == UINT32_MAX || node_mapping[to] == UINT32_MAX)
        continue;


      if ( to != UINT32_MAX && node_mapping[from] < node_mapping[to] ) {
        tmp  = from;
        from = to;
        to   = tmp;
      }
     from = node_mapping[from];
     if (to != UINT32_MAX) 
       to = node_mapping[to];

    for (uint32_t j = 0; j< ptr->model->num_nodes; j++) {
        if (node_mapping[j] == from)
          node_mapping[j] = to;
        else if (node_mapping[j] > from && node_mapping[j] != UINT32_MAX) 
          node_mapping[j] -= 1;
    }
  }

  freeosdiparam( tmparam );
  freeosdiparam( tiparam );

  if( ! V_BOOL_TAB[ __AVT_USE_CACHE_OSDI ].VALUE ) {
    freeosdiparam( mparam );
    freeosdiparam( iparam );
  }
}

int osdi_initialize( osdi_trs             *ptr,
                mcc_modellist   *mccmodel, 
                elp_lotrs_param *lotrsparam,
                double           L,
                double           W,
                double           temp,
                osditunedparam   *tuned 
              )
{
  int         i ;
  ptype_list *ptl ;
  osdicachemodel *cache = NULL;
  OsdiInitInfo  *info_inst, *info_model;
  double *charge, *solve, *info_state;
  uint32_t from, to, tmp;

  ptr->mdata = NULL ;
  ptr->idata = NULL ;
  
  if( mccmodel->MODELTYPE == MCC_PSPVA )
    ptr->model = dyn_osdi_psp103va ;
  else
  if( mccmodel->MODELTYPE == MCC_PSPTVA )
    ptr->model = dyn_osdi_psp103tva ;
  else
    ptr->model = dyn_osdi_pspnqs103va ;

  ptr->mccmodel = mccmodel ;

  info_inst = (OsdiInitInfo*)calloc(1,sizeof(OsdiInitInfo));
  info_model = (OsdiInitInfo*)calloc(1,sizeof(OsdiInitInfo));
  charge = (double*)calloc(ptr->model->num_nodes, sizeof(double));
  solve  = (double*)calloc(ptr->model->num_nodes, sizeof(double));
  info_state = (double*)calloc(ptr->model->num_states, sizeof(double));

  ptr->simdata.prev_solve = solve ;
  ptr->simdata.prev_state = ptr->simdata.next_state = info_state ;
  ptr->simdata.flags = CALC_REACT_RESIDUAL | CALC_RESIST_RESIDUAL | CALC_OP;
 
  if( V_BOOL_TAB[ __AVT_USE_CACHE_OSDI ].VALUE && !tuned ) {
    ptl = getptype( mccmodel->USER, OSDICACHETRS ) ;
    if( ptl ) {
      cache = (osdicachemodel*)ptl->DATA ;
      ptr->mdata     = cache->mdata ;
      ptr->idata     = cache->idata ;
      ptr->cleanmidata = 0 ;
    }
  }
  if( !ptr->mdata ) 
    ptr->mdata = calloc( 1, ptr->model->model_size );

  osdi_loadmodel( ptr, L, W, tuned, lotrsparam );

    ptr->model->setup_model( NULL, ptr->mdata, NULL, info_model );
    if (info_model->num_errors)
    for(int i=0; i<info_model->num_errors; i++)
    switch( info_model->errors[i].code ) {
    case INIT_ERR_OUT_OF_BOUNDS :
      printf( "  ->Model Init out of bounds error %s:%d\n\n" , ptr->mccmodel->NAME, info_model->errors[i].payload.parameter_id);
      return 0 ;
      break ;
    default :
      printf( "  ->unknown error\n\n" );
      return 0 ;
      break ;
    }

  if( !ptr->idata ) 
  ptr->idata = calloc( 1, ptr->model->instance_size );

  uint32_t *node_mapping = 
       (uint32_t *)((char*)ptr->idata + ptr->model->node_mapping_offset);
  for(i=0;i<ptr->model->num_nodes;i++)
       node_mapping[i] = i;

    
  ptr->model->setup_instance( NULL,
                                    ptr->idata, 
                                    ptr->mdata, 
                                    temp,
                                    ptr->model->num_nodes,
                                    NULL,
                                    info_inst
                                  );
  if (info_inst->num_errors)
    for(int i=0; i<info_inst->num_errors; i++)
    switch( info_inst->errors[i].code ) {
    case INIT_ERR_OUT_OF_BOUNDS :
      printf( "  ->Inst Init out of bounds error %s:%d\n\n" , ptr->mccmodel->NAME, info_inst->errors[i].payload.parameter_id);
      return 0 ;
      break ;
    default :
      printf( "  ->unknown error\n\n" );
      return 0 ;
      break ;
    }
  for(i=0; i < ptr->model->num_collapsible; i++) {
      if (0&& ! (bool*)((char*) ptr->idata + ptr->model->collapsed_offset)[i] )
        continue;
      from = ptr->model->collapsible[i].node_1;
      to   = ptr->model->collapsible[i].node_2;

      if (to == UINT32_MAX || node_mapping[to] == UINT32_MAX)
        continue;


      if ( to != UINT32_MAX && node_mapping[from] < node_mapping[to] ) {
        tmp  = from;
        from = to;
        to   = tmp;
      }
     from = node_mapping[from];
     if (to != UINT32_MAX) 
       to = node_mapping[to];

    for (uint32_t j = 0; j< ptr->model->num_nodes; j++) {
        if (node_mapping[j] == from)
          node_mapping[j] = to;
        else if (node_mapping[j] > from && node_mapping[j] != UINT32_MAX) 
          node_mapping[j] -= 1;
    }
  }


    uint32_t code = ptr->model->eval (  NULL,
                        ptr->idata, 
                        ptr->mdata, 
                        &(ptr->simdata)
                     );
    if (code & EVAL_RET_FLAG_LIM)
      printf( "  eval LIM\n\n" );
    if (code & EVAL_RET_FLAG_FATAL)
      printf( "  eval FATAL\n\n" );
    if (code & EVAL_RET_FLAG_FINISH)
      printf( "  eval FINISH\n\n" );
    if (code & EVAL_RET_FLAG_STOP)
      printf( "  eval STOP\n\n" );

    if( V_BOOL_TAB[ __AVT_USE_CACHE_OSDI ].VALUE && !tuned ) {
      ptr->cleanmidata = 0 ;
      if(!cache) {
         cache = (osdicachemodel*)mbkalloc( sizeof( osdicachemodel ) );
         cache->mdata     = ptr->mdata ;
         cache->idata     = ptr->idata ;
         mccmodel->USER = addptype( mccmodel->USER, OSDICACHETRS, cache );
      }
    }
    else
      ptr->cleanmidata = 1 ;

  for( i=0 ; i<OSDI_OP_PARAM_NB ; i++ )
    ptr->tabid[i] = OSDI_UNDEF ;
      
  return 1 ;
}

void osdi_terminate( osdi_trs *ptr )
{
  if( ptr->cleanmidata ) {
      mbkfree( ptr->mdata );
      mbkfree( ptr->idata );
  }
}

uint32_t osdi_get_id_param( osdi_trs *ptr, int param )
{
  uint32_t      base ;
  uint32_t      i ;

  if( ptr->tabid[ param ] == OSDI_UNDEF ) {
    base = ptr->model->num_params ;
    for( i = 0 ; i < ptr->model->num_opvars ; i++ ) {
      if( !strcasecmp( *ptr->model->param_opvar[i+base].name, 
                       osdi_op_param_label[ param ] ) )
        break ;
    }
    if( i < ptr->model->num_opvars )
      ptr->tabid[ param ] = i ;
    else {
      printf( "can't find operating point #%d\n", param );
      exit(1);
    }
  }

  return ptr->tabid[ param ] ;
}

void osdi_set_polarization( osdi_trs *ptr, double vgs, double vds, double vbs )
{
  uint32_t       i ;
  int            j ;
  double         v[2];
  double         vg, vd, vs, vb ;
  int            flag ;

  vs = 0.0 ;
  vg = vgs + vs ;
  vd = vds + vs ;
  vb = vbs + vs ;

  ptr->simdata.prev_solve[0] = vd;
  ptr->simdata.prev_solve[1] = vg;
  ptr->simdata.prev_solve[2] = vs;
  ptr->simdata.prev_solve[3] = vb;
  
  ptr->simdata.flags = CALC_REACT_RESIDUAL | CALC_RESIST_RESIDUAL | CALC_OP ;

  uint32_t code = ptr->model->eval( NULL, ptr->idata, ptr->mdata, &ptr->simdata );
    if (code & EVAL_RET_FLAG_LIM)
      printf( "  eval LIM\n\n" );
    if (code & EVAL_RET_FLAG_FATAL)
      printf( "  eval FATAL\n\n" );
    if (code & EVAL_RET_FLAG_FINISH)
      printf( "  eval FINISH\n\n" );
    if (code & EVAL_RET_FLAG_STOP)
      printf( "  eval STOP\n\n" );


}

int osdi_loaddynamiclibrary( void )
{
  OsdiDescriptor *handle ;
  chain_list *curr;
  uint32_t num, i;

  Dl_info m ;

  printf( "Scanning dynamic libraries .\n" );
  curr = OSDI_HANDLE;
  if( !curr ) {
    printf( "can't load dynamic library \n" );
    return 0 ;
  }
  while ( curr ) {
     num = *(uint32_t*)dlsym( curr->DATA, "OSDI_NUM_DESCRIPTORS" );
     handle = dlsym( curr->DATA, "OSDI_DESCRIPTORS" );
     for (i = 0; i < num; i++ ) {
      if( !handle ) {
        printf( "can't find dynamic function !\n" );
        return 0 ;
      }
      if( dladdr( handle->setup_model, &m ) != -1 )
        printf( "  -> %s\n", m.dli_fname );
      else
        printf( "  -> can't get information about the loaded library.\n" );
      if( strcmp(handle->name, "PSP103VA") == 0) dyn_osdi_psp103va = handle;
      if( strcmp(handle->name, "PSP103TVA") == 0) dyn_osdi_psp103tva = handle;
      if( strcmp(handle->name, "PSPNQS103VA") == 0) dyn_osdi_pspnqs103va = handle;
   
       addchain(OSDI_HEAD, handle);
       handle++;
     }
     curr = curr->NEXT;
  }
  
  return 1 ;
}

void mcc_clean_osdi_interface( mcc_modellist *mccmodel, int check )
{
  ptype_list *ptl ;
  osdicachemodel *cache ;

  ptl = getptype( mccmodel->USER, OSDICACHEINSTANCE );
  if( ptl ) {
    if( check )
      printf( "warning : non empty instance cache.\n" );
    freeosdiparam( (osdimodelparam*)ptl->DATA );
    mccmodel->USER = delptype( mccmodel->USER, OSDICACHEINSTANCE );
  }

  ptl = getptype( mccmodel->USER, OSDICACHETRS );
  if( ptl ) {
    if( check )
      printf( "warning : non empty model cache.\n" );
    cache = (osdicachemodel*)ptl->DATA ;
    mbkfree( cache->mdata );
    mbkfree( cache->idata );
    mbkfree( cache );
    mccmodel->USER = delptype( mccmodel->USER, OSDICACHETRS );
  }
}

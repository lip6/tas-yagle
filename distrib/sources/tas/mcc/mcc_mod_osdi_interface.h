
#define OSDI_UNDEF 314159265

#define OSDI_OP_PARAM_IDS  0
#define OSDI_OP_PARAM_VTH  1
#define OSDI_OP_PARAM_WEFF 2
#define OSDI_OP_PARAM_CGD  3
#define OSDI_OP_PARAM_CGDOL  4
#define OSDI_OP_PARAM_CJD  5
#define OSDI_OP_PARAM_CJDSTI  6
#define OSDI_OP_PARAM_CJDGAT  7
#define OSDI_OP_PARAM_SWGEO  8
#define OSDI_OP_PARAM_NB   9

#define OSDI_FIND_MPARAM 0x01
#define OSDI_FIND_IPARAM 0x02
#define OSDI_FIND_OPARAM 0x04

/* USER field of mcc_modellist* */
#define OSDIHASHINSTANCEPARAM 0x00001
/* index des parametres d'instance. global à un modele                  */
#define OSDIHASHMODELPARAM    0x00002
/* index des parametres de modele. global à un modele                   */
#define OSDIHASHOPVARPARAM    0x00007
/* index des parametres d'operational. global à un modele                   */
#define OSDICACHECHARGE       0x00003
/* charges. depend de chaque instance                                   */
#define OSDICACHEMODEL        0x00004
/* tableau des couples idx-valeur des parametres de modele. global à un 
   modèle                                                               */
#define OSDICACHEOPVARS      0x00008
/* tableau des couples idx-valeur des parametres d'operations. 
   depend de chaque instance                                            */
#define OSDICACHEINSTANCE     0x00005
/* tableau des couples idx-valeur des parametres d'instance. 
   depend de chaque instance                                            */
#define OSDICACHETRS          0x00006
/* cache du model osdi integrale. depend de chaque instance              */

typedef struct {
  void              *mdata ;
  void              *idata ;
} osdicachemodel ;

typedef struct {
  OsdiDescriptor    *model ;
  mcc_modellist     *mccmodel ;
  void              *mdata ;
  void              *idata ;
  OsdiSimInfo        simdata ;
  uint32_t           tabid[OSDI_OP_PARAM_NB] ;
  char               cleanmidata ;
} osdi_trs ;

typedef struct {
  uint32_t         *index ;
  double          *value ;
  uint32_t          n ;
  uint32_t          max ;
} osdimodelparam ;

typedef struct {
  int               n ;
  char            **param ;
  double          *value ;
} osditunedparam ;


int osdi_initialize( osdi_trs             *ptr, 
                mcc_modellist   *mccmodel, 
                elp_lotrs_param *lotrsparam,
                double           L,
                double           W,
                double           temp,
                osditunedparam   *tuned
              );
extern char* osdi_op_param_label[];
void osdi_terminate( osdi_trs *ptr );
uint32_t osdi_get_id_param( osdi_trs *ptr, int param );
void osdi_set_polarization( osdi_trs *ptr, double vgs, double vds, double vbs );
int osdi_loaddynamiclibrary( void );
void mcc_clean_osdi_interface( mcc_modellist *mccmodel, int check );
void *osdi_access_ptr(osdi_trs *ptr, int index, int *aflag, int write) ;
uint32_t osdi_getindexparam( osdi_trs *ptr, char *name, int paramtype );
osdimodelparam* duposdimodelparam( osdimodelparam *p );

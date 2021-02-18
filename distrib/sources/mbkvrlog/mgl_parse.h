/****************************************************************************/
/*                                                                          */
/*                      Chaine de CAO & VLSI   AVERTEC                      */
/*                                                                          */
/*    Produit : Grammar for Verilog                                         */
/*    Fichier : mgl_parse.h                                                 */
/*                                                                          */
/*    (c) copyright 2000 AVERTEC                                            */
/*    Tous droits reserves                                                  */
/*                                                                          */
/*    Auteur(s) : Anthony LESTER                                            */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

extern int           MGL_PARSE_ESC_VECTORS;

extern int           MGL_USE_LIBRARY;

extern char         *MGL_VDD;
extern char         *MGL_VSS;

extern FILE    *mgl_scompin;

extern int      mgl_scompdebug;

#ifdef __cplusplus
extern "C" {
#endif
void             *mgl_initparser(FILE *ptinbuf);
void              mgl_delparser(void *parm);
mgl_scompcontext *mgl_getcontext(void *parm);
int               mgl_scompparse(void *parm);
void 	          mgl_scompclean(mgl_scompcontext *context);
#ifdef __cplusplus
}
#endif


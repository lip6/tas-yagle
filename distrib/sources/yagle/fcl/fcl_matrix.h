/****************************************************************************/
/*                                                                          */
/*                      Chaine de CAO & VLSI   Alliance                     */
/*                                                                          */
/*    Produit : FCL  v1.02                                                  */
/*    Fichier : fcl_matrix.h                                                */
/*                                                                          */
/*    (c) copyright 1996 Laboratoire MASI equipe CAO & VLSI                 */
/*    Tous droits reserves                                                  */
/*    Support : e-mail alliance-support@asim.lip6.fr                        */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#define BUFFER  100 

typedef struct matrice {
     struct matrice *NEXT;
     struct matrice *PREV;
     void           **tab;
 } matrice_list;     

int           fclFindCorrespondingLosig(matrice_list *ptmatrice, int niveau, int connecteur, fcl_label prevcon, int precniveau, int loop, int niveau_donedrainsource);
int           fclFindCorrespondingLotrs(matrice_list *ptmatrice, int niveau, int connecteur, int loop);
int           fclFindCorrespondingLoins(matrice_list *ptmatrice, int niveau, int connecteur, int loop);

chain_list   *fils_sig( losig_list *ptlosig, locon_list *ptlocon );

void          tab_ajoute( matrice_list *ptmatrice, int niveau, void *ptdata );
matrice_list *tab_del( matrice_list *ptmatrice_begin, matrice_list *ptmatrice );
matrice_list *tab_init();
void          tab_free( matrice_list *ptmatrice );
matrice_list *tab_copy( matrice_list *ptmatrice_begin, matrice_list *ptmatrice, int niveau );
chain_list   *fils_sig( losig_list *ptlosig, locon_list *ptlocon );
void          printmatrice( matrice_list *ptmatrice, int niveau );
void          displaymatrice( matrice_list *ptmatrice, int niveau, int full );

extern int           niveau_last;


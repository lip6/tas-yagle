/****************************************************************************/
/*                                                                          */
/*                      Chaine de CAO & VLSI   Alliance                     */
/*                                                                          */
/*    Produit : FCL  v1.02                                                  */
/*    Fichier : fcl_util.h                                                  */
/*                                                                          */
/*    (c) copyright 1996 Laboratoire MASI equipe CAO & VLSI                 */
/*    Tous droits reserves                                                  */
/*    Support : e-mail alliance-support@asim.lip6.fr                        */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

void        fclExit();
fcl_label   CountChain( chain_list *pthead );
fcl_label   getlabeltrs( lotrs_list *ptlotrs );
fcl_label   getlabellocon( locon_list *ptlocon );
fcl_label   getlabelsig( losig_list *ptlosig );
fcl_label   getlabelsigtrs( lotrs_list *ptlotrs, losig_list *ptlosig );
void        fclCorruptLotrs( lotrs_list *ptlotrs );
void        fclCorruptLosig( losig_list *ptlosig );
lotrs_list *whichtrs( lofig_list *ptlofig, locon_list *ptlocon );
void        fclPrintTime( struct rusage *start, struct rusage *end, time_t rstart, time_t rend );
void        fclChrono( struct rusage *t, time_t *rt );

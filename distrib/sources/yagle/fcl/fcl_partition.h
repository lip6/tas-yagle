/****************************************************************************/
/*                                                                          */
/*                      Chaine de CAO & VLSI   Alliance                     */
/*                                                                          */
/*    Produit : FCL  v1.02                                                  */
/*    Fichier : fcl_partition.h                                             */
/*                                                                          */
/*    (c) copyright 1996 Laboratoire MASI equipe CAO & VLSI                 */
/*    Tous droits reserves                                                  */
/*    Support : e-mail alliance-support@asim.lip6.fr                        */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

typedef int fcl_label;

typedef struct partition {
     struct partition *NEXT;
     char              TYPE;    /* soit T (transistor), soit S (signal) */
     fcl_label         LABEL;
     int               NBELEM;
     void             *DATA;
 } partition_list;     

extern partition_list *addpartition( partition_list *ptpartition, void *ptdata, fcl_label valeur, char type );
extern partition_list *getpartition( partition_list *ptpartition, fcl_label valeur, char type );
extern partition_list *delpartition( partition_list *ptpartition, fcl_label valeur, char type );
extern partition_list *smallestpartition( partition_list *ptpartition );
extern partition_list *comparepartition( partition_list *ptpartition_c, partition_list *ptpartition_m );
extern partition_list *removepartitionelement( partition_list *ptpartition, fcl_label valeur, char type );

extern void printpartition( partition_list *ptpartition );
extern void displaypartition( partition_list *ptpartition );


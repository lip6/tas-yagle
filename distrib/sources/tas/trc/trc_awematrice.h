/****************************************************************************/
/*                                                                          */
/*                      Chaine de CAO & VLSI   AVERTEC                      */
/*                                                                          */
/*    Produit : RCX - AWE support.                                          */
/*    Fichier : trc_awematrice.h                                            */
/*                                                                          */
/*    (c) copyright 2000 AVERTEC                                            */
/*    Tous droits reserves                                                  */
/*                                                                          */
/*    Auteur(s) : Gr�goire Avot                                             */
/*                                                                          */
/****************************************************************************/

/* CVS informations :

Revision : $Revision: 1.6 $
Author   : $Author: gregoire $
Date     : $Date: 2002/12/06 14:52:38 $

*/

/* Le champs data est d�fini par :

        [ 0 1 2 ]
        [ 3 4 5 ]
        [ 6 7 8 ]

  on acc�de � l'�l�ment (ligne,colonne) par la macro 
  MATELEM( matrice, ligne, colonne )
*/

#define MATELEM(a,l,c) (a->data[c+l*a->col]) 

/*
 * La r�servation m�moire est unique, elle correspond � une matrice de dimension
 * MAT_ALLOCLINE et MAT_ALLOCCOL. Seule une partie correspondant � la taille 
 * r��elle de la matrice est utilis�e. */

#define MAT_ALLOCLINE 10
#define MAT_ALLOCCOL  10

typedef struct smatrice {
  struct smatrice *NEXT;
  int lin, col ;
  RCXFLOAT *data ;
} matrice ;

matrice*  mat_create( int, int );
void      mat_free( matrice* );
int       mat_solve( matrice*, matrice*, matrice* );
void      mat_mult(matrice*, matrice*, matrice* );
matrice*  mat_dup( matrice* );
void      mat_sub( matrice*, matrice*, matrice* );
RCXFLOAT mat_sq( matrice* );

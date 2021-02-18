/****************************************************************************/
/*                                                                          */
/*                      Chaine de CAO & VLSI   AVERTEC                      */
/*                                                                          */
/*    Produit : RCX - AWE support.                                          */
/*    Fichier : trc_awepolynomes.h                                          */
/*                                                                          */
/*    (c) copyright 2000 AVERTEC                                            */
/*    Tous droits reserves                                                  */
/*                                                                          */
/*    Auteur(s) : Gr�goire Avot                                             */
/*                                                                          */
/****************************************************************************/

/* CVS informations :

Revision : $Revision: 1.5 $
Author   : $Author: gregoire $
Date     : $Date: 2002/12/06 14:52:38 $

*/


/*
 * Manipulation de polynomes.
 *
 * Les coefficients d'un polynome de degr�s n sont stock� dans un tableau
 * de n+1 �lements : a0, a1, a2, ..., an
 */

int poly_findroot ( RCXFLOAT*, int, RCXFLOAT* );

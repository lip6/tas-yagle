/* 
 * This file is part of the Alliance CAD System
 * Copyright (C) Laboratoire LIP6 - D�partement ASIM
 * Universite Pierre et Marie Curie
 * 
 * Home page          : http://www-asim.lip6.fr/alliance/
 * E-mail support     : mailto:alliance-support@asim.lip6.fr
 * 
 * This library is free software; you  can redistribute it and/or modify it
 * under the terms  of the GNU Library General Public  License as published
 * by the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * Alliance VLSI  CAD System  is distributed  in the hope  that it  will be
 * useful, but WITHOUT  ANY WARRANTY; without even the  implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 * 
 * You should have received a copy  of the GNU General Public License along
 * with the GNU C Library; see the  file COPYING. If not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


/*******************************************************************************
*                                                                              *
*  Tool        : Spice parser / driver v 7.00                                  *
*  Author(s)   : Gregoire AVOT                                                 *
*  Updates     : March, 18th 1998                                              *
*                                                                              *
*******************************************************************************/

#ifndef mbk_int
#define mbk_int


/* Lors qu'on ajoute quelque chose dans la liste, le pointeur DATA ne doit pas
   etre NULL, cette valeur etant utilis�e en interne sur ce champs */

#define SPI_MBK_TABLEINTMAX  256
#define SPI_MBK_TABLEINTMASK 0x000000FF

mbk_tableint*	creatmbk_tableint __P(());
/* Cr�e une nouvelle table						*/

int		setmbk_tableint __P(( mbk_tableint *table, int value, int data ));
/* Ajoute un �l�ment dans la table. Renvoie 1 si l'�l�ment est nouveau,
   0 sinon (l'�l�ment pr�c�dent est alors remplac�)			*/

void            removembk_tableint __P(( mbk_tableint*, int ));

int		tstmbk_tableint __P(( mbk_tableint *table, int value ));
/* Renvoie l'�l�ment recherch� par la valeur				*/

void		freembk_tableint __P(( mbk_tableint *table ));
/* Lib�re la table							*/

int             mbk_scanint __P(( mbk_tableint *table, int n ));
/* Parcour des �l�ments m�moris�s dans la table. On commence par mettre
  l'entier d'entr�e � 0, puis pour les elements suivant on met la valeur
  du noeud retourn� pr�c�dement. La fin est indiqu� par un retour 0	*/


#endif //spi_int

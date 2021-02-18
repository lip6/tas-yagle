/*
 * This file is part of the Alliance CAD System
 * Copyright (C) Laboratoire LIP6 - D�partement ASIM
 * Universite Pierre et Marie Curie
 *
 * Home page          : http://www-asim.lip6.fr/alliance/
 * E-mail support     : mailto:alliance-support@asim.lip6.fr
 *
 * This progam is  free software; you can redistribute it  and/or modify it
 * under the  terms of the GNU  General Public License as  published by the
 * Free Software Foundation;  either version 2 of the License,  or (at your
 * option) any later version.
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

/*
 * Tool    : ABL, BDD, HT Librarie
 * Date    : 1991,92
 * Author  : Luc Burgun
 * Modified by Czo <Olivier.Sirol@lip6.fr> 1996,97
 */



#ident "$Id: log_thashloc.c,v 1.4 2005/05/03 08:26:33 antony Exp $"

/*--------------------------------------------------------------------------
   la table de hachage local 
   la version du 14.12.90 
  -------------------------------------------------------------------------- */
#include <stdlib.h>
#include  MUT_H
#include LOG_H
#include AVT_H


/* les fonction de l'utilisateurs :
   ------------------------------

   initialisation au depart des pointeurs de vertex a NULL.
   Les fonctions qui utilisent le hachage peuvent renvoye DDBTABLE_PLEINE
   si la table est trop remplie.
   Il n'y a pas de destruction possible d'un element.

   a. creation de table 

   pTableLoc createTabLoc(len)
   int len;

   b. destruction de la table

   void destroyTabLoc(pTab)
   pTableLoc pTab;

   c. re-allocation d'une table de hachage

   reAllocTabLoc(pTab)
   pTableLoc pTab;

   d. recherche d'un element

   pNode searchTabLoc(pTab,high,low)
   pTableLoc pTab;
   pVertexLoc high,low;

   e. ajout d'un element

   int addTabLoc(pTab,pLoc)
   pTableLoc pTab;
   pVertexLoc pLoc;

   f. affichage total d'une table 

   void displayLoc(pTab)
   pTableLoc pTab;


 */

/*-------------------- la fonction de hachage ---------------------------- */

long 
hashLoc (high, low)
     pNode high, low;
{
  return (abs (high->index + low->index +
	       (((long) high) >> 4) + (((long) low) >> 5) +
	       (long) high + (long) low));
}

/*--------------- La table de hachage pour des LOC ------------ */

/* la fonction de creation de table de hachage pour les LOC .
   On alloue en premier lieu une structure TABLE, puis une table
   qui n'est rien d'autre qu'un tableau de pointeurs de LOC. Il est
   donc possible de travailler avec plusieurs table a la fois. */

pTableLoc 
createTabLoc (len)
     int len;
{
  pTableLoc pTab;
  pVertexLoc pLoc;
  int i;

  if (!(pTab = (pTableLoc) mbkalloc (sizeof (struct tableLoc))))
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"158");
//      printf ("allocation impossible\n");
//      EXIT (-1);
    }
  pTab->lenTabLoc = len;

  if (!(pLoc = (pVertexLoc) mbkalloc (len * sizeof (struct vertexLoc))))
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"159");
//      printf ("allocation impossible\n");
//      EXIT (-1);
    }
  pTab->pLoc = pLoc;
  for (i = 0; i < len; i++)
    {
      pLoc->oper = EMPTYTH;	/* -1 */
      pLoc++;
    }
  return (pTab);
}

   /* destruction d'une table de hachage */

void 
destroyTabLoc (pTab)
     pTableLoc pTab;
{
  mbkfree (pTab->pLoc);
  mbkfree (pTab);
}

  /* recherche d'un element dans la table
     renvoie NULL si la recherche echoue. */

pNode 
searchTabLoc (pTab, high, low, oper)
     pTableLoc pTab;
     pNode high, low;
     short oper;
{
  pVertexLoc pLoc;
  int indice;

  /* un seul acces permis */

  indice = hashLoc (high, low) % pTab->lenTabLoc;
  pLoc = pTab->pLoc;
  if (pLoc[indice].oper == oper)
    if (high == pLoc[indice].high && low == pLoc[indice].low)
      return (pLoc[indice].father);
  return (NULL);
}

 /* ajout d'un element a la table */


long 
addTabLoc (pTab, high, low, father, oper)
     pTableLoc pTab;
     pNode high, low, father;
     short oper;
{
  pVertexLoc pLoc;
  long indice;

  /* un seul acces permis */

  indice = hashLoc (high, low) % pTab->lenTabLoc;
  pLoc = pTab->pLoc;
  pLoc[indice].oper = oper;
  pLoc[indice].high = high;
  pLoc[indice].low = low;
  pLoc[indice].father = father;
  return (indice);		/* retourne la place utilisee */
}


/* affichage des elements de la table */

void 
displayLoc (pTab)
     pTableLoc pTab;
{
  int i;
  int co = 0;
  pVertexLoc pLoc;

  pLoc = pTab->pLoc;

  printf ("\n---------------------------------------------------------------");
  printf ("---------\n               AFFICHAGE DE LA TABLE DE HACHAGE\n\n");

  for (i = 0; i < pTab->lenTabLoc; i++)
    {
      if (pLoc[i].oper != EMPTYTH)
	{
	  co++;
	  printf ("****** indice %d ****** \n", i);
	  printf ("HIGH %ld LOW %ld FATHER %ld\n", (long) pLoc[i].high, (long) pLoc[i].low, (long) pLoc[i].father);
	  printf ("\n");
	}
    }
  printf ("\n****** Nombre de noeuds dans la table  = %d\n", co);
}



/*------------------------------------------------------------------------------
videTabLoc       :vide une table de hachage locale.
-------------------------------------------------------
parametres 	 :une table de hashage locale.
-------------------------------------------------------
return 		 :rien.
------------------------------------------------------------------------------*/

void 
videTabLoc (pTab)
     pTableLoc pTab;
{
  pVertexLoc pLoc;
  int i;

  pLoc = pTab->pLoc;


  for (i = 0; i < pTab->lenTabLoc; i++)
    {
      pLoc[i].oper = EMPTYTH;
    }
}

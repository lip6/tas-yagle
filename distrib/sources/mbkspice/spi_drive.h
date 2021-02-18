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

#ifndef SPI_DRIVER
#define SPI_DRIVER

#define SPI_DRIVER_PTYPE 6809
#define SPI_SIG_INODE    6810
#define SPI_LOCON_INODE  6811
#define SPI_EXTLOCON     6812
#define SPI_CHOOSEN_LOCON 6813


typedef struct sconvindex
{
  losig_list    *sig;           /* Pointeur sur le signal */
  int            premier;       /* Index du premier noeud sur le signal */
} convindex;

/* NOTE : 
  Le champs premier contient le numerot du premier noeud dans le fichier
  Spice. Dans les vues RCN, le premier num�rot de noeud est le 1. Le num�rot
  de noeud Spice est donc donn� par la relation : 
      noeud_spice = noeud_rcn + premier - 1
*/
void   spi_vect             __P(( char* ));
void   cherche_alim         __P(( lofig_list *ptfig, char **vdd, char **vss ));
void   sortrcn              __P(( lofig_list *ptfig, FILE *df, char *vss ));
void   signalnoeud          __P(( lofig_list *ptfig ));
void   sortconnecteur       __P(( FILE *df, locon_list *c ));
void   sortnet              __P(( lofig_list *ptfig, FILE *df ));
void   sortinstance         __P(( lofig_list *ptfig, FILE *df ));
void   sorttransistormos    __P(( lofig_list *ptfig,
                                  FILE *df,
                                  char *vss,
                                  char *vdd
                               ));
void   sortcircuit          __P(( lofig_list *ptfig, FILE *df ));
void   tooutput             __P((FILE *,...)); /* va_list. */
void   spi_print            __P((FILE *, ...));
void   spi_env              __P(());
void   sortconnecteur_ordre __P(( FILE *df, chain_list*, locon_list* ));
void (*spi_getfuncinode(void))( FILE*, lofig_list*, void* );
void *spi_getdatainode(void);
void spi_cleanextlocon      __P(( lofig_list *lofig ));
locon_list* spi_getextlocon __P(( losig_list *losig ));

/*** Fonctions utilisateurs pour configurer le driver ***/
extern void spi_setfuncinode( void (*fn)( FILE*, lofig_list*, void* ), void *data );
extern num_list* spi_getinode( locon_list *locon );
extern void spi_setinode( locon_list *locon, num_list *head );
extern void spi_clearinode( locon_list *locon );
locon_list* spichooseonelocon( losig_list *losig );

#endif

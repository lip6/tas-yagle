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
 
/* ###--------------------------------------------------------------### */
/*									*/
/* file		: mvl_parse.c						*/
/* date		: Feb 15 1995						*/
/* author	: L.A TABUSSE & H.N. VUONG & P. BAZARGAN-SABET		*/
/* description	: Parser VHDL --> MBK					*/
/*									*/
/* ###--------------------------------------------------------------### */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MACOS
#include <sys/types.h>
#endif
#include AVT_H
#include MUT_H
#include MLO_H
#include "mvl_parse.h"
#include "mvl_utdef.h"

char          MVL_MBKMOD;               /* The mode of getlofig         */
char          MVL_ERRFLG = 0;           /* if = 1 no structure is made  */
int           MVL_LINNUM = 1;           /* file's line number           */
char          MVL_CURFIL[200];          /* current file's name          */

struct dct_entry **MVL_HSHTAB;          /* dictionnary's entry points   */

extern void mvl_scomprestart (FILE*);

static int   MVL_PRELOAD = 0;
void
parsevhdlnetlist(libname)
    char           *libname;
{
    MVL_PRELOAD = 1;
    vhdlloadlofig(NULL, libname, 'A');
    MVL_PRELOAD = 0;
}

void vhdlloadlofig (pt_lofig, figname, mode)

struct lofig *pt_lofig;
char         *figname ;
char          mode    ;

  {
  struct lofig *pt_lofig_tmp;
  struct loins *pt_loins   ;
  struct loins *pt_loinsaux;
  struct locon *pt_locon   ;
  struct locon *pt_loconaux;
  struct losig *pt_losig   ;
  struct losig *pt_losigaux;
  char          filename[120];
  char          local_mbk_fast_mode;
  static int    call_nbr = 0;
  char         *suffix = NULL;
  char         *str;

  if (!MVL_PRELOAD) {
    suffix = V_STR_TAB[__MVL_FILE_SUFFIX].VALUE;
    if (suffix == NULL) suffix = IN_LO;
    sprintf(filename,"%s.%s",figname,suffix);
  }
  else sprintf(filename,"%s",figname);

  /* Initialization of some variables */
  MVL_LINNUM = 1;
  MVL_MBKMOD = mode;
  strcpy (MVL_CURFIL, figname);

  /* FAST_MODE asked for MBK */
  local_mbk_fast_mode = FAST_MODE;
  FAST_MODE = 'Y';

  if ((mode != 'A') && (mode != 'P') && (mode != 'C'))
    {
    printf("vhdlloadfig : Bad mode '%c' asked\n", mode);
    FAST_MODE = local_mbk_fast_mode;
    EXIT(1);
    }

  /* Opening file */
  mvl_scompin = (FILE *) mbkfopen(figname, suffix, READ_TEXT);

  if(mvl_scompin == NULL)
    {
    (void)fprintf(stderr,"\n*** mbk error *** can't open file : %s\n",
						filename);
    FAST_MODE = local_mbk_fast_mode;
    EXIT(1);
    }

  /* TRACE_MODE asked for MBK */
  if(TRACE_MODE == 'Y')
    {
    (void)printf("\n--- mbk --- parsing file : %s in mode : %c\n",
                filename, mode);
    }

  MVL_LOFPNT = pt_lofig; /* passing main parameter */
  MVL_TOPPNT = pt_lofig; /* passing main parameter */

  if (call_nbr != 0)
    mvl_scomprestart (mvl_scompin);

  call_nbr ++;

  /* -------------------------------------------------------------------*/
  /* Parsing  : If mode is P or A, then normal parsing, if mode is C	*/
  /* then parsing of a new figure, then from the new one, we fill the   */
  /* old one.								*/
  /* -------------------------------------------------------------------*/

  if((mode == 'P') || (mode == 'A'))
    {
    if(mvl_scompparse() != 0)
      {
      (void)fprintf(stderr,"\n*** mbk error *** abnormal parsing for : %s\n",filename);
      FAST_MODE = local_mbk_fast_mode;
      EXIT(1);
      }
    }

  if(mode == 'C')
    {
    /* ----------------------------------------------------------------	*/
    /* Saving the lofig pointer, creating a new one to allow the 	*/
    /* parsing of the figure in 'A' mode.				*/
    /* ----------------------------------------------------------------	*/
    pt_lofig_tmp            = pt_lofig; 
    MVL_LOFPNT              = (lofig_list *)mbkalloc(sizeof(lofig_list));
    MVL_LOFPNT->MODE        = 'A';
    MVL_LOFPNT->NAME        = namealloc(figname);
    MVL_LOFPNT->MODELCHAIN  = NULL;
    MVL_LOFPNT->LOINS       = NULL;
    MVL_LOFPNT->LOTRS       = NULL;
    MVL_LOFPNT->LOCON       = NULL;
    MVL_LOFPNT->LOSIG       = NULL;
    mbk_init_NewBKSIG(&MVL_LOFPNT->BKSIG);
    MVL_LOFPNT->USER        = NULL;
    MVL_LOFPNT->NEXT        = NULL;

    MVL_MBKMOD = 'A';

    if(mvl_scompparse() != 0)
      {
      (void)fprintf(stderr,"\n*** mbk error *** abnormal parsing for : %s\n",filename);
      FAST_MODE = local_mbk_fast_mode;
      EXIT(1);
      }
    /* ----------------------------------------------------------------	*/
    /* Now, with the new figure, we duplicate the new informations	*/
    /* to fill the old one.						*/
    /* ----------------------------------------------------------------	*/
    pt_lofig = mvl_fill(pt_lofig_tmp, MVL_LOFPNT);
    }

  MVL_MBKMOD = mode;
  /* Closing file */
  if(fclose(mvl_scompin) != 0)
    {
    (void)fprintf(stderr,"\n*** mbk error *** can't close file : %s\n",filename);
    FAST_MODE = local_mbk_fast_mode;
    EXIT(1);
    }

  if(strcmp(IN_LO,"vbe") == 0)
    {
    strcpy(IN_LO,"vst");
    return;
    }

 if (mode == 'P' )
   {
   pt_locon = MVL_LOFPNT->LOCON;
   while (pt_locon != NULL)
     {
     if (pt_locon->TYPE == 'I')
       {
       pt_loconaux = pt_locon;
       pt_locon = pt_locon->NEXT;
       dellocon(MVL_LOFPNT, pt_loconaux->NAME);
       }
     else
       {
       /* pt_locon->SIG = NULL; */
       pt_locon = pt_locon->NEXT;
       }
     }
   pt_losig = MVL_LOFPNT->LOSIG;
   while (pt_losig != NULL)  
     {
     if (pt_losig->TYPE == 'E') 
       {
       pt_losig = pt_losig->NEXT;
       continue;
       }
     pt_losigaux = pt_losig;
     pt_losig = pt_losig->NEXT;
     dellosig(MVL_LOFPNT, pt_losigaux->INDEX);
     }
   pt_loins = MVL_LOFPNT->LOINS;
   while (pt_loins != NULL)
     {
     pt_loinsaux = pt_loins;
     pt_loins = pt_loins->NEXT;
     delloins(MVL_LOFPNT, pt_loinsaux->INSNAME);
     }
   }
 FAST_MODE = local_mbk_fast_mode;
 }

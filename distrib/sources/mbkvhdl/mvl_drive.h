/* 
 * This file is part of the Alliance CAD System
 * Copyright (C) Laboratoire LIP6 - Département ASIM
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
/* file		: mvl_drive.h						*/
/* date		: Nov 12 91						*/
/* author	: L.A. TABUSSE & P. BAZARGAN-SABET & VUONG H.N.		*/
/*									*/
/* content	: declaration of functions and global variables used by	*/
/*		  the structural MBK-> VHDL driver			*/
/*									*/
/* ###--------------------------------------------------------------### */


extern void  mvl_addtab ( struct dct_entry **head, char *key_str, char *ctx_str, int field, long valu );
extern long  mvl_chktab ( struct dct_entry **head, char *key_str, char *ctx_str, int field );
extern void  mvl_toolbug ( int code, char *str1, char *str2, int nbr1 );
extern void  mvl_message ( int code, char *str1, int nmb1 );
extern void  mvl_error ( int code, char *str1 );
extern char *mvl_vhdlname ( char* name );
extern void  mvl_name ( char *name, char *new_name  );
extern char *mvl_abl2str ();
extern void *mvl_vectnam( void *pt_list, int *left, int *right, char **name, char type );
extern char     mvl_chkconsistency();

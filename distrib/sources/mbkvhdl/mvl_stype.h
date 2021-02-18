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
/* file		: mvl_stype.h						*/
/* date		: Oct 15 1991						*/
/* author	: P. BAZARGAN, L.A. TABUSSE, VUONG H.N.			*/
/*									*/
/* contents	: This file contains defines and structure definitions	*/
/*		  for the structural compiler				*/
/*									*/
/* ###--------------------------------------------------------------### */

typedef struct
  {
  char            *NAME;          /* identifier name */
  short            LEFT;          /* vector's left index */
  short            RIGHT;         /* vector's right index */
  short            ERR_FLG;
  }
mvl_name;

typedef struct
  {
  short         WIDTH;			/* expression's width		*/
  struct chain *LIST;			/* list of losig pointers	*/
  }
mvl_expr;

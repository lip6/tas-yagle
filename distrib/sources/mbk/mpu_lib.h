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

/* 
 * Purpose : physical utilities declaration
 * Date    : 05/08/93
 * Author  : Frederic Petrot <Frederic.Petrot@lip6.fr>
 * Modified by Czo <Olivier.Sirol@lip6.fr> 1997,98
 * $Id: mpu_lib.h,v 1.1 2002/10/15 14:27:25 gregoire Exp $
 */

#ifndef _MPU_H_
#define _MPU_H_

#ifndef __P
# if defined(__STDC__) ||  defined(__GNUC__)
#  define __P(x) x
# else
#  define __P(x) ()
# endif
#endif

  extern           void  flattenphfig __P((phfig_list *ptfig, char *insname, char concat));
  extern           void  rflattenphfig __P((phfig_list *ptfig, char concat, char catal));
  extern           char  instanceface __P((char face, char sym));
  extern           void  bigvia __P((phfig_list *f, char via, long x, long y, long dx, long dy));
  extern           void  loadphfig __P((phfig_list *ptfig, char *figname, char mode));
  extern           void  savephfig __P((phfig_list *ptfig));
  extern           void  mphdebug __P((void *head_pnt, char *stru_name));
  
#endif /* !MPU */


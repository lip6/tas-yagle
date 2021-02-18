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




#ident "$Id: log_prefbib.c,v 1.57 2009/04/07 10:24:58 fabrice Exp $"

/****************************************************************************/
/*    Produit :  librairie ABL - Gestion de representations prefixees       */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include  MUT_H
#include LOG_H
#include AVT_H

#undef MIN_OPER
#undef MAX_OPER
#define MIN_OPER 0
#define MAX_OPER 1000

/*-------------------------------------------------------------------------
gensym_abl         : calcule un identificateur avec concatenation d'un numero.
---------------------------------------------------------------------------
retour                : un pointeur de chaine de caracteres.
---------------------------------------------------------------------------*/
char *
gensym_abl (name, num)
     char *name;
     int num;
{
  char *name1;
  char number[3];
  char *ret;

  name1 = (char *) mbkalloc (strlen (name) + 5);
  (void *) strcpy (name1, name);
  sprintf (number, "%d", num);
  strcat (name1, number);
  ret = namealloc (name1);
  mbkfree (name1);

  return ret;
}

/*================ Gestion des erreurs ====================*/

/*------------------------------------------------------------------
  IDENT : ablError
  ------------------------------------------------------------------
  FUNCTION : verifie que expr est bien une expression. func est le nom
             la fonction appelante.
  ------------------------------------------------------------------ */
void 
ablError (expr, func)
     chain_list *expr;
     char *func;
{
  int oper = OPER (expr);

  if (oper < MIN_OPER || oper > MAX_OPER)
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"162");
//      printf ("\nablError :\n%s: error - bad operator\n", func);
//      EXIT (-1);
    }
  func = NULL;
}

/*================ Creation et destruction d'ABL ====================*/

/*------------------------------------------------------------------
  IDENT : createAtom
  ------------------------------------------------------------------
  FUNCTION : cree un atome (un doublet) avec son CAR qui pointe sur
             name et son CDR egal a NULL. 
  ------------------------------------------------------------------ */
chain_list *
createAtom_no_NA (name)
     char *name;
{
  return addchain (NULL, name);
}

chain_list *
createAtom (name)
     char *name;
{
  return createAtom_no_NA (namealloc (name));
}

/*------------------------------------------------------------------
  IDENT : createExpr
  ------------------------------------------------------------------
  FUNCTION : cree une tete d'expression ayant oper comme operateur. 
  ------------------------------------------------------------------ */
chain_list *
createExpr (oper)
     short oper;
{
  /* POUR GNU : le cast (void *)oper a ete enleve' */
  if (oper >= MIN_OPER && oper <= MAX_OPER)
    return addchain (NULL, (void *) addchain (NULL, (void*)((long)oper)));
  else
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"163");
//      printf ("\ncreateExpr: error - unknown operator %d\n", oper);
//      EXIT (-1);
    }
#ifndef __ALL__WARNING__
  return 0; // ..anto..
#endif
}

/*------------------------------------------------------------------
  IDENT : notExpr
  ------------------------------------------------------------------
  FUNCTION : ajoute un NOT sur une expression.
  ------------------------------------------------------------------ */
chain_list *
notExpr (expr)
     chain_list *expr;
{
  chain_list *expr1;

  if (!ATOM (expr) && OPER (expr) == NOT)        /* simplification */
    {
      expr1 = CAR (CDR (expr));
      delchain (CAR (expr), CAR (expr));        /* desalloue l'operateur */
      freechain (expr);
      return expr1;
    }
  expr1 = createExpr (NOT);
  addHExpr (expr1, expr);
  return expr1;
}

/*------------------------------------------------------------------
  IDENT : createBinExpr
  ------------------------------------------------------------------
  FUNCTION : creer une expression a partir de deux expressions 
             avec fusion possible des operateurs.
  ------------------------------------------------------------------ */
chain_list *
createBinExpr (oper, expr1, expr2)
     short oper;
     chain_list *expr1, *expr2;
{
  chain_list *expr, *exprp;

  if (oper < MIN_OPER || oper > MAX_OPER || oper == NOT)
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"164");
//      printf ("\ncreateBinExpr: error - bad operator %d\n", oper);
//      EXIT (-1);
    }

  if (ATOM (expr1) || (OPER (expr1) != oper) ||
      ((oper == NOR) || (oper == NAND) || (oper == NXOR)))
    {
      /* non-compactage sur l'operateur */
      expr = createExpr (oper);
      addHExpr (expr, expr1);
    }
  else
    {
      expr = expr1;
    }

  if (ATOM (expr2) || (OPER (expr2) != oper) ||
      ((oper == NOR) || (oper == NAND) || (oper == NXOR)))
    {
      /* non-compactage sur l'operateur */
      addHExpr (expr, expr2);
    }
  else
    {
      /* Fusion de expr2 dans expr */
      exprp = expr2;
      while (CDR (exprp) != NULL)
        {
          exprp = CDR (exprp);
        }
      exprp->NEXT = expr->NEXT;
      expr->NEXT = expr2->NEXT;
      expr2->NEXT = NULL;        /* desa. de la tete de expr2 */
      freechain(expr2->DATA); freechain(expr2);
//      freeExpr (expr2);
    }
  return (expr);
}

/*------------------------------------------------------------------
  IDENT : addQExpr
  ------------------------------------------------------------------
  FUNCTION : ajoute a la fin d'une expression un nouvel argument. 
  ------------------------------------------------------------------ */
void 
addQExpr (expr1, expr2)
     chain_list *expr1, *expr2;
{
  ablError (expr1, "addQExpr");

  while (CDR (expr1))
    expr1 = CDR (expr1);
  expr1->NEXT = addchain (NULL, (void *) expr2);
}

/*------------------------------------------------------------------
  IDENT : addHExpr
  ------------------------------------------------------------------
  FUNCTION : ajoute au debut d'une expression un nouvel argument. 
  ------------------------------------------------------------------ */
void 
addHExpr (expr1, expr2)
     chain_list *expr1, *expr2;
{
  ablError (expr1, "addHExpr");

  expr1->NEXT = addchain (expr1->NEXT, (void *) expr2);
}

/*------------------------------------------------------------------
  IDENT : freeExpr
  ------------------------------------------------------------------
  FUNCTION : Desalloue l'expression jusqu'aux atomes. 
  ------------------------------------------------------------------ */
void 
freeExpr (expr)
     chain_list *expr;
{
  chain_list *pt;
 
 if (expr)
  {
      if (ATOM (expr))
      {
          freechain (expr);
      }
      else
      {
          pt = expr;
          delchain (CAR (expr), CAR (expr));        /* desalloue l'operateur */
          while ((expr = CDR (expr)))        /* puis les ss-expr...   */
              freeExpr (CAR (expr));
          freechain (pt);
      }
  }
}

/*============================ Affichage ============================*/

/*------------------------------------------------------------------
  IDENT : operToChar
  ------------------------------------------------------------------
  FUNCTION : renvoie la chaine de car. correspondant a l'operateur
             oper.
  ------------------------------------------------------------------ */
char *
operToChar (oper)
     short oper;
{
  if (oper >= MIN_OPER && oper <= MAX_OPER)
    {
      switch (oper)
        {
        case NOT:
          return (namealloc ("not"));
        case AND:
          return (namealloc ("and"));
        case OR:
          return (namealloc ("or"));
        case NAND:
          return (namealloc ("nand"));
        case NOR:
          return (namealloc ("nor"));
        case XOR:
          return (namealloc ("xor"));
        case NXOR:
          return (namealloc ("nxor"));
        case STABLE:
          return (namealloc ("stable"));
    case CAT:
      return (namealloc ("&"));
    case REPLICATE:
      return (namealloc ("rep"));
        default:
          return (namealloc ("???"));
        }
    }
  else
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"165");
//      printf ("\nablError - operToChar: unknown operator %d\n", oper);
//      EXIT (-1);
    }
#ifndef __ALL__WARNING__
  return 0; // ..anto..
#endif
}

/*------------------------------------------------------------------
  IDENT : charToOper
  ------------------------------------------------------------------
  FUNCTION : renvoie le numero d'oper. correspondant a l'ident name.
  ------------------------------------------------------------------ */
short 
charToOper (name)
     char *name;
{
  if (name == namealloc ("not"))
    return (NOT);
  if (name == namealloc ("and"))
    return (AND);
  if (name == namealloc ("or"))
    return (OR);
  if (name == namealloc ("xor"))
    return (XOR);
  if (name == namealloc ("nand"))
    return (NAND);
  if (name == namealloc ("nor"))
    return (NOR);
  if (name == namealloc ("nxor"))
    return (NXOR);
  if (name == namealloc ("stable"))
    return (STABLE);
  return (-1);
}

/*------------------------------------------------------------------
  IDENT : displayExprInt
  ------------------------------------------------------------------
  FUNCTION : affiche l'expression a la sauce LISP --> (oper list_arg) 
  ------------------------------------------------------------------ */
static int LOG_LOGLIB;
static int LOG_LOGLEVEL;

static void 
displayExprIntLog(chain_list *expr)
{
  int oper;

  if (ATOM (expr)) 
    avt_log (LOG_LOGLIB,LOG_LOGLEVEL," %s", VALUE_ATOM (expr));
  else
  {
    oper = OPER (expr);
    if (oper >= MIN_OPER && oper <= MAX_OPER)
    {
      avt_log (LOG_LOGLIB,LOG_LOGLEVEL," (%s", operToChar (OPER (expr)));
      if (oper==REPLICATE) {
        avt_log(LOG_LOGLIB,LOG_LOGLEVEL," %ld",REP_NB(expr));
      }
      mapExpr(displayExprIntLog, expr);
        avt_log(LOG_LOGLIB,LOG_LOGLEVEL,")");
    }
    else
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"166");
    }
  }
}

static void 
displayExprInt(chain_list *expr)
{
  int oper;

  if (ATOM (expr)) 
    printf(" %s", VALUE_ATOM (expr));
  else
  {
    oper = OPER (expr);
    if (oper >= MIN_OPER && oper <= MAX_OPER)
    {
      printf(" (%s", operToChar (OPER (expr)));
      if (oper==REPLICATE) {
        printf(" %ld",REP_NB(expr));
      }
      mapExpr(displayExprInt, expr);
        printf(")");
    }
    else
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"166");
    }
  }
}

/*------------------------------------------------------------------
  IDENT : displayExpr
  ------------------------------------------------------------------
  FUNCTION : affiche l'expression a la sauce LISP --> (oper list_arg) 
             + saut a la ligne
  ------------------------------------------------------------------ */

void 
displayExprLog (int lib, int loglevel, chain_list *expr)
{
  LOG_LOGLIB = lib;
  LOG_LOGLEVEL = loglevel;
  displayExprIntLog(expr);
  avt_log(lib, loglevel, "\n");
}

void 
displayExpr (chain_list *expr)
{
  displayExprInt(expr);
  fflush(stdout);
}

/*------------------------------------------------------------------
  IDENT : displayInfExpr
  ------------------------------------------------------------------
  FUNCTION : affiche l'ABL "expr" en infixe
  ------------------------------------------------------------------ */
void 
displayInfExprLog (int lib, int loglevel, chain_list *expr)
{
  int oper;

  if (ATOM (expr))		/* Traitement atomique */
    avt_log(lib,loglevel,"%s", VALUE_ATOM (expr));
  else
    /* Traitement d'expression */
    {
      oper = OPER (expr);	/* Sauvegarde de l'operateur */
      if (oper == NOT)		/* Traitement du NOT */
        {
          avt_log(lib,loglevel,"not ");
          displayInfExprLog(lib, loglevel, CAR (CDR (expr)));
        }
      else
        /* Traitement des autres operateurs */
        {
          avt_log(lib,loglevel,"(");
          while ((expr = CDR (expr)))
            {
              displayInfExprLog (lib, loglevel, CAR (expr));
              if (CDR (expr))	/* Un operateur en moins que  */
                avt_log(lib,loglevel," %s ", operToChar (oper));	/* le nombre d'arguments      */
            }
          avt_log(lib,loglevel,")");
        }
    }
}

void 
displayInfExpr (chain_list *expr)
{
  int oper;

  if (ATOM (expr))		/* Traitement atomique */
    printf("%s", VALUE_ATOM (expr));
  else
    /* Traitement d'expression */
    {
      oper = OPER (expr);	/* Sauvegarde de l'operateur */
      if (oper == NOT)		/* Traitement du NOT */
        {
          printf("not ");
          displayInfExpr (CAR (CDR (expr)));
        }
      else
        /* Traitement des autres operateurs */
        {
          printf("(");
          while ((expr = CDR (expr)))
            {
              displayInfExpr (CAR (expr));
              if (CDR (expr))	/* Un operateur en moins que  */
                printf(" %s ", operToChar (oper));	/* le nombre d'arguments      */
            }
          printf(")");
        }
    }
    fflush(stdout);
}

/*------------------------------------------------------------------
  IDENT : exprToCharInt
  ------------------------------------------------------------------
  FUNCTION : renvoie une chaine de caracteres correspondant a la fonction
             en mode infixe (mode = 1) ou en mode prefixe (mode = 0)
  ------------------------------------------------------------------ */
char *
exprToCharInt (expr, mode, chaine, taille)
     chain_list *expr;
     int mode;
     char *chaine;
     int *taille;
{
  char *oper;
  int lenVA;

  if (ATOM (expr))
    {
      if ( (unsigned int)*taille < strlen (chaine) + 50)   /* re-allocation */
        {
          *taille = *taille + 100;
          chaine = (char *) mbkrealloc (chaine, *taille);
        }
      lenVA = strlen (VALUE_ATOM (expr));

      strncat (chaine, VALUE_ATOM (expr), lenVA);
      return (chaine);
    }
  else
    {
      oper = operToChar (OPER (expr));

      if (mode)
/*----  infixe ----*/
        {
          if (OPER (expr) == NOT)
            {
              strncat (chaine, "not ", 4);
              chaine = exprToCharInt (CAR (CDR (expr)), mode, chaine, taille);
              return (chaine);
            }
          else if (OPER (expr) == STABLE)
            {
              chaine = exprToCharInt (CAR (CDR (expr)), mode, chaine, taille);
              strncat (chaine, "'stable ", 8);
              return (chaine);
            }
          else
            {
              strncat (chaine, "(", 1);
              while ((expr = CDR (expr)))
                {
                  chaine = exprToCharInt (CAR (expr), mode, chaine, taille);
                  if (CDR (expr))
                    {
                      strncat (chaine, " ", 1);
                      strcat (chaine, oper);
                      strncat (chaine, " ", 1);
                    }
                }
              strncat (chaine, ")", 1);
              return (chaine);
            }
        }
      else
/*----- prefixe -----*/
        {
          int lenOPER = strlen (oper);
          strncat (chaine, "(", 1);
          strncat (chaine, oper, lenOPER);
          while ((expr = CDR (expr)))
            {
              strncat (chaine, " ", 1);
              chaine = exprToCharInt (CAR (expr), mode, chaine, taille);
            }
          strncat (chaine, ")", 1);
          return (chaine);
        }
    }
}

/*------------------------------------------------------------------
  IDENT : exprToChar
  ------------------------------------------------------------------
  FUNCTION : renvoie une chaine de caracteres correspondant a la fonction
             en mode infixe (mode = 1) ou en mode prefixe (mode = 0)
  ------------------------------------------------------------------ */
char *
exprToChar (expr, mode)
     chain_list *expr;
     int mode;
{
  char *chaine;
  int taille;

  taille = 100;
  chaine = (char *) mbkalloc (100);
  chaine[0] = '\0';

  return exprToCharInt (expr, mode, chaine, &taille);
}

/*------------------------------------------------------------------
  IDENT : identExprInt
  ------------------------------------------------------------------
  FUNCTION : renvoie une chaine de caracteres correspondant a un nom
             d'expression aux variables pres. 
  ------------------------------------------------------------------ */
char *
identExprInt (expr, chaine, taille)
     chain_list *expr;
     char *chaine;
     int *taille;
{
  if (ATOM (expr))
    return chaine;
  else
    {
      chain_list *expr1;
      char arite[3];
      int oper = OPER (expr);

      switch (oper)
        {
        case NOT:
          strncat (chaine, "n", 1);
          break;
        case AND:
          strncat (chaine, "a", 1);
          break;
        case OR:
          strncat (chaine, "o", 1);
          break;
        case NAND:
          strncat (chaine, "na", 2);
          break;
        case NOR:
          strncat (chaine, "no", 2);
          break;
        case XOR:
          strncat (chaine, "x", 1);
          break;
        case NXOR:
          strncat (chaine, "nx", 2);
          break;
        default:
          strncat (chaine, "???", 3);
        }
      expr1 = expr;
      while ((expr = CDR (expr)))
        {
          chaine = identExprInt (CAR (expr), chaine, taille);
        }
      if (oper != NOT)
        {
          sprintf (arite, "%d", lengthExpr (expr1));
          strcat (chaine, arite);
        }
      return (chaine);
    }
}


/*------------------------------------------------------------------
  IDENT : identExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie une chaine de caracteres correspondant a un nom
             d'expression aux variables pres. 
  ------------------------------------------------------------------ */
char *
identExpr (expr)
     chain_list *expr;
{
  char *chaine, *bidon;
  int taille, i;
  chain_list *expr1;
  int flag_alloc = 0;

  taille = 100;

  /* normalisation */
  flatArityExpr (expr);
  if (searchOperExpr (expr, NOT))
    {
      if (!searchOperExpr (expr, XOR))
        {
          flag_alloc = 1;
          expr1 = flatPolarityExpr (expr, 1);
          if (searchOperExpr (expr1, NOT))
            {
              expr = notExpr (flatPolarityExpr (expr1, 0));
              freeExpr (expr1);
            }
          else
            {
              expr = expr1;
            }
        }
    }

  if (profExpr (expr) == 1)
    {
      chaine = gensym_abl (operToChar (OPER (expr)), lengthExpr (expr));
    }
  else if (profExpr (expr) == 2 && OPER (expr) == NOT)
    {
      bidon = (char *) mbkalloc (taille);
      bidon[0] = 'n';
      bidon[1] = '\0';
      strcat (bidon, operToChar (OPER (CADR (expr))));
      chaine = gensym_abl (bidon, lengthExpr (CADR (expr)));
      mbkfree (bidon);
    }
  else
    {
      bidon = (char *) mbkalloc (taille);
      bidon[0] = '\0';
      bidon = identExprInt (expr, bidon, &taille);
      for (i = strlen (bidon) - 1; bidon[i] == '_'; i--)
        bidon[i] = '\0';
      chaine = namealloc (bidon);
      mbkfree (bidon);
    }
  if (flag_alloc)
    freeExpr (expr);

  return chaine;
}

/*------------------------------------------------------------------
  IDENT : profExpr
  ------------------------------------------------------------------
  FUNCTION : calcule la profondeur d'une expression.
  ------------------------------------------------------------------*/
int 
profExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return 0;
  else
    {
      int max = 0, profCar;

      while ((expr = CDR (expr)))
        {
          profCar = profExpr (CAR (expr));
          if (profCar > max)
            max = profCar;
        }
      return ++max;
    }
}

/*------------------------------------------------------------------
  IDENT : profAOExpr
  ------------------------------------------------------------------
  FUNCTION : calcule la profondeur d'une expression sans tenir compte
  des inverseurs.
  ------------------------------------------------------------------*/
int 
profAOExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return 0;
  else
    {
      int max = 0, profCar, oper = OPER (expr);

      while ((expr = CDR (expr)))
        {
          profCar = profAOExpr (CAR (expr));
          if (profCar > max)
            max = profCar;
        }
      if (oper == NOT)
        return max;

      return ++max;
    }
}

/*======================= Les fonctionnelles ========================*/

/*------------------------------------------------------------------
  IDENT : mapCarExpr
  ------------------------------------------------------------------
  FUNCTION : applique la fonction func sur l'expression expr et
             renvoie une expression ayant oper comme operateur et 
             ayant pour ss-expr. le resultat de l'application de 
             func sur les ss-expr. de expr.
  ------------------------------------------------------------------*/
chain_list *
mapCarExpr (func, oper, expr)
     chain_list *(*func) ();
     short oper;
     chain_list *expr;
{
  chain_list *pt;

  pt = createExpr (oper);
  while ((expr = CDR (expr)))
    addQExpr (pt, (*func) (CAR (expr)));
  return pt;
}

/*------------------------------------------------------------------
  IDENT : mapExpr
  ------------------------------------------------------------------
  FUNCTION : applique la fonction func sur l'expression expr.
             func doit renvoye un void et prendre comme argument un 
             pointeur de chain_list.
  ------------------------------------------------------------------*/
void 
mapExpr (func, expr)
     void (*func) ();
     chain_list *expr;
{
  while ((expr = CDR (expr)))
    (*func) (CAR (expr));
}

/*------------------------------------------------------------------
  IDENT : anyExpr
  ------------------------------------------------------------------
  FUNCTION : applique la fonction func sur l'expression expr.
             renvoie 1 s'il y a eu au moins un appel de func qui a
             renvoye 1. (OU  logique)
  ------------------------------------------------------------------*/
int 
anyExpr (func, expr)
     int (*func) ();
     chain_list *expr;
{
  while ((expr = CDR (expr)))
    if ((*func) (CAR (expr)))
      return 1;

  return 0;
}

/*------------------------------------------------------------------
  IDENT : EveryExpr
  ------------------------------------------------------------------
  FUNCTION : applique la fonction func sur l'expression expr.
             renvoie 1 si tout les appels de func ont renvoye 1.
            (ET logique)
  ------------------------------------------------------------------*/
int 
everyExpr (func, expr)
     int (*func) ();
     chain_list *expr;
{
  while ((expr = CDR (expr)))
    if (!((*func) (CAR (expr))))
      return 0;

  return 1;
}

/*====================== Manipulations d'ABL =========================*/

/*------------------------------------------------------------------
  IDENT : searchOperExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie 1 si oper existe dans expr, 0 sinon.
  ------------------------------------------------------------------*/
int 
searchOperExpr (expr, oper)
     chain_list *expr;
     short oper;
{
  if (ATOM (expr))
    return 0;
  else if (OPER (expr) == oper)
    return 1;
  else
    {
      while ((expr = CDR (expr)))
        if (searchOperExpr (CAR (expr), oper) == 1)
          return 1;
      return 0;
    }
}

/*------------------------------------------------------------------
  IDENT : searchExprLow
  ------------------------------------------------------------------
  FUNCTION : renvoie 1 si name existe dans expr, 0 sinon.
  ------------------------------------------------------------------*/
short 
searchExprLow (expr, name)
     chain_list *expr;
     char *name;
{
  if (ATOM (expr))
    if (name == VALUE_ATOM (expr))
      return 1;
    else
      return 0;
  else
    {
      while ((expr = CDR (expr)))
        if (searchExprLow (CAR (expr), name))
          return 1;
      return 0;
    }
}

/*------------------------------------------------------------------
  IDENT : searchExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie 1 si namealloc(name) existe dans expr, 0 sinon.
  ------------------------------------------------------------------*/
int 
searchExpr (expr, name)
     chain_list *expr;
     char *name;
{

  return searchExprLow (expr, namealloc (name));
}

/*------------------------------------------------------------------
  IDENT : equalExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie 1 si expr1 = expr2, 0 sinon.
  ------------------------------------------------------------------*/
int 
equalExpr (expr1, expr2)
     chain_list *expr1, *expr2;
{
  if (ATOM (expr1)) {
    if (ATOM (expr2) && (CAR (expr1) == CAR (expr2)))
      return 1;
    else
      return 0;
  }
  if (ATOM (expr2))
    return 0;
  else
    {
      if (OPER (expr1) != OPER (expr2))                /* egalite des operateurs */
        return 0;
      expr1 = CDR (expr1);
      expr2 = CDR (expr2);
      while (expr1 && expr2)
        {
          if (!(equalExpr (CAR (expr1), CAR (expr2))))
            return 0;
          expr1 = CDR (expr1);
          expr2 = CDR (expr2);
        }
      if (expr2 || expr1)
        return 0;
      else
        return 1;
    }
}

/*------------------------------------------------------------------
  IDENT : equalVarExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie 1 si expr1 = expr2 aux variables pres, 0 sinon.
  ------------------------------------------------------------------*/
int 
equalVarExpr (expr1, expr2)
     chain_list *expr1, *expr2;
{
  if (ATOM (expr1))
    {
      if (ATOM (expr2))
        return 1;
      else
        return 0;
    }
  if (ATOM (expr2))
    return 0;
  else
    {
      if (OPER (expr1) != OPER (expr2))                /* egalite des operateurs */
        return 0;
      expr1 = CDR (expr1);
      expr2 = CDR (expr2);
      while (expr1 && expr2)
        {
          if (!(equalVarExpr (CAR (expr1), CAR (expr2))))
            return 0;
          expr1 = CDR (expr1);
          expr2 = CDR (expr2);
        }
      if (expr2 || expr1)
        return 0;
      else
        return 1;
    }
}

/*------------------------------------------------------------------
  IDENT : lengthExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie la taille de l'expression au premier niveau.
             (en fait l'arite de l'operateur)
  ------------------------------------------------------------------*/
int 
lengthExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return 1;
  else
  {
    int cpt = 0;

    while ((expr = CDR (expr)))
      cpt++;
    return cpt;
  }
}

/*------------------------------------------------------------------
  IDENT : numberOperBinExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le nombre d'operateurs binaires (equivalent) 
             d'une expression.
  ------------------------------------------------------------------*/
int 
numberOperBinExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return 0;
  else
    {
      int cpt = -1;

      if (OPER (expr) == NOT)
        return numberOperBinExpr (CADR (expr));

      while ((expr = CDR (expr)))
        cpt = 1 + cpt + numberOperBinExpr (CAR (expr));

      return cpt;
    }
}

/*------------------------------------------------------------------
  IDENT : numberAtomExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le nombre d'atomes de l'expression.
  ------------------------------------------------------------------*/
int 
numberAtomExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return 1;
  else
    {
      int cpt = 0;

      while ((expr = CDR (expr)))
        cpt += numberAtomExpr (CAR (expr));

      return cpt;
    }
}

/*------------------------------------------------------------------
  IDENT : copyExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie une expression qui correspond a une recopie de
             expr (avec l'arite).
  ------------------------------------------------------------------*/
chain_list *
copyExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return addchain (NULL, expr->DATA);
  else
    {
      chain_list *auxExpr;
      chain_list *oldExpr = expr;
      chain_list *expr1;

      auxExpr = createExpr (OPER (expr));
      expr1 = auxExpr;

      while ((expr = CDR (expr)))
        {
          /* on optimise a fond la recopie de la dorsale */
          expr1->NEXT = addchain (NULL, (void *) copyExpr (CAR (expr)));
          expr1 = CDR (expr1);

          /* remplace : addQExpr(auxExpr,copyExpr(CAR(expr))); */
        }
      //
      // warning CDR(CAR(expr)) can contain a *pNode or *chain_list
      // to treat the two cases i make a test on the value pointed
      // by CDR(CAR(expr))
      //                                                anto

      if (OPER(oldExpr)==REPLICATE)
      {
        CDAR(auxExpr) = addchain(NULL,(void *)REP_NB(oldExpr));
        CDAR(auxExpr) = addchain(CDAR(auxExpr),(void *)VECT_S(oldExpr));
      }
      else if ((CDAR(oldExpr))&&(!CDR(CADR(oldExpr))))
        CDAR(auxExpr) = addchain(NULL,(void *)VECT_S(oldExpr));
      else  /* on sauvegarde la poubelle */
        ((chain_list *) auxExpr->DATA)->NEXT = CDR (CAR (oldExpr));
      
      return (auxExpr);
    }
}
/*------------------------------------------------------------------
  IDENT : substPhyExpr
  ------------------------------------------------------------------
  FUNCTION : substitue toutes les occurences de "name" dans expr1 par
             l'expression expr2.
  ------------------------------------------------------------------*/
void 
substPhyExpr (expr1, name, expr2)
     chain_list *expr1;
     char *name;
     chain_list *expr2;
{
  if (!ATOM (expr1))
    {
      while ((expr1 = CDR (expr1)))
        if (ATOM (CAR (expr1)))
          {
            if (!strcmp (name, VALUE_ATOM (CAR (expr1))))
              {
                freechain (CAR (expr1));
                expr1->DATA = (void *) copyExpr (expr2);
              }
          }
        else
          substPhyExpr (CAR (expr1), name, expr2);
    }
}

/*------------------------------------------------------------------
  IDENT : substExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie une expression dans laquelle on a substitue 
             toutes les occurences de "name" par l'expression expr2.
  ------------------------------------------------------------------*/
chain_list *
substExpr (expr1, name, expr2)
     chain_list *expr1;
     char *name;
     chain_list *expr2;
{
  if (ATOM (expr1))
    if (!strcmp (name, VALUE_ATOM (expr1)))
      return copyExpr (expr2);
    else
      return copyExpr (expr1);
  else
    {
      chain_list *res = createExpr (OPER (expr1));

      while ((expr1 = CDR (expr1)))
        addQExpr (res, substExpr (CAR (expr1), name, expr2));

      return res;
    }
}

/*------------------------------------------------------------------
  IDENT : devXorExpr
  ------------------------------------------------------------------
  FUNCTION : developpe recursivement les XOR dans une expression. 
             a xor b = ((not a) and b) or (a and (not b))
  ------------------------------------------------------------------*/
chain_list *
devXorExpr (expr)
     chain_list *expr;
{
  chain_list    *tmp;
  if (ATOM (expr))
    return copyExpr (expr);
  else
    {
      if (OPER (expr) != XOR)
      {// modified by anto to compute vectabl
        tmp = mapCarExpr (devXorExpr, OPER (expr), expr);
        if (OPER (expr) == REPLICATE)
        {
          CDAR(tmp) = addchain(NULL     ,(void *)REP_NB(expr));
          CDAR(tmp) = addchain(CDAR(tmp),(void *)VECT_S(expr));
        }
        else if (OPER (expr) == CAT)
          CDAR(tmp) = addchain(NULL     ,(void *)VECT_S(expr));
        return tmp;
      }
      else if (lengthExpr (expr) == 2)
        {
          chain_list *dev, *ret;

          dev = devXor2Expr (expr);
          ret = devXorExpr (dev);
          freeExpr (dev);
          return (ret);
        }
      else
        {
          /* XOR N>2 */
          chain_list *xor, *dev, *inter;

          xor = createExpr (XOR);
          addHExpr (xor, copyExpr (CADR (expr)));

          /* elimination du premier arg dans expr */

          inter = createExpr (XOR);
          expr = CDR (expr);
          while ((expr = CDR (expr)))
            addHExpr (inter, copyExpr (CAR (expr)));

          dev = devXorExpr (inter);        /* on developpe expr */
          freeExpr (inter);

          addQExpr (xor, dev);
          dev = devXor2Expr (xor);
          freeExpr (xor);
          return (dev);
        }
    }
}

/*------------------------------------------------------------------
  IDENT : devXor2Expr
  ------------------------------------------------------------------
  FUNCTION : developpe un XOR 2
             a xor b = ((not a) and b) or (a and (not b))
  ------------------------------------------------------------------*/
chain_list *
devXor2Expr (expr)
     chain_list *expr;
{
  if (ATOM (expr) || OPER (expr) != XOR || lengthExpr (expr) > 2)
    return (copyExpr (expr));
  else
    {
      chain_list *ret;
      chain_list *and1, *and2;
      chain_list *a, *b;

      a = CADR (expr);                /* premier argument */
      b = CADR (CDR (expr));        /* deuxieme argument */
      ret = createExpr (OR);
      and1 = createExpr (AND);
      and2 = createExpr (AND);
      addQExpr (and1, copyExpr (a));
      addQExpr (and2, copyExpr (b));

      addQExpr (and1, notExpr (copyExpr (b)));
      addQExpr (and2, notExpr (copyExpr (a)));
      addQExpr (ret, and2);
      addQExpr (ret, and1);
      return ret;
    }
}


/*------------------------------------------------------------------
  IDENT : flatPolarityExpr
  ------------------------------------------------------------------
  FUNCTION : Reporte sur les entrees primaires (atomes) la negation
             par applications successives des lois de Morgan et de
             simplifications de l'algebre de Boole.
             Si le signe est egal a 0 a l'appel, on obtient un
             report sur la fonction inverse.
  ------------------------------------------------------------------ */
chain_list *
flatPolarityExpr (expr, signe)
     chain_list *expr;
     int signe;
{
  short oper;
  chain_list *pt;

  if (ATOM (expr))
    if (signe)
      return (createAtom (VALUE_ATOM (expr)));
    else
    {
      pt = createExpr (NOT);
      addQExpr (pt, createAtom (VALUE_ATOM (expr)));
      return pt;
    }
  else
  {
    oper = OPER (expr);
    if (oper == NOT)                /* operateur NOT */
    {
      expr = CADR (expr);
      if (signe)
        return (flatPolarityExpr (expr, 0));
      else
        return (flatPolarityExpr (expr, 1));
    }
    else
    {
      if (oper == OR || oper == AND)
      {
        if (signe)        /* aucun changement a ce niveau */
          pt = createExpr (oper);
        else
          /* application des lois de Morgan */ if (oper == AND)
            pt = createExpr (OR);
          else
            pt = createExpr (AND);
        while ((expr = CDR (expr)))
          addQExpr (pt, flatPolarityExpr (CAR (expr), signe));
        return (pt);
      } 
      else if (oper == NOR || oper == NAND)
      {
        pt = copyExpr (expr);
        if (OPER (expr) == NOR)
          changeOperExpr (pt, OR);
        else
          changeOperExpr (pt, AND);
        if (signe)
          expr = flatPolarityExpr (pt, 0);
        else
          expr = flatPolarityExpr (pt, 1);
        freeExpr (pt);
        return (expr);
      }
      // added by anto to compute vectabl
      else if (oper == CAT || oper == REPLICATE)
      {
        pt = createExpr(oper);
        switch (oper)
        {
          case REPLICATE : CDAR(pt) = addchain(NULL,(void *)REP_NB(expr));
          case CAT       : CDAR(pt) = addchain(CDAR(pt),(void *)VECT_S(expr));
        }
        while ((expr = CDR (expr)))
          addQExpr (pt, flatPolarityExpr (CAR (expr), signe));
        return (pt);
      }
      ///////////////////////////////////
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"167");
//      printf ("flatPolarityExpr: error - unknown operator %s\n", operToChar (oper));
//      EXIT (-1);
    }
  }
#ifndef __ALL__WARNING__
  return 0; // ..anto..
#endif
}

/*------------------------------------------------------------------
  IDENT : flatArityExpr
  ------------------------------------------------------------------
  FUNCTION : Ecrase une expression par elimination des operateurs redondants.  
  ------------------------------------------------------------------ */
void 
flatArityExpr (expr)
     chain_list *expr;
{
  int oper;
  chain_list *expr1, *expr2, *exprBefore;

  if (!ATOM (expr))
    {
      oper = OPER (expr);
      if (oper == NOT || oper == STABLE)
        flatArityExpr (CADR (expr));
      else if (oper == NOR || oper == NAND || oper == NXOR)
        {
          expr1 = expr;
          while ((expr1 = CDR (expr1)))
            flatArityExpr (CAR (expr1));
        }
      else
        {
          expr1 = expr;
          exprBefore = expr;
          while ((expr1 = CDR (expr1)))
            {
              expr2 = CAR (expr1);
              if (!ATOM (expr2))
                {
                  if (oper == OPER (expr2))        /* detection OK */
                    {
                      exprBefore->NEXT = CDR (expr2);
                      while (CDR (expr2))
                        expr2 = CDR (expr2);
                      CDR (expr2) = CDR (expr1);

                      CDR (expr1) = NULL;
                      expr2 = CAR (expr1);
                      CDR (expr2) = NULL;
                      freechain (expr1);
                      freechain (expr2);
                      expr1 = exprBefore;
                    }
                  else
                    {
                      flatArityExpr (expr2);
                      exprBefore = CDR (exprBefore);
                    }
                }
              else
                exprBefore = CDR (exprBefore);
            }
        }
    }
}

/*------------------------------------------------------------------
  IDENT : supportChain_listExprInt
  ------------------------------------------------------------------
  FUNCTION : renvoie le support d'une expression dans une liste chainee.
             (1 seule occurence par variable).
  ------------------------------------------------------------------ */
void 
supportChain_listExprInt (expr, ppCL)
     chain_list *expr;
     chain_list **ppCL;
{
  if (ATOM (expr))
    {
      chain_list *pCL = *ppCL;
      if (strcmp (VALUE_ATOM (expr), "'d'") &&
          strcmp (VALUE_ATOM (expr), "'u'") &&
          strcmp (VALUE_ATOM (expr), "'z'") &&
          strcmp (VALUE_ATOM (expr), "'1'") &&
          strcmp (VALUE_ATOM (expr), "'0'") &&
        // added by anto for bitstring
      !isBitStr(VALUE_ATOM (expr))
      )
        {
          while (pCL)                /* recherche occurence */
            {
              if ((char *) pCL->DATA == VALUE_ATOM (expr))
                return;
              pCL = pCL->NEXT;
            }
          *ppCL = addchain (*ppCL, (void *) VALUE_ATOM (expr));
        }
    }
  else
    {
      while ((expr = CDR (expr)))
        supportChain_listExprInt (CAR (expr), ppCL);
    }
}

/*------------------------------------------------------------------
  IDENT : supportChain_listExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le support d'une expression dans une liste chainee.
             (1 seule occurence par variable).
  ------------------------------------------------------------------ */
chain_list *
supportChain_listExpr (expr)
     chain_list *expr;
{
  chain_list *ppCL;
  chain_list *res;

//  ppCL = (chain_list **) mbkalloc (sizeof (chain_list *));
//  *
  ppCL = NULL;
  supportChain_listExprInt (expr, &ppCL);
  res = reverse (ppCL);
//  mbkfree (ppCL);

  return res;
}

/*------------------------------------------------------------------
  IDENT : fullSupportChain_listExprInt
  ------------------------------------------------------------------
  FUNCTION : renvoie le support d'une expression dans une liste chainee.
  ------------------------------------------------------------------ */
void 
fullSupportChain_listExprInt (expr, ppCL)
     chain_list *expr;
     chain_list **ppCL;
{
  if (ATOM (expr))
    *ppCL = addchain (*ppCL, (void *) VALUE_ATOM (expr));
  else
    while ((expr = CDR (expr)))
          fullSupportChain_listExprInt (CAR (expr), ppCL);
}

/*------------------------------------------------------------------
  IDENT : fullSupportChain_listExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le support d'une expression dans une liste chainee.
  ------------------------------------------------------------------ */
chain_list *
fullSupportChain_listExpr (expr)
     chain_list *expr;
{
  chain_list *ppCL;
  chain_list *res;

//  ppCL = (chain_list **) mbkalloc (sizeof (chain_list *));
//  *
  ppCL = NULL;
  fullSupportChain_listExprInt (expr, &ppCL);
  res = reverse (ppCL);
//  res = reverse (*ppCL);
//  mbkfree (ppCL);

  return res;
}

/*------------------------------------------------------------------
  IDENT : supportPtype_listExprInt
  ------------------------------------------------------------------
  FUNCTION : renvoie le support d'une expression dans une liste typee.
             (1 seule occurence par variable).
  ------------------------------------------------------------------ */
void 
supportPtype_listExprInt (expr, ppCL)
     chain_list *expr;
     ptype_list **ppCL;
{
  if (ATOM (expr))
    {
      ptype_list *pCL = *ppCL;
      if (strcmp (VALUE_ATOM (expr), "'d'") &&
          strcmp (VALUE_ATOM (expr), "'u'") &&
          strcmp (VALUE_ATOM (expr), "'z'") &&
          strcmp (VALUE_ATOM (expr), "'1'") &&
          strcmp (VALUE_ATOM (expr), "'0'"))
        {
          while (pCL)                /* recherche occurence */
            {
              if ((char *) pCL->DATA == VALUE_ATOM (expr))
                return;
              pCL = pCL->NEXT;
            }
          *ppCL = addptype (*ppCL, 0, (void *) VALUE_ATOM (expr));
        }
    }
  else
    {
      while ((expr = CDR (expr)))
        supportPtype_listExprInt (CAR (expr), ppCL);
    }
}

/*------------------------------------------------------------------
  IDENT : supportPtype_listExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le support d'une expression dans une liste typee.
             (1 seule occurence par variable).
  ------------------------------------------------------------------ */
ptype_list *
supportPtype_listExpr (expr)
     chain_list *expr;
{
  ptype_list **ppCL;
  ptype_list *res;

  ppCL = (ptype_list **) mbkalloc (sizeof (ptype_list *));
  *ppCL = NULL;
  supportPtype_listExprInt (expr, ppCL);
  res = (ptype_list *) reverse ((chain_list *) * ppCL);
  mbkfree (ppCL);
  return (res);
}

/*------------------------------------------------------------------
  IDENT : maxExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le max des arguments d'une expression avec une
  fonction de comparaison func renvoyant un int.
  ------------------------------------------------------------------ */
chain_list *
maxExpr (expr, func)
     chain_list *expr;
     int (*func) ();
{
  chain_list *expr1;

  expr = CDR (expr);
  expr1 = expr;
  while ((expr = CDR (expr)))
    if ((*func) (CAR (expr)) > (*func) (CAR (expr1)))
      expr1 = expr;

  return CAR (expr1);
}

/*------------------------------------------------------------------
  IDENT : minExpr
  ------------------------------------------------------------------
  FUNCTION : renvoie le min des arguments d'une expression avec une
  fonction de comparaison func renvoyant un int.
  ------------------------------------------------------------------ */

chain_list *
minExpr (expr, func)
     chain_list *expr;
     int (*func) ();
{
  chain_list *expr1;

  expr = CDR (expr);
  expr1 = expr;
  while ((expr = CDR (expr)))
    if ((*func) (CAR (expr)) < (*func) (CAR (expr1)))
      expr1 = expr;

  return CAR (expr1);
}

/*------------------------------------------------------------------
  IDENT : sortExpr
  ------------------------------------------------------------------
  FUNCTION : trie les arguments d'une expression a partir d'une fonction
  func renvoyant un int. 
            direction = 1 croissant  ;  direction = 0 decroissant
  ------------------------------------------------------------------ */
void 
sortExpr (expr, func, direction)
     chain_list *expr;
     long (*func) ();
     int direction;
{
  chain_list *expr1;
  int okPermu = 1;

  if (!ATOM (expr))                /* ce n'est pas un atome */
    {
      while (okPermu)
        {
          expr1 = CDR (expr);
          okPermu = 0;
          while (CDR (expr1))
            {
              if ((direction == 1 &&
                   ((*func) CAR (expr1) > (*func) CADR (expr1))) ||
                  (direction == 0 &&
                   ((*func) CAR (expr1) < (*func) CADR (expr1))))
                {
                  chain_list *permu;

                  permu = CAR (expr1);
                  expr1->DATA = (void *) CADR (expr1);
                  (expr1->NEXT)->DATA = (void *) permu;
                  okPermu = 1;
                }
              expr1 = CDR (expr1);
            }
        }
    }
}

/*------------------------------------------------------------------
  IDENT : funcNormExpr
  ------------------------------------------------------------------
  FUNCTION : normalise une expression a tous les niveaux
                variable = pointeur d'allocation
                expression = numero d'operateur. 
             ex : (and a (xor w g) b (or x f) (or e g o) (xor t h k))
                --> (and (or x f) (or e o g) (xor g w) (xor t h k) b a)
  ------------------------------------------------------------------ */
long 
funcNormExpr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
    return ((long) VALUE_ATOM (expr));
  else if (OPER (expr) == NOT) {
    if (ATOM (CADR (expr)))
      return (NOT * 100 + 10);
    else
      return (NOT * 100 + OPER (CADR (expr)));
  }
  return OPER (expr) * 100 + lengthExpr (expr);
}

/*------------------------------------------------------------------
  IDENT : normExpr
  ------------------------------------------------------------------
  FUNCTION : normalise une expression a tous les niveaux
                variable = pointeur d'allocation
                expression = numero d'operateur. 
             ex : (and a (xor w g) b (or x f) (or e g o) (xor t h k))
                --> (and (or x f) (or e o g) (xor g w) (xor t h k) b a)
  ------------------------------------------------------------------ */
void 
normExpr (expr)
     chain_list *expr;
{
  if (!ATOM (expr))
    {
      sortExpr (expr, funcNormExpr, 1);
      while ((expr = CDR (expr)))
        normExpr (CAR (expr));
    }
}

/*----------------------------------------------------------------------------
deleteNumExpr         : elimine le i-eme element de expr.
                  le premier argument commence a la place 0.
------------------------------------------------------
parametres          : un ABL et une place. 
-------------------------------------------------------
return                 : rien. 
----------------------------------------------------------------------------*/
void 
deleteNumExpr (expr, i)
     chain_list *expr;
     int i;
{
  int co = 0;
  chain_list /*e = expr,*/ *pt;

  if (lengthExpr (expr) <= i)
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"168");
//      printf ("deleteNumExpr : error - index %d out of bound for \n", i);
//      displayExpr (e);
//      EXIT (-1);
    }

  pt = expr;
  while ((i >= co) && (expr = CDR (expr)))
    {
      if (i == co)
        {
          CDR (pt) = CDR (expr);
          CDR (expr) = NULL;        /*desallocation du i-eme arg */
          freeExpr (CAR (expr));
          freechain (expr);
          return;
        }
      else
        pt = CDR (pt);
      co++;
    }
}

/*----------------------------------------------------------------------------
searchNumExpr         : renvoie le i-eme element de expr.
                  le premier argument commence a la place 0.
------------------------------------------------------
parametres          : un ABL et une place. 
-------------------------------------------------------
return                 : un chain_list * 
----------------------------------------------------------------------------*/
chain_list *
searchNumExpr (expr, i)
     chain_list *expr;
     int i;
{
  int co = 0;
//  chain_list *e = expr;

  while ((expr = CDR (expr)))
    {
      if (i == co)
        return CAR (expr);
      co++;
    }

  avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"169");
//  printf ("\nsearchNumExpr : error - index %d out of bound for \n", i);
//  displayExpr (e);
//  EXIT (-1);
#ifndef __ALL__WARNING__
  return 0; // ..anto..
#endif
}

/*----------------------------------------------------------------------------
numberOccExpr         : renvoie le nombre d'occurences d'une variable dans expr. 
------------------------------------------------------
parametres          : un ABL et une variable. 
-------------------------------------------------------
return                 : un int. 
----------------------------------------------------------------------------*/
int 
numberOccExpr (exp, name)
     chain_list *exp;
     char *name;
{
  if (ATOM (exp))
    {
      if ((long) name == (long) VALUE_ATOM (exp))
        return 1;
      else
        return 0;
    }
  else
    {
      int cpt = 0;

      while ((exp = CDR (exp)))
        cpt += numberOccExpr (CAR (exp), name);
      return cpt;
    }
}

/*----------------------------------------------------------------------------
changeOperExpr         : change l'operateur d'une expression 
------------------------------------------------------
parametres          : un ABL
-------------------------------------------------------
return                 : un void. 
----------------------------------------------------------------------------*/
void 
changeOperExpr (expr, oper)
     chain_list *expr;
     short oper;
{
  int op = (int) oper;

  if (ATOM (expr))
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"170");
//      printf ("changeOperExpr : error - atomic expression\n");
//      EXIT (-1);
    }
  else
    {
      expr = CAR (expr);
      expr->DATA = (void *)((long) op);
    }
}

/*----------------------------------------------------------------------------
simplif10Expr         : simplifie une expression en eliminant les constantes 0 et 1. 
------------------------------------------------------
parametres          : un ABL
-------------------------------------------------------
return                 : un ABL. 
----------------------------------------------------------------------------*/
chain_list *
simplif10Expr (expr)
     chain_list *expr;
{
  if (ATOM (expr))
  {
    return copyExpr (expr);
  }
  else
  {
    short oper = OPER (expr);
    chain_list *auxExpr = createExpr (oper), *args, *res;
    int cpt1 = 0;
    int i;

    switch (oper)
    {
      case REPLICATE :
           CDAR(auxExpr) = addchain(NULL     ,(void *)REP_NB(expr));
      case CAT       : 
           CDAR(auxExpr) = addchain(CDAR(auxExpr),(void *)VECT_S(expr));
      default        : ;
    }

    for (args = CDR (expr); args; args = CDR (args))
    {
      res = simplif10Expr (CAR (args));
      if (!ATOM (res))
      {
        if (oper == STABLE)
        {
          addQExpr (auxExpr, res);
          return auxExpr;
        }
        /* simplif (not (not a)) -> a */

        if ((oper == NOT) && (OPER (res) == NOT))
        {
          freeExpr (auxExpr);
          return notExpr (res);
        }
        else
          addQExpr (auxExpr, res);
      }
      else
      {
        switch (oper)
        {
          case NOT:
               if (!strcmp (VALUE_ATOM (res), "'0'"))
               {
                 freeExpr (res);
                 freeExpr (auxExpr);
                 return createAtom ("'1'");
               }
               if (!strcmp (VALUE_ATOM (res), "'1'"))
               {
                 freeExpr (res);
                 freeExpr (auxExpr);
                 return createAtom ("'0'");
               }
               addQExpr (auxExpr, res);
               return auxExpr;
               break;
               
              // added by anto to compute vectabl 
          case CAT:
               addQExpr (auxExpr, res);
               break;

          case REPLICATE:
               OPER(auxExpr) = CAT;
               for (i=0;i<REP_NB(expr);i++) 
                 addQExpr (auxExpr, copyExpr(res));
               freeExpr(res);
               break;
              ////////////////////////////////////
               
          case STABLE:
               addQExpr (auxExpr, res);
               return auxExpr;
               break;

          case AND:
               if (!strcmp (VALUE_ATOM (res), "'0'"))
               {
                 freeExpr (res);
                 freeExpr (auxExpr);
                 return createAtom ("'0'");
               }
               if (strcmp (VALUE_ATOM (res), "'1'"))
                 addQExpr (auxExpr, res);
               else
                 freeExpr(res);
               break;
          case NAND:
               if (!strcmp (VALUE_ATOM (res), "'0'"))
               {
                 freeExpr (res);
                 freeExpr (auxExpr);
                 return createAtom ("'1'");
               }
               if (strcmp (VALUE_ATOM (res), "'1'"))
                 addQExpr (auxExpr, res);
               else
                 freeExpr(res);
               break;
          case OR:
               if (!strcmp (VALUE_ATOM (res), "'1'"))
               {
                 freeExpr (res);
                 freeExpr (auxExpr);
                 return createAtom ("'1'");
               }
               if (strcmp (VALUE_ATOM (res), "'0'"))
                 addQExpr (auxExpr, res);
               else
                 freeExpr(res);
               break;
          case NOR:
               if (!strcmp (VALUE_ATOM (res), "'1'"))
               {
                 freeExpr (res);
                 freeExpr (auxExpr);
                 return createAtom ("'0'");
               }
               if (strcmp (VALUE_ATOM (res), "'0'"))
                 addQExpr (auxExpr, res);
               else
                 freeExpr(res);
               break;
          case XOR:
               if (!strcmp (VALUE_ATOM (res), "'1'"))
               {
                 freeExpr (res);
                 cpt1++;
               }
               else if (strcmp (VALUE_ATOM (res), "'0'"))
                 addQExpr (auxExpr, res);
               else
                 freeExpr(res);
               break;
          case NXOR:
               if (!strcmp (VALUE_ATOM (res), "'1'"))
               {
                 freeExpr (res);
                 cpt1++;
               }
               else if (strcmp (VALUE_ATOM (res), "'0'"))
                 addQExpr (auxExpr, res);
               else
                 freeExpr(res);
               break;
        }
      }
    }
    /* cas des XOR */

    if (((cpt1 % 2) == 1) && CDR (auxExpr))
    {
      if (lengthExpr (auxExpr) < 2)
      {
        if (!CDR (auxExpr))
        {
          freeExpr (auxExpr);
          return createAtom ("'1'");
        }
        else
          changeOperExpr (auxExpr, NOT);
      }
      else
      {
        chain_list *aux = createExpr (NOT);
        addQExpr (aux, auxExpr);
        auxExpr = aux;
      }
    }
    /* expression normale */

    if ((lengthExpr (auxExpr) > 1) || (OPER (auxExpr) == NOT) ||
        (OPER (auxExpr) == REPLICATE) )
      return auxExpr;

    /* traitement pas d'argument */
    if (!CDR (auxExpr))
    {
      freeExpr (auxExpr);
      if ((oper == NAND) || (oper == OR))
        return createAtom ("'0'");
      if ((oper == AND) || (oper == NOR))
        return createAtom ("'1'");
      if (oper == XOR) {
        if ((cpt1 % 2) == 1)
          return createAtom ("'1'");
        else
          return createAtom ("'0'");
      }
      if (oper == NXOR) {
        if ((cpt1 % 2) == 1)
          return createAtom ("'0'");
        else
          return createAtom ("'1'");
      }
    }
    /* traitement 1 seul argument */

    if ((oper == AND) || (oper == OR))
    {
      res=CADR (auxExpr);
      CADR (auxExpr)=createAtom ("'1'"); freeExpr(auxExpr);
      return res;
    }
    if ((oper == NAND) || (oper == NOR))
    {
      changeOperExpr (auxExpr, NOT);
      return auxExpr;
    }
    if (oper == XOR) {
      if ((cpt1 % 2) == 1)
      {
        changeOperExpr (auxExpr, NOT);
        return auxExpr;
      }
      else
      {
        res=CADR (auxExpr);
        CADR (auxExpr)=createAtom ("'1'"); freeExpr(auxExpr);
        return res;
      }
    }
    if (oper == NXOR) {
      if ((cpt1 % 2) == 1)
      {
        res=CADR (auxExpr);
        CADR (auxExpr)=createAtom ("'1'"); freeExpr(auxExpr);
        return res;
      }
      else
      {
        changeOperExpr (auxExpr, NOT);
        return auxExpr;
      }
    }
  }

  avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"171");
//  printf ("## ERROR in simplif10Expr, %s.%d\n", __FILE__, __LINE__);
//  EXIT (-1);

#ifndef __ALL__WARNING__
  return 0; // ..anto..
#endif
}

/*----------------------------------------------------------------------------
simplifNotExpr         : simplifie une expression en appliquant (not (not a)) --> a 
------------------------------------------------------
parametres          : un ABL
-------------------------------------------------------
return                 : un ABL. 
----------------------------------------------------------------------------*/
chain_list *
simplifNotExpr (exp)
     chain_list *exp;
{
  if (ATOM (exp))
    return copyExpr (exp);
  if (OPER (exp) == NOT)
    {
      if (!ATOM (CADR (exp)) && OPER (CADR (exp)) == NOT)        /* (not (not a)) */
        return (simplifNotExpr (CADR (CADR (exp))));
      else
        {
          chain_list *expr = createExpr (NOT);

          addQExpr (expr, simplifNotExpr (CADR (exp)));
          return expr;
        }
    }
  else
    {
      return (mapCarExpr (simplifNotExpr, OPER (exp), exp));
    }
}

/*----------------------------------------------------------------------------
charToExprInt :        analyse une chaine de caracteres decrivant une
                equation logique et remplit l'expression de retour. 
-------------------------------------------------------
parametres :        char.
-------------------------------------------------------
return :        expr.
----------------------------------------------------------------------------*/
chain_list *
charToExprInt (stringExpr, cptCar)
     char *stringExpr;
     int *cptCar;
{
  char *tokOper;
  char *tok;
//  char *tokenExpr ();
  chain_list *expr;

  tok = tokenExpr (stringExpr, cptCar);
  if (!tok)
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"172");
//      printf ("charToExpr : error - bad string to convert\n");
//      printf ("%s\n", stringExpr);
//      EXIT (-1);
    }

  if (strcmp (tok, "("))        /* noeud terminal direct */
    {
      return createAtom (tok);
    }

  mbkfree (tok);

  tokOper = tokenExpr (stringExpr, cptCar);
  if (!tokOper)
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"173");
//      printf ("charToExpr : error - bad string to convert\n");
//      printf ("%s\n", stringExpr);
//      EXIT (-1);
    }
  switch (*tokOper)
    {
    case 'o':
      expr = createExpr (OR);
      break;
    case 'a':
      expr = createExpr (AND);
      break;
    case 'x':
      expr = createExpr (XOR);
      break;
    case 's':
      expr = createExpr (STABLE);
      break;
    case 'n':
      switch (*(tokOper + 2))
        {
        case 't':
          expr = createExpr (NOT);
          break;
        case 'r':
          expr = createExpr (NOR);
          break;
        case 'n':
          expr = createExpr (NAND);
          break;
        case 'o':
          expr = createExpr (NXOR);
          break;
        default:
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"174");
//          printf ("charToExpr : error - bad operator %s\n", tokOper);
//          EXIT (-1);
        }
      break;
    default:
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"175");
//      printf ("charToExpr : error - bad operator %s\n", tokOper);
//      EXIT (-1);
    }
  mbkfree (tokOper);

  tok = tokenExpr (stringExpr, cptCar);                /* lecture des operandes */
  if (!tok)
    {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"176");
//      printf ("charToExpr : error - bad string to convert\n");
//      printf ("%s\n", stringExpr);
//      EXIT (-1);
    }

  do
    {
      /* traitement des elements terminaux */

      if (strcmp (tok, "("))
        {
          addQExpr (expr, createAtom (tok));
          mbkfree (tok);
        }
      else
        {
          /* traitement des sous-expressions */

          mbkfree (tok);
          *cptCar = *cptCar - 1;        /*compteur revient sur "(" */

          addQExpr (expr, charToExprInt (stringExpr, cptCar));
        }

      tok = tokenExpr (stringExpr, cptCar);
      if (!tok)
        {
      avt_errmsg(LOG_ERRMSG,"000",AVT_FATAL,"177");
//          printf ("charToExpr : error - bad string to convert\n");
//          printf ("%s\n", stringExpr);
//          EXIT (-1);
        }
    }
  while (strcmp (tok, ")"));
  mbkfree (tok);

  /* on applique maintenant l'operateur */

  return (expr);
}

/*----------------------------------------------------------------------------
charToExpr :        analyse une chaine de caracteres decrivant une
                equation logique et remplit l'expression de retour. 
-------------------------------------------------------
parametres :        char.
-------------------------------------------------------
return :        expr.
----------------------------------------------------------------------------*/
chain_list *
charToExpr (stringExpr)
     char *stringExpr;
{
  int cptCar = 0;

  return charToExprInt (stringExpr, &cptCar);
}

/*----------------------------------------------------------------------------
tokenExpr :        cherche les tokenExprs.
-------------------------------------------------------
parametres :
-------------------------------------------------------
return :        (, ), ou une chaine de caractere.
----------------------------------------------------------------------------*/
char *
tokenExpr (stringExpr, cptCar)
     char *stringExpr;
     int *cptCar;
{
  char *tok;
  int i;
  int compteur = *cptCar;

  for (; isspace ((int)stringExpr[compteur]); compteur++);        /* suppression ds blancs */

  if (stringExpr[compteur] == '\0')        /* probleme ... */
    return (NULL);

  i = compteur;

  /* on avance jusqu'au prochain tokenExpr */

  for (; isspace ((int)stringExpr[i]) == 0 && strncmp (stringExpr + i, ")", 1) != 0 && strncmp (stringExpr + i, "(", 1) != 0; i++);

  /* parenthese ouvrante ou fermante */

  if (!strncmp (stringExpr + compteur, "(", 1) || !strncmp (stringExpr + compteur, ")", 1))
    {
      tok = (char *) mbkalloc (2);
      strncpy (tok, stringExpr + i, 1);
      tok[1] = '\0';
      compteur++;
      *cptCar = compteur;
      return (tok);
    }

  /* operateur ou atome */

  tok = (char *) mbkalloc (i - compteur + 1);
  strncpy (tok, stringExpr + compteur, i - compteur + 1);
  tok[i - compteur] = '\0';
  *cptCar = i;
  return (tok);
}

/*----------------------------------------------------------------------------
PMExprInt   : Pattern Matching entre expr et pattern
           (and (not a) (not a)) avec le pattern (and x x) renvoie 1
------------------------------------------------------------------------------
retour   : 0 ou 1
----------------------------------------------------------------------------*/ int 
PMExprInt (expr, pattern, bind)
     chain_list *expr, *pattern;
     ptype_list **bind;
{
  if (ATOM (pattern))
    {
      ptype_list *find = getptype (*bind, (long)VALUE_ATOM (pattern));

      /* le couple existe deja */

      if (find)
        {
          if (!equalExpr (expr, (chain_list *) find->DATA))
            return (0);
        }
      else
        *bind = addptype (*bind, (long)VALUE_ATOM (pattern), expr);

      return (1);
    }
  else
    /* pattern n'est pas un atome */
    {
      if (ATOM (expr))
        return (0);
      else
        /* expr n'est pas un atome */
        {
          if (OPER (expr) != OPER (pattern))
            return (0);

          if (lengthExpr (expr) != lengthExpr (pattern))
            return (0);

          while ((pattern = CDR (pattern)))
            {
              expr = CDR (expr);
              if (PMExprInt (CAR (expr), CAR (pattern), bind) == 0)
                return (0);
            }
          return (1);
        }
    }
}

/*----------------------------------------------------------------------------
PMExpr   : Pattern Matching entre expr et pattern
           (and (not a) (not a)) avec le pattern (and x x) renvoie 1
------------------------------------------------------------------------------
retour   : 0 ou 1
----------------------------------------------------------------------------*/
int 
PMExpr (expr, pattern)
     chain_list *expr, *pattern;
{
  int ret;
  ptype_list *bind = NULL;

  /* on mettra dans "bind" le couple (char * pattern, expr) */

  ret = PMExprInt (expr, pattern, &bind);
  freeptype (bind);

  return ret;
}

/*}}}************************************************************************/
/*{{{                    replaceAtomByExpr()                                */
/*                                                                          */
/* replace an atom by a expression                                          */
/****************************************************************************/
void replaceAtomByExpr(chain_list *expr, char *name, chain_list *newexpr)
{
  if (VALUE_ATOM(expr) == name)
  {
    chain_list  *copy;

    copy = copyExpr(newexpr);
    expr->DATA = copy->DATA;
    expr->NEXT = delchain(copy,copy);
  }
  else
  {
    while ((expr = CDR (expr)))
    {
      replaceAtomByExpr (CAR (expr), name, newexpr);
    }
  }
}

static chain_list *do_not_AND(ptype_list *andlist)
{
  chain_list *cl=NULL;
  ptype_list *next;
  while (andlist!=NULL)
  {
    next=andlist->NEXT;
    andlist->TYPE=(andlist->TYPE+1)&1;
    andlist->NEXT=NULL;
    cl=addchain(cl, andlist);
    andlist=next;
  }
  return cl;
}

static chain_list *do_AND(chain_list *neworlist, chain_list *orlist)
{
   chain_list *tmporlist, *cl, *ch;
   ptype_list *andlist;
   if (neworlist==NULL) return orlist;
   else if (orlist==NULL) return neworlist;
   else
    {
      tmporlist=NULL;
      for (cl=neworlist; cl!=NULL; cl=cl->NEXT)
        for (ch=orlist; ch!=NULL; ch=ch->NEXT)
        {
          andlist=(ptype_list *)append((chain_list *)dupptypelst((ptype_list *)cl->DATA), (chain_list *)dupptypelst((ptype_list *)ch->DATA));
          tmporlist=addchain(tmporlist, andlist);
        }
      for (cl=neworlist; cl!=NULL; cl=delchain(cl, cl)) freeptype((ptype_list *)cl->DATA);
      for (cl=orlist; cl!=NULL; cl=delchain(cl, cl)) freeptype((ptype_list *)cl->DATA);
      neworlist=tmporlist;
    }
   return neworlist;
}


static chain_list *expandabl(chain_list *abl)
{
  short oper;
  chain_list *orlist, *neworlist=NULL, *tmporlist, *expr;
  
  if (ATOM (abl)) 
    return addchain(NULL, addptype(NULL, 0, VALUE_ATOM (abl)));
  else
  {
    oper = OPER (abl);
    if (oper >= MIN_OPER && oper <= MAX_OPER)
    {
      if (oper==AND || oper==OR || oper==NOT)
      {
         expr=abl;
         while ((expr = CDR (expr)))
         {
           orlist=expandabl(CAR (expr));
           switch (oper)
           {
             case OR:
                neworlist=append(neworlist, orlist);
                break;
             case AND:
                neworlist=do_AND(neworlist, orlist);
                break;
             case NOT:
                while (orlist!=NULL)
                {
                  tmporlist=do_not_AND((ptype_list *)orlist->DATA);
                  neworlist=do_AND(neworlist, tmporlist);
                  orlist=delchain(orlist, orlist);
                }
                break;
           }
         }
      }
    }
  }
 return neworlist;
}

static chain_list *getexpr(ptype_list *pt)
{
  chain_list *expr;
  expr=createAtom((char *)pt->DATA);
  if (pt->TYPE) expr=notExpr(expr);
  return expr;
}

static void replace_with_expressions(chain_list *mterm)
{
  chain_list *cl, *expr;
  ptype_list *pt;
  for (cl=mterm; cl!=NULL; cl=cl->NEXT)
  {
     pt=(ptype_list *)cl->DATA;
     if (pt->NEXT==NULL) expr=getexpr(pt);
     else
     {
       expr=createExpr (AND);
       while (pt!=NULL)
       {
         addQExpr(expr, getexpr(pt));
         pt=pt->NEXT;
       }
     }
    freeptype((ptype_list *)cl->DATA);
    cl->DATA=expr;
  }
  
}

chain_list *ablCreateMintermList(chain_list *expr)
{
  chain_list *lst;
  lst=expandabl(expr);
  replace_with_expressions(lst);
  return lst;
}


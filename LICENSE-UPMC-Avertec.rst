

.. -*- Mode: rst -*-

.. role:: raw-latex(raw)
   :format: latex

.. role:: ul
.. role:: cb
.. role:: sc

.. |VHDL|                             replace:: :sc:`vhdl`
.. |Verilog|                          replace:: :sc:`Verilog`
.. |HiTas|                            replace:: :sc:`HiTas`
.. |Tas|                              replace:: :sc:`Tas`
.. |Yagle|                            replace:: :sc:`Yagle`


.. raw:: latex

   \tableofcontents
   \newpage


=========================
SOFWARE LICENSE AGREEMENT
=========================

:raw-latex:`\noindent`
Between

  **Université Pierre et Marie Curie  (Paris 6)**, a not for profit corporation
  under  the laws of  France, N°  SIRET :  19751722000012 -  Code APE  : 8542Z,
  having  its place  of business  at  4, place  Jussieu, 75252  Paris cedex  5,
  France,  represented by  its  President Mr  Jean :sc:`Chambaz`,  (hereinafter
  called "UPMC"),

  **CNRS**

:raw-latex:`\medskip\noindent`
On the one hand, 

:raw-latex:`\medskip\noindent`
And

  **The RECIPIENT**
  
:raw-latex:`\medskip\noindent`
On the other hand,

:raw-latex:`\medskip\noindent`
Hereinafter solely or collectively designed "Party/Parties"

:raw-latex:`\medskip`
Whereas the  Laboratoire LIP6 (UMR 7606)  at UPMC has developed  a new software
called  HITAS/YAGLE, hereafter  referred  to as  the  "SOFTWARE", which  allows
hierarchical static timing analysis of VLSI designs.

|HiTas| is a static timing analysis tool.  Its strength lies in its transparent
hierarchical  approach combined  with the  ability to  perform analysis  at the
transistor-level, cell-level or a mixture of the two.

The  transistor-level analysis  brings the  possiblity of  handling full-custom
circuits not individually but also as  blocks within the hierachy of a complete
chip.  In addition, working at transistor level removes the need for costly and
time consuming  re-characterization when  performing the analysis  at different
corners.  The fact that delays are dynamically calculated throughout the design
means  that the  differences  in  local context  are  automatically taken  into
account  such as  power  supply variations  due  to IR  drop,  or simply  using
different voltages for low-power applications.

|Yagle|  provides automatic  generation  of  a behaviral  model  (in |VHDL|  or
|Verilog|)  directly from  a transistor  netlist.   Its major  strength is  the
ability  to take  into account  functional correlation  between  signals.  This
optimizes the partitioning and functional characterization as well as providing
the means for automatically identifying  and characterizing most kind of memory
elements.   Yagle  is  able  to  mix  the totally  automatic  approach  with  a
pattern-matching  approach.   This allows  functional  abstraction of  circuits
containing a  mixture of analog and digital.   This can be used  to handle RAMs
which contain sense amplifiers for example.

:raw-latex:`\medskip`
Whereas UPMC desire to make the SOFTWARE available for public use and benefit.

:raw-latex:`\medskip`
Whereas the RECIPIENT, wishes to use the SOFTWARE **for research purposes** and
asked UPMC for a copy of the SOFTWARE.

:raw-latex:`\medskip`
**NOW THEREFORE, IT IS HEREBY AGREED BETWEEN THE PARTIES:**


1. GRANT OF RIGHTS:
===================

Subject  to the  provisions contained  herein, UPMC  hereby grants  RECIPIENT a
non-exclusive royalty-free non-transferable rights  to use the SOFTWARE and all
relating  documentations  for  research  purposes  for a  period  of  2  years,
effective from the date of download of the SOFTWARE and the present agreement.

Upon execution of  the present agreement by the  RECIPIENT, UPMC authorizes the
download of the SOFTWARE by the RECIPIENT.

In  the event  that the  RECIPIENT wishes  to use  the SOFTWARE  for commercial
and/or industrialization  purpose, it recognizes that this  action requires the
express  and prior authorization  of UPMC,  in accordance  with article  6. The
terms  of obtaining  a commercial  license have  to be  negotiated  between the
Parties, on a case by case base.

All  rights, title, interest  and copyright  to the  SOFTWARE, to  all portions
thereof,  and to any  associated documentation  shall at  all times  remain the
property of  UPMC. RECIPIENT  agrees that  to use the  SOFTWARE solely  for non
commercial purposes and in full  compliance with the disposition of the present
agreement.

No  part  of the  SOFTWARE  or  of the  accompanying  written  material may  be
reproduced, transmitted to another location or to any other person, stored in a
retrieval system, or  translated into any language or  computer language in any
form other than  granted above by any means,  electronic, mechanical, magnetic,
optical, chemical, manual, or  otherwise without the express written permission
of UPMC.

The RECIPIENT  may not modify,  reverse engineer, decompile or  disassemble the
SOFTWARE or related documentations. The RECIPIENT may not use, copy, modify, or
transfer the SOFTWARE or documentation or any copy except as expressly provided
in this agreement.


2. FDA AND OTHER APPROVALS
==========================

The RECIPIENT  agrees that  this SOFTWARE has  not been reviewed,  nor received
clearance for marketing from any health- regulation-agency such as the Food and
Drug Administration,  Health Canada or  Agence Française de  Sécurité Sanitaire
des Produits de Santé, in any country.


3. DISCLAIMER OF WARRANTY:
==========================

This SOFTWARE is ©  copyright UPMC – 2011. UPMC holds all  the ownership on the
SOFTWARE.

UPMC and the authors of the SOFTWARE are hereinafter called the “DISCLOSERS”.

RECIPIENT acknowledges that  the SOFTWARE is a research tool,  that it is being
supplied "as is" and that DISCLOSERS are not committed to provide any services,
improvements or updates.

DISCLOSERS make no representation or  warranties, express or implied. By way of
example, but not limitation to, DISCLOSERS make no representation or warranties
of merchantability or fitness for any particular purpose or that the use of the
SOFTWARE  will  not infringe  any  patents,  copyrights,  trademarks, or  other
rights.  DISCLOSERS shall  not  be liable  for  any liability  or damages  with
respect to any claim by RECIPIENT or  any third party on account of, or arising
from, this licence or use of the SOFTWARE.

DISCLOSERS  shall not  be held  liable for  any liability  nor for  any direct,
indirect, or  consequential damages with respect  to any claim  by RECIPIENT or
any third party on account of or arising from the use of the SOFTWARE.

DISCLOSERS are not liable for  any hardware components used in conjunction with
this SOFTWARE. DISCLOSERS are not liable for any failure of hardware components
used  with this  SOFTWARE. If  failure of  the disk  or hardware  component has
resulted from  accident, abuse, or  misapplication of the  SOFTWARE, DISCLOSERS
shall have  no responsibility to replace  the disk or  hardware component under
this limited warranty.

The entire risk as to the results and performance of the SOFTWARE is assumed by
the RECIPIENT. Should  the SOFTWARE prove defective, the  RECIPIENT will assume
all costs of  necessary service, repair, or correction.  Further, DISCLOSERS do
not warrant,  guarantee, or make any  representations regarding the  use of the
SOFTWARE  in  terms  of  correctness, accuracy,  reliability,  currentness,  or
otherwise ; and the RECIPIENT relies on  the SOFTWARE and the results solely at
his own risk.


4. LIMITATION OF LIABILITY - INDEMNITY:
=======================================

Under no circumstances and under no  legal theory, whether in tort, contract or
otherwise, shall  UPMC or anyone  else who has  been involved in  the creation,
production, or  delivery of this SOFTWARE  be liable to RECIPIENT  or any other
person for any direct,  indirect, special, incidental, or consequential damages
of any character  including, without limitation, damages for  loss of goodwill,
work stoppage, computer failure or malfunction, or any and all other damages or
losses, arising out  of the use, the  results of use, or inability  to use such
product,  even if  UPMC shall  have been  informed of  the possibility  of such
damages, or for any claim by any other party.

To the  extent allowed  by law, RECIPIENT  shall indemnify, hold  harmless, and
defend UPMC, its officers, employees,  students, and agents against any and all
claims  arising  out  of the  exercise  of  any  rights under  this  agreement,
including,  without  limiting the  generality  of  the  foregoing, against  any
damages, losses, or  liabilities whatsoever with respect to  death or injury to
person or  damage to property  arising from or  out of the possession,  use, or
operation of the SOFTWARE by the RECIPIENT.


5. PUBLICATION - ACKNOWLEDGMENT OF CONTRIBUTION - USE OF NAME
=============================================================

For  any publication  or  communication of  results,  information or  knowledge
obtained with the utilization of the SOFTWARE, and for any published work based
on the  SOFTWARE, the RECIPIENT commits  itself to indicating  that the results
are obtained using **"HITAS/YAGLE SOFTWARE, UPMC/LIP6"**.

5.2 Nothing however  in this agreement shall be  construed as conferring rights
to use in advertising, publicity, or otherwise  the name of UPMC, of any of its
employees or any of its marks.


6. NOTICES
==========

Any notices or disclosures required or  provided by the terms of this agreement
shall be in writing, and shall  be delivered personally or sent by certified or
registered   mail,   return   receipt   requested,  postage   prepaid   or   by
internationally-recognized   express  mail   service   providing  evidence   of
delivery. The effective  date of any notice shall be the  date of first receipt
by the receiving Party or the date of refusal of receipt. Notices shall be sent
to the addresses/addressees given below:

- Technical Contact:

    | Laboratoire LIP6
    | Département SoC
    | Équipe CIAN
    | 4, place Jussieu
    | F-75252 Paris cedex 05
    | c/o: Mr. Jean-Paul Chaput

- Administrative Contact:

    | Université Pierre et Marie Curie (Paris 6)
    | Direction de la Recherche et du Transfert de Technologies
    | Tour Zamansky
    | 4, Place Jussieu
    | 75252 Paris cedex 05
    | c/o: Mr Laurent BUISSON - Ref. UPMC : X11xxx


7. TERMINATION
==============

This  Agreement will  be terminated  by  UPMC two  (2) years  after a  SOFTWARE
download by  the RECIPIENT. Upon such termination,  RECIPIENT shall immediately
cease all uses of the Software and destroy the SOFTWARE and all copies thereof.

The provisions of Articles 3, 4, 5, 8 and 9 shall survive any termination.


8. GOVERNING LAW
================

This agreement shall  be construed, interpreted and applied  in accordance with
the laws of France. If any claims or lawsuits concerning this license agreement
of the  SOFTWARE are brought  against UPMC, the  Parties agree to  endeavour to
seek  an amicable  solution to  any disagreements  or disputes  that  may arise
during the  performance of the  agreement. Failing an amicable  solution within
three (3) months as from their occurrence, and unless emergency proceedings are
necessary, the disagreements or disputes  shall be referred to the jurisdiction
of Paris which will be exclusively competent.

9. MISCELLANEOUS:
=================

This Agreement  constitutes the complete  and exclusive agreement  between UPMC
and the RECIPIENT with respect to the subject matter hereof, and supersedes all
prior  oral  or  written   understandings,  communications  or  agreements  not
specifically incorporated  herein. This agreement  may not be modified.  If any
provision of  this agreement is held  to be unenforceable for  any reason, such
provision  shall  be  reformed  only   to  the  extent  necessary  to  make  it
enforceable, and such decision shall  not affect the enforceability (i) of such
provision under other circumstances, or (ii) of the remaining provisions hereof
under all circumstances.

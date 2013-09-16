/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2010-2012                                                 |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 |      version: 0.2                                                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/*
 *    This file is part of GF2 Toolkit.
 *
 *    GF2 toolkit -- A small library for Factorization in Galois Field F2.
 *    Copyright (C) 2010-2012 by Enrico Bertolazzi.
 *    All rights reserved.
 *
 *    GF2 toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    GF2 toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_MatrixLower.hh"
#include "GF2toolkit_4Russian.hh"
#include "GF2toolkit_Strassen.hh"

namespace GF2toolkit {

  /*
  //  Moltiplico a destra per le matrici
  //
  //  +---+---+-...--+---+  +---+---+---+---+-----+---+ (-1)
  //  |   |   |      |   |  | U0|   |   |   |     |   |
  //  |   |   |      |   |  +---+ U1|   |   |     |   |
  //  |   |   |      |   |  |   |   | U2|   |     |   |
  //  |   |   |      |   |  |   +---|   | U3|     |   |
  //  | A0| A1|      | An|  |       |   |   |     |   |
  //  |   |   |      |   |  |       +---|   |     |   |
  //  |   |   |      |   |  |           |   |     | Un|
  //  |   |   |      |   |  |           +---+     |   |
  //  |   |   |      |   |  |                 .   |   |
  //  |   |   |      |   |  |                   . |   |
  //  +---+---+-...--+---+  |                     |   |
  //                        +---------------------+---+
  //
  //  +---+   +---+ +---+
  //  |   |   |   | | U0|
  //  |   |   |   | +---+
  //  |   |   |   |
  //  |   |   |   |
  //  | A0| = | X0|
  //  |   |   |   |
  //  |   |   |   |
  //  |   |   |   |
  //  |   |   |   |
  //  |   |   |   |
  //  +---+   +---+
  //
  //  +---+   +---+---+ +---+          +---+   +---+ +---+   +---+ +...+
  //  |   |   |   |   | |   |          |   |   |   | | U1|   |   | :   :
  //  |   |   |   |   | | U1|          |   |   |   | +---+   |   | +---+
  //  |   |   |   |   | |   |          |   |   |   | :   :   |   | | U1|
  //  |   |   |   |   | +---+          |   |   |   | +...+   |   | +---+  
  //  | A1| = | X0| X1|                | A1| - | X0|       = | X1|        
  //  |   |   |   |   |                |   |   |   |         |   |        
  //  |   |   |   |   |                |   |   |   |         |   |        
  //  |   |   |   |   |                |   |   |   |         |   |        
  //  |   |   |   |   |                |   |   |   |         |   |        
  //  |   |   |   |   |                |   |   |   |         |   |        
  //  +---+   +---+---+                +---+   +---+         +---+          
  //
  //  +---+   +---+---+---+ +---+      +---+   +---+---+ +---+   +---+ +...+
  //  |   |   |   |   |   | |   |      |   |   |   |   | |   |   |   | :   :
  //  |   |   |   |   |   | |   |      |   |   |   |   | | U2|   |   | :   :
  //  |   |   |   |   |   | | U2|      |   |   |   |   | |   |   |   | :   :
  //  |   |   |   |   |   | |   |      |   |   |   |   | +---+   |   | +---+  
  //  | A2| = | X0| X1| X2| |   |      | A2| - | X0| X1| :   : = | X2| | U1|  
  //  |   |   |   |   |   | +---+      |   |   |   |   | +...+   |   | +---+      
  //  |   |   |   |   |   |            |   |   |   |   |         |   |        
  //  |   |   |   |   |   |            |   |   |   |   |         |   |        
  //  |   |   |   |   |   |            |   |   |   |   |         |   |        
  //  |   |   |   |   |   |            |   |   |   |   |         |   |        
  //  +---+   +---+---+---+            +---+   +---+---+         +---+          
  //
  */
  template <typename UNSIGNED>
  void
  UinvApply( UNSIGNED const * U, unsigned dimRowsU, unsigned nbStart, unsigned nColsBlockU,
             UNSIGNED       * A, unsigned dimRowsA, unsigned nRowsA ) {

    MM4R<UNSIGNED,8> MMop ;
    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;

    for ( unsigned nb = nbStart ; nb < nColsBlockU ; ++nb ) {
      UNSIGNED const * Uk = U + nb * dimRowsU ;
      UNSIGNED       * Ak = A + nb * dimRowsA ;
      for ( unsigned j = 0 ; j < nb ; ++j ) {
        MMop.makeTable( Uk + j * numBits ) ;
        MMop.multRightAdd( A + j * dimRowsA, nRowsA, Ak ) ;
      }
      MMop . makeTable( Uk + nb * numBits ) ; // Uk sulla diagonale è l'inversa del blocco diagonale
      // MMop . mul( Ak, nRowsA ) ; vecchia instruzione
      MMop . multRight( Ak, nRowsA ) ;
    }
  }
  
  template
  void
  UinvApply( uint32_t const * U, unsigned dimRowsU, unsigned nbStart, unsigned nColsBlockU,
             uint32_t       * A, unsigned dimRowsA, unsigned nRowsA ) ;

  template
  void
  UinvApply( uint64_t const * U, unsigned dimRowsU, unsigned nbStart, unsigned nColsBlockU,
             uint64_t       * A, unsigned dimRowsA, unsigned nRowsA ) ;

  template
  void
  UinvApply( uint128_t const * U, unsigned dimRowsU, unsigned nbStart, unsigned nColsBlockU,
             uint128_t       * A, unsigned dimRowsA, unsigned nRowsA ) ;

}

// EOF GF2toolkit_MatrixSolve.cc

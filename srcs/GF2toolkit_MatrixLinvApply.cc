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

#include "GF2toolkit_Matrix.hh"
#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_MatrixLower.hh"
#include "GF2toolkit_4Russian.hh"
#include "GF2toolkit_Strassen.hh"

namespace GF2toolkit {

  /*
  //  Moltiplico a sinistra per le matrici
  //                                (-1)                             (-1)
  //  +---+---+---------------------+  +---+---+---------------------+  +---+
  //  | I | 0 | 0 . . . . . . . . 0 |  | L | 0 | 0 . . . . . . . . 0 |  | M0|
  //  +---O---+---------------------+  +---O---+---------------------+  +---+
  //  | 0 | L | 0 . . . . . . . . 0 |  |   | I | 0 . . . . . . . . 0 |  | M1|
  //  | : +---+---------------------+  |   +---+---------------------+  +---+
  //  | : |   | I                   |  |   | 0 | I                   |  | M2|
  //  | : |   |   .                 |  |   | : |   .                 |  | : |
  //  | : |   |     .               |  |   | : |     .               |  | : |
  //  | : |   |       .             |  | A0| : |       .             |  | : |
  //  | : | A1|         .           |  |   | : |         .           |  | : |
  //  | : |   |           .         |  |   | : |           .         |  | : |
  //  | : |   |             .       |  |   | : |             .       |  | : |
  //  | : |   |               .     |  |   | : |               .     |  | : |
  //  | : |   |                 .   |  |   | : |                 .   |  | : |
  //  | 0 |   |                   I |  |   | 0 |                   I |  | Mn|
  //  +---+---+---------------------+  +---+---+---------------------+  +---+ 
  //
  //  +---+
  //  | X0|
  //  +---+   +---+ +---+    
  //  | M1|   |   | | X0|    
  //  +---+   |   | +---+   +---+ +---+
  //  | M2|   |   |         |   | | X1|
  //  | : |   |   |         |   | +---+                   
  //  | : |   |   |         |   |                           
  //  | : |   | A0|         |   |                           
  //  | : | + |   |       + | A1|       + .... + +---+ +---+ 
  //  | : |   |   |         |   |                |   | | Xk|  
  //  | : |   |   |         |   |                |   | +---+  
  //  | : |   |   |         |   |                | Ak|           
  //  | : |   |   |         |   |                |   |           
  //  | Mn|   |   |         |   |                |   |           
  //  +---+   +---+         +---+                +---+                  
  //
  //  X0 = L^(-1) M0
  //  X1 = L^(-1) (M1-A0_1 * X0)
  //  ....
  //
  */
  template <typename UNSIGNED>
  void
  LinvApply( UNSIGNED const * L, unsigned dimRows, unsigned nRows, unsigned nColsBlock,
             UNSIGNED       * M ) {
 
    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;
    UNSIGNED Linv[numBits] ;

    for ( unsigned nb = 0 ; nb < nColsBlock ; ++nb ) {
      GF2toolkit::Linverse<UNSIGNED>( L, Linv ) ;
      GF2toolkit::MMassRight<UNSIGNED>( Linv, M ) ; // M = Linv * M
      L     += numBits ;
      nRows -= numBits ;
      GF2toolkit::MMadd<UNSIGNED>( L, nRows, M, M+numBits ) ;
      // passo a colonna successiva
      M     += numBits ;
      L     += dimRows ;
    }
  }

  /// versione ricorsiva per ottimizzare l'accesso in memoria
  template <typename UNSIGNED>
  void
  LinvApply( UNSIGNED       * W, unsigned Wsize,
             UNSIGNED const * L, unsigned dimRows,  unsigned nRows, unsigned nColsBlock,
             UNSIGNED       * M, unsigned dimRowsM,                 unsigned nColsBlockM ) {
 
    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;

    if ( numBits*nColsBlock > L_INV_APPLY_MIN_SIZE ) {
      /*
      //  +----+----+ +----+   +----+
      //  | L1 |  0 | | M1 |   | R1 |
      //  +----+----+ +----+ = +----+  
      //  | L2 | L3 | | M2 |   | R2 |
      //  +----+----+ +----+   +----+
      //
      //  M1 = L1^(-1) R1
      //  
      //  R2' = R2 - L2 * M1
      //
      //  M2 = L3^(-1) R2'
      //
      */
      unsigned nb = nColsBlock>>1 ;
      unsigned nr = nb*numBits ;
      
      UNSIGNED const * L1 = L ;
      UNSIGNED const * L2 = L  + nr ;
      UNSIGNED const * L3 = L2 + nb * dimRows ;

      UNSIGNED       * M1 = M ;
      UNSIGNED       * M2 = M + nr ;
      
      GF2toolkit::LinvApply<UNSIGNED>( W, Wsize,
                                       L1, dimRows, nr, nb,
                                       M1, dimRowsM, nColsBlockM ) ;

      GF2toolkit::MMaddStrassen<UNSIGNED>( W, Wsize,
                                           L2, dimRows, nRows - nr, nb,
                                           M1, dimRowsM, nr, nColsBlockM,
                                           M2, dimRowsM ) ;

      GF2toolkit::LinvApply<UNSIGNED>( W, Wsize, 
                                       L3, dimRows, nRows - nr, nColsBlock - nb,
                                       M2, dimRowsM, nColsBlockM ) ;

    } else {

      UNSIGNED Linv[numBits] ;
      for ( unsigned nb = 0 ; nb < nColsBlock ; ++nb ) {
        GF2toolkit::Linverse<UNSIGNED>( L, Linv ) ;
        GF2toolkit::MMassRight<UNSIGNED>( Linv, M, dimRowsM, nColsBlockM ) ; // M = Linv * M
        L     += numBits ;
        nRows -= numBits ;

        GF2toolkit::MMadd<UNSIGNED>( L, dimRows,  nRows,   1,
                                     M, dimRowsM, numBits, nColsBlockM,
                                     M+numBits, dimRowsM ) ;

        // passo a colonna successiva
        M     += numBits ;
        L     += dimRows ;
      }

    }

  }

  //////////////////////////////////////////////////////////////////////////////
  // INSTANZIAZIONE ESPLICITA
  //////////////////////////////////////////////////////////////////////////////

  template
  void
  LinvApply( uint32_t const *, unsigned, unsigned, unsigned, uint32_t * );

  template
  void
  LinvApply( uint64_t const *, unsigned, unsigned, unsigned, uint64_t * ) ;

  template
  void
  LinvApply( uint128_t const *, unsigned, unsigned, unsigned, uint128_t * ) ;

  template
  void
  LinvApply( uint32_t       * W, unsigned Wsize,
             uint32_t const * L, unsigned dimRows,  unsigned nRows, unsigned nColsBlock,
             uint32_t       * M, unsigned dimRowsM,                 unsigned nColsBlockM ) ;

  template
  void
  LinvApply( uint64_t       * W, unsigned Wsize,
             uint64_t const * L, unsigned dimRows,  unsigned nRows, unsigned nColsBlock,
             uint64_t       * M, unsigned dimRowsM,                 unsigned nColsBlockM ) ;

  template
  void
  LinvApply( uint128_t       * W, unsigned Wsize,
             uint128_t const * L, unsigned dimRows,  unsigned nRows, unsigned nColsBlock,
             uint128_t       * M, unsigned dimRowsM,                 unsigned nColsBlockM ) ;


}

// EOF GF2toolkit_MatrixLinvApply.cc

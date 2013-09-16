/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2010                                                      |
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
 |      version: 0.2 21-01-2010                                             |
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

#ifndef GF2TOOLKIT_MATRIX_MULT_HH
#define GF2TOOLKIT_MATRIX_MULT_HH

namespace GF2toolkit {

  /*
  //  +---+ +---+    
  //  |   | | G |    
  //  |   | +---+
  //  |   |
  //  |   |
  //  |   |
  //  | A |
  //  |   |
  //  |   |
  //  |   |
  //  |   |
  //  |   |
  //  |   |
  //  +---+
  //
  //  moltiplico matrice colonna per matrice quadrata
  */

  //! A = A*G
  template <typename UNSIGNED>
  void
  MMassLeft( UNSIGNED * A, unsigned nRowsA, UNSIGNED const * G ) ;

  //! B = A*G
  template <typename UNSIGNED>
  void
  MMass( UNSIGNED const * A, unsigned nRowsA, UNSIGNED const * G, UNSIGNED * B ) ;

  //! B ^= A*G
  template <typename UNSIGNED>
  void
  MMadd( UNSIGNED const * A, unsigned nRowsA, UNSIGNED const * G, UNSIGNED * B ) ;
  
  /*
  //  +---+ +--------------------+    
  //  | G | |         A          |
  //  +---+ +--------------------+
  //
  //  moltiplico matrice quadrata per matrice riga
  */
  //! A = G*A
  template <typename UNSIGNED>
  void
  MMassRight( UNSIGNED const * G, UNSIGNED * A ) ;

  //! A = G*A
  template <typename UNSIGNED>
  void
  MMassRight( UNSIGNED const * G, UNSIGNED * A, unsigned dimRowsA, unsigned nColsBlockA ) ;

  /* 
  //                   _        _                           _        _      
  //   _ __ ___   __ _| |_ _ __(_)_  __     _ __ ___   __ _| |_ _ __(_)_  __
  //  | '_ ` _ \ / _` | __| '__| \ \/ /____| '_ ` _ \ / _` | __| '__| \ \/ /
  //  | | | | | | (_| | |_| |  | |>  <_____| | | | | | (_| | |_| |  | |>  < 
  //  |_| |_| |_|\__,_|\__|_|  |_/_/\_\    |_| |_| |_|\__,_|\__|_|  |_/_/\_\
  */

  //! C = A*B
  template <typename UNSIGNED>
  void
  MMass( UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         UNSIGNED       * C, unsigned dimRowsC ) ;


  //! C += A*B
  template <typename UNSIGNED>
  void
  MMadd( UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         UNSIGNED       * C, unsigned dimRowsC ) ;

  //////////////////////////////////////////////////////////////////////////////

  //! C += A*B
  template <typename UNSIGNED>
  void
  MMaddRecurr( UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
               UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
               UNSIGNED       * C, unsigned dimRowsC ) ;

  //////////////////////////////////////////////////////////////////////////////

  template <typename UNSIGNED>
  void
  UinvApply( UNSIGNED const * U, unsigned dimRowsU, unsigned nbStart, unsigned nColsBlockU,
             UNSIGNED       * A, unsigned dimRowsA, unsigned nRowsA ) ;

}

#endif

// EOF GF2toolkit_MatrixMult.hh

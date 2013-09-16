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

#ifndef GF2TOOLKIT_MATRIX_LOWER_HH
#define GF2TOOLKIT_MATRIX_LOWER_HH

#include "GF2toolkit_common.hh"
#include "GF2toolkit_4Russian.hh"

namespace GF2toolkit {

  template <typename UNSIGNED>
  void
  Linverse ( UNSIGNED const * L, UNSIGNED * Linv ) ;

  template <typename UNSIGNED>
  void
  Linverse ( UNSIGNED const * L, UNSIGNED * Linv, unsigned nrows ) ;

  template <typename UNSIGNED>
  inline
  void
  Lapply ( UNSIGNED const * L, UNSIGNED * B, unsigned nrows ) {
    MM4R<UNSIGNED,4> mm4r ; // ottimo con 4
    mm4r.makeTable( B ) ;
    mm4r.multRightAss( L, nrows, B ) ;
  }

  //////////////////////////////////////////////////////////////////////////////

  template <typename UNSIGNED>
  void
  LinvApply( UNSIGNED const * L, unsigned dimRows, unsigned nRows, unsigned nColsBlock,
             UNSIGNED       * M ) ;

  template <typename UNSIGNED>
  void
  LinvApply( UNSIGNED       * W, unsigned Wsize, // work for Strassen 
             UNSIGNED const * L, unsigned dimRows,  unsigned nRows, unsigned nColsBlock,
             UNSIGNED       * M, unsigned dimRowsM, unsigned nColsBlockM ) ;

}

#endif

// EOF GF2toolkit_MatrixLower.hh

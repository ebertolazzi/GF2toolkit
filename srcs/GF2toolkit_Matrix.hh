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

#ifndef GF2TOOLKIT_MATRIX_HH
#define GF2TOOLKIT_MATRIX_HH

#include "GF2toolkit_common.hh"

namespace GF2toolkit {

  using namespace std ;

  /*
  //   ____
  //  / ___| _ __  _   _
  //  \___ \| '_ \| | | |
  //   ___) | |_) | |_| |
  //  |____/| .__/ \__, |
  //        |_|    |___/
  */
  template <typename UNSIGNED>
  void
  spy ( char const fname[],
        UNSIGNED const * A, unsigned dimRowsA,
        unsigned nRows,
        unsigned nColsBlock ) ;
        
  /*
  //   _     _____              
  //  (_)___|__  /___ _ __ ___  
  //  | / __| / // _ \ '__/ _ \ 
  //  | \__ \/ /|  __/ | | (_) |
  //  |_|___/____\___|_|  \___/ 
  */

  //! check if A is 0
  template <typename UNSIGNED>
  bool
  isZero( UNSIGNED * A, unsigned nRows ) ;

  template <typename UNSIGNED>
  bool
  isZero ( UNSIGNED * A, unsigned dimRowsA,
           unsigned nRows,
           unsigned nColsBlock ) ;

  /*
  //   _______ _ __ ___  
  //  |_  / _ \ '__/ _ \ 
  //   / /  __/ | | (_) |
  //  /___\___|_|  \___/ 
  */                   
  //! A = 0
  template <typename UNSIGNED>
  void
  zero ( UNSIGNED * A, unsigned nRows ) ;

  template <typename UNSIGNED>
  void
  zero ( UNSIGNED * A, unsigned dimRowsA,
         unsigned nRows,
         unsigned nColsBlock ) ;

  /*
  //    ___ ___  _ __  _   _ 
  //   / __/ _ \| '_ \| | | |
  //  | (_| (_) | |_) | |_| |
  //   \___\___/| .__/ \__, |
  //            |_|    |___/
  */
  //! B = A
  template <typename UNSIGNED>
  void
  copy ( UNSIGNED const * From, unsigned nRows, UNSIGNED * To ) ;

  template <typename UNSIGNED>
  void
  copy ( UNSIGNED const * From, unsigned dimRowsA,
         UNSIGNED       * To,   unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) ;
  /*
  //                         ____  _     _  __ _   
  //    ___ ___  _ __  _   _/ ___|| |__ (_)/ _| |_ 
  //   / __/ _ \| '_ \| | | \___ \| '_ \| | |_| __|
  //  | (_| (_) | |_) | |_| |___) | | | | |  _| |_ 
  //   \___\___/| .__/ \__, |____/|_| |_|_|_|  \__|
  //            |_|    |___/                       
  */
  template <typename UNSIGNED>
  void
  copyLeftShift( UNSIGNED const * From, unsigned nrows, UNSIGNED * To, unsigned lshift ) ;

  template <typename UNSIGNED>
  void
  copyRightShift( UNSIGNED const * From, unsigned nrows, UNSIGNED * To, unsigned rshift ) ;

  /*
  //             _     _ 
  //    __ _  __| | __| |
  //   / _` |/ _` |/ _` |
  //  | (_| | (_| | (_| |
  //   \__,_|\__,_|\__,_|
  */
  //! C = A+B
  template <typename UNSIGNED>
  void
  add ( UNSIGNED const * A, UNSIGNED const * B, UNSIGNED * C, unsigned nRows ) ;

  template <typename UNSIGNED>
  void
  add ( UNSIGNED const * A, unsigned dimRowsA,
        UNSIGNED const * B, unsigned dimRowsB,
        UNSIGNED       * C, unsigned dimRowsC,
        unsigned nRows,
        unsigned nColsBlock ) ;

  /*
  //             _     _ _____     
  //    __ _  __| | __| |_   _|__  
  //   / _` |/ _` |/ _` | | |/ _ \ 
  //  | (_| | (_| | (_| | | | (_) |
  //   \__,_|\__,_|\__,_| |_|\___/
  */
  //! B += A
  template <typename UNSIGNED>
  void
  addTo ( UNSIGNED const * A, UNSIGNED * B, unsigned nRows ) ;

  template <typename UNSIGNED>
  void
  addTo ( UNSIGNED const * A, unsigned dimRowsA,
          UNSIGNED       * B, unsigned dimRowsB,
          unsigned nRows,
          unsigned nColsBlock ) ;

  /*
  //                                   _       
  //   _ __   ___ _ __ _ __ ___  _   _| |_ ___ 
  //  | '_ \ / _ \ '__| '_ ` _ \| | | | __/ _ \
  //  | |_) |  __/ |  | | | | | | |_| | ||  __/
  //  | .__/ \___|_|  |_| |_| |_|\__,_|\__\___|
  //  |_|
  */
  //! swap row i with row j
  template <typename UNSIGNED>
  void
  permute ( UNSIGNED * A, unsigned nRows, int const iperm[] ) ;

  template <typename UNSIGNED>
  void
  permute ( UNSIGNED * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock,
            int32_t const iperm[] ) ;

}

#endif

// EOF GF2toolkit_Matrix.hh

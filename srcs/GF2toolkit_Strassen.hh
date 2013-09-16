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

#ifndef GF2TOOLKIT_STRASSEN_HH
#define GF2TOOLKIT_STRASSEN_HH

#include "GF2toolkit_common.hh"

namespace GF2toolkit {

  /*
  //   __  __ __  __     _
  //  |  \/  |  \/  |___| |_ _ __ __ _ ___ ___  ___ _ __
  //  | |\/| | |\/| / __| __| '__/ _` / __/ __|/ _ \ '_ \
  //  | |  | | |  | \__ \ |_| | | (_| \__ \__ \  __/ | | |
  //  |_|  |_|_|  |_|___/\__|_|  \__,_|___/___/\___|_| |_|
  */
  template <typename UNSIGNED>
  void
  MMaddStrassen( UNSIGNED * work, unsigned nwork,
                 unsigned nBlock,
                 UNSIGNED const * Amat, unsigned dimRowA,
                 UNSIGNED const * Bmat, unsigned dimRowB,
                 UNSIGNED       * Cmat, unsigned dimRowC ) ;

  //! C += A*B
  template <typename UNSIGNED>
  void
  MMaddStrassen( UNSIGNED * work, unsigned nwork,
                 UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
                 UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
                 UNSIGNED       * C, unsigned dimRowsC ) ;

}

#endif

// EOF GF2toolkit_Strassen.hh

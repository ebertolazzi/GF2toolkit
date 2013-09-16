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

#include "GF2toolkit_m4ri.hh"

namespace GF2toolkit {

  using namespace std ;

  template <>
  mzd_t *
  toM4RI( uint32_t const * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {
    ASSERT( (nColsBlock & 0x01) == 0, "nColsBlock = " << nColsBlock << " must be even" ) ;
    unsigned nCols = nColsBlock*sizeof(uint32_t)*CHAR_BIT ;
    mzd_t * mat = mzd_init( nRows, nCols );
    for ( unsigned i = 0 ; i < nRows ; ++i )
      for ( unsigned jb = 0 ; jb < nColsBlock ; jb += 2 ) {
        uint64_t lo = A[i+jb*dimRowsA]     & 0xFFFFFFFFU ;
        uint64_t hi = A[i+(jb+1)*dimRowsA] & 0xFFFFFFFFU ;
        mat -> rows[i][jb>>1] = lo | (hi<<32) ;
      }
    return mat ;
  }

  template <>
  void
  fromM4RI( mzd_t * mat, uint32_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {
    unsigned nCols = nColsBlock*sizeof(uint32_t)*CHAR_BIT ;
    ASSERT( (unsigned)mat -> nrows <= nRows, "mat -> nrows = " << mat -> nrows << " greather that " << nRows ) ;
    ASSERT( (unsigned)mat -> ncols <= nCols, "mat -> ncols = " << mat -> ncols << " greather that " << nCols ) ;
    unsigned nblk = (mat -> ncols)>>6 ;
    ASSERT( nColsBlock >= 2*nblk, "fromM4RI, not enoght block in out matrix" ) ;
    for ( unsigned i = 0 ; i < mat -> nrows ; ++i )
      for ( unsigned jb = 0 ; jb < nblk ; ++jb ) {
        uint64_t tmp = mat -> rows[i][jb] ;
        A[2*jb*dimRowsA+i]     = (uint32_t)(tmp) ;
        A[(2*jb+1)*dimRowsA+i] = (uint32_t)(tmp>>32);
      }
  }

  // ---------------------------------------------------------------------------

  template <>
  mzd_t *
  toM4RI( uint64_t const * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {
    unsigned nCols = nColsBlock*64 ;
    mzd_t * mat = mzd_init( nRows, nCols );
    for ( unsigned i = 0 ; i < nRows ; ++i )
      for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb )
        mat -> rows[i][jb] = A[jb*dimRowsA+i] ;
    return mat ;
  }

  template <>
  void
  fromM4RI( mzd_t * mat, uint64_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {
    unsigned nCols = nColsBlock*sizeof(uint64_t)*CHAR_BIT ;
    ASSERT( (unsigned)mat -> nrows <= nRows, "mat -> nrows = " << mat -> nrows << " greather that " << nRows ) ;
    ASSERT( (unsigned)mat -> ncols <= nCols, "mat -> ncols = " << mat -> ncols << " greather that " << nCols ) ;
    unsigned nblk = mat -> ncols>>6 ;
    for ( unsigned i = 0 ; i < mat -> nrows ; ++i )
      for ( unsigned jb = 0 ; jb < nblk ; ++jb )
        A[jb*dimRowsA+i] = mat -> rows[i][jb] ;
  }

  // ---------------------------------------------------------------------------

  template <>
  mzd_t *
  toM4RI( uint128_t const * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {
    unsigned nCols = nColsBlock*sizeof(uint128_t)*CHAR_BIT ;
    mzd_t * mat = mzd_init( nRows, nCols );
    for ( unsigned i = 0 ; i < nRows ; ++i )
      for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb ) {
        uint128_t tmp = A[jb*dimRowsA+i] ;
        mat -> rows[i][2*jb+0] = tmp.i64.lo ;
        mat -> rows[i][2*jb+1] = tmp.i64.hi ;
      }
    return mat ;
  }

  template <>
  void
  fromM4RI( mzd_t * mat, uint128_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {
    unsigned nCols = nColsBlock*sizeof(uint128_t)*CHAR_BIT ;
    ASSERT( (unsigned)mat -> nrows <= nRows, "mat -> nrows = " << mat -> nrows << " greather that " << nRows ) ;
    ASSERT( (unsigned)mat -> ncols <= nCols, "mat -> ncols = " << mat -> ncols << " greather that " << nCols ) ;
    ASSERT( (mat -> nrows & 0x7F) == 0, "mat -> nrows = " << mat -> nrows << " must be a multiple of 128" ) ;
    unsigned nblk = mat -> ncols>>7 ;
    for ( unsigned i = 0 ; i < mat -> nrows ; ++i )
      for ( unsigned jb = 0 ; jb < nblk ; ++jb ) {
        A[jb*dimRowsA+i].i64.lo = mat -> rows[i][2*jb+0] ;
        A[jb*dimRowsA+i].i64.hi = mat -> rows[i][2*jb+1] ;
      }
  }

}

// EOF GF2toolkit_m4ri.cc

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
#include <fstream>

namespace GF2toolkit {

  using namespace std ;

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
  isZero( UNSIGNED * A, unsigned nRows ) {
    bool ok = true ;
    GF2_UNROLL_BY_16( nRows, ok = *A++ == 0 ; if ( !ok ) break ) ;
    return ok ;
  }

  template <>
  bool
  isZero( uint128_t * A, unsigned nRows ) {
    bool ok = true ;
    GF2_UNROLL_BY_16( nRows, 
                      ok = (A->i64.lo == 0 && A->i64.hi == 0) ;
                      ++A ;
                      if ( !ok ) break ) ;
    return ok ;
  }

  template bool isZero( uint32_t  * A, unsigned nRows ) ;
  template bool isZero( uint64_t  * A, unsigned nRows ) ;

  template <typename UNSIGNED>
  bool
  isZero ( UNSIGNED * A, unsigned dimRowsA,
           unsigned nRows,
           unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return true ;

    ASSERT( dimRowsA >= nRows, "GF2::zero nRows = " << nRows << " dimRows = " << dimRowsA ) ;

    bool ok = true ;
    for ( unsigned jb = 0 ; jb < nColsBlock && ok ; ++jb, A += dimRowsA )
      ok = isZero( A, nRows ) ;

    return ok ;
  }

  template bool isZero( uint32_t  * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) ;
  template bool isZero( uint64_t  * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) ;
  template bool isZero( uint128_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) ;

  /*
  //   _______ _ __ ___
  //  |_  / _ \ '__/ _ \
  //   / /  __/ | | (_) |
  //  /___\___|_|  \___/
  */
  template <typename UNSIGNED>
  void
  zero( UNSIGNED * A, unsigned nRows ) 
  { GF2_UNROLL_BY_16( nRows, *A++ = UNSIGNED(0) ) ; }

  template <>
  void
  zero( uint128_t * A, unsigned nRows ) {
    GF2_UNROLL_BY_16( nRows, SSE_ZERO( A ) ; ++A ) ;
  }

  template void zero( uint32_t * A, unsigned nRows ) ;
  template void zero( uint64_t * A, unsigned nRows ) ;

  template <typename UNSIGNED>
  void
  zero( UNSIGNED * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::zero nRows = " << nRows << " dimRows = " << dimRowsA ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA )
      zero( A, nRows ) ;
  }

  template <>
  void
  zero( uint128_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::zero nRows = " << nRows << " dimRows = " << dimRowsA ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA )
      zero( A, nRows ) ;
  }
  
  template void zero( uint32_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) ;
  template void zero( uint64_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock ) ;

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
  copy ( UNSIGNED const * A, unsigned nRows, UNSIGNED * B ) 
  { GF2_UNROLL_BY_16( nRows, *B++ = *A++ ) ; }
  
  template <>
  void
  copy ( uint128_t const * A, unsigned nRows, uint128_t * B ) 
  { GF2_UNROLL_BY_16( nRows, B++->i128 = A++->i128 ) ; }

  template void copy ( uint32_t const * A, unsigned nRow, uint32_t * Bs ) ;
  template void copy ( uint64_t const * A, unsigned nRows, uint64_t * B ) ;

  //! B = A
  template <typename UNSIGNED>
  void
  copy ( UNSIGNED const * A, unsigned dimRowsA,
         UNSIGNED       * B, unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::copy nRows = " << nRows << " dimRowsA = " << dimRowsA ) ;
    ASSERT( dimRowsB >= nRows, "GF2::copy nRows = " << nRows << " dimRowsB = " << dimRowsB ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA, B += dimRowsB )
      copy<UNSIGNED>( A, nRows, B ) ;
  }

  template <>
  void
  copy ( uint128_t const * A, unsigned dimRowsA,
         uint128_t       * B, unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::copy nRows = " << nRows << " dimRowsA = " << dimRowsA ) ;
    ASSERT( dimRowsB >= nRows, "GF2::copy nRows = " << nRows << " dimRowsB = " << dimRowsB ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA, B += dimRowsB )
      copy<uint128_t>( A, nRows, B ) ;
  }
  
  template
  void
  copy ( uint32_t const * A, unsigned dimRowsA,
         uint32_t       * B, unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) ;

  template
  void
  copy ( uint64_t const * A, unsigned dimRowsA,
         uint64_t       * B, unsigned dimRowsB,
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
  template <>
  void
  copyLeftShift( uint32_t const * From, unsigned nrows, uint32_t * To, unsigned lshift ) {
    GF2_UNROLL_BY_16( nrows, *To++ = (*From++)<<lshift ) ;
  }

  template <>
  void
  copyLeftShift( uint64_t const * From, unsigned nrows, uint64_t * To, unsigned lshift ) {
    GF2_UNROLL_BY_16( nrows, *To++ = (*From++)<<lshift ) ;
  }

  template <>
  void
  copyLeftShift( uint128_t const * From, unsigned nrows, uint128_t * To, unsigned lshift  ) {
    if ( lshift == 64 ) {
      GF2_UNROLL_BY_16( nrows, To->i64.lo = 0 ; To++->i64.hi = From++->i64.lo ) ;
    } else if ( lshift > 64 ) {
      GF2_UNROLL_BY_16( nrows, To->i64.lo = 0 ; To++->i64.hi = ((From++->i64.lo)<<(lshift-64)) ) ;
    } else { // n < 64
      GF2_UNROLL_BY_16( nrows, To->i64.hi   = ((From->i64.lo)>>(64-lshift)) | ((From->i64.hi)<<lshift) ;
                               To++->i64.lo = ((From++->i64.lo)<<lshift) ) ;
    }
  }

  template <>
  void
  copyRightShift( uint32_t const * From, unsigned nrows, uint32_t * To, unsigned rshift ) {
    GF2_UNROLL_BY_16( nrows, *To++ = (*From++)>>rshift ) ;  
  }

  template <>
  void
  copyRightShift( uint64_t const * From, unsigned nrows, uint64_t * To, unsigned rshift ) {
    GF2_UNROLL_BY_16( nrows, *To++ = (*From++)>>rshift ) ;  
  }

  template <>
  void
  copyRightShift( uint128_t const * From, unsigned nrows, uint128_t * To, unsigned rshift ) {
    if ( rshift == 64 ) {
      GF2_UNROLL_BY_16( nrows, To->i64.hi = 0 ; To++->i64.lo = From++->i64.hi ) ;
    } else if ( rshift > 64 ) {
      GF2_UNROLL_BY_16( nrows, To->i64.hi = 0 ; To++->i64.lo = ((From++->i64.hi)>>(rshift-64)) ) ;
    } else { // n < 64
      GF2_UNROLL_BY_16( nrows, To->i64.lo   = ((From->i64.lo)>>rshift) | ((From->i64.hi)<<(64-rshift)) ;
                               To++->i64.hi = ((From++->i64.hi)>>rshift) ) ;
    }
  }

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
  add( UNSIGNED const * A,
       UNSIGNED const * B,
       UNSIGNED       * C,
       unsigned         nRows ) {
    GF2_UNROLL_BY_16( nRows, *C++ = *A++ ^ *B++ ) ;
  }
  
  template
  void
  add( uint32_t const * A,
       uint32_t const * B,
       uint32_t       * C,
       unsigned         nRows ) ;

  template
  void
  add( uint64_t const * A,
       uint64_t const * B,
       uint64_t       * C,
       unsigned         nRows ) ;

  template <>
  void
  add( uint128_t const * A,
       uint128_t const * B,
       uint128_t       * C,
       unsigned          nRows ) {
    GF2_UNROLL_BY_16( nRows, SSE_XOR( C, A, B ) ; ++C ; ++A ; ++B ) ;
  }

  //! C = A+B
  template <typename UNSIGNED>
  void
  add( UNSIGNED const * A, unsigned dimRowsA,
       UNSIGNED const * B, unsigned dimRowsB,
       UNSIGNED       * C, unsigned dimRowsC,
       unsigned nRows,
       unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::add nRows = " << nRows << " dimRowsA = " << dimRowsA ) ;
    ASSERT( dimRowsB >= nRows, "GF2::add nRows = " << nRows << " dimRowsB = " << dimRowsB ) ;
    ASSERT( dimRowsC >= nRows, "GF2::add nRows = " << nRows << " dimRowsC = " << dimRowsC ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA, B += dimRowsB, C += dimRowsC ) {
      UNSIGNED const * pA = A ;
      UNSIGNED const * pB = B ;
      UNSIGNED       * pC = C ;
      GF2_UNROLL_BY_16( nRows, *pC++ = *pA++ ^ *pB++ ) ;
    }
  }

  template <>
  void
  add( uint128_t const * A, unsigned dimRowsA,
       uint128_t const * B, unsigned dimRowsB,
       uint128_t       * C, unsigned dimRowsC,
       unsigned nRows,
       unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::add nRows = " << nRows << " dimRowsA = " << dimRowsA ) ;
    ASSERT( dimRowsB >= nRows, "GF2::add nRows = " << nRows << " dimRowsB = " << dimRowsB ) ;
    ASSERT( dimRowsC >= nRows, "GF2::add nRows = " << nRows << " dimRowsC = " << dimRowsC ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA, B += dimRowsB, C += dimRowsC ) {
      uint128_t const * pA = A ;
      uint128_t const * pB = B ;
      uint128_t       * pC = C ;
      GF2_UNROLL_BY_16( nRows, SSE_XOR( pC, pA, pB ) ; ++pC ; ++pA ; ++pB ) ;
      //GF2_UNROLL_BY_16( nRows, pC->i64.lo = pA->i64.lo^pB->i64.lo ;
      //                         pC++->i64.hi = pA++->i64.hi^pB++->i64.hi ) ;
    }
  }

  template
  void
  add( uint32_t const * A, unsigned dimRowsA,
       uint32_t const * B, unsigned dimRowsB,
       uint32_t       * C, unsigned dimRowsC,
       unsigned nRows,
       unsigned nColsBlock ) ;

  template
  void
  add( uint64_t const * A, unsigned dimRowsA,
       uint64_t const * B, unsigned dimRowsB,
       uint64_t       * C, unsigned dimRowsC,
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
  addTo( UNSIGNED const * A,
         UNSIGNED       * B,
         unsigned         nRows ) {
    GF2_UNROLL_BY_16( nRows, *B++ ^= *A++ ) ;
  }

  template
  void
  addTo( uint32_t const * A,
         uint32_t       * B,
         unsigned         nRows ) ;
  
  template
  void
  addTo( uint64_t const * A,
         uint64_t       * B,
         unsigned         nRows ) ;

  template <>
  void
  addTo( uint128_t const * A,
         uint128_t       * B,
         unsigned          nRows ) {
    GF2_UNROLL_BY_16( nRows, SSE_ASS_XOR( B, A ) ; ++A ; ++B ) ;
  }

  //! B += A
  template <typename UNSIGNED>
  void
  addTo( UNSIGNED const * A, unsigned dimRowsA,
         UNSIGNED       * B, unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::add nRows = " << nRows << " dimRowsA = " << dimRowsA ) ;
    ASSERT( dimRowsB >= nRows, "GF2::add nRows = " << nRows << " dimRowsB = " << dimRowsB ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA, B += dimRowsB ) {
      UNSIGNED const * pA = A ;
      UNSIGNED       * pB = B ;
      GF2_UNROLL_BY_16( nRows, *pB++ ^= *pA++ ) ;
    }
  }
  
  template <>
  void
  addTo( uint128_t const * A, unsigned dimRowsA,
         uint128_t       * B, unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) {

    if ( nRows == 0 || nColsBlock == 0 ) return ;

    ASSERT( dimRowsA >= nRows, "GF2::add nRows = " << nRows << " dimRowsA = " << dimRowsA ) ;
    ASSERT( dimRowsB >= nRows, "GF2::add nRows = " << nRows << " dimRowsB = " << dimRowsB ) ;

    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA, B += dimRowsB ) {
      uint128_t const * pA = A ;
      uint128_t       * pB = B ;
      GF2_UNROLL_BY_16( nRows, SSE_ASS_XOR( pB, pA ) ; ++pA ; ++pB ) ;
    }
  }

  template
  void
  addTo( uint32_t const * A, unsigned dimRowsA,
         uint32_t       * B, unsigned dimRowsB,
         unsigned nRows,
         unsigned nColsBlock ) ;

  template
  void
  addTo( uint64_t const * A, unsigned dimRowsA,
         uint64_t       * B, unsigned dimRowsB,
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
  //! apply permutatuion
  template <typename UNSIGNED>
  void
  permute ( UNSIGNED * A, unsigned nRows, int32_t const iperm[] ) {
    int32_t j = 0 ;
    for ( unsigned i = 0 ; i < nRows ; ++i ) {
      int32_t ip = iperm[i] ;
      if ( ip < 0  ) continue ;
      if ( ip != j ) std::swap( A[j], A[ip] ) ; 
      ++j ;
    }
  }

  template void permute ( uint32_t  * A, unsigned nRows, int32_t const iperm[] ) ;
  template void permute ( uint64_t  * A, unsigned nRows, int32_t const iperm[] ) ;
  template void permute ( uint128_t * A, unsigned nRows, int32_t const iperm[] ) ;

  template <typename UNSIGNED>
  void
  permute ( UNSIGNED * A, unsigned dimRowsA,
            unsigned nRows, unsigned nColsBlock,
            int32_t const iperm[] ) {
    for ( unsigned jb = 0 ; jb < nColsBlock ; ++jb, A += dimRowsA )
      permute<UNSIGNED>( A, nRows, iperm ) ;
  }

  template void permute ( uint32_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock, int32_t const iperm[] ) ;
  template void permute ( uint64_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock, int32_t const iperm[] ) ;
  template void permute ( uint128_t * A, unsigned dimRowsA, unsigned nRows, unsigned nColsBlock, int32_t const iperm[] ) ;

  /*
  //   ____
  //  / ___| _ __  _   _
  //  \___ \| '_ \| | | |
  //   ___) | |_) | |_| |
  //  |____/| .__/ \__, |
  //        |_|    |___/
  */
  template <typename UNSIGNED>
  bool
  isSet( UNSIGNED const aa, unsigned nbit ) {
    return (aa >> nbit) & UNSIGNED(1) ;
  }

  template <>
  bool
  isSet( uint128_t const aa, unsigned nbit ) {
    if ( nbit < 64 ) return (aa.i64.lo >> nbit) & uint64_t(1) ;
    else             return (aa.i64.hi >> (nbit-64)) & uint64_t(1) ;
  }

  template <typename UNSIGNED>
  void
  spy ( char const fname[],
        UNSIGNED const * A, unsigned dimRowsA,
        unsigned numRows,
        unsigned numColsBlock ) {
    /*
     // based on:
     //
     // PSPLTM - PostScript PLoTer of a (sparse) Matrix
     // By Loris Renggli (renggli@masg1.epfl.ch) and Youcef Saad
     */

    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;
    unsigned const numCols = numColsBlock * numBits ;

    // units (cm or in) to dot conversion factor and paper size
    double u2dot  = 72.0/2.54 ;
    double paperx = 21.0 ;
    double xsize  = 18.0 ;
    double m      = std::max( numRows, numCols ) ;

    // left and right margins (drawing is centered)
    double leftRightMargin = (paperx-xsize)/2.0 ;

    // bottom margin : 2 cm
    double bottomMargin = 2.0 ;

    // scaling factor
    double scalingFactor = xsize/m ;

    // matrix frame line witdh
    double sepLaneWidth = 0.01 ;

    // matrix frame line witdh
    double frameLaneWidth  = std::max(sepLaneWidth,0.001*m) ;
    double frameLaneWidth2 = 4*frameLaneWidth ;

    // almost exact bounding box
    double wd = scalingFactor*frameLaneWidth/2 ;
    double xl = leftRightMargin                      - wd ;
    double xr = (leftRightMargin+xsize)              + wd ;
    double yb = bottomMargin                         - wd ;
    double yt = (bottomMargin+scalingFactor*numRows) + wd ;

    // add some room to bounding box
    double delt = 15.0 ;
    xl = xl * u2dot - delt ;
    xr = xr * u2dot + delt ;
    yb = yb * u2dot - delt ;
    yt = yt * u2dot + delt ;

    ofstream fd( fname ) ;

    double xmin = 0.5-frameLaneWidth2/2-frameLaneWidth ;
    double xmax = numCols+0.5+frameLaneWidth2/2+frameLaneWidth ;
    double ymin = 0.5-frameLaneWidth2/2-frameLaneWidth ;
    double ymax = numRows+0.5+frameLaneWidth2/2+frameLaneWidth ;

    // begin of output
    fd << "%!PS-Adobe-3.0 EPSF-3.0\n"
       << "%%Creator: SPARSELIB by Enrico Bertolazzi\n"
       << "%%BoundingBox: " << int(xl+0.5) << ' '
       << int(yb+0.5) << ' '
       << int(xr+0.5) << ' '
       << int(yt+0.5) << '\n'
       << "%%EndComments\n"
       << "/cm {72 mul 2.54 div} def\n"
       << "gsave\n"
       << leftRightMargin << " cm "
       << bottomMargin << " cm translate\n"
       << xsize << " cm " << m << " div dup scale\n"
       // draw a frame around the matrix
       << frameLaneWidth2 << " setlinewidth\n"
       << "newpath\n"
       << xmin << ' ' << ymin << " moveto\n"
       << xmax << ' ' << ymin << " lineto\n"
       << xmax << ' ' << ymax << " lineto\n"
       << xmin << ' ' << ymax << " lineto\n"
       << "closepath stroke\n" ;

    // drawing the separation lines (if required)

    unsigned nsep = numBits ;
    fd << sepLaneWidth << " setlinewidth\n" ;
    for ( unsigned i = nsep ; i < numRows ; i += nsep ) {
      double y = numRows - i + 0.5 ;
      fd << xmin << ' ' << y << " moveto " << xmax << ' ' << y << " lineto stroke\n" ;
    }
    fd << sepLaneWidth << " setlinewidth\n" ;
    for ( unsigned i = nsep ; i < numCols ; i += nsep ) {
      double x = i+0.5 ;
      fd << x << ' ' << ymin << " moveto " << x << ' ' << ymax << " lineto stroke\n" ;
    }

    // ----------- plotting loop
    fd << "1 1 translate\n"
       << "0.8 setlinewidth\n"
       << "/p {moveto 0 -.40 rmoveto 0 .80 rlineto stroke} def\n" ;

    for ( unsigned colBock = 0 ; colBock < numColsBlock ; ++colBock ) {
      for ( unsigned row = 0 ; row < numRows ; ++row ) {
        UNSIGNED aa = A[dimRowsA*colBock+row] ;
        for ( unsigned nbit = 0 ; nbit < numBits ; ++nbit )
          if ( isSet(aa,nbit) )
            fd << colBock*numBits+nbit << ' ' << (numRows - row - 1) << " p\n" ;
      }
    }

    fd << "showpage\n" ;
    fd . close() ;
  }

  template
  void
  spy ( char const fname[],
        uint32_t const * A, unsigned dimRowsA,
        unsigned numRows,
        unsigned numColsBlock ) ;

  template
  void
  spy ( char const fname[],
        uint64_t const * A, unsigned dimRowsA,
        unsigned numRows,
        unsigned numColsBlock ) ;

  template
  void
  spy ( char const fname[],
        uint128_t const * A, unsigned dimRowsA,
        unsigned numRows,
        unsigned numColsBlock ) ;


}

// EOF GF2toolkit_Matrix.cc

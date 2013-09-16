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

#include "GF2toolkit_common.hh"
#include "GF2toolkit_LU.hh"
#include "GF2toolkit_Matrix.hh"

#include <iostream>
#include <iomanip>

namespace GF2toolkit {

  #define POS(r,jb,dim) ((r)+(jb)*(dim))

  /*
  //             _                  _      _    
  //    _____  _| |_ _ __ __ _  ___| |_   / \   
  //   / _ \ \/ / __| '__/ _` |/ __| __| / _ \  
  //  |  __/>  <| |_| | | (_| | (__| |_ / ___ \ 
  //   \___/_/\_\\__|_|  \__,_|\___|\__/_/   \_\
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::extractA( UNSIGNED * A, unsigned dimRowsA, unsigned numBlocksA,
                          unsigned & nr, unsigned & nc ) {

    ASSERT( dimRowsA >= numRows,
            "not enought rows, required = " << numRows <<
            " availabes = " << dimRowsA ) ;

    ASSERT( numBlocksA >= numColsBlock,
            "not enought column block, required = " << numColsBlock <<
            " availabes = " << numBlocksA ) ;

    nr = numRows ;
    nc = numColsBlock<<numShift ;

    copy<UNSIGNED>( &block(0,0), dimRows, A, dimRowsA, numRows, numColsBlock ) ;

  }

  /*
  //             _                  _   _     
  //    _____  _| |_ _ __ __ _  ___| |_| |    
  //   / _ \ \/ / __| '__/ _` |/ __| __| |    
  //  |  __/>  <| |_| | | (_| | (__| |_| |___ 
  //   \___/_/\_\\__|_|  \__,_|\___|\__|_____|
  //  
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::extractL( UNSIGNED * L, unsigned dimRowsL, unsigned numBlocksL,
                          unsigned & nr, unsigned & nc ) {

    ASSERT( (numBlocksL<<numShift) >= rowStart[numColsBlock],
            "not enought column, required = " << rowStart[numColsBlock] <<
            " availabes = " << (numBlocksL<<numShift) ) ;

    ASSERT( numRows <= dimRowsL,
            "not enought rows, required = " << numRows <<
            " availabes = " << dimRowsL ) ;

    nr = numRows ;
    nc = rowStart[numColsBlock] ;

    GF2toolkit::zero<UNSIGNED>( L, dimRowsL, dimRowsL, numBlocksL ) ;
    unsigned nblk = (nc+numBits-1)>>numShift ;
    for ( unsigned jb = 0 ; jb < nblk ; ++jb ) {
      unsigned nrows = jb*numBits ;
      GF2toolkit::copy<UNSIGNED>( &block(nrows,jb),
                                  numRows - nrows,
                                  L + POS(nrows,jb,dimRowsL) ) ;
    }
  }

  /*
  //             _                  _   _   _ 
  //    _____  _| |_ _ __ __ _  ___| |_| | | |
  //   / _ \ \/ / __| '__/ _` |/ __| __| | | |
  //  |  __/>  <| |_| | | (_| | (__| |_| |_| |
  //   \___/_/\_\\__|_|  \__,_|\___|\__|\___/
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::extractU( UNSIGNED * U, unsigned dimRowsU, unsigned numBlocksU,
                          unsigned & nr, unsigned & nc ) {

    ASSERT( dimRowsU >= rowStart[numColsBlock],
            "not enought rows, required = " << rowStart[numColsBlock] <<
            " availabes = " << dimRows ) ;

    ASSERT( numBlocksU >= numColsBlock,
            "not enought column block, required = " << numColsBlock <<
            " availabes = " << numBlocksU ) ;

    nr = rowStart[numColsBlock] ;
    nc = numColsBlock<<numShift ;

    GF2toolkit::zero<UNSIGNED>( U, dimRowsU, dimRowsU, numColsBlock ) ;

    for ( unsigned jb = 0 ; jb < numColsBlock ; ++jb ) {
      GF2toolkit::copy<UNSIGNED>( &block(0,jb),
                                  rowStart[jb],
                                  U + POS(0,jb,dimRowsU) ) ;

      GF2toolkit::copy<UNSIGNED>( &Ublock(jb),
                                  rowStart[jb+1]-rowStart[jb], 
                                  U + POS(rowStart[jb],jb,dimRowsU) ) ;
    }
  }

  /*
  //                     _       ____                           _       
  //    __ _ _ __  _ __ | |_   _|  _ \ ___ _ __ _ __ ___  _   _| |_ ___ 
  //   / _` | '_ \| '_ \| | | | | |_) / _ \ '__| '_ ` _ \| | | | __/ _ \
  //  | (_| | |_) | |_) | | |_| |  __/  __/ |  | | | | | | |_| | ||  __/
  //   \__,_| .__/| .__/|_|\__, |_|   \___|_|  |_| |_| |_|\__,_|\__\___|
  //        |_|   |_|      |___/
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::applyPermute( UNSIGNED * M, unsigned dimRowsM, unsigned numBlocksM ) {
  
    ASSERT( dimRowsM >= numRows,
            "not enought rows, required = " << numRows <<
            " availabes = " << dimRowsM ) ;

    ASSERT( numBlocksM >= numColsBlock,
            "not enought column block, required = " << numColsBlock <<
            " availabes = " << numBlocksM ) ;
    
    int32_t const * P = Perm ;
    for ( unsigned ib = 0 ; ib < numColsBlock ; ++ib, P += numBits )
      if ( rowStart[ib+1] > rowStart[ib] ) // ci sono righe da permutare ?
        for ( unsigned jb = 0 ; jb < numBlocksM ; ++jb )
          GF2toolkit::permute<UNSIGNED>( M + POS(rowStart[ib],jb,dimRowsM), numBits, P ) ;
  }

  #define EXPLICT_INTANTIATION(TYPE) \
  template void LU<TYPE>::extractA( TYPE * A, unsigned dimRowsA, unsigned numBlocksA, unsigned & nr, unsigned & nc ) ; \
  template void LU<TYPE>::extractL( TYPE * L, unsigned dimRowsL, unsigned numBlocksL, unsigned & nr, unsigned & nc ) ; \
  template void LU<TYPE>::extractU( TYPE * U, unsigned dimRowsU, unsigned numBlocksU, unsigned & nr, unsigned & nc ) ; \
  template void LU<TYPE>::applyPermute( TYPE * M, unsigned dimRowsM, unsigned numBlocksM )

  EXPLICT_INTANTIATION(uint32_t) ;
  EXPLICT_INTANTIATION(uint64_t) ;
  EXPLICT_INTANTIATION(uint128_t) ;

}

// EOF GF2toolkit_LUextract.cc

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

#include "GF2toolkit_LU.hh"
#include "GF2toolkit_popCounts.hh"

#include <iostream>
#include <iomanip>
#include <vector>

#include "GF2toolkit_Matrix.hh"
#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_MatrixLower.hh"
#include "GF2toolkit_4Russian.hh"
#include "GF2toolkit_Strassen.hh"

using namespace std ;

namespace GF2toolkit {

  template <> unsigned const LU<uint32_t>::numBits  = 32   ;
  template <> unsigned const LU<uint32_t>::numShift = 5    ;
  template <> unsigned const LU<uint32_t>::numMask  = 0x1F ;

  template <> unsigned const LU<uint64_t>::numBits  = 64   ;
  template <> unsigned const LU<uint64_t>::numShift = 6    ;
  template <> unsigned const LU<uint64_t>::numMask  = 0x3F ;

  template <> unsigned const LU<uint128_t>::numBits  = 128  ;
  template <> unsigned const LU<uint128_t>::numShift = 7    ;
  template <> unsigned const LU<uint128_t>::numMask  = 0x7F ;

  /*
  //   _    _   _ 
  //  | |  | | | |
  //  | |  | | | |
  //  | |__| |_| |
  //  |_____\___/ 
  */
  template <typename UNSIGNED>
  LU<UNSIGNED>::LU() 
  : dataMatrix(0)
  , dataWork(0)
  , dataUmat(0)
  , dataLmat(0)
  , sizeWork(0)
  , rowStart(0)
  , Perm(0)
  , numRows(0)
  , numColsBlock(0)
  , dimRows(0)
  {
  }

  template <typename UNSIGNED>
  LU<UNSIGNED>::LU( unsigned nr, unsigned nc ) 
  : dataMatrix(0)
  , dataWork(0)
  , sizeWork(0)
  , numRows(0)
  , numColsBlock(0)
  , dimRows(0)  
  {
    resize( nr, nc ) ;
  }

  /*
  //   /\/|    _   _ 
  //  |/\/ |  | | | |
  //     | |  | | | |
  //     | |__| |_| |
  //     |_____\___/ 
  */
  template <typename UNSIGNED>
  LU<UNSIGNED>::~LU() {
    if ( dataMatrix != 0 ) delete [] dataMatrix ; dataMatrix = 0 ;
  }

  /*                _         
  //  _ __ ___  ___(_)_______ 
  // | '__/ _ \/ __| |_  / _ \
  // | | |  __/\__ \ |/ /  __/
  // |_|  \___||___/_/___\___|
  */
  template <typename UNSIGNED>
  void 
  LU<UNSIGNED>::resize( unsigned nr, unsigned nc ) {

    numRows      = nr ;
    dimRows      = nr ;
    numColsBlock = nc>>numShift ;

    if ( nc & numMask ) ++numColsBlock ; // se non divisibile aggiunge blocco colonna

    sizeWork = numRows * numColsBlock ;

    if ( dataMatrix != 0 ) {
      delete [] dataMatrix ;
      delete [] Perm ;
      delete [] rowStart ;
    }

    dataMatrix = new UNSIGNED[ sizeWork + 
                               numColsBlock * ( dimRows + 2 * numBits ) ] ;

    dataLmat = dataMatrix + numColsBlock * dimRows ;
    dataUmat = dataLmat   + numColsBlock * numBits ;
    dataWork = dataUmat   + numColsBlock * numBits ;

    Perm     = new int32_t  [ numColsBlock * numBits ] ;
    rowStart = new unsigned [ numColsBlock + 1 ] ;

    this -> clear() ;
  }

  /*
  //  _ __  _ __  ____
  // | '_ \| '_ \|_  /
  // | | | | | | |/ / 
  // |_| |_|_| |_/___|
  */
  template <typename UNSIGNED>
  unsigned
  LU<UNSIGNED>::nnz() const {
    unsigned tot = 0 ;
    for ( unsigned cb = 0 ; cb < numColsBlock ; ++cb )
      for ( unsigned r = 0 ; r < numRows ; ++r )
        tot += popCounts( block( r, cb ) ) ;
    return tot ;
  }

  /*
  //    __ _ _ _ ____        ____  _ _       
  //   / _(_) | | __ ) _   _| __ )(_) |_ ___ 
  //  | |_| | | |  _ \| | | |  _ \| | __/ __|
  //  |  _| | | | |_) | |_| | |_) | | |_\__ \
  //  |_| |_|_|_|____/ \__, |____/|_|\__|___/
  //                   |___/                 
  */
  //! fill a matrix with user values
  template <>
  void
  LU<uint32_t>::fillByBits( uint64_t const vec[], unsigned sizeVec ) {
    unsigned nblk = numRows*numColsBlock ;
    ASSERT( 2*sizeVec == nblk, "GF2Rank::fillByBits sizeVec = " << sizeVec << " expected = " << nblk/2.0 ) ;
    uint64_t const *pvec = vec ;
    for ( unsigned r = 0 ; r < numRows ; ++r ) {
      for ( unsigned cb = 0 ; cb < numColsBlock ; cb += 2, ++pvec ) {
        block( r, cb   ) = (*pvec)&0xFFFFFFFF ;
        block( r, cb+1 ) = (*pvec)>>32 ;
      }
    }
  }

  template <>
  void
  LU<uint64_t>::fillByBits( uint64_t const vec[], unsigned sizeVec ) {
    unsigned nblk = numRows*numColsBlock ;
    ASSERT( sizeVec == nblk, "GF2Rank::fillByBits sizeVec = " << sizeVec << " expected = " << nblk ) ;
    uint64_t const *pvec = vec ;
    for ( unsigned r = 0 ; r < numRows ; ++r )
      for ( unsigned cb = 0 ; cb < numColsBlock ; ++cb )
        block( r, cb ) = *pvec++ ;
  }

  template <>
  void
  LU<uint128_t>::fillByBits( uint64_t const vec[], unsigned sizeVec ) {
    unsigned nblk = numRows*numColsBlock ;
    ASSERT( sizeVec == 2*nblk, "GF2Rank::fillByBits sizeVec = " << sizeVec << " expected = " << 2*nblk ) ;
    uint64_t const *pvec = vec ;
    for ( unsigned r = 0 ; r < numRows ; ++r ) {
      for ( unsigned cb = 0 ; cb < numColsBlock ; ++cb ) {
        uint128_t & bk = block( r, cb ) ;
        bk.i64.lo = *pvec++ ;
        bk.i64.hi = *pvec++ ;
      }
    }
  }

  /*
  //    __ _ _ _ ____        ____        _       
  //   / _(_) | | __ ) _   _| __ ) _   _| |_ ___ 
  //  | |_| | | |  _ \| | | |  _ \| | | | __/ _ \
  //  |  _| | | | |_) | |_| | |_) | |_| | ||  __/
  //  |_| |_|_|_|____/ \__, |____/ \__, |\__\___|
  //                   |___/       |___/         
  */
  //! fill a matrix with user values
  template <>
  void
  LU<uint32_t>::fillByByte( uint8_t const vec[], unsigned sizeVec ) {
    unsigned nbb = (numBits*numColsBlock)/256 ;
    ASSERT( sizeVec == numRows*nbb, "GF2Rank::fillByByte sizeVec = " << sizeVec << " expected = " << numRows*nbb ) ;

    GF2toolkit::zero( dataMatrix, dimRows, numRows, numColsBlock );

    uint8_t const *pvec = vec ;
    for ( unsigned r = 0 ; r < numRows ; ++r ) {
      for ( unsigned cb = 0 ; cb < nbb ; ++cb ) {
        uint8_t v = *pvec++ ;
        block( r, 8*cb+(v>>5) ) |= uint32_t(1)<<(v&0x1F) ;
      }
    }
  }

  template <>
  void
  LU<uint64_t>::fillByByte( uint8_t const vec[], unsigned sizeVec ) {
    unsigned nbb = (numBits*numColsBlock)/256 ;
    ASSERT( sizeVec == numRows*nbb, "GF2Rank::fillByByte sizeVec = " << sizeVec << " expected = " << numRows*nbb ) ;

    GF2toolkit::zero( dataMatrix, dimRows, numRows, numColsBlock );

    uint8_t const *pvec = vec ;
    for ( unsigned r = 0 ; r < numRows ; ++r ) {
      for ( unsigned cb = 0 ; cb < nbb ; ++cb ) {
        uint8_t v = *pvec++ ;
        block( r, 4*cb+(v>>6) ) |= uint64_t(1)<<(v&0x3F) ;
      }
    }
  }

  template <>
  void
  LU<uint128_t>::fillByByte( uint8_t const vec[], unsigned sizeVec ) {
    unsigned nbb = (numBits*numColsBlock)/256 ;
    ASSERT( sizeVec == numRows*nbb, "GF2Rank::fillByByte sizeVec = " << sizeVec << " expected = " << numRows*nbb ) ;

    GF2toolkit::zero( dataMatrix, dimRows, numRows, numColsBlock );

    uint8_t const *pvec = vec ;
    for ( unsigned r = 0 ; r < numRows ; ++r ) {
      for ( unsigned cb = 0 ; cb < nbb ; ++cb ) {
        uint8_t v = *pvec++ ;
        uint128_t & bk = block( r, 2*cb+(v>>7) ) ;
        if ( v & 0x40 ) bk.i64.hi |= uint64_t(1)<<(v&0x3F) ;
        else            bk.i64.lo |= uint64_t(1)<<(v&0x3F) ;
      }
    }
  }

  /*
  //             _     _ _____     
  //    __ _  __| | __| |_   _|__  
  //   / _` |/ _` |/ _` | | |/ _ \ 
  //  | (_| | (_| | (_| | | | (_) |
  //   \__,_|\__,_|\__,_| |_|\___/ 
  */
  template <>
  void
  LU<uint32_t>::addTo( unsigned i, unsigned j ) {
    for ( unsigned cb = 0 ; cb < numColsBlock ; ++cb )
      block(j,cb) ^= block(i,cb) ;
  }

  template <>
  void
  LU<uint64_t>::addTo( unsigned i, unsigned j ) {
    for ( unsigned cb = 0 ; cb < numColsBlock ; ++cb )
      block(j,cb) ^= block(i,cb) ;
  }

  template <>
  void
  LU<uint128_t>::addTo( unsigned i, unsigned j ) {
    for ( unsigned cb = 0 ; cb < numColsBlock ; ++cb ) {
      uint128_t & res = block(j,cb) ;
      uint128_t & val = block(i,cb) ;
      res.i128 = _mm_xor_si128( res.i128, val.i128 ) ;
    }
  }
  
  /*
  //   _     _ _                                 
  //  | |__ (_) |_    __ _  ___ ___ ___  ___ ___ 
  //  | '_ \| | __|  / _` |/ __/ __/ _ \/ __/ __|
  //  | |_) | | |_  | (_| | (_| (_|  __/\__ \__ \
  //  |_.__/|_|\__|  \__,_|\___\___\___||___/___/
  */
  //! access to the (i,j) element of the martrix
  template <>
  unsigned
  LU<uint32_t>::operator() ( unsigned i, unsigned j ) const {
    uint32_t tmp = block( i, (j>>5) ) >> (j&0x1F) ;
    return tmp & 0x01 ;
  }

  template <>
  unsigned
  LU<uint64_t>::operator() ( unsigned i, unsigned j ) const {
    uint64_t tmp = block( i, (j>>6) ) >> (j&0x3F) ;
    return tmp & 0x01 ;
  }

  template <>
  unsigned
  LU<uint128_t>::operator() ( unsigned i, unsigned j ) const {
    uint128_t tmp = block( i, (j>>7) ) ;
    if ( j & 0x40 ) return (tmp.i64.hi>>(j&0x3F)) & 0x01 ;
    else            return (tmp.i64.lo>>(j&0x3F)) & 0x01 ;
  }

  //----------------------------------------------------------------------------

  template <>
  void
  LU<uint32_t>::set( unsigned i, unsigned j ) {
    uint32_t & bk = block( i, (j>>5) ) ;
    bk |= uint32_t(1) << (j&0x1F) ;
  }

  template <>
  void
  LU<uint64_t>::set( unsigned i, unsigned j ) {
    uint64_t & bk = block( i, (j>>6) ) ;
    bk |= uint64_t(1) << (j&0x3F) ;
  }

  template <>
  void
  LU<uint128_t>::set( unsigned i, unsigned j ) {
    uint128_t & bk = block( i, (j>>7) ) ;
    if ( j & 0x40 ) bk.i64.hi |= uint64_t(1) << (j&0x3F) ;
    else            bk.i64.lo |= uint64_t(1) << (j&0x3F) ;
  }

  //----------------------------------------------------------------------------

  template <>
  void
  LU<uint32_t>::clear( unsigned i, unsigned j ) {
    uint32_t & bk = block( i, (j>>5) ) ;
    bk &= ~(uint32_t(1) << (j&0x1F)) ;
  }

  template <>
  void
  LU<uint64_t>::clear( unsigned i, unsigned j ) {
    uint64_t & bk = block( i, (j>>6) ) ;
    bk &= ~(uint64_t(1) << (j&0x3F)) ;
  }

  template <>
  void
  LU<uint128_t>::clear( unsigned i, unsigned j ) {
    uint128_t & bk = block( i, (j>>7) ) ;
    if ( j & 0x40 ) bk.i64.hi &= ~(uint64_t(1) << (j&0x3F)) ;
    else            bk.i64.lo &= ~(uint64_t(1) << (j&0x3F)) ;
  }

  /*              _       _
  //   _ __  _ __(_)_ __ | |_ 
  //  | '_ \| '__| | '_ \| __|
  //  | |_) | |  | | | | | |_ 
  //  | .__/|_|  |_|_| |_|\__|
  //  |_|                     
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::print( ostream & s ) const {
    for ( unsigned row = 0 ; row < numRows ; ++row ) {
      s << setw(3) << row << "   " ;
      for ( unsigned col = 0 ; col < numColsBlock*numBits ; ++col )
        s << ( (*this)(row,col) ? '*' : '.' ) ;
      s << '\n' ;
    }
  }
  
  /*
  //   ___ _ __  _   _ 
  //  / __| '_ \| | | |
  //  \__ \ |_) | |_| |
  //  |___/ .__/ \__, |
  //      |_|    |___/ 
  */
  //! Plot a sihouette of the sparse pattern to a file
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::spy( char const fname[] ) const {
    GF2toolkit::spy<UNSIGNED>( fname, dataMatrix, dimRows, numRows, numColsBlock ) ;
  }

  /*
  //            _                
  //   ___  ___| |__  _   _ _ __ 
  //  / __|/ __| '_ \| | | | '__|
  //  \__ \ (__| | | | |_| | |   
  //  |___/\___|_| |_|\__,_|_|   
  //
  // complemento di shur interno alla matrice
  //          jb      je
  //  +---+ ~ +-------+ rcBegin
  //  |   |   |   C   |
  //  +---+ ~ +-------+ rcEnd
  //  |   |   |       |
  //  | D |   |   E   |
  //  |   |   |       |
  //  +---+ ~ +-------+ rowEnd
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::schur( unsigned rcBegin, // riga colonna top
                       unsigned rcEnd, 
                       unsigned rowEnd,
                       unsigned jBlockBegin,
                       unsigned jBlockEnd ) {

    if ( rcBegin == rcEnd ) return ; // niente da fare, esco

    unsigned jbBegin = rcBegin>>numShift ;
    unsigned jbEnd   = rcEnd>>numShift ;
    unsigned nrows   = rowEnd - rcEnd ;
    unsigned rrb     = rcBegin & numMask ;
    unsigned nblk    = jbEnd - jbBegin ;

    if ( nblk == 0 ) { // caso degenere sullo stesso blocco colonna
      UNSIGNED const * C = &block(rcBegin,jBlockBegin) ;
      UNSIGNED const * D = &block(rcEnd,jbBegin) ;
      UNSIGNED       * E = &block(rcEnd,jBlockBegin) ;
      unsigned nr  = rcEnd - rcBegin ;
      for ( unsigned jb = jBlockBegin ; jb < jBlockEnd ; ++jb, C += dimRows, E += dimRows ) {
        mm4r_reduced.makeTable( C, nr ) ;
        if ( rrb ) mm4r_reduced.multRightAddShift( D, nrows, E, rrb ) ;
        else       mm4r_reduced.multRightAdd( D, nrows, E ) ;
      }
    } else {
      unsigned rre = rcEnd & numMask ;
      if ( rre ) {
        UNSIGNED const * C = &block(rcEnd-rre,jBlockBegin) ;
        UNSIGNED const * D = &block(rcEnd,jbEnd) ; // ultimo blocco colonna
        UNSIGNED       * E = &block(rcEnd,jBlockBegin) ;
        for ( unsigned jb = jBlockBegin ; jb < jBlockEnd ; ++jb, C += dimRows, E += dimRows ) {
          mm4r_reduced.makeTable( C, rre ) ;
          mm4r_reduced.multRightAdd( D, nrows, E ) ;
        }
      }
      if ( rrb ) {
        UNSIGNED const * C = &block(rcBegin,jBlockBegin) ;
        UNSIGNED const * D = &block(rcEnd,jbBegin) ; // primo blocco colonna
        UNSIGNED       * E = &block(rcEnd,jBlockBegin) ;
        unsigned nr = numBits - rrb ;
        for ( unsigned jb = jBlockBegin ; jb < jBlockEnd ; ++jb, C += dimRows, E += dimRows ) {
          mm4r_reduced.makeTable( C, nr ) ;
          mm4r_reduced.multRightAddShift( D, nrows, E, rrb ) ;
        }
        --nblk ; // colonna gia usata
        ++jbBegin ; // passo colonna sucessiva
        rcBegin += nr ;
      }
      ASSERT( nblk * numBits == (rcEnd-rre) - rcBegin, 
              "Errore in dimensioni, nblk = " << nblk <<
              " rcBegin = " << rcBegin <<
              " rcEnd = " << rcEnd <<
              " rre = " << rre ) ;

      if ( nblk > 0 ) { // solo se necessario
        UNSIGNED const * C = &block(rcBegin,jBlockBegin) ;
        UNSIGNED const * D = &block(rcEnd,jbBegin) ; // primo blocco colonna
        UNSIGNED       * E = &block(rcEnd,jBlockBegin) ;
        GF2toolkit::MMaddStrassen<UNSIGNED>( dataWork, sizeWork,
                                             D, dimRows, nrows, nblk,
                                             C, dimRows, nblk*numBits, jBlockEnd-jBlockBegin,
                                             E, dimRows ) ;
      }
    }
  }
  
  /*
  //////////////////////////////////////////////////////////////////////////////
  // applicazione parziale
  //
  //  +--|--+-----+-----+
  //  |     |     |     |
  //  |  +~~|~~~~~|~~~+ |
  //  |  :  |     |   : |
  //  |  :  |     |   : |
  //  |  :  |     |   : |   L 
  //  |  :  |     |   : |
  //  |  :  |     |   : |
  //  |  +~~|~~~~~|~~~+ |
  //  +-----+-----+-----+
  //
  */
  /*
  //   _     _             _                _       
  //  | |   (_)_ ____   __/ \   _ __  _ __ | |_   _ 
  //  | |   | | '_ \ \ / / _ \ | '_ \| '_ \| | | | |
  //  | |___| | | | \ V / ___ \| |_) | |_) | | |_| |
  //  |_____|_|_| |_|\_/_/   \_\ .__/| .__/|_|\__, |
  //                           |_|   |_|      |___/ 
  */
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::LinvApply( unsigned rcBegin,     // riga colonna top
                           unsigned rcEnd,       // riga colonna bot
                           unsigned jBlockBegin,
                           unsigned jBlockEnd ) {

    UNSIGNED Linv[numBits], Lsmall[numBits] ;

    if ( rcBegin == rcEnd ) return ; // niente da fare, esco

    unsigned jbBegin = rcBegin>>numShift ;
    unsigned jbEnd   = rcEnd>>numShift ;
    unsigned rrb     = rcBegin & numMask ;
    unsigned nblk    = jbEnd - jbBegin ;
    //MM4Rreduced<UNSIGNED> mm4r_reduced ;

    if ( nblk == 0 ) { // caso degenere sullo stesso blocco colonna
      UNSIGNED const * L = &block(rcBegin,jbBegin) ;
      UNSIGNED       * C = &block(rcBegin,jBlockBegin) ;
      unsigned nr = rcEnd - rcBegin ;

      // estrae L
      if ( rrb == 0 ) GF2toolkit::copy<UNSIGNED>( L, nr, Lsmall ) ;
      else            GF2toolkit::copyRightShift<UNSIGNED>( L, nr, Lsmall, rrb ) ;
      GF2toolkit::Linverse<UNSIGNED>( Lsmall, Linv, nr ) ;

      for ( unsigned jb = jBlockBegin ; jb < jBlockEnd ; ++jb, C += dimRows ) {
        mm4r_reduced.makeTable( C, nr ) ;
        mm4r_reduced.multRightAss( Linv, nr, C ) ; // L^(-1)*C --> C
      }
    } else {
      if ( rrb ) {
        UNSIGNED const * L = &block(rcBegin,jbBegin) ;
        UNSIGNED       * C = &block(rcBegin,jBlockBegin) ;
        unsigned nr    = numBits - rrb ;
        unsigned nrows = rcEnd - rcBegin - nr ; // righe sotto la L
        // estrae L
        GF2toolkit::copyRightShift<UNSIGNED>( L, nr, Lsmall, rrb ) ;
        GF2toolkit::Linverse<UNSIGNED>( Lsmall, Linv, nr ) ;

        for ( unsigned jb = jBlockBegin ; jb < jBlockEnd ; ++jb, C += dimRows ) {
          mm4r_reduced.makeTable( C, nr ) ;
          mm4r_reduced.multRightAss( Linv, nr, C ) ; // L^(-1)*C --> C
          mm4r_reduced.makeTable( C, nr ) ;
          mm4r_reduced.multRightAddShift( L+nr, nrows, C+nr, rrb ) ;
        }
        --nblk ; // colonna gia usata
        ++jbBegin ; // passo colonna sucessiva
        rcBegin += nr ; // non considero righe gia usare
      }
      if ( nblk > 0 ) { // solo se necessario
        GF2toolkit::LinvApply<UNSIGNED>( dataWork, sizeWork,
                                         &block(rcBegin,jbBegin),     dimRows, rcEnd - rcBegin, nblk,
                                         &block(rcBegin,jBlockBegin), dimRows, jBlockEnd - jBlockBegin ) ;
        jbBegin += nblk ; // passo colonna sucessiva
        rcBegin += nblk * numBits ;
      }
      unsigned rre = rcEnd & numMask ;
      if ( rre ) {
        UNSIGNED const * L = &block(rcBegin,jbBegin) ;
        UNSIGNED       * C = &block(rcBegin,jBlockBegin) ;
        // estrae L
        GF2toolkit::copy<UNSIGNED>( L, rre, Lsmall ) ;
        GF2toolkit::Linverse<UNSIGNED>( Lsmall, Linv, rre ) ;

        for ( unsigned jb = jBlockBegin ; jb < jBlockEnd ; ++jb, C += dimRows ) {
          mm4r_reduced.makeTable( C, rre ) ;
          mm4r_reduced.multRightAss( Linv, rre, C ) ; // L^(-1)*C --> C
        }
      }
    }
  }

  template class LU<uint32_t> ;
  template class LU<uint64_t> ;
  template class LU<uint128_t> ;

}

// EOF GF2toolkit_LU.cc

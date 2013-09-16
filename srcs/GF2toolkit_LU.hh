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

#ifndef GF2TOOLKIT_RANK_HH
#define GF2TOOLKIT_RANK_HH

#include "TimeMeter.hh"
#include "GF2toolkit_common.hh"
#include "GF2toolkit_Matrix.hh"
#include "GF2toolkit_m4ri.hh"
#include "GF2toolkit_4Russian.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

///////////////////////////////////////////////////////////////////
//                                                               //
//  Structure of the n x m  matrix B(k) is a block of 256 bits   //
//  N = m/256                                                    //
//                                                               //
//  B(0)   B(n+0)   B(2*n+0)  ...   B(n*(N-1)+0)                 //
//  B(1)   B(n+1)   B(2*n+1)  ...   B(n*(N-1)+1)                 //
//  B(2)   B(n+2)   B(2*n+2)  ...   B(n*(N-1)+2)                 //
//                                                               //
//  ........                                                     //
//  ........                                                     //
//                                                               //
//  B(n-1) B(2*n-1) B(3*n-1)  ...   B(n*N-1)                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#define USE_TIMER

namespace GF2toolkit {

  using namespace std ;

  #ifdef USE_TIMER
    #define TIMER_BEGIN(A) tm[A].start()
    #define TIMER_END(A)   tm[A].stop()
  #else
    #define TIMER_BEGIN(A)
    #define TIMER_END(A)  
  #endif

  /*
  //   _    _   _ 
  //  | |  | | | |
  //  | |  | | | |
  //  | |__| |_| |
  //  |_____\___/ 
  */
  
  #define N_TIMERS 15

  template <typename UNSIGNED>
  class LU {

    static unsigned const numBits  ;
    static unsigned const numShift ;
    static unsigned const numMask  ;

    #ifdef USE_TIMER
    TimeMeter tm[N_TIMERS] ; // 15 timer
    #endif
    
    MM4Rreduced<UNSIGNED> mm4r_reduced ;

    UNSIGNED * dataMatrix ;
    UNSIGNED * dataWork ;
    UNSIGNED * dataUmat ;
    UNSIGNED * dataLmat ;
    unsigned   sizeWork ;

    unsigned * rowStart ;
    int32_t  * Perm ;

    unsigned   numRows ;
    unsigned   numColsBlock ;
    unsigned   dimRows ;

    void
    start() {
      #ifdef USE_TIMER
      for ( unsigned i = 0 ; i < N_TIMERS ; ++i ) tm[i].reset() ;
      #endif
    }

    /* applico L^(-1) al blocco C alla destra
    //          jb      je
    //  +---+ ~ +-------+ rcBegin
    //  | L |   |   C   |
    //  +---+ ~ +-------+ rcEnd
    //  |   |   |       |
    //  :   :   :       :
    //  |   |   |       |
    //  +---+ ~ +-------+
    */
    void
    LinvApply( unsigned rcBegin,     // riga colonna top
               unsigned rcEnd,       // riga colonna bot
               unsigned jBlockBegin,
               unsigned jBlockEnd ) ;

    /* complemento di shur interno alla matrice
    //          jb      je
    //  +---+ ~ +-------+ rcBegin
    //  |   |   |   C   |
    //  +---+ ~ +-------+ rcEnd
    //  |   |   |       |
    //  | D |   |   E   |
    //  |   |   |       |
    //  +---+ ~ +-------+ rowEnd
    */
    void
    schur( unsigned rcBegin,     // riga colonna top
           unsigned rcEnd,       // riga colonna bot
           unsigned rowEnd,      // ultima riga
           unsigned jBlockBegin,
           unsigned jBlockEnd ) ;

    void computeRankRecurr( unsigned jbBegin, unsigned jbEnd ) ;

  public:

    //! initialize an empty matrix
    LU() ;
    
    //! initialize a matrix with \c nr rows and \c nc columns, \nc must be a multiple of 256  
    LU( unsigned nr, unsigned nc ) ;

    //! matrix destroier (automatically called)
    ~LU() ;

    //! allocate memory for a \c nr x \c nc matrix
    void resize( unsigned nr, unsigned nc ) ;

    //! access to the (i,j) element of the martrix
    unsigned 
    operator() ( unsigned i, unsigned j ) const ;

    void set( unsigned i, unsigned j ) ;
    void clear( unsigned i, unsigned j ) ;
    void addTo( unsigned i, unsigned j ) ;

    unsigned nnz() const ;

    //! fill a matrix with zeros
    void clear() { GF2toolkit::zero( dataMatrix, dimRows, numRows, numColsBlock ) ; }

    //! fill a matrix with user values
    void fillByBits( uint64_t const vec[], unsigned sizeVec ) ;
    void fillByByte( uint8_t const vec[], unsigned sizeVec ) ;

    mzd_t * 
    toM4RI() const 
    { return GF2toolkit::toM4RI<UNSIGNED>( dataMatrix, dimRows, numRows, numColsBlock ) ; }
  
    void
    fromM4RI( mzd_t * mat ) const 
    { GF2toolkit::fromM4RI<UNSIGNED>( mat, dataMatrix, dimRows, numRows, numColsBlock ) ; }

    // organization, block of columns
    UNSIGNED & block( unsigned row, unsigned jb ) const 
    { return dataMatrix[row+jb*dimRows] ; }

    // organization, block of columns
    UNSIGNED & Ublock( unsigned jb ) const 
    { return dataUmat[rowStart[jb]] ; }

    // organization, block of columns
    UNSIGNED & Lblock( unsigned jb ) const 
    { return dataLmat[jb*numBits] ; }

    unsigned computeRank() ;
    unsigned computeRankRecurr() ;

    //! Plot a sihouette of the sparse pattern to a file
    void spy( char const fname[] ) const ;
    
    //! print the matrix in human readable form to the stream \c s
    void print( std::ostream & s ) const ;

    void
    elapsed( ostream & stream ) {
      #ifdef USE_TIMER
      for ( unsigned i = 0 ; i < N_TIMERS ; ++i )
        if ( tm[i].totalElapsedMilliseconds() > 0 ) 
          stream << "Timer" << i << " = " << tm[i].totalElapsedMilliseconds() << "ms\n" ;
      #endif
    }

    void extractA( UNSIGNED * A, unsigned dimRows, unsigned numBlocks,
                   unsigned & nr, unsigned & nc ) ;

    void extractL( UNSIGNED * L, unsigned dimRows, unsigned numBlocks,
                   unsigned & nr, unsigned & nc ) ;

    void extractU( UNSIGNED * U, unsigned dimRows, unsigned numBlocks,
                   unsigned & nr, unsigned & nc ) ;

    void applyPermute( UNSIGNED * M, unsigned dimRows, unsigned numBlocks ) ;

  } ;
}

template <typename UNSIGNED>
inline 
std::ostream & 
operator << ( std::ostream & s, GF2toolkit::LU<UNSIGNED> const & zm )
{ zm.print( s ) ; return s ; }

#endif

// EOF GF2toolkit_LU.hh

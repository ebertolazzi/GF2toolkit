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
#include "GF2toolkit_InvertBlock.hh"
#include "GF2toolkit_Matrix.hh"
#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_MatrixLower.hh"
#include "GF2toolkit_4Russian.hh"

#include <iostream>
#include <iomanip>

namespace GF2toolkit {
/*
//    __           _   _                ____  
//   / _| __ _ ___| |_| |    ___   __ _|___ \ 
//  | |_ / _` / __| __| |   / _ \ / _` | __) |
//  |  _| (_| \__ \ |_| |__| (_) | (_| |/ __/ 
//  |_|  \__,_|___/\__|_____\___/ \__, |_____|
//                                |___/       
*/

////////////////////////////////////////////////////////////

static
inline
unsigned
fastLog2( uint16_t n ) {
  unsigned r, shift ;
  r     = unsigned( n > 0xFFu ) << 3 ; n >>= r;
  shift = unsigned( n > 0xFu  ) << 2 ; n >>= shift; r |= shift;
  shift = unsigned( n > 0x3u  ) << 1 ; n >>= shift; r |= shift;
  r |= (n >> 1);
  return r ;
}

////////////////////////////////////////////////////////////
  
static
inline
unsigned
fastLog2( uint32_t n ) {
  unsigned r, shift ;
  r     = unsigned( n > 0xFFFFu ) << 4 ; n >>= r;
  shift = unsigned( n > 0xFFu   ) << 3 ; n >>= shift; r |= shift;
  shift = unsigned( n > 0xFu    ) << 2 ; n >>= shift; r |= shift;
  shift = unsigned( n > 0x3u    ) << 1 ; n >>= shift; r |= shift;
  r |= (n >> 1);
  return r ;
}

////////////////////////////////////////////////////////////

static
inline
unsigned
fastLog2( uint64_t n ) {
  unsigned r, shift ;
  r     = unsigned( n > 0xFFFFFFFFu ) << 5 ; n >>= r     ;
  shift = unsigned( n > 0xFFFFu     ) << 4 ; n >>= shift ; r |= shift ;
  shift = unsigned( n > 0xFFu       ) << 3 ; n >>= shift ; r |= shift ;
  shift = unsigned( n > 0xFu        ) << 2 ; n >>= shift ; r |= shift ;
  shift = unsigned( n > 0x3u        ) << 1 ; n >>= shift ; r |= shift ;
  r |= unsigned(n >> 1);
  return r ;
}

  /*
  //                                   _       ____             _    
  //    ___ ___  _ __ ___  _ __  _   _| |_ ___|  _ \ __ _ _ __ | | __
  //   / __/ _ \| '_ ` _ \| '_ \| | | | __/ _ \ |_) / _` | '_ \| |/ /
  //  | (_| (_) | | | | | | |_) | |_| | ||  __/  _ < (_| | | | |   < 
  //   \___\___/|_| |_| |_| .__/ \__,_|\__\___|_| \_\__,_|_| |_|_|\_\
  //                      |_|
  */
  
  template <typename UNSIGNED>
  unsigned
  LU<UNSIGNED>::computeRankRecurr() {

    #ifdef DEBUG
    cout << "SIZE = " << sizeof(UNSIGNED) 
         << " numBits = "  << numBits
         << " numShift = " << numShift 
         << " numColsBlock = " << numColsBlock 
         << '\n' ;
    #endif

    start() ;
    rowStart[0] = 0 ;
    LU<UNSIGNED>::computeRankRecurr( 0, numColsBlock ) ;

    //elapsed( cout ) ;
    //cout << "NCOL = " << numColsBlock*numBits 
    //     << " NROW = " << rowStart[numColsBlock]
    //     << '\n' ;

    return rowStart[numColsBlock] ;
  }
  
  template <typename UNSIGNED>
  void
  LU<UNSIGNED>::computeRankRecurr( unsigned jBlockBegin, unsigned jBlockEnd ) {

    //cout << jBlockBegin << " - " << jBlockEnd << '\n' ;
    
    unsigned lg2     = fastLog2((jBlockEnd - jBlockBegin)>>1) ;
    unsigned jjBlock = 1<<lg2 ;
    
    bool minsizeC = (jjBlock<<numShift) >= STRASSEN_MIN_SIZE ;
    bool minsizeR = numRows - rowStart[jBlockBegin] >= STRASSEN_MIN_SIZE  ;

    if ( minsizeC && minsizeR ) {
      /*
      // +----+------------+
      // | LU |     C      |
      // +----+------------+
      // |    |            |
      // |    |            |
      // | D  |     E      |
      // |    |            |
      // |    |            |
      // +----+------------+
      */
      
      /*
      //
      //  +---+---+ +---+---+
      //  | L | 0 | | U |iLC|
      //  +---+---+ +---+---+
      //  | D | I | | 0 | S |
      //  +---+---+ +---+---+
      //
      //  S + D L^(-1) C = E
      //
      //
      */
            
      unsigned jBlockMid = jBlockBegin+jjBlock ;
      LU<UNSIGNED>::computeRankRecurr( jBlockBegin, jBlockMid ) ;
      
      /*
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
      // applico L^(-1) a C
      #if 1
      unsigned rcBegin = rowStart[jBlockBegin] ;
      unsigned rcMid   = rowStart[jBlockMid] ;
      TIMER_BEGIN(7) ;
      LinvApply( rcBegin, rcMid, jBlockMid, jBlockEnd ) ;
      TIMER_END(7) ;
      TIMER_BEGIN(8) ;
      schur( rcBegin, rcMid, numRows, jBlockMid, jBlockEnd ) ;
      TIMER_END(8) ;
      #else
      for ( unsigned jb = jBlockBegin ; jb < jBlockMid ; ++jb )
        schur( rowStart[jb], rowStart[jb+1], numRows, jBlockMid, jBlockEnd ) ;
      #endif
      LU<UNSIGNED>::computeRankRecurr( jBlockMid, jBlockEnd ) ;

    } else {

      UNSIGNED Z[numBits] ;
      UNSIGNED *B, *C, *D ;

      int32_t * P = Perm + jBlockBegin * numBits ;
      unsigned  n = jBlockBegin ; 
      for ( ; n < jBlockEnd ; ++n, P += numBits ) {

        /*
        // +---+------------+
        // | B |     C      |
        // +---+------------+
        // |   |            |
        // |   |            |
        // | D |     E      |
        // |   |            |
        // |   |            |
        // +---+------------+
        */

        unsigned rstart = rowStart[n] ;
        unsigned nrows  = numRows-rstart ;

        if ( nrows == 0 ) break ; // non ci sono piu righe da elaborare!

        B = &block(rstart,n) ;

        // estrazione blocco non singolare
        TIMER_BEGIN(0);
        unsigned NR = GF2toolkit::invertBlock<UNSIGNED>( B, nrows, P, Z ) ;
        rowStart[n+1] = rstart + NR ; // NR = righe linearmente indipendenti
        TIMER_END(0) ;

        if ( NR == 0 ) continue ; // niente da fare per questo blocco

        // applico permutazione (SU TUTTA LA MATRICE!)
        TIMER_BEGIN(1) ;
        C = &block(rstart,0) ;
        for ( unsigned jb = 0 ; jb < numColsBlock ; ++jb, C += dimRows )
          if ( jb != n )
            GF2toolkit::permute<UNSIGNED>( C, numBits, P ) ;
        TIMER_END(1) ;

        // estraggo U e salvo in matrice U e metto pezzo indentità al suo posto
        TIMER_BEGIN(2) ;
        memcpy( &Ublock(n), B, sizeof(UNSIGNED) * NR ) ; // salvo U
        setIdentity<UNSIGNED>( B, NR ) ; // copia pezzo identità
        TIMER_END(2) ;

        /*
        //  +---+ +---+    
        //  |   | | Z |    
        //  | D | +---+
        //  |   |
        //  +---+
        //  moltiplico matrice colonna per matrice quadrata
        */
        TIMER_BEGIN(3);
        C      = B + dimRows ; // rimetto B alla sua posizione
        D      = B + NR ;
        nrows -= NR ;
        GF2toolkit::MMassLeft( D, nrows, Z ) ;
        TIMER_END(3);

        // complemento di schur (SU MATRICE PARZIALE)
        for ( unsigned jb = n+1 ; jb < jBlockEnd ; ++jb, C += dimRows ) { // <----------------
          TIMER_BEGIN(4);
          mm4r_reduced.makeTable( C, NR ) ;
          TIMER_END(4);
          TIMER_BEGIN(5);
          mm4r_reduced.multRightAdd( D, nrows, C+NR ) ; // E += D*C
          TIMER_END(5);
        }

        // compatto colonne di L
        TIMER_BEGIN(6);
        if ( rstart < (n<<numShift) ) { // posso compattare
          unsigned nblk     = rstart>>numShift ;
          unsigned startcol = rstart&numMask ;

          if ( nblk+1 == n )
            compactColumn<UNSIGNED>( &block(rstart,n),      // blocco colonna da compattare
                                     NR,                    // quante colonne contiene
                                     &block(rstart,nblk),   // blocco colonna di arrivo
                                     startcol, // quante colonne sono gia occupate in arrivo
                                     numRows-rstart ) ;        
          else
            compactColumn<UNSIGNED>( &block(rstart,n),      // blocco colonna da compattare
                                     NR,                    // quante colonne contiene
                                     &block(rstart,nblk),   // blocco colonna di arrivo
                                     &block(rstart,nblk+1), // blocco colonna successivo
                                     startcol, // quante colonne sono gia occupate in arrivo
                                     numRows-rstart ) ;
        }
        TIMER_END(6);
      }

      while ( n < jBlockEnd ) { // <----------------
        for ( unsigned i = 0 ; i < numBits ; ++i ) *P++ = -1 ; // nessuna permutazione
        rowStart[n+1] = rowStart[n] ; ++n ;
      }

    }
  }

  template unsigned LU<uint32_t>::computeRankRecurr() ;
  template unsigned LU<uint64_t>::computeRankRecurr() ;
  template unsigned LU<uint128_t>::computeRankRecurr() ;

  template void LU<uint32_t>::computeRankRecurr(unsigned,unsigned) ;
  template void LU<uint64_t>::computeRankRecurr(unsigned,unsigned) ;
  template void LU<uint128_t>::computeRankRecurr(unsigned,unsigned) ;

}

// EOF GF2toolkit_RankComputeRecurr.cc

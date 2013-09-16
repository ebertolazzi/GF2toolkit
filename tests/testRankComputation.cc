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
 |      Dipartimento di Ingegneria Meccanica e Strutturale                  |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 |      version: 0.2 21-01-2010                                             |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "TimeMeter.hh"
#include "Random.hh"
#include "GF2toolkit_Matrix.hh"
#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_Strassen.hh"
#include "GF2toolkit_LU.hh"
#include "GF2toolkit_m4ri.hh"
#include "GF2toolkit_popCounts.hh"
#include <iostream>
#include <iomanip>

using namespace std ;

//RNG::MersenneTwister mt(912130) ;
RNG::CMRG mt(21232) ;

template <typename UNSIGNED>
UNSIGNED 
random() ;

template <>
uint32_t 
random()
{ return mt.rand() ; }

template <>
uint64_t 
random() {
  uint64_t r1 = mt.rand() ;
  uint64_t r2 = mt.rand() ;
  return r1|(r2<<32) ;
}

template <>
uint128_t 
random() {
  uint128_t tmp ;
  tmp.i64.lo = random<uint64_t>() ;
  tmp.i64.hi = random<uint64_t>() ;
  return tmp ;
}

template <typename UNSIGNED>
unsigned
doTest( unsigned nr, unsigned nc, uint64_t const TMP[] ) {

  unsigned numBits = CHAR_BIT * sizeof(UNSIGNED) ;
  unsigned nblk    = nc/64 ;
  unsigned num     = nr*nblk ;

  TimeMeter                tm ;
  GF2toolkit::LU<UNSIGNED> mat ;

  mat . resize(nr,nc) ;
  mat . fillByBits( TMP, num ) ;
  
  unsigned dimRowsA   = nr ;
  unsigned numBlocksA = (nc+numBits-1)/numBits ;

  UNSIGNED * Amat = new UNSIGNED[dimRowsA*numBlocksA] ;
  mat .extractA( Amat, dimRowsA, numBlocksA, nr, nc ) ;

  cout << "Start..." << flush ;
  tm . reset() ;
  tm . start() ;
  //unsigned rk = mat.computeRank() ;
  unsigned rk = mat.computeRankRecurr() ;
  tm . stop() ;
  cout << "done\nRank (" << setw(3) << numBits << ")= "
       << rk << " Time = " << tm . partialElapsedMilliseconds()
       << "ms\n" ;

  // chek solution
  cout << "Check solution\n" << flush ;
  tm . start() ;

  unsigned dimRowsL   = dimRowsA ;
  unsigned numBlocksL = numBlocksA ;
  unsigned dimRowsU   = dimRowsA ;
  unsigned numBlocksU = numBlocksA ;
  unsigned nrL, ncL, nrU, ncU ;

  UNSIGNED * L = new UNSIGNED[dimRowsL*numBlocksL] ;
  UNSIGNED * U = new UNSIGNED[dimRowsU*numBlocksU] ;
    
  unsigned sizeW = nr*(nc/numBits) ;
  UNSIGNED * W = new UNSIGNED[sizeW] ;

  mat.extractL( L, dimRowsL, numBlocksL, nrL, ncL ) ;
  mat.extractU( U, dimRowsU, numBlocksU, nrU, ncU ) ;

  //cout << "nrL = " << nrL << '\n' ;
  //cout << "ncL = " << ncL << '\n' ;
  //cout << "nrU = " << nrU << '\n' ;
  //cout << "ncU = " << ncU << '\n' ;
    
  nrU = ((nrU+numBits-1)/numBits)*numBits ;
    
  unsigned nColBlockL = (ncL+numBits-1)/numBits ;
  unsigned nColBlockU = (ncU+numBits-1)/numBits ;

  mat.applyPermute( Amat, dimRowsA, numBlocksA ) ;
  tm . stop() ;
  cout << "permute,  Time = " << tm . partialElapsedMilliseconds() << "ms\n" << flush  ;

  tm . start() ;
  GF2toolkit::MMaddStrassen<UNSIGNED>( W, sizeW,
                                       L, dimRowsL, nrL, nColBlockL,
                                       U, dimRowsU, nrU, nColBlockU,
                                       Amat, dimRowsA ) ;
  tm . stop() ;
  cout << "strassen, Time = " << tm . partialElapsedMilliseconds() << "ms\n" << flush   ;

  tm . start() ;
  bool ok = GF2toolkit::isZero<UNSIGNED>( Amat, dimRowsA, nrL, nColBlockU ) ;
  tm . stop() ;
  cout << "check,    Time = " << tm . partialElapsedMilliseconds()
       << (ok?"ms --- OK\n":"ms --- ERROR!\n") ;

  delete [] Amat ;
  delete [] W ;
  delete [] L ;
  delete [] U ;
  
  return rk ;
}


int
main() {

  TimeMeter          tm ;

  //unsigned nr = 64000 ; // 32768 ;
  //unsigned nc = 64000 ; // 32768 ;

  unsigned nr = 32768 ;
  unsigned nc = 32768 ;

  //unsigned nr = 4*128 ;
  //unsigned nc = 6*128 ;

  //unsigned nr = 50000 ; // 32768 ;
  //unsigned nc = 200*128 ; // 32768 ;

  //unsigned nr = 400*128 ;
  //unsigned nc = 20480 ;

  //unsigned nr = 19968 ;
  //unsigned nc = 19968 ;

  //unsigned nr = 128 ;
  //unsigned nc = 128 ;

  //unsigned nr = 16128 ;
  //unsigned nc = 16128 ;

  //unsigned nr = 3*10240 ;
  //unsigned nc = 3*10240 ;

  //unsigned nr = 31876 ; // 32768 ;
  //unsigned nc = 32768 ; // 32768 ;
  //unsigned nr = 256 ; // 32768 ;
  //unsigned nc = 10*256 ; // 32768 ;
  //unsigned nr = 2*64 ;
  //unsigned nc = 2*64 ;

  cout << "nrow = " << nr << " ncol = " << nc << '\n' ;

  unsigned nblk = nc/64 ;
  unsigned num  = nr*nblk ;
  uint64_t * TMP = new uint64_t[num] ;

  memset( TMP, 0, num*sizeof(uint64_t) ) ;

  for ( unsigned r = 0 ; r < nr ; ++r )
    for ( unsigned jb = 0 ; jb < nblk ; ++jb )
      TMP[jb+r*nblk] = random<uint64_t>() & 0xF37F1966FFFF6519ULL ;

  // riduco artificialmente rango
  #if 0
  for ( unsigned r = 0 ; r < nr ; ++r ) {
    unsigned r1 = random<unsigned>() % nr ;
    for ( unsigned jb = 0 ; jb < nblk ; ++jb )
      TMP[jb+r*nblk] = TMP[jb+r1*nblk] ;
  }
  #endif
  #if 0
  for ( unsigned jb = 1 ; jb < nblk ; ++jb )
    if ( (random<unsigned>() % 1) == 0 )
      for ( unsigned r = 0 ; r < nr ; ++r )
        TMP[jb+r*nblk] = TMP[(jb-1)+r*nblk] ;
  #endif
  #if 0
  TMP[1]    = 0xFFFFFFFFFFFFFFFFu ;
  TMP[nblk] = 0xFFFFFFFFFFFFFFFFu ;
  for ( unsigned r = 1 ; r < nr ; ++r ) {
    TMP[0+r*nblk] ^= r+(uint64_t(r)<<32) ;
    TMP[1+r*nblk] = TMP[0+r*nblk] ;
  }
  #endif
  
  unsigned rk, rk1, rk2, rk3 ;

  rk = doTest<uint32_t>(nr,nc,TMP) ;
  rk = doTest<uint64_t>(nr,nc,TMP) ;
  rk = doTest<uint128_t>(nr,nc,TMP) ;
  
  #if 1

  mzd_t *A = mzd_init(nr, nc);
  for ( unsigned i = 0 ; i < nr ; ++i )
    for ( unsigned jb = 0 ; jb < nblk ; ++jb )
      A->rows[i][jb] = TMP[jb+i*nblk] ;

  mzd_t *B = mzd_copy(NULL, A);
  mzd_t *C = mzd_copy(NULL, A);
  mzd_t *D = mzd_copy(NULL, A);
  mzp_t* P = mzp_init(nr);
  mzp_t* Q = mzp_init(nc);

  tm . start() ;
  rk1 = mzd_ple(B, P, Q, 0);
  tm . stop() ;
  cout << "Rank (ple)  = " << rk1 << " Time = " << tm.partialElapsedMilliseconds() << "ms\n" ;

  ASSERT( rk == rk1, "Errore nel rango" ) ;

  tm . start() ;
  rk2 = mzd_pluq(C, P, Q, 0);
  tm . stop() ;
  cout << "Rank (pluq) = " << rk2 << " Time = " << tm.partialElapsedMilliseconds() << "ms\n" ;

  tm . start() ;
  rk3 = mzd_echelonize_pluq(D, 1);
  tm . stop() ;
  cout << "Rank (m4ri) = " << rk3 << " Time = " << tm.partialElapsedMilliseconds() << "ms\n" ;

  ASSERT( rk == rk1 && rk == rk2 && rk == rk3, "Errore nel rango" ) ;

  #endif

  cout << "\n\n\n----------------- DONE -------------------\n\n\n" ;

  return 0 ;
}

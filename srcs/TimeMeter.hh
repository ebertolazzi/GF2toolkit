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
 |      version: 0.2 07-03-2011                                             |
 |                                                                          |
 \*--------------------------------------------------------------------------*/

#ifndef TIME_METER_HH
#define TIME_METER_HH

bool getTime( long & sec, long & usec ) ;

class TimeMeter {
  
  double elapsedTotal, elapsedPartial ;
  long   sec, usec ;

  TimeMeter( TimeMeter const & ) ;
  TimeMeter const & operator = ( TimeMeter const & ) const ;

public:

  TimeMeter() : elapsedTotal(0) { start() ; }
  ~TimeMeter() {} ;

  void reset() { elapsedTotal = 0 ; }
  void start() { getTime( sec, usec ) ; }
  void stop() {
    long new_sec, new_usec ;
    getTime( new_sec, new_usec ) ;
    elapsedPartial = (new_sec-sec)+(new_usec-usec)*1E-6 ; 
    elapsedTotal  += elapsedPartial ;
  }
  double totalElapsedSeconds()      const { return elapsedTotal ; }
  double totalElapsedMilliseconds() const { return elapsedTotal*1000 ; }

  double partialElapsedSeconds()      const { return elapsedPartial ; }
  double partialElapsedMilliseconds() const { return elapsedPartial*1000 ; }

} ;

#endif

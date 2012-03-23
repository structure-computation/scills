#ifndef CROUT_H
#define CROUT_H
//
// C++ Interface: crout
//
// Description: 
//
//
// Author: Alain CAIGNOT <caignot@lmt.ens-cachan.fr>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
using namespace LMT;


struct Crout {
  void open(unsigned mpi_rank) {
    f.open( ("/tmp/res"+to_string(mpi_rank)).c_str() );
  }
  std::ofstream f;
};
template<class T>
    Crout &operator<<(Crout &croute, const T &toto) {
    //std::cout << toto;
  croute.f << toto << std::flush;
  return croute;
    }

#endif //CROUT_H

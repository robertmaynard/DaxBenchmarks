#ifndef __ZCURVE_H
#define __ZCURVE_H

#include <cstdint>
#include <array>

namespace zcurve
{
  /** Converts the given point to a zcurve index, assuming uniform size dims */
  // L_max: length per dimension (must be power of 2)
  template <class Index,int_fast8_t D>
  Index convertPointToZ(std::array<Index,D> p,int_fast8_t L_max)
  {
    Index I=0;
    for(int_fast8_t L=0;L<=L_max;)
      for(int_fast8_t i=0;i<D;++i,++L)
      {
	I|=(p[i]&1)<<L;
	p[i]>>=1;
      }
    return I;
  }
 
  //! Converts the given index to a point
  template <class Index,int_fast8_t D>
  std::array<Index,D> convertZToPoint(Index I)
  {
    std::array<Index,D> P;
    std::fill(&P[0],&P[D],0);
    for(int_fast8_t L=0;I;++L,I>>=1)
      P[L%D]|=(I&1LL)<<(L/D);
    return P;
  }
}

void gen_zcurve(int dim[4], std::vector<int> &mappingID);

#endif

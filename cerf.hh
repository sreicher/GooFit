#ifndef CERF_HH
#define CERF_HH

#include <cmath>
#include "devcomplex.hh"

namespace Cerf{

  // For documentation, see RooMath.h or ask @author Manuel Schiller <manuel.schiller@nikhef.nl>
  static devcomplex faddeeva(devcomplex z);
  static devcomplex faddeeva_fast(devcomplex z);
  static devcomplex erf(const devcomplex z);
  static devcomplex erf_fast(const devcomplex z);
  static devcomplex erfc(const devcomplex z);
  static devcomplex erfc_fast(const devcomplex z);

}; 

#endif 

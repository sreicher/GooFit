#ifndef RESONANCE_PDF_HH
#define RESONANCE_PDF_HH

#include "GooPdf.hh" 
#include "devcomplex.hh"
enum DecayChannels {pipi = 0, KK =1, four_pi =2, hh =3, hhprime =4};
 
typedef devcomplex<fptype> (*resonance_function_ptr) (fptype, fptype, fptype, unsigned int*); 

class ResonancePdf : public GooPdf {
  // Service class intended to hold parametrisations of
  // resonances on Dalitz plots. Don't try to use this
  // as a standalone PDF! It should only be used as a
  // component in one of the friend classes. It extends
  // GooPdf so as to take advantage of the 
  // infrastructure, but will crash if used on its own. 

  friend class TddpPdf;
  friend class DalitzPlotPdf; 
  friend class IncoherentSumPdf; 
public:
  // Constructor for the K Matrix formalism
  ResonancePdf (string name,
                            Variable* ar,             
                            Variable* ai,             
                            Variable* Spr0,
                            unsigned int term,
                            unsigned int sp,         
                            unsigned int cyc);


  // Constructor for regular BW 
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int sp, 
			  unsigned int cyc); 

  // Gounaris-Sakurai
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  unsigned int sp, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int cyc); 
  // LASS 
  ResonancePdf (string name,
                          Variable* ar,
                          Variable* ai,
			  Variable* mass,
			  unsigned int sp,
                          Variable* width,
                          unsigned int cyc);

  // Nonresonant constructor
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai);  

  // Gaussian constructor
  ResonancePdf (string name,
			  Variable* ar, 
			  Variable* ai,
			  Variable* mean, 
			  Variable* sigma,
			  unsigned int cyc);

private:
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}

  Variable* amp_real;
  Variable* amp_imag;
/*
  Variable* mass;
  Variable* width;
  unsigned int spin;
  unsigned int cyclic_index;
  unsigned int eval_type;
  unsigned int resonance_type; 
*/
};

#endif

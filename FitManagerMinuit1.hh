#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "TMinuit.hh" 
#include <sys/time.h>
#include <sys/times.h>
extern PdfBase* pdfPointer; 
extern int numPars; 
#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag); 

class FitManager { 
public:
  FitManager (PdfBase* dat);
  ~FitManager ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void setupMinuit ();
  void runMigrad (); 
  void fit (); 
  TMinuit* getMinuitObject () {return minuit;} 
  void getMinuitValues () const;
  TMinuit* minuit; 
private:
  double overrideCallLimit; 
};

#endif 

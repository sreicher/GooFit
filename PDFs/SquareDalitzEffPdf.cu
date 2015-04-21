#include "SquareDalitzEffPdf.hh"

EXEC_TARGET fptype inPS(fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) {

  if (m12 < pow(mKS0 + mh1, 2)) return 0;
  if (m12 > pow(mD - mh2, 2)) return 0;

  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - mh1*mh1 + mKS0*mKS0) / sqrt(m12);
  fptype e3star = 0.5 * (mD*mD - m12 - mh2*mh2) / sqrt(m12);

  fptype minimum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - mKS0*mKS0) + sqrt(e3star*e3star - mh2*mh2), 2);
  if (m13 < minimum) return 0;
  fptype maximum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - mKS0*mKS0) - sqrt(e3star*e3star - mh2*mh2), 2);
  if (m13 > maximum) return 0;

  return 1;
}

EXEC_TARGET fptype mprime (fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) {
  // Helper function to calculate m'^2
  fptype m23 = mD*mD + mKS0*mKS0 + mh1*mh1 + mh2*mh2 - m12 - m13; 
  fptype rootPi = -2.*ATAN2(-1.0,0.0); // Pi

  if (m23 < 0) return -99;
  fptype tmp = ((2.0*(SQRT(m23) - (mh1 + mh2))/(mD - mKS0 - (mh1 + mh2))) - 1.0);
  if (isnan(tmp)) tmp = -99;
  return tmp;
}

EXEC_TARGET fptype thetaprime (fptype m12, fptype m13, fptype mD, fptype mKS0, fptype mh1, fptype mh2) {
  // Helper function to calculate theta'
  fptype m23 = mD*mD + mKS0*mKS0 + mh1*mh1 + mh2*mh2 - m12 - m13; 
  if (m23 < 0) return -99;

  fptype num = m23*( m12 - m13) + (mh2*mh2 - mh1*mh1)*(mD*mD - mKS0*mKS0);
  fptype denum = SQRT(((m23 - mh1*mh1 + mh2*mh2)*(m23 - mh1*mh1 + mh2*mh2) - 4*m23*mh2*mh2))*SQRT(((mD*mD - mKS0*mKS0 - m23)*(mD*mD - mKS0*mKS0 -m23) - 4*m23*mKS0*mKS0));
  fptype theta = -99 ;
  if (isnan(denum)) return -99;

  if (denum != 0.){
    theta = num/denum;
  }

  return theta;
}

EXEC_TARGET fptype device_SquareDalitzEff (fptype* evt, fptype* p, unsigned int* indices) {
  // Implementation of m23, mprime, thetaprime and also of efficiency tested
  // eff =  [0]*pow(theta,2) + [1]*m*pow(theta,2) +  exp([2] + [3]*dtime + ([4] + [5]*dtime)*m)

  // Define observables 
  fptype x = evt[indices[2 + indices[0] + 0]]; // m12   
  fptype y = evt[indices[2 + indices[0] + 1]]; // m13   
  //fptype z = evt[indices[2 + indices[0] + 2]]; // dtime   

  // Define constvals
  fptype mD   = p[indices[1]];
  fptype mKS0 = p[indices[2]];
  fptype mh1  = p[indices[3]];
  fptype mh2  = p[indices[4]];

  // Define coefficients
  fptype c0 = p[indices[5]];   
  fptype c1 = p[indices[6]];   
  fptype c2 = p[indices[7]];   
  fptype c3 = p[indices[8]];   
  fptype c4 = p[indices[9]];   
  fptype c5 = p[indices[10]];   
  // Check phase space
  if (inPS == 0) return 0;
  
  // Call helper functions
  fptype thetap = thetaprime(x,y,mD,mKS0,mh1,mh2); 
  if (thetap > 1. || thetap < -1.) return 0; 

  fptype m23 = mD*mD + mKS0*mKS0 + mh1*mh1 + mh2*mh2 - x - y; 
  if (m23 < 0) return 0;

  fptype ret = c0*m23*m23 + c1*m23 + c2*m23*thetap*thetap + c3*thetap*thetap + c4*thetap + c5; 

  //fptype mp     = mprime(x,y,mD,mKS0,mh1,mh2); 
  //if (mp > 1. || mp < -1.) return 0;
  //fptype rootPi = -2.*ATAN2(-1.0,0.0); // Pi
  // Calculate acutual m'^2 and theta'
  //mp     = POW(ACOS( mp )/rootPi,2);
  //thetap = ACOS( thetap )/rootPi;

  //fptype tmp = c2 + c3*mp + c4*POW(mp,2);
  //fptype ret = c0*POW(thetap,2) + c1*mp*POW(thetap,2) + tmp;
  //printf("Efficiency %f and m12 %f and m13 %f and dtime %f\n", ret, x, y, z);
  //printf("Efficiency %f and mp %f and thetap %f and dtime %f\n", ret, mp, thetap, z);
  
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_SquareDalitzEff = device_SquareDalitzEff; 

__host__ SquareDalitzEffPdf::SquareDalitzEffPdf (std::string n, vector<Variable*> obses, vector<Variable*> coeffs, vector<Variable*> constvals) 
  : GooPdf(0, n) 
{
  // Register observables - here m12, m13 and dtime
  for (unsigned int i = 0; i < obses.size(); ++i) {
    registerObservable(obses[i]);
  }

  std::vector<unsigned int> pindices;
  // Register constvals
  for (vector<Variable*>::iterator v = constvals.begin(); v != constvals.end(); ++v) {
    pindices.push_back(registerParameter(*v));
  }
  // Register coefficients
  for (vector<Variable*>::iterator c = coeffs.begin(); c != coeffs.end(); ++c) {
    pindices.push_back(registerParameter(*c));
  }

  GET_FUNCTION_ADDR(ptr_to_SquareDalitzEff);
  initialise(pindices);
}

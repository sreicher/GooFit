#include "ResonancePdf.hh" 

EXEC_TARGET fptype twoBodyCMmom (fptype rMassSq, fptype d1m, fptype d2m) {
  // For A -> B + C, calculate momentum of B and C in rest frame of A. 
  // PDG 38.16.

  fptype kin1 = 1 - POW(d1m+d2m, 2) / rMassSq;
  if (kin1 >= 0) kin1 = SQRT(kin1);
  else kin1 = 1;
  fptype kin2 = 1 - POW(d1m-d2m, 2) / rMassSq;
  if (kin2 >= 0) kin2 = SQRT(kin2);
  else kin2 = 1; 

  return 0.5*SQRT(rMassSq)*kin1*kin2; 
}


EXEC_TARGET fptype dampingFactorSquare (fptype cmmom, int spin, fptype mRadius) {
  fptype square = mRadius*mRadius*cmmom*cmmom;
  fptype dfsq = 1 + square; // This accounts for spin 1
  if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.   

  // Spin 3 and up not accounted for. 
  return dfsq; 
}

EXEC_TARGET fptype spinFactor (unsigned int spin, fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, fptype m12, fptype m13, fptype m23, unsigned int cyclic_index) {
  if (0 == spin) return 1; // Should not cause branching since every thread evaluates the same resonance at the same time. 
  /*
  // Copied from BdkDMixDalitzAmp
   
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass)); 
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass)); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 
    
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 

  // The above, collapsed into single tests where possible. 
  fptype _mA = (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass);
  fptype _mB = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 

  fptype _mAC = (PAIR_23 == cyclic_index ? m13 : m23);
  fptype _mBC = (PAIR_12 == cyclic_index ? m13 : m12);
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 
  */

  // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
  // Again, all threads should get the same branch. 
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug2Mass));
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass));
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12));
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13));
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

 /*
  // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
  // Again, all threads should get the same branch. 
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 
*/

  fptype massFactor = 1.0/_mAB;
  fptype sFactor = -1; 
  sFactor *= ((_mBC - _mAC) + (massFactor*(motherMass*motherMass - _mC*_mC)*(_mA*_mA-_mB*_mB)));
  if (2 == spin) {
    sFactor *= sFactor; 
    fptype extraterm = ((_mAB-(2*motherMass*motherMass)-(2*_mC*_mC))+massFactor*pow((motherMass*motherMass-_mC*_mC),2));
    extraterm *= ((_mAB-(2*_mA*_mA)-(2*_mB*_mB))+massFactor*pow((_mA*_mA-_mB*_mB),2));
    extraterm /= 3;
    sFactor -= extraterm;
  }
  return sFactor; 
}

EXEC_TARGET devcomplex<fptype> rhoF( fptype mCh, fptype m23 )
{

  fptype rhoSq = 1. - POW( mCh, 2 ) / m23;

  if ( rhoSq >= 0. )
    return devcomplex<fptype>(SQRT( rhoSq ),0);

  devcomplex<fptype> I;
  I.real = 0.0;
  I.imag = SQRT(- rhoSq);

  return I;
}

EXEC_TARGET devcomplex<fptype> rhoFourPiF( fptype mCh, fptype m23)
{
  if ( m23 > 1. )
    return rhoF( 4. * mCh, m23 );

  fptype m4 = POW( m23, 2 );
  fptype m6 = POW( m23, 3 );
  fptype m8 = POW( m23, 4 );

  fptype term = 0.;
  term += 0.00370909 / m4;
  term -= 0.111203 / m23;
  term += 1.2274;
  term -= 6.39017 * m23;
  term += 16.8358 * m4;
  term -= 21.8845 * m6;
  term += 11.3153 * m8;

  return rhoF( 4. * mCh, 1. ) * devcomplex<fptype>(term,0);
}

EXEC_TARGET devcomplex <fptype> Get_kMatrix (fptype m12, fptype m13, fptype m23, unsigned int* indices) {

  fptype Spr0                   = cudaArray[indices[2]];
  unsigned int term             = indices[3]; 
  unsigned int spin             = indices[4];
  unsigned int cyclic_index     = indices[5]; 

  assert(term >= 1 && term <= 6);

  // Particle Masses (notation: h = \eta)
  fptype pionMass   = 0.13957018; // PDG: (139.57018 \pm 0.00035) MeV
  fptype pionMassSq = pionMass*pionMass;
  fptype kMass      = 0.497614;   // PDG: (493.677 \pm 0.016) MeV
  fptype hMass      = 0.547853;   // PDG: (547.853 \pm 0.024) MeV
  fptype hprimeMass = 0.95778;    // PDG: (957.78 \pm 0.06) MeV
  
  // Invariant mass squared of resonant particles (pi+ pi- = m23)
  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));     

  // Pole Masses:           pipi    KK    4pi    hh   hhprime
  fptype poleMassesSq[5] = {0.651,1.2036,1.55817,1.21,1.82206};  // GeV
  for (int m=0; m<=4; m++){ poleMassesSq[m] *= poleMassesSq[m]; } // GeV^2

  // Other Fixed Single Variables
  fptype Ssc_0 = -3.92637; // GeV^2
  fptype sA_0 = -0.15;     // GeV^2
  fptype sA = 1.0;

  // Defining symmetric f^scattering Matrix of elements
  fptype Fsc_matrix[25];
  // non-zero terms;
  Fsc_matrix[0] = 0.23399;
  Fsc_matrix[1] = 0.15044;   Fsc_matrix[5] = 0.15044;
  Fsc_matrix[2] = -0.20545;  Fsc_matrix[10] = -0.20545;
  Fsc_matrix[3] = 0.32825;   Fsc_matrix[15] = 0.32825;
  Fsc_matrix[4] = 0.35412;   Fsc_matrix[20] = 0.35412;

  // all other terms are zero
  for (int i=1; i<=4; i++){
      for (int j=1; j<=4; j++){
          Fsc_matrix[5*i+j]=0.0;
      }
  }

  // Defining g^0 Matrix of elements
  fptype g0_matrix[25];

  g0_matrix[0]  = 0.22889;
  g0_matrix[1]  = -0.55377;
  g0_matrix[2]  = 0.0;
  g0_matrix[3]  = -0.39899;
  g0_matrix[4]  = -0.34639;

  g0_matrix[5]  = 0.94128;
  g0_matrix[6]  = 0.55095;
  g0_matrix[7]  = 0.0;
  g0_matrix[8]  = 0.39065;
  g0_matrix[9]  = 0.31503;

  g0_matrix[10] = 0.36856;
  g0_matrix[11] = 0.23888;
  g0_matrix[12] = 0.55639;
  g0_matrix[13] = 0.1834;
  g0_matrix[14] = 0.18681;

  g0_matrix[15] = 0.3365;
  g0_matrix[16] = 0.40907;
  g0_matrix[17] = 0.85679;
  g0_matrix[18] = 0.19906;
  g0_matrix[19] = -0.00984;

  g0_matrix[20] = 0.18171;
  g0_matrix[21] = -0.17558;
  g0_matrix[22] = -0.79658;
  g0_matrix[23] = -0.00355;
  g0_matrix[24] = 0.22358;

  // Computing the (real and symmetric) K-Matrix
  fptype AdlerZero = (1.0-sA_0)/(rMassSq-sA_0);
         AdlerZero *= (rMassSq - 0.5*sA*pionMassSq);

  devcomplex<fptype> iK_matrix[25]; // i*K-matrix
  for (int i=0; i<=24; i++){ iK_matrix[i] = 0.0;}

  fptype sum_term = 0.0, nonRes_term = 0.0; int irow,jcol,alpha; 
 
  for (irow=0; irow<=4; irow++){
      for (jcol=irow; jcol<=4; jcol++){

	  nonRes_term = Fsc_matrix[5*irow+jcol]*(1.0 - Ssc_0)/(rMassSq - Ssc_0); // Background non-resonant contribution

	  sum_term = 0.0;
	  for (alpha=0; alpha<=4; alpha++){		
	      sum_term += g0_matrix[5*alpha+irow]*g0_matrix[5*alpha+jcol] / (poleMassesSq[alpha]-rMassSq); 
	  }

	 iK_matrix[5*irow + jcol].imag = (sum_term + nonRes_term)*AdlerZero;
	 iK_matrix[5*jcol + irow].imag = (sum_term + nonRes_term)*AdlerZero; // K Matrix is symmetric

      } 
  }
 
  // Calculating pseudo propagator (I - iKp)^-1
  // Calculating Phase Spaces (p = \rho)
  devcomplex<fptype> rho[5];
  for (int d=0; d<=4; d++){rho[d] = 0.0;}
  
  rho[0]    = rhoF(2*pionMass, rMassSq);
  rho[1]      = rhoF(2*kMass, rMassSq);
  rho[2] = rhoFourPiF(pionMass, rMassSq);
  rho[3]      = rhoF(2*hMass, rMassSq);
  rho[4] = rhoF(hMass + hprimeMass, rMassSq);

  // Multiplying iK by diagonal -p matrix
  for (irow=0; irow<=4; irow++){
      for (jcol=0; jcol<=4; jcol++){
          iK_matrix[irow*5+jcol] *= rho[jcol];
          iK_matrix[irow*5+jcol] *= -1;
      }
  }

  // Adding Identity matrix to obtain I-iKp
  for (int d=0; d<=4; d++){
      iK_matrix[d*5+d].real = 1.0 + iK_matrix[d*5+d].real;
  }


  //Explicitly Defining Matrix to be Inverted
  devcomplex <fptype> n11 = iK_matrix[0];
  devcomplex <fptype> n12 = iK_matrix[1];
  devcomplex <fptype> n13 = iK_matrix[2];
  devcomplex <fptype> n14 = iK_matrix[3];
  devcomplex <fptype> n15 = iK_matrix[4];

  devcomplex <fptype> n21 = iK_matrix[5];
  devcomplex <fptype> n22 = iK_matrix[6];
  devcomplex <fptype> n23 = iK_matrix[7];
  devcomplex <fptype> n24 = iK_matrix[8];
  devcomplex <fptype> n25 = iK_matrix[9];

  devcomplex <fptype> n31 = iK_matrix[10];
  devcomplex <fptype> n32 = iK_matrix[11];
  devcomplex <fptype> n33 = iK_matrix[12];
  devcomplex <fptype> n34 = iK_matrix[13];
  devcomplex <fptype> n35 = iK_matrix[14];

  devcomplex <fptype> n41 = iK_matrix[15];
  devcomplex <fptype> n42 = iK_matrix[16];
  devcomplex <fptype> n43 = iK_matrix[17];
  devcomplex <fptype> n44 = iK_matrix[18];
  devcomplex <fptype> n45 = iK_matrix[19];

  devcomplex <fptype> n51 = iK_matrix[20];
  devcomplex <fptype> n52 = iK_matrix[21];
  devcomplex <fptype> n53 = iK_matrix[22];
  devcomplex <fptype> n54 = iK_matrix[23];
  devcomplex <fptype> n55 = iK_matrix[25];

// Computing elements of the first row of (I-iKp)^-1
// Formulae for inverted matrix elements obtained from Maple 

  devcomplex <fptype> inv11(0.0,0.0);
  inv11 += (n22*n33*n44*n55-n22*n33*n45*n54-n22*n34*n43*n55+n22*n34*n45*n53+n22*n35*n43*n54-n22*n35*n44*n53-n23*n32*n44*n55+n23*n32*n45*n54+n23*n34*n42*n55-n23*n34*n45*n52-n23*n35*n42*n54+n23*n35*n44*n52+n24*n32*n43*n55-n24*n32*n45*n53-n24*n33*n42*n55+n24*n33*n45*n52+n24*n35*n42*n53-n24*n35*n43*n52-n25*n32*n43*n54+n25*n32*n44*n53+n25*n33*n42*n54-n25*n33*n44*n52-n25*n34*n42*n53+n25*n34*n43*n52)/(n11*n22*n33*n44*n55-n11*n22*n33*n45*n54-n11*n22*n34*n43*n55+n11*n22*n34*n45*n53+n11*n22*n35*n43*n54-n11*n22*n35*n44*n53-n11*n23*n32*n44*n55+n11*n23*n32*n45*n54+n11*n23*n34*n42*n55-n11*n23*n34*n45*n52-n11*n23*n35*n42*n54+n11*n23*n35*n44*n52+n11*n24*n32*n43*n55-n11*n24*n32*n45*n53-n11*n24*n33*n42*n55+n11*n24*n33*n45*n52+n11*n24*n35*n42*n53-n11*n24*n35*n43*n52-n11*n25*n32*n43*n54+n11*n25*n32*n44*n53+n11*n25*n33*n42*n54-n11*n25*n33*n44*n52-n11*n25*n34*n42*n53+n11*n25*n34*n43*n52-n12*n21*n33*n44*n55+n12*n21*n33*n45*n54+n12*n21*n34*n43*n55-n12*n21*n34*n45*n53-n12*n21*n35*n43*n54+n12*n21*n35*n44*n53+n12*n23*n31*n44*n55-n12*n23*n31*n45*n54-n12*n23*n34*n41*n55+n12*n23*n34*n45*n51+n12*n23*n35*n41*n54-n12*n23*n35*n44*n51-n12*n24*n31*n43*n55+n12*n24*n31*n45*n53+n12*n24*n33*n41*n55-n12*n24*n33*n45*n51-n12*n24*n35*n41*n53+n12*n24*n35*n43*n51+n12*n25*n31*n43*n54-n12*n25*n31*n44*n53-n12*n25*n33*n41*n54+n12*n25*n33*n44*n51+n12*n25*n34*n41*n53-n12*n25*n34*n43*n51+n13*n21*n32*n44*n55-n13*n21*n32*n45*n54-n13*n21*n34*n42*n55+n13*n21*n34*n45*n52+n13*n21*n35*n42*n54-n13*n21*n35*n44*n52-n13*n22*n31*n44*n55+n13*n22*n31*n45*n54+n13*n22*n34*n41*n55-n13*n22*n34*n45*n51-n13*n22*n35*n41*n54+n13*n22*n35*n44*n51+n13*n24*n31*n42*n55-n13*n24*n31*n45*n52-n13*n24*n32*n41*n55+n13*n24*n32*n45*n51+n13*n24*n35*n41*n52-n13*n24*n35*n42*n51-n13*n25*n31*n42*n54+n13*n25*n31*n44*n52+n13*n25*n32*n41*n54-n13*n25*n32*n44*n51-n13*n25*n34*n41*n52+n13*n25*n34*n42*n51-n14*n21*n32*n43*n55+n14*n21*n32*n45*n53+n14*n21*n33*n42*n55-n14*n21*n33*n45*n52-n14*n21*n35*n42*n53+n14*n21*n35*n43*n52+n14*n22*n31*n43*n55-n14*n22*n31*n45*n53-n14*n22*n33*n41*n55+n14*n22*n33*n45*n51+n14*n22*n35*n41*n53-n14*n22*n35*n43*n51-n14*n23*n31*n42*n55+n14*n23*n31*n45*n52+n14*n23*n32*n41*n55-n14*n23*n32*n45*n51-n14*n23*n35*n41*n52+n14*n23*n35*n42*n51+n14*n25*n31*n42*n53-n14*n25*n31*n43*n52-n14*n25*n32*n41*n53+n14*n25*n32*n43*n51+n14*n25*n33*n41*n52-n14*n25*n33*n42*n51+n15*n21*n32*n43*n54-n15*n21*n32*n44*n53-n15*n21*n33*n42*n54+n15*n21*n33*n44*n52+n15*n21*n34*n42*n53-n15*n21*n34*n43*n52-n15*n22*n31*n43*n54+n15*n22*n31*n44*n53+n15*n22*n33*n41*n54-n15*n22*n33*n44*n51-n15*n22*n34*n41*n53+n15*n22*n34*n43*n51+n15*n23*n31*n42*n54-n15*n23*n31*n44*n52-n15*n23*n32*n41*n54+n15*n23*n32*n44*n51+n15*n23*n34*n41*n52-n15*n23*n34*n42*n51-n15*n24*n31*n42*n53+n15*n24*n31*n43*n52+n15*n24*n32*n41*n53-n15*n24*n32*n43*n51-n15*n24*n33*n41*n52+n15*n24*n33*n42*n51); 

  devcomplex <fptype> inv12(0.0,0.0);
  inv12 += -1*(n12*n33*n44*n55-n12*n33*n45*n54-n12*n34*n43*n55+n12*n34*n45*n53+n12*n35*n43*n54-n12*n35*n44*n53-n13*n32*n44*n55+n13*n32*n45*n54+n13*n34*n42*n55-n13*n34*n45*n52-n13*n35*n42*n54+n13*n35*n44*n52+n14*n32*n43*n55-n14*n32*n45*n53-n14*n33*n42*n55+n14*n33*n45*n52+n14*n35*n42*n53-n14*n35*n43*n52-n15*n32*n43*n54+n15*n32*n44*n53+n15*n33*n42*n54-n15*n33*n44*n52-n15*n34*n42*n53+n15*n34*n43*n52)/(n11*n22*n33*n44*n55-n11*n22*n33*n45*n54-n11*n22*n34*n43*n55+n11*n22*n34*n45*n53+n11*n22*n35*n43*n54-n11*n22*n35*n44*n53-n11*n23*n32*n44*n55+n11*n23*n32*n45*n54+n11*n23*n34*n42*n55-n11*n23*n34*n45*n52-n11*n23*n35*n42*n54+n11*n23*n35*n44*n52+n11*n24*n32*n43*n55-n11*n24*n32*n45*n53-n11*n24*n33*n42*n55+n11*n24*n33*n45*n52+n11*n24*n35*n42*n53-n11*n24*n35*n43*n52-n11*n25*n32*n43*n54+n11*n25*n32*n44*n53+n11*n25*n33*n42*n54-n11*n25*n33*n44*n52-n11*n25*n34*n42*n53+n11*n25*n34*n43*n52-n12*n21*n33*n44*n55+n12*n21*n33*n45*n54+n12*n21*n34*n43*n55-n12*n21*n34*n45*n53-n12*n21*n35*n43*n54+n12*n21*n35*n44*n53+n12*n23*n31*n44*n55-n12*n23*n31*n45*n54-n12*n23*n34*n41*n55+n12*n23*n34*n45*n51+n12*n23*n35*n41*n54-n12*n23*n35*n44*n51-n12*n24*n31*n43*n55+n12*n24*n31*n45*n53+n12*n24*n33*n41*n55-n12*n24*n33*n45*n51-n12*n24*n35*n41*n53+n12*n24*n35*n43*n51+n12*n25*n31*n43*n54-n12*n25*n31*n44*n53-n12*n25*n33*n41*n54+n12*n25*n33*n44*n51+n12*n25*n34*n41*n53-n12*n25*n34*n43*n51+n13*n21*n32*n44*n55-n13*n21*n32*n45*n54-n13*n21*n34*n42*n55+n13*n21*n34*n45*n52+n13*n21*n35*n42*n54-n13*n21*n35*n44*n52-n13*n22*n31*n44*n55+n13*n22*n31*n45*n54+n13*n22*n34*n41*n55-n13*n22*n34*n45*n51-n13*n22*n35*n41*n54+n13*n22*n35*n44*n51+n13*n24*n31*n42*n55-n13*n24*n31*n45*n52-n13*n24*n32*n41*n55+n13*n24*n32*n45*n51+n13*n24*n35*n41*n52-n13*n24*n35*n42*n51-n13*n25*n31*n42*n54+n13*n25*n31*n44*n52+n13*n25*n32*n41*n54-n13*n25*n32*n44*n51-n13*n25*n34*n41*n52+n13*n25*n34*n42*n51-n14*n21*n32*n43*n55+n14*n21*n32*n45*n53+n14*n21*n33*n42*n55-n14*n21*n33*n45*n52-n14*n21*n35*n42*n53+n14*n21*n35*n43*n52+n14*n22*n31*n43*n55-n14*n22*n31*n45*n53-n14*n22*n33*n41*n55+n14*n22*n33*n45*n51+n14*n22*n35*n41*n53-n14*n22*n35*n43*n51-n14*n23*n31*n42*n55+n14*n23*n31*n45*n52+n14*n23*n32*n41*n55-n14*n23*n32*n45*n51-n14*n23*n35*n41*n52+n14*n23*n35*n42*n51+n14*n25*n31*n42*n53-n14*n25*n31*n43*n52-n14*n25*n32*n41*n53+n14*n25*n32*n43*n51+n14*n25*n33*n41*n52-n14*n25*n33*n42*n51+n15*n21*n32*n43*n54-n15*n21*n32*n44*n53-n15*n21*n33*n42*n54+n15*n21*n33*n44*n52+n15*n21*n34*n42*n53-n15*n21*n34*n43*n52-n15*n22*n31*n43*n54+n15*n22*n31*n44*n53+n15*n22*n33*n41*n54-n15*n22*n33*n44*n51-n15*n22*n34*n41*n53+n15*n22*n34*n43*n51+n15*n23*n31*n42*n54-n15*n23*n31*n44*n52-n15*n23*n32*n41*n54+n15*n23*n32*n44*n51+n15*n23*n34*n41*n52-n15*n23*n34*n42*n51-n15*n24*n31*n42*n53+n15*n24*n31*n43*n52+n15*n24*n32*n41*n53-n15*n24*n32*n43*n51-n15*n24*n33*n41*n52+n15*n24*n33*n42*n51); 

  devcomplex <fptype> inv13(0.0,0.0);
  inv13 += (n12*n23*n44*n55-n12*n23*n45*n54-n12*n24*n43*n55+n12*n24*n45*n53+n12*n25*n43*n54-n12*n25*n44*n53-n13*n22*n44*n55+n13*n22*n45*n54+n13*n24*n42*n55-n13*n24*n45*n52-n13*n25*n42*n54+n13*n25*n44*n52+n14*n22*n43*n55-n14*n22*n45*n53-n14*n23*n42*n55+n14*n23*n45*n52+n14*n25*n42*n53-n14*n25*n43*n52-n15*n22*n43*n54+n15*n22*n44*n53+n15*n23*n42*n54-n15*n23*n44*n52-n15*n24*n42*n53+n15*n24*n43*n52)/(n11*n22*n33*n44*n55-n11*n22*n33*n45*n54-n11*n22*n34*n43*n55+n11*n22*n34*n45*n53+n11*n22*n35*n43*n54-n11*n22*n35*n44*n53-n11*n23*n32*n44*n55+n11*n23*n32*n45*n54+n11*n23*n34*n42*n55-n11*n23*n34*n45*n52-n11*n23*n35*n42*n54+n11*n23*n35*n44*n52+n11*n24*n32*n43*n55-n11*n24*n32*n45*n53-n11*n24*n33*n42*n55+n11*n24*n33*n45*n52+n11*n24*n35*n42*n53-n11*n24*n35*n43*n52-n11*n25*n32*n43*n54+n11*n25*n32*n44*n53+n11*n25*n33*n42*n54-n11*n25*n33*n44*n52-n11*n25*n34*n42*n53+n11*n25*n34*n43*n52-n12*n21*n33*n44*n55+n12*n21*n33*n45*n54+n12*n21*n34*n43*n55-n12*n21*n34*n45*n53-n12*n21*n35*n43*n54+n12*n21*n35*n44*n53+n12*n23*n31*n44*n55-n12*n23*n31*n45*n54-n12*n23*n34*n41*n55+n12*n23*n34*n45*n51+n12*n23*n35*n41*n54-n12*n23*n35*n44*n51-n12*n24*n31*n43*n55+n12*n24*n31*n45*n53+n12*n24*n33*n41*n55-n12*n24*n33*n45*n51-n12*n24*n35*n41*n53+n12*n24*n35*n43*n51+n12*n25*n31*n43*n54-n12*n25*n31*n44*n53-n12*n25*n33*n41*n54+n12*n25*n33*n44*n51+n12*n25*n34*n41*n53-n12*n25*n34*n43*n51+n13*n21*n32*n44*n55-n13*n21*n32*n45*n54-n13*n21*n34*n42*n55+n13*n21*n34*n45*n52+n13*n21*n35*n42*n54-n13*n21*n35*n44*n52-n13*n22*n31*n44*n55+n13*n22*n31*n45*n54+n13*n22*n34*n41*n55-n13*n22*n34*n45*n51-n13*n22*n35*n41*n54+n13*n22*n35*n44*n51+n13*n24*n31*n42*n55-n13*n24*n31*n45*n52-n13*n24*n32*n41*n55+n13*n24*n32*n45*n51+n13*n24*n35*n41*n52-n13*n24*n35*n42*n51-n13*n25*n31*n42*n54+n13*n25*n31*n44*n52+n13*n25*n32*n41*n54-n13*n25*n32*n44*n51-n13*n25*n34*n41*n52+n13*n25*n34*n42*n51-n14*n21*n32*n43*n55+n14*n21*n32*n45*n53+n14*n21*n33*n42*n55-n14*n21*n33*n45*n52-n14*n21*n35*n42*n53+n14*n21*n35*n43*n52+n14*n22*n31*n43*n55-n14*n22*n31*n45*n53-n14*n22*n33*n41*n55+n14*n22*n33*n45*n51+n14*n22*n35*n41*n53-n14*n22*n35*n43*n51-n14*n23*n31*n42*n55+n14*n23*n31*n45*n52+n14*n23*n32*n41*n55-n14*n23*n32*n45*n51-n14*n23*n35*n41*n52+n14*n23*n35*n42*n51+n14*n25*n31*n42*n53-n14*n25*n31*n43*n52-n14*n25*n32*n41*n53+n14*n25*n32*n43*n51+n14*n25*n33*n41*n52-n14*n25*n33*n42*n51+n15*n21*n32*n43*n54-n15*n21*n32*n44*n53-n15*n21*n33*n42*n54+n15*n21*n33*n44*n52+n15*n21*n34*n42*n53-n15*n21*n34*n43*n52-n15*n22*n31*n43*n54+n15*n22*n31*n44*n53+n15*n22*n33*n41*n54-n15*n22*n33*n44*n51-n15*n22*n34*n41*n53+n15*n22*n34*n43*n51+n15*n23*n31*n42*n54-n15*n23*n31*n44*n52-n15*n23*n32*n41*n54+n15*n23*n32*n44*n51+n15*n23*n34*n41*n52-n15*n23*n34*n42*n51-n15*n24*n31*n42*n53+n15*n24*n31*n43*n52+n15*n24*n32*n41*n53-n15*n24*n32*n43*n51-n15*n24*n33*n41*n52+n15*n24*n33*n42*n51);

  devcomplex <fptype> inv14(0.0,0.0);
  inv14 += -1*(n12*n23*n34*n55-n12*n23*n35*n54-n12*n24*n33*n55+n12*n24*n35*n53+n12*n25*n33*n54-n12*n25*n34*n53-n13*n22*n34*n55+n13*n22*n35*n54+n13*n24*n32*n55-n13*n24*n35*n52-n13*n25*n32*n54+n13*n25*n34*n52+n14*n22*n33*n55-n14*n22*n35*n53-n14*n23*n32*n55+n14*n23*n35*n52+n14*n25*n32*n53-n14*n25*n33*n52-n15*n22*n33*n54+n15*n22*n34*n53+n15*n23*n32*n54-n15*n23*n34*n52-n15*n24*n32*n53+n15*n24*n33*n52)/(n11*n22*n33*n44*n55-n11*n22*n33*n45*n54-n11*n22*n34*n43*n55+n11*n22*n34*n45*n53+n11*n22*n35*n43*n54-n11*n22*n35*n44*n53-n11*n23*n32*n44*n55+n11*n23*n32*n45*n54+n11*n23*n34*n42*n55-n11*n23*n34*n45*n52-n11*n23*n35*n42*n54+n11*n23*n35*n44*n52+n11*n24*n32*n43*n55-n11*n24*n32*n45*n53-n11*n24*n33*n42*n55+n11*n24*n33*n45*n52+n11*n24*n35*n42*n53-n11*n24*n35*n43*n52-n11*n25*n32*n43*n54+n11*n25*n32*n44*n53+n11*n25*n33*n42*n54-n11*n25*n33*n44*n52-n11*n25*n34*n42*n53+n11*n25*n34*n43*n52-n12*n21*n33*n44*n55+n12*n21*n33*n45*n54+n12*n21*n34*n43*n55-n12*n21*n34*n45*n53-n12*n21*n35*n43*n54+n12*n21*n35*n44*n53+n12*n23*n31*n44*n55-n12*n23*n31*n45*n54-n12*n23*n34*n41*n55+n12*n23*n34*n45*n51+n12*n23*n35*n41*n54-n12*n23*n35*n44*n51-n12*n24*n31*n43*n55+n12*n24*n31*n45*n53+n12*n24*n33*n41*n55-n12*n24*n33*n45*n51-n12*n24*n35*n41*n53+n12*n24*n35*n43*n51+n12*n25*n31*n43*n54-n12*n25*n31*n44*n53-n12*n25*n33*n41*n54+n12*n25*n33*n44*n51+n12*n25*n34*n41*n53-n12*n25*n34*n43*n51+n13*n21*n32*n44*n55-n13*n21*n32*n45*n54-n13*n21*n34*n42*n55+n13*n21*n34*n45*n52+n13*n21*n35*n42*n54-n13*n21*n35*n44*n52-n13*n22*n31*n44*n55+n13*n22*n31*n45*n54+n13*n22*n34*n41*n55-n13*n22*n34*n45*n51-n13*n22*n35*n41*n54+n13*n22*n35*n44*n51+n13*n24*n31*n42*n55-n13*n24*n31*n45*n52-n13*n24*n32*n41*n55+n13*n24*n32*n45*n51+n13*n24*n35*n41*n52-n13*n24*n35*n42*n51-n13*n25*n31*n42*n54+n13*n25*n31*n44*n52+n13*n25*n32*n41*n54-n13*n25*n32*n44*n51-n13*n25*n34*n41*n52+n13*n25*n34*n42*n51-n14*n21*n32*n43*n55+n14*n21*n32*n45*n53+n14*n21*n33*n42*n55-n14*n21*n33*n45*n52-n14*n21*n35*n42*n53+n14*n21*n35*n43*n52+n14*n22*n31*n43*n55-n14*n22*n31*n45*n53-n14*n22*n33*n41*n55+n14*n22*n33*n45*n51+n14*n22*n35*n41*n53-n14*n22*n35*n43*n51-n14*n23*n31*n42*n55+n14*n23*n31*n45*n52+n14*n23*n32*n41*n55-n14*n23*n32*n45*n51-n14*n23*n35*n41*n52+n14*n23*n35*n42*n51+n14*n25*n31*n42*n53-n14*n25*n31*n43*n52-n14*n25*n32*n41*n53+n14*n25*n32*n43*n51+n14*n25*n33*n41*n52-n14*n25*n33*n42*n51+n15*n21*n32*n43*n54-n15*n21*n32*n44*n53-n15*n21*n33*n42*n54+n15*n21*n33*n44*n52+n15*n21*n34*n42*n53-n15*n21*n34*n43*n52-n15*n22*n31*n43*n54+n15*n22*n31*n44*n53+n15*n22*n33*n41*n54-n15*n22*n33*n44*n51-n15*n22*n34*n41*n53+n15*n22*n34*n43*n51+n15*n23*n31*n42*n54-n15*n23*n31*n44*n52-n15*n23*n32*n41*n54+n15*n23*n32*n44*n51+n15*n23*n34*n41*n52-n15*n23*n34*n42*n51-n15*n24*n31*n42*n53+n15*n24*n31*n43*n52+n15*n24*n32*n41*n53-n15*n24*n32*n43*n51-n15*n24*n33*n41*n52+n15*n24*n33*n42*n51);

  devcomplex <fptype> inv15(0.0,0.0);
  inv15 += (n12*n23*n34*n45-n12*n23*n35*n44-n12*n24*n33*n45+n12*n24*n35*n43+n12*n25*n33*n44-n12*n25*n34*n43-n13*n22*n34*n45+n13*n22*n35*n44+n13*n24*n32*n45-n13*n24*n35*n42-n13*n25*n32*n44+n13*n25*n34*n42+n14*n22*n33*n45-n14*n22*n35*n43-n14*n23*n32*n45+n14*n23*n35*n42+n14*n25*n32*n43-n14*n25*n33*n42-n15*n22*n33*n44+n15*n22*n34*n43+n15*n23*n32*n44-n15*n23*n34*n42-n15*n24*n32*n43+n15*n24*n33*n42)/(n11*n22*n33*n44*n55-n11*n22*n33*n45*n54-n11*n22*n34*n43*n55+n11*n22*n34*n45*n53+n11*n22*n35*n43*n54-n11*n22*n35*n44*n53-n11*n23*n32*n44*n55+n11*n23*n32*n45*n54+n11*n23*n34*n42*n55-n11*n23*n34*n45*n52-n11*n23*n35*n42*n54+n11*n23*n35*n44*n52+n11*n24*n32*n43*n55-n11*n24*n32*n45*n53-n11*n24*n33*n42*n55+n11*n24*n33*n45*n52+n11*n24*n35*n42*n53-n11*n24*n35*n43*n52-n11*n25*n32*n43*n54+n11*n25*n32*n44*n53+n11*n25*n33*n42*n54-n11*n25*n33*n44*n52-n11*n25*n34*n42*n53+n11*n25*n34*n43*n52-n12*n21*n33*n44*n55+n12*n21*n33*n45*n54+n12*n21*n34*n43*n55-n12*n21*n34*n45*n53-n12*n21*n35*n43*n54+n12*n21*n35*n44*n53+n12*n23*n31*n44*n55-n12*n23*n31*n45*n54-n12*n23*n34*n41*n55+n12*n23*n34*n45*n51+n12*n23*n35*n41*n54-n12*n23*n35*n44*n51-n12*n24*n31*n43*n55+n12*n24*n31*n45*n53+n12*n24*n33*n41*n55-n12*n24*n33*n45*n51-n12*n24*n35*n41*n53+n12*n24*n35*n43*n51+n12*n25*n31*n43*n54-n12*n25*n31*n44*n53-n12*n25*n33*n41*n54+n12*n25*n33*n44*n51+n12*n25*n34*n41*n53-n12*n25*n34*n43*n51+n13*n21*n32*n44*n55-n13*n21*n32*n45*n54-n13*n21*n34*n42*n55+n13*n21*n34*n45*n52+n13*n21*n35*n42*n54-n13*n21*n35*n44*n52-n13*n22*n31*n44*n55+n13*n22*n31*n45*n54+n13*n22*n34*n41*n55-n13*n22*n34*n45*n51-n13*n22*n35*n41*n54+n13*n22*n35*n44*n51+n13*n24*n31*n42*n55-n13*n24*n31*n45*n52-n13*n24*n32*n41*n55+n13*n24*n32*n45*n51+n13*n24*n35*n41*n52-n13*n24*n35*n42*n51-n13*n25*n31*n42*n54+n13*n25*n31*n44*n52+n13*n25*n32*n41*n54-n13*n25*n32*n44*n51-n13*n25*n34*n41*n52+n13*n25*n34*n42*n51-n14*n21*n32*n43*n55+n14*n21*n32*n45*n53+n14*n21*n33*n42*n55-n14*n21*n33*n45*n52-n14*n21*n35*n42*n53+n14*n21*n35*n43*n52+n14*n22*n31*n43*n55-n14*n22*n31*n45*n53-n14*n22*n33*n41*n55+n14*n22*n33*n45*n51+n14*n22*n35*n41*n53-n14*n22*n35*n43*n51-n14*n23*n31*n42*n55+n14*n23*n31*n45*n52+n14*n23*n32*n41*n55-n14*n23*n32*n45*n51-n14*n23*n35*n41*n52+n14*n23*n35*n42*n51+n14*n25*n31*n42*n53-n14*n25*n31*n43*n52-n14*n25*n32*n41*n53+n14*n25*n32*n43*n51+n14*n25*n33*n41*n52-n14*n25*n33*n42*n51+n15*n21*n32*n43*n54-n15*n21*n32*n44*n53-n15*n21*n33*n42*n54+n15*n21*n33*n44*n52+n15*n21*n34*n42*n53-n15*n21*n34*n43*n52-n15*n22*n31*n43*n54+n15*n22*n31*n44*n53+n15*n22*n33*n41*n54-n15*n22*n33*n44*n51-n15*n22*n34*n41*n53+n15*n22*n34*n43*n51+n15*n23*n31*n42*n54-n15*n23*n31*n44*n52-n15*n23*n32*n41*n54+n15*n23*n32*n44*n51+n15*n23*n34*n41*n52-n15*n23*n34*n42*n51-n15*n24*n31*n42*n53+n15*n24*n31*n43*n52+n15*n24*n32*n41*n53-n15*n24*n32*n43*n51-n15*n24*n33*n41*n52+n15*n24*n33*n42*n51); 

  // Computing the F0 element
  fptype temp = (1.0-Spr0)/(rMassSq-Spr0);
  devcomplex <fptype> F0(0.0,0.0);
  devcomplex <fptype> _fr12prod(1.87981, -0.628378);
  devcomplex <fptype> _fr13prod(4.3242, 2.75019);
  devcomplex <fptype> _fr14prod(3.22336, 0.271048);
  devcomplex <fptype> _fr15prod(0., 0.);

  if (term <= 5){ // beta factors
    for (int j = 0; j < 5; j++) {
      fptype tmp = 1./(poleMassesSq[term-1]-rMassSq);
      if (j == 0) F0 += tmp*inv11*g0_matrix[5*(term-1)+j]; 
      else if (j == 1) F0 += tmp*inv12*g0_matrix[5*(term-1)+j]; 
      else if (j == 2) F0 += tmp*inv13*g0_matrix[5*(term-1)+j]; 
      else if (j == 3) F0 += tmp*inv14*g0_matrix[5*(term-1)+j]; 
      else if (j == 4) F0 += tmp*inv15*g0_matrix[5*(term-1)+j]; 
    }
  }
  else { // fprod factors wrt to fprod11
    F0 += inv11;
    F0 += inv12*_fr12prod;
    F0 += inv13*_fr13prod;
    F0 += inv14*_fr14prod;
    F0 += inv15*_fr15prod;
    F0 *= temp;
  }

 return F0;

}

EXEC_TARGET devcomplex<fptype> plainBW (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int spin             = indices[4];
  unsigned int cyclic_index     = indices[5]; 

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, 
					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, 
					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius); 
  }  
 
  // RBW evaluation
  fptype A = (resmass - rMassSq); 
  fptype B = resmass*reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor / SQRT(rMassSq);
  //fptype C = (SQRT(resmass)*reswidth) / (A*A + B*B); 
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> ret(A*C, B*C); // Dropping F_D=1

  ret *= SQRT(frFactor); 
  fptype spinF = spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index); 
  ret *= spinF; 

  return ret; 
}

EXEC_TARGET devcomplex<fptype> gaussian (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  // indices[1] is unused constant index, for consistency with other function types. 
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int cyclic_index     = indices[4]; 

  // Notice sqrt - this function uses mass, not mass-squared like the other resonance types. 
  fptype massToUse = SQRT(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  massToUse -= resmass;
  massToUse /= reswidth;
  massToUse *= massToUse;
  fptype ret = EXP(-0.5*massToUse); 

  // Ignore factor 1/sqrt(2pi). 
  ret /= reswidth;

  return devcomplex<fptype>(ret, 0); 
}

EXEC_TARGET fptype hFun (fptype s, fptype daug2Mass, fptype daug3Mass) {
  // Last helper function
  const fptype _pi = 3.14159265359;
  fptype sm   = daug2Mass + daug3Mass;
  fptype SQRTs = sqrt(s);
  fptype k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);

  fptype val = ((2/_pi) * (k_s/SQRTs) * log( (SQRTs + 2*k_s)/(sm)));

  return val;
}

EXEC_TARGET fptype dh_dsFun (fptype s, fptype daug2Mass, fptype daug3Mass) {
  // Yet another helper function
  const fptype _pi = 3.14159265359;
  fptype k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);
  
  fptype val = (hFun(s, daug2Mass, daug3Mass) * (1.0/(8.0*pow(k_s, 2)) - 1.0/(2.0 * s)) + 1.0/(2.0* _pi*s));
  return val;
}


EXEC_TARGET fptype dFun (fptype s, fptype daug2Mass, fptype daug3Mass) {
  // Helper function used in Gronau-Sakurai
  const fptype _pi = 3.14159265359;
  fptype sm   = daug2Mass + daug3Mass;
  fptype sm24 = sm*sm/4.0;
  fptype m    = sqrt(s);
  fptype k_m2 = twoBodyCMmom(s, daug2Mass, daug3Mass);
 
  fptype val = 3.0/_pi * sm24/pow(k_m2, 2) * log((m + 2*k_m2)/sm) + m/(2*_pi*k_m2) - sm24*m/(_pi * pow(k_m2, 3));
  return val;
}

EXEC_TARGET fptype fsFun (fptype s, fptype m2, fptype gam, fptype daug2Mass, fptype daug3Mass) {
  // Another G-S helper function
   
  fptype k_s   = twoBodyCMmom(s,  daug2Mass, daug3Mass);
  fptype k_Am2 = twoBodyCMmom(m2, daug2Mass, daug3Mass);
   
  fptype f     = gam * m2 / POW(k_Am2, 3);
  f           *= (POW(k_s, 2) * (hFun(s, daug2Mass, daug3Mass) - hFun(m2, daug2Mass, daug3Mass)) + (m2 - s) * pow(k_Am2, 2) * dh_dsFun(m2, daug2Mass, daug3Mass));
 
  return f;
}

EXEC_TARGET devcomplex<fptype> gouSak (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int spin             = indices[4];
  unsigned int cyclic_index     = indices[5]; 

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  resmass *= resmass; 
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius); 
  }
  
  // Implement Gou-Sak:

  //fptype D = resmass*(1.0 + dFun(resmass, daug2Mass, daug3Mass) * reswidth/SQRT(resmass));
  fptype D = (1.0 + dFun(resmass, daug2Mass, daug3Mass) * reswidth/SQRT(resmass));
  fptype E = resmass - rMassSq + fsFun(rMassSq, resmass, reswidth, daug2Mass, daug3Mass);
  fptype F = SQRT(resmass) * reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor;

  D       /= (E*E + F*F);
  devcomplex<fptype> retur(D*E, D*F); // Dropping F_D=1
  retur *= SQRT(frFactor);
  retur *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);

  return retur; 
}


  EXEC_TARGET devcomplex<fptype> lass (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
//  fptype lass_a                 = cudaArray[indices[4]];
//  fptype lass_r                 = cudaArray[indices[5]];
//  fptype lass_B                 = cudaArray[indices[6]];
//  fptype lass_phiB              = cudaArray[indices[7]];
//  fptype lass_R                 = cudaArray[indices[8]];
//  fptype lass_phiR              = cudaArray[indices[9]];
  unsigned int spin             = indices[10];
  unsigned int cyclic_index     = indices[11];

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  resmass *= resmass;
  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2).
  
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_23 == cyclic_index ? daug3Mass : daug2Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_23 == cyclic_index ? daug3Mass : daug2Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius);
  }

  //Implement LASS:

  fptype q = measureDaughterMoms;
  fptype g = reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor / SQRT(rMassSq);
  //fptype g = reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor * ( SQRT(resmass) / SQRT(rMassSq)); // as in PDG
  fptype lass_a    = 0.22357;
  fptype lass_r    = -15.042;
  fptype lass_R    = 1; // ?
  fptype lass_phiR = 1.10644;
  fptype lass_B    = 0.614463;
  fptype lass_phiB = -0.0981907;

  // background phase motion
  fptype cot_deltaB = (1.0 / (lass_a*q)) + 0.5*lass_r*q;
  fptype qcot_deltaB = (1.0 / lass_a) + 0.5*lass_r*q*q;

  // calculate resonant part
  devcomplex<fptype> expi2deltaB = devcomplex<fptype>(qcot_deltaB,q)/devcomplex<fptype>(qcot_deltaB,-q);
  devcomplex<fptype>  resT = devcomplex<fptype>(cos(lass_phiR+2*lass_phiB),sin(lass_phiR+2*lass_phiB))*lass_R;

  devcomplex<fptype> prop = devcomplex<fptype>(1, 0)/devcomplex<fptype>(resmass-rMassSq, SQRT(resmass)*g);
  //devcomplex<fptype> prop = devcomplex<fptype>(1, 0)/devcomplex<fptype>(resmass-rMassSq, -SQRT(resmass)*g); // as in EvtGen and Papers
  //resT *= prop*(SQRT(resmass)*reswidth/nominalDaughterMoms)*expi2deltaB;
  resT *= prop*(resmass*reswidth/nominalDaughterMoms)*expi2deltaB;

  // calculate bkg part
  resT += devcomplex<fptype>(cos(lass_phiB),sin(lass_phiB))*lass_B*q*(cos(lass_phiB)+cot_deltaB*sin(lass_phiB))/devcomplex<fptype>(qcot_deltaB,-q);
  //resT += devcomplex<fptype>(cos(_phiB),sin(_phiB))*_B*(cos(_phiB)+cot_deltaB*sin(_phiB))*SQRT(rMassSq)/devcomplex<fptype>(qcot_deltaB,-q);

  resT *= SQRT(frFactor);
  resT *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);

  return resT;
}

EXEC_TARGET devcomplex<fptype> nonres (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  return devcomplex<fptype>(1, 0); 
}


EXEC_TARGET void getAmplitudeCoefficients (devcomplex<fptype> a1, devcomplex<fptype> a2, fptype& a1sq, fptype& a2sq, fptype& a1a2real, fptype& a1a2imag) {
  // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
  a1sq = a1.abs2();
  a2sq = a2.abs2();
  a1 *= conj(a2);
  a1a2real = a1.real;
  a1a2imag = a1.imag; 
}

 MEM_DEVICE resonance_function_ptr ptr_to_RBW = plainBW;
 MEM_DEVICE resonance_function_ptr ptr_to_GOUSAK = gouSak; 
 MEM_DEVICE resonance_function_ptr ptr_to_GAUSSIAN = gaussian;
 MEM_DEVICE resonance_function_ptr ptr_to_NONRES = nonres;
 MEM_DEVICE resonance_function_ptr ptr_to_LASS = lass;
 MEM_DEVICE resonance_function_ptr ptr_to_kMatrix = Get_kMatrix; 

 ResonancePdf::ResonancePdf (string name,
						Variable* ar,             
                                                Variable* ai,             
						Variable* Spr0,
                                                unsigned int term,
						unsigned int sp,         	
						unsigned int cyc)               
  :GooPdf(0,name)
  ,amp_real(ar)                                                         
  ,amp_imag(ai)
{
  vector <unsigned int> pindices; 
  pindices.push_back(0); 
  pindices.push_back(registerParameter(Spr0));
  pindices.push_back(term);
  pindices.push_back(sp);
  pindices.push_back(cyc); 


 
  GET_FUNCTION_ADDR(ptr_to_kMatrix);
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name, 
						Variable* ar, 
						Variable* ai, 
						Variable* mass, 
						Variable* width, 
						unsigned int sp, 
						unsigned int cyc) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex. 

  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc); 

  GET_FUNCTION_ADDR(ptr_to_RBW);
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name, 
						Variable* ar, 
						Variable* ai, 
						unsigned int sp, 
						Variable* mass, 
						Variable* width, 
						unsigned int cyc) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
{
  // Same as BW except for function pointed to. 
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc); 

  GET_FUNCTION_ADDR(ptr_to_GOUSAK);
  initialise(pindices); 
} 


ResonancePdf::ResonancePdf (string name,
                                                Variable* ar,
                                                Variable* ai,
						Variable* mass,
//						Variable* a,
//						Variable* r,
//						Variable* B,
//						Variable* phiB,
//						Variable* R,
//						Variable* phiR,
                                                unsigned int sp,
                                                Variable* width,
                                                unsigned int cyc)
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
{
  // Same as BW except for function pointed to.
  vector<unsigned int> pindices;
  pindices.push_back(0);
  pindices.push_back(registerParameter(mass));
  //pindices.push_back(registerParameter(a));
  //pindices.push_back(registerParameter(r));
  //pindices.push_back(registerParameter(B));
  //pindices.push_back(registerParameter(phiB));
  //pindices.push_back(registerParameter(R));
  //pindices.push_back(registerParameter(phiR));
  pindices.push_back(sp);
  pindices.push_back(registerParameter(width));
  pindices.push_back(cyc);

  GET_FUNCTION_ADDR(ptr_to_LASS);
  initialise(pindices);
}

ResonancePdf::ResonancePdf (string name, 
						Variable* ar, 
						Variable* ai) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  GET_FUNCTION_ADDR(ptr_to_NONRES);
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name,
						Variable* ar, 
						Variable* ai,
						Variable* mean, 
						Variable* sigma,
						unsigned int cyc) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma)); 
  pindices.push_back(cyc); 

  GET_FUNCTION_ADDR(ptr_to_GAUSSIAN);
  initialise(pindices); 

}



#include "StepPdf.hh"

EXEC_TARGET fptype device_Step (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]]; 
  fptype x0 = p[indices[1]];
  return (x > x0 ? 1 : 0);
}

MEM_DEVICE device_function_ptr ptr_to_Step = device_Step; 
device_function_ptr hptr_to_Step = device_Step; 

__host__ StepPdf::StepPdf (std::string n, Variable* _x, Variable* x0)
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(x0));
  MEMCPY_FROM_SYMBOL((void**) &host_fcn_ptr, ptr_to_Step, sizeof(void*), 0, cudaMemcpyDeviceToHost);
  initialise(pindices); 
}

__host__ fptype StepPdf::integrate (fptype lo, fptype hi) const {
  unsigned int* indices = host_indices+parameters; 
  fptype x0 = host_params[indices[1]];
  return (hi - x0);
}


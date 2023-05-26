#include "rns-gsw.h"

// We follow the methods defined in https://eprint.iacr.org/2017/996.pdf#page=27

void clwe_sample_uniform(uint64_t * out, uint64_t size, uint64_t p){
  generate_random_bytes(sizeof(uint64_t)*size, (uint8_t *) out);
  array_mod_switch_from_2k(out, out, p, p, size); 
  out[0] = 0;
  for (size_t i = 1; i < size; i++){
    out[0] = intel::hexl::SubUIntMod(out[0], out[i], p);
  }
}

void clwe_sample_gaussian(uint64_t * out, uint64_t size, double sigma){
  memset(out, 0, sizeof(uint64_t)*size);
  clwe_sample_gaussian_addto(out, size, sigma);
}

void clwe_sample_gaussian_addto(uint64_t * out, uint64_t size, double sigma){
  const uint64_t buffer_size = pow(sigma, 2)*size/2;
  uint64_t a[buffer_size], b[buffer_size];
  // Generate uniform dist. from 0 to 2^64
  generate_random_bytes(sizeof(uint64_t)*buffer_size, (uint8_t *) a); 
  generate_random_bytes(sizeof(uint64_t)*buffer_size, (uint8_t *) b); 
  // Mod switch to get an uniform dist. from 0 to size.
  array_mod_switch_from_2k(a, a, size, size, buffer_size); 
  array_mod_switch_from_2k(b, b, size, size, buffer_size); 

  // e = sum (X^ai - X^bi)
  for (size_t i = 0; i < buffer_size; i++){
    out[a[i]] += 1;
    out[b[i]] -= 1;
  }
}
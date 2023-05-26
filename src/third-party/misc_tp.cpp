#include "rns-gsw.h"

// Function adapted from OpenFHE
// BSD 2-Clause License https://www.openfhe.org/
// the main rounding operation used in ModSwitch (as described in Section 3 of
// https://eprint.iacr.org/2014/816) The idea is that Round(x) = 0.5 + Floor(x)
uint64_t RoundqQ(uint64_t v, uint64_t q, uint64_t Q) {
  return ((uint64_t) floor(0.5 + ((double) v) * ((double) q) / ((double) Q)) % q);
}

// Functions adapted from MOSFHET

#ifndef PORTABLE_BUILD
// TODO: add code src.
void generate_rnd_seed(uint64_t * p){
  if(0 == _rdrand64_step ((unsigned long long *) p) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[1])) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[2])) ||
     0 == _rdrand64_step ((unsigned long long *) &(p[3]))){
    printf("Random Generation Failed\n");
    return;
  }
}
#else 
void generate_rnd_seed(uint64_t * p){
  FILE *fp;
  fp = fopen("/dev/urandom", "r");
  fread(p, 1, 32, fp);
  fclose(fp);
}
#endif

#ifdef PORTABLE_BUILD
#include "sha3/fips202.c"
#else
#include "aes_rng.c"
// void aes_prng(uint8_t *output, uint64_t outlen, const uint8_t *input,  uint64_t inlen);
#endif

void get_rnd_from_hash(uint64_t amount, uint8_t * pointer){
  uint64_t rnd[4];
  generate_rnd_seed(rnd);
  #ifdef USE_SHAKE
  shake256(pointer, amount, (uint8_t *) rnd, 32);
  #else
  aes_prng(pointer, amount, (uint8_t *) rnd, 32);
  #endif
}

void get_rnd_from_buffer(uint64_t amount, uint8_t * pointer){
  static uint8_t buffer[1024];
  static uint64_t idx = 1024;
  if(amount > (1024 - idx)){
    idx = 0;
    get_rnd_from_hash(1024, buffer);
  }
  memcpy(pointer, buffer+idx, amount);
  idx += amount;
}

void generate_random_bytes(uint64_t amount, uint8_t * pointer){
  if(amount < 512) get_rnd_from_buffer(amount, pointer);
  else get_rnd_from_hash(amount, pointer);
}


double int2double(uint64_t x){
  return ((double) x)/18446744073709551616.0;
}


#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
double generate_normal_random(double sigma){
  uint64_t rnd[2];
  generate_random_bytes(16, (uint8_t *) rnd);
  return cos(2.*M_PI*int2double(rnd[0]))*sqrt(-2.*log(int2double(rnd[1])))*sigma;
}


// Mem alloc

uint64_t _glb_mem_count = 0;
// safe_malloc
// https://stackoverflow.com/questions/48043811/creating-a-function-to-check-if-malloc-succeeded
void * safe_malloc(size_t size){
  void *ptr = malloc(size);
  if (!ptr && (size > 0)) {
    perror("malloc failed!");
    exit(EXIT_FAILURE);
  }
  // memset(ptr, 0, size); 
  return ptr;
}

void * safe_aligned_malloc(size_t size){
  void * ptr;
  #ifdef AVX512_OPT
  int err = posix_memalign(&ptr, 64, size);
  #else
  int err = posix_memalign(&ptr, 32, size);
  #endif
  if (err || (!ptr && (size > 0))) {
    perror("aligned malloc failed!");
    exit(EXIT_FAILURE);
  }
  // memset(ptr, 0, size);
  return ptr;
}
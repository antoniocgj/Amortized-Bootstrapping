#include "rns-gsw.h"

void bit_rev_list(uint64_t * out, uint64_t * in, uint64_t n){
  assert(n == 512 || n == 1024 || n == 2048);
  for (size_t i = 0; i < n; i++){
    out[i] = in[bit_rev[n>>10][i]];
  }
}

void ntt_forward(uint64_t * out, uint64_t * in, uint64_t * ws, uint64_t Q, uint64_t n){
  assert(in != out);
  bit_rev_list(out, in, n);
  uint64_t m = 2;
  while(m <= n){
    const uint64_t w_m = ws[n/m];
    uint64_t w = 1;
    for (size_t j = 0; j < m/2; j++){
      // const uint64_t w = ws[j*n/m];
      for (size_t k = 0; k < n-1; k+=m){
        const uint64_t t = intel::hexl::MultiplyMod(w, out[k + j + m/2], Q);
        const uint64_t u = out[k + j];
        out[k + j] = intel::hexl::AddUIntMod(u, t, Q);
        out[k + j + m/2] = intel::hexl::SubUIntMod(u, t, Q);
      }
      w = intel::hexl::MultiplyMod(w, w_m, Q);
    }
    m *= 2;
  }
}

void calcRootsOfUnit(uint64_t * out, uint64_t min_RoU, uint64_t q, uint64_t n){
  uint64_t r = 1;
  for (size_t i = 0; i < n; i++){
    out[i] = r;
    r = intel::hexl::MultiplyMod(r, min_RoU, q);
  }
}

void entry_wise_prod(uint64_t * out, uint64_t * in1, uint64_t * in2, uint64_t n){
  for (size_t i = 0; i < n; i++){
    out[i] = in1[i] * in2[i];
  }
}

// negacyclic forward ntt
void nc_ntt_forward(uint64_t * out, uint64_t * in, uint64_t rou_2nth, uint64_t q, uint64_t n){
  uint64_t rou_vec[n], in2[n]; 
  // memcpy(in2, in, sizeof(*in)*n);
  calcRootsOfUnit(rou_vec, rou_2nth, q, n);
  entry_wise_prod(in2, in, rou_vec, n);
  const uint64_t rou = intel::hexl::MultiplyMod(rou_2nth, rou_2nth, q);
  calcRootsOfUnit(rou_vec, rou, q, n);
  ntt_forward(out, in2, rou_vec, q, n);
}



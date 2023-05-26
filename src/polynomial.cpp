#include "rns-gsw.h"

RNS_Polynomial polynomial_new_RNS_polynomial(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_Polynomial res;
  const uint64_t N = next_power_of_2(Q) << 1;
  res = (RNS_Polynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t **) safe_malloc(sizeof(uint64_t*) * l);
  for (size_t i = 0; i < l; i++){
    res->coeffs[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t) * N);
  }
  res->N = N;
  res->Q = Q;
  res->l = l;
  res->ntt = ntt;
  return res;
}

void polynomial_copy_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in){
  for (size_t i = 0; i < in->l; i++){
    memcpy(out->coeffs[i], in->coeffs[i], sizeof(uint64_t)*out->N);
  }
}

void polynomial_RNS_zero(RNS_Polynomial p){
  for (size_t i = 0; i < p->l; i++){
    memset(p->coeffs[i], 0, sizeof(uint64_t)*p->N);
  }
}

void free_RNS_polynomial(void * p){
  RNS_Polynomial pp = (RNS_Polynomial) p;
  for (size_t i = 0; i < pp->l; i++){
    free(pp->coeffs[i]);
  }
  free(pp->coeffs);
  free(pp);
}

RNS_Polynomial * polynomial_new_array_of_RNS_polynomials(uint64_t Q, uint64_t l, uint64_t size, intel::hexl::NTT ** ntt){
  RNS_Polynomial * res = (RNS_Polynomial *) safe_malloc(sizeof(RNS_Polynomial)*size);
  for (size_t i = 0; i < size; i++) res[i] = polynomial_new_RNS_polynomial(Q, l, ntt);
  return res;
}

// out = RNS(in)
// Assumes ||in||_inf < min(p)^2
void polynomial_to_RNS(RNS_Polynomial out, IntPolynomial in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->Q; j++){
      out->coeffs[i][j] = (p + in->coeffs[j])%p;
    }    
    memset(&out->coeffs[i][out->Q], 0, sizeof(uint64_t)*(out->N - out->Q));
    out->ntt[i]->ComputeForward(out->coeffs[i], out->coeffs[i], 1, 1);
  }
}

// todo: fix uniform distribution
void polynomial_gen_random_RNSc_polynomial(RNSc_Polynomial out){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    memset(out->coeffs[i], 0, sizeof(uint64_t)*out->N);
    clwe_sample_uniform(out->coeffs[i], out->Q, p);
  }
}

// reduce p mod (X^d - 1) and its coefficients mod q.
void polynomial_reduce_mod_Xd_minus_1_mod_q(uint64_t * p, int64_t d, uint64_t q, size_t size){
  for (int64_t i = size - 1; i >= d; i--){
    p[i - d] = intel::hexl::AddUIntMod(p[i - d], p[i], q);
    p[i] = 0;
  }
}

// reduce p mod (X^d + 1) and its coefficients mod q.
void polynomial_reduce_mod_Xd_plus_1_mod_q(uint64_t * p, int64_t d, uint64_t q, size_t size){
  for (int64_t i = size - 1; i >= d; i--){
    p[i - d] = intel::hexl::SubUIntMod(p[i - d], p[i], q);
    p[i] = 0;
  }
}

void polynomial_reduce_mod_Xd_minus_1(uint64_t * p, int64_t d, size_t size){
  for (int64_t i = size - 1; i >= d; i--){
    p[i - d] += p[i];
    p[i] = 0;
  }
}

// reduce the polynomial  mod (X^Q - 1) and its coefficients mod q.
void polynomial_RNSc_reduce_mod_XQ_minus_1(RNSc_Polynomial out){
  for (size_t i = 0; i < out->l; i++){
    polynomial_reduce_mod_Xd_minus_1_mod_q(out->coeffs[i], out->Q, out->ntt[i]->GetModulus(), out->N);
  }
}

/* out = in1*in2 mod (X^Q - 1) */
void polynomial_mul_RNS_polynomial(RNSc_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseMultMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q, 1);
    out->ntt[i]->ComputeInverse(out->coeffs[i], out->coeffs[i], 1, 1);
    polynomial_reduce_mod_Xd_minus_1_mod_q(out->coeffs[i], out->Q, q, out->N);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

/* out = in1 - in2 mod (X^Q - 1) */
void polynomial_sub_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseSubMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

/* out = in1 - in2 mod (X^Q - 1) */
void polynomial_sub_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseSubMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}


void polynomial_add_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    intel::hexl::EltwiseAddMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

void polynomial_RNSc_negate(RNSc_Polynomial out, RNSc_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->Q; j++){
      out->coeffs[i][j] = intel::hexl::SubUIntMod(0, in->coeffs[i][j], q);
    }
  }
}


void polynomial_RNSc_to_RNS(RNS_Polynomial out, RNSc_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    out->ntt[i]->ComputeForward(out->coeffs[i], in->coeffs[i], 1, 1);
  }
}

void polynomial_RNS_to_RNSc(RNSc_Polynomial out, RNS_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    out->ntt[i]->ComputeInverse(out->coeffs[i], in->coeffs[i], 1, 1);
  }
}

// out = in + e
// e <- Gaussian(0, sigma)
// Assumes sigma <= min(p)^2
void polynomial_RNSc_add_noise(RNSc_Polynomial out, RNSc_Polynomial in, double sigma){
  uint64_t * e = (uint64_t*) safe_aligned_malloc(sizeof(uint64_t)*out->Q);
  clwe_sample_gaussian(e, out->Q, sigma);
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->Q; j++){
      out->coeffs[i][j] = intel::hexl::AddUIntMod(in->coeffs[i][j], (p + e[j])%p, p);
    }
  }
  free(e);
}

/* Assumes q_i/q_j < 2 for all i, j*/
void polynomial_base_extend_RNSc(RNSc_Polynomial out, uint64_t * in, uint64_t p){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    memset(&out->coeffs[i][out->Q], 0, sizeof(uint64_t)*(out->N - out->Q));
    if(q >= p){ // copy
      memcpy(out->coeffs[i], in, sizeof(uint64_t)*out->Q);
    }else{ // reduce
      intel::hexl::EltwiseReduceMod(out->coeffs[i], in, out->Q, q, 2, 1);
    }
  }
}

// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc_wo_free(RNSc_Polynomial out){
  const uint64_t l = out->l, p = out->ntt[l-1]->GetModulus();
  const uint64_t Q = out->Q;
  for (size_t i = 0; i < l - 1; i++){
    const uint64_t q = out->ntt[i]->GetModulus();
    const uint64_t inv_p = intel::hexl::InverseMod(p, q);
    if(q >= p){ // copy
      for (size_t j = 0; j < Q; j++){
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }else{ // reduce
      for (size_t j = 0; j < Q; j++){
        const uint64_t in_j = intel::hexl::ReduceMod<2>(out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], in_j, q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }
  } 
  out->l -= 1;
}


// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc(RNSc_Polynomial out){
  polynomial_base_reduce_RNSc_wo_free(out);
  free(out->coeffs[out->l]);
}


void polynomial_RNSc_permute(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t gen){
  const uint64_t Q = out->Q;
  assert(gen < Q);
  assert(gen > 0);
  uint64_t idx = 0;
  polynomial_RNS_zero((RNS_Polynomial) out);
  for (size_t i = 0; i < Q; i++){
    for (size_t j = 0; j < out->l; j++){
      out->coeffs[j][idx] = in->coeffs[j][i];
    }
    idx = intel::hexl::AddUIntMod(idx, gen, Q);
  }
}

void polynomial_int_permute_mod_Q(IntPolynomial out, IntPolynomial in, uint64_t gen){
  const uint64_t Q = in->N;
  uint64_t idx = 0;
  for (size_t i = 0; i < Q; i++){
    out->coeffs[idx] = in->coeffs[i];
    idx = intel::hexl::AddUIntMod(idx, gen, Q);
  }
}

void polynomial_RNSc_mul_by_xai(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t a){
  assert(in != out);
  const uint64_t Q = out->Q;
  assert(a < Q);
  polynomial_RNS_zero((RNS_Polynomial) out);
  for (size_t j = 0; j < out->l; j++){
    for (size_t i = 0; i < Q-a; i++){
      out->coeffs[j][i+a] = in->coeffs[j][i];
    }
    for (size_t i = Q-a; i < Q; i++){
      out->coeffs[j][i-Q+a] = in->coeffs[j][i];
    }
  }
}

void polynomial_int_decompose_i(IntPolynomial out, IntPolynomial in, uint64_t Bg_bit, uint64_t l, uint64_t q, uint64_t bit_size, uint64_t i){
  const uint64_t N = in->N;
  // const uint64_t half_Bg = (1UL << (Bg_bit - 1));
  const uint64_t h_mask = (1UL << Bg_bit) - 1;
  const uint64_t h_bit = bit_size - (i + 1) * Bg_bit;

  uint64_t offset = 1ULL << (bit_size - l * Bg_bit - 1);
  // for (size_t i = 1; i < l; i++){
  //   offset += (1UL << (bit_size - i * Bg_bit - 1));
  // }
  // offset = 0;
  
  for (size_t c = 0; c < N; c++){
    const uint64_t coeff_off = in->coeffs[c] + offset;
    out->coeffs[c] = (coeff_off>>h_bit) & h_mask;
  }
}

IntPolynomial polynomial_new_int_polynomial(uint64_t N){
  IntPolynomial res;
  res = (IntPolynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t) * N);
  res->N = N;
  return res;
}

void free_polynomial(void * p){
  free(((IntPolynomial) p)->coeffs);
  free(p);
}
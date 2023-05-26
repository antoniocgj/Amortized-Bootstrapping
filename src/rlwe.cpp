#include "rns-gsw.h"

// standard (cyclotomic power of 2) RLWE basic functions

RLWE_Key rlwe_alloc_key(uint64_t N, uint64_t k, intel::hexl::NTT * ntt){
  RLWE_Key key;
  key = (RLWE_Key) safe_malloc(sizeof(*key));
  key->s = (IntPolynomial *) safe_malloc(sizeof(*key->s)*k);
  key->s_DFT = (IntPolynomial *) safe_malloc(sizeof(*key->s_DFT)*k);
  for (size_t i = 0; i < k; i++){
    key->s[i] = polynomial_new_int_polynomial(N);
    key->s_DFT[i] = polynomial_new_int_polynomial(N);
  }
  key->N = N;
  key->k = k;
  key->ntt = ntt;
  return key;
}

RLWE rlwe_alloc_sample(uint64_t N, uint64_t k, intel::hexl::NTT * ntt){
  RLWE c;
  c = (RLWE) safe_malloc(sizeof(*c));
  c->a = (IntPolynomial *) safe_malloc(sizeof(*c->a)*k);
  for (size_t i = 0; i < k; i++){
    c->a[i] = polynomial_new_int_polynomial(N);
  }  
  c->b = polynomial_new_int_polynomial(N);
  c->k = k;
  c->ntt = ntt;
  return c;
}

RLWE_Key rlwe_new_key(uint64_t N, uint64_t q, uint64_t k, double sec_sigma, double err_sigma, intel::hexl::NTT * ntt){
  assert(q == ntt->GetModulus());
  RLWE_Key key = rlwe_alloc_key(N, k, ntt);
  for (size_t j = 0; j < k; j++){
    for (size_t i = 0; i < N; i++){
      key->s[j]->coeffs[i] = (q + double2int(generate_normal_random(sec_sigma)))%q;
    }
    ntt->ComputeForward(key->s_DFT[j]->coeffs, key->s[j]->coeffs, 1, 1);
  }
  key->N = N;
  key->q = q;
  key->sigma = err_sigma;
  return key;
}

// Generates a sparse ternary key with Hamming Weight h, balanced
RLWE_Key rlwe_new_sparse_ternary_key(uint64_t N, uint64_t q, uint64_t k, uint64_t h, double err_sigma, intel::hexl::NTT * ntt){
  assert(q == ntt->GetModulus());
  RLWE_Key key = rlwe_alloc_key(N, k, ntt);
  for (size_t j = 0; j < k; j++){
    gen_sparse_ternary_array(key->s[j]->coeffs, N, h, q);
    ntt->ComputeForward(key->s_DFT[j]->coeffs, key->s[j]->coeffs, 1, 1);
  }
  key->N = N;
  key->q = q;
  key->sigma = err_sigma;
  return key;
}

void rlwe_key_mod_switch(RLWE_Key key, intel::hexl::NTT * ntt){
  const uint64_t p = key->q;
  const uint64_t q = ntt->GetModulus();
  key->ntt = ntt;
  key->q = q;
  for (size_t i = 0; i < key->k; i++){
    array_additive_inverse_mod_switch(key->s[i]->coeffs, key->s[i]->coeffs, p, q, key->N);
    ntt->ComputeForward(key->s_DFT[i]->coeffs, key->s[i]->coeffs, 1, 1);
  }
}

void rlwe_sample(RLWE c, IntPolynomial m, RLWE_Key key){
  generate_random_bytes(key->N*sizeof(*c->a[0]->coeffs), (uint8_t *) c->a[0]->coeffs);
  array_mod_switch_from_2k(c->a[0]->coeffs, c->a[0]->coeffs, key->q, key->q, key->N);
  key->ntt->ComputeForward(c->a[0]->coeffs, c->a[0]->coeffs, 1, 1);
  intel::hexl::EltwiseMultMod(c->b->coeffs, c->a[0]->coeffs, key->s_DFT[0]->coeffs, key->N, key->q, 1);
  key->ntt->ComputeInverse(c->a[0]->coeffs, c->a[0]->coeffs, 1, 1);
  for (size_t i = 1; i < key->k; i++){
    generate_random_bytes(key->N*sizeof(*c->a[i]->coeffs), (uint8_t *) c->a[i]->coeffs);
    array_mod_switch_from_2k(c->a[i]->coeffs, c->a[i]->coeffs, key->q, key->q, key->N);
    key->ntt->ComputeForward(c->a[i]->coeffs, c->a[i]->coeffs, 1, 1);
    for (size_t j = 0; j < key->N; j++){
      const uint64_t as_i = intel::hexl::MultiplyMod(c->a[i]->coeffs[j], key->s_DFT[i]->coeffs[j], key->q);
      c->b->coeffs[j] = intel::hexl::AddUIntMod(c->b->coeffs[j], as_i, key->q);
    }
    key->ntt->ComputeInverse(c->a[i]->coeffs, c->a[i]->coeffs, 1, 1);
  }
  key->ntt->ComputeInverse(c->b->coeffs, c->b->coeffs, 1, 1);
  for (size_t i = 0; i < key->N; i++){
    c->b->coeffs[i] = intel::hexl::AddUIntMod(c->b->coeffs[i],  (key->q + double2int(generate_normal_random(key->sigma)))%key->q, key->q);
  }
  if(m){
    intel::hexl::EltwiseAddMod(c->b->coeffs, c->b->coeffs, m->coeffs, m->N, key->q);
  }
}

RLWE rlwe_new_sample(IntPolynomial m, RLWE_Key key){
  RLWE c = rlwe_alloc_sample(key->N, key->k, key->ntt);
  rlwe_sample(c, m, key);
  return c;
}

void rlwe_trivial_sample(RLWE c, IntPolynomial m){
  const uint64_t N = c->a[0]->N, k = c->k;
  for (size_t i = 0; i < k; i++){
    memset(c->a[i]->coeffs, 0, sizeof(*c->a[i]->coeffs)*N);
  }
  if(m){
    memcpy(c->b->coeffs, m->coeffs, sizeof(*m->coeffs)*m->N);
  }else{
    memset(c->b->coeffs, 0, sizeof(*c->b->coeffs)*N);
  }
}

RLWE rlwe_new_trivial_sample(IntPolynomial m, uint64_t N, uint64_t k, intel::hexl::NTT * ntt){
  RLWE c = rlwe_alloc_sample(N, k, ntt);
  rlwe_trivial_sample(c, m);
  return c;
}

void rlwe_phase(IntPolynomial out, RLWE c, RLWE_Key key){
  assert(c->k == key->k);
  key->ntt->ComputeForward(c->a[0]->coeffs, c->a[0]->coeffs, 1, 1);
  intel::hexl::EltwiseMultMod(out->coeffs, c->a[0]->coeffs, key->s_DFT[0]->coeffs, key->N, key->q, 1);
  key->ntt->ComputeInverse(c->a[0]->coeffs, c->a[0]->coeffs, 1, 1);
  for (size_t i = 1; i < key->k; i++){
    key->ntt->ComputeForward(c->a[i]->coeffs, c->a[i]->coeffs, 1, 1);
    for (size_t j = 0; j < key->N; j++){
      const uint64_t as_i = intel::hexl::MultiplyMod(c->a[i]->coeffs[j], key->s_DFT[i]->coeffs[j], key->q);
      out->coeffs[j] = intel::hexl::AddUIntMod(out->coeffs[j], as_i, key->q);
    }
    key->ntt->ComputeInverse(c->a[i]->coeffs, c->a[i]->coeffs, 1, 1);
  }
  key->ntt->ComputeInverse(out->coeffs, out->coeffs, 1, 1);
  intel::hexl::EltwiseSubMod(out->coeffs, c->b->coeffs, out->coeffs, out->N, key->q);
}

void rlwe_to_DFT(RLWE c){
  for (size_t i = 0; i < c->k; i++){
    c->ntt->ComputeForward(c->a[i]->coeffs, c->a[i]->coeffs, 1, 1);    
  }
  c->ntt->ComputeForward(c->b->coeffs, c->b->coeffs, 1, 1);    
}

void rlwe_from_DFT(RLWE c){
  for (size_t i = 0; i < c->k; i++){
    c->ntt->ComputeInverse(c->a[i]->coeffs, c->a[i]->coeffs, 1, 1);    
  }
  c->ntt->ComputeInverse(c->b->coeffs, c->b->coeffs, 1, 1);    
}

RLWE_KS_Key rlwe_new_full_packing_KS_key(RLWE_Key out_key, LWE_Key in_key, uint64_t t, uint64_t base_bit){
  assert((in_key->n%512) == 0);
  const uint64_t bit_size = ((uint64_t)log2(in_key->q)) + 1;
  RLWE_KS_Key res;
  res = (RLWE_KS_Key) safe_malloc(sizeof(*res));

  IntPolynomial dec_poly = polynomial_new_int_polynomial(out_key->N);
  memset(dec_poly->coeffs, 0, sizeof(*dec_poly->coeffs)*dec_poly->N);

  res->t = t;
  res->base_bit = base_bit;
  res->s = (RLWE **) safe_malloc(sizeof(RLWE *)*in_key->n);
  for (size_t i = 0; i < in_key->n; i++){
    res->s[i] = (RLWE *) safe_malloc(sizeof(RLWE)*t);
    for (size_t j = 0; j < t; j++){
      dec_poly->coeffs[0] = intel::hexl::MultiplyMod(in_key->s[i], (1UL << (bit_size - (j + 1) * base_bit)), in_key->q);
      res->s[i][j] = rlwe_new_sample(dec_poly, out_key);
      rlwe_to_DFT(res->s[i][j]);
    }
  }
  free_polynomial(dec_poly);
  return res;
}


void rlwe_full_packing_keyswitch(RLWE out, LWE * in, uint64_t size, RLWE_KS_Key ks_key){
  const uint64_t N = out->a[0]->N, q = out->ntt->GetModulus(), k = out->k;
  const uint64_t bit_size = ((uint64_t)log2(in[0]->q)) + 1;

  IntPolynomial dec_in_a = polynomial_new_int_polynomial(N);
  IntPolynomial a_i = polynomial_new_int_polynomial(N);
  IntPolynomial tmp = polynomial_new_int_polynomial(N);


  memset(a_i->coeffs, 0, sizeof(*a_i->coeffs)*N);
  rlwe_trivial_sample(out, 0);

  for (size_t i = 0; i < in[0]->n; i++){
    for (size_t j = 0; j < size; j++){
      a_i->coeffs[j] = in[j]->a[i];
    }
    for (size_t j = 0; j < ks_key->t; j++){
      polynomial_int_decompose_i(dec_in_a, a_i, ks_key->base_bit, ks_key->t, q, bit_size, j);
      out->ntt->ComputeForward(dec_in_a->coeffs, dec_in_a->coeffs, 1, 1);
      // a
      for (size_t m = 0; m < k; m++){
        intel::hexl::EltwiseMultMod(tmp->coeffs, dec_in_a->coeffs, ks_key->s[i][j]->a[m]->coeffs, N, q, 1);
        intel::hexl::EltwiseSubMod(out->a[m]->coeffs, out->a[m]->coeffs, tmp->coeffs, N, q);
      }  
      // b
      intel::hexl::EltwiseMultMod(tmp->coeffs, dec_in_a->coeffs, ks_key->s[i][j]->b->coeffs, N, q, 1);
      intel::hexl::EltwiseSubMod(out->b->coeffs, out->b->coeffs, tmp->coeffs, N, q);
    }
  }
  rlwe_from_DFT(out);
  for (size_t j = 0; j < size; j++){
    a_i->coeffs[j] = in[j]->b;
  }
  intel::hexl::EltwiseAddMod(out->b->coeffs, out->b->coeffs,  a_i->coeffs, N, q);

  free_polynomial(a_i);
  free_polynomial(tmp);
  free_polynomial(dec_in_a);
}

RLWE_KS_Key rlwe_new_KS_key(RLWE_Key out_key, RLWE_Key in_key, uint64_t t, uint64_t base_bit){
  const uint64_t bit_size = ((uint64_t)log2(in_key->q)) + 1, N = in_key->N;
  assert(N == out_key->N);
  RLWE_KS_Key res;
  res = (RLWE_KS_Key) safe_malloc(sizeof(*res));

  IntPolynomial dec_poly = polynomial_new_int_polynomial(in_key->N);
  memset(dec_poly->coeffs, 0, sizeof(*dec_poly->coeffs)*dec_poly->N);

  res->t = t;
  res->base_bit = base_bit;
  res->s = (RLWE **) safe_malloc(sizeof(RLWE *)*in_key->k);
  for (size_t i = 0; i < in_key->k; i++){
    res->s[i] = (RLWE *) safe_malloc(sizeof(RLWE)*t);
    for (size_t j = 0; j < t; j++){
      for (size_t i2 = 0; i2 < N; i2++){
        dec_poly->coeffs[i2] = in_key->s[i]->coeffs[i2] * (1UL << (bit_size - (j + 1) * base_bit));
      }
      res->s[i][j] = rlwe_new_sample(dec_poly, out_key);
      rlwe_to_DFT(res->s[i][j]);
    }
  }
  free_polynomial(dec_poly);
  return res;
}

void rlwe_keyswitch(RLWE out, RLWE in, RLWE_KS_Key ks_key){
  const uint64_t N = out->a[0]->N, q = out->ntt->GetModulus();
  const uint64_t bit_size = ((uint64_t)log2(in->ntt->GetModulus())) + 1;

  IntPolynomial dec_in_a = polynomial_new_int_polynomial(N);
  IntPolynomial a_i = polynomial_new_int_polynomial(N);
  IntPolynomial tmp = polynomial_new_int_polynomial(N);


  memset(a_i->coeffs, 0, sizeof(*a_i->coeffs)*N);
  rlwe_trivial_sample(out, 0);

  for (size_t i = 0; i < in->k; i++){
    for (size_t j = 0; j < ks_key->t; j++){
      polynomial_int_decompose_i(dec_in_a, in->a[i], ks_key->base_bit, ks_key->t, q, bit_size, j);
      out->ntt->ComputeForward(dec_in_a->coeffs, dec_in_a->coeffs, 1, 1);
      // a
      for (size_t m = 0; m < out->k; m++){
        intel::hexl::EltwiseMultMod(tmp->coeffs, dec_in_a->coeffs, ks_key->s[i][j]->a[m]->coeffs, N, q, 1);
        intel::hexl::EltwiseSubMod(out->a[m]->coeffs, out->a[m]->coeffs, tmp->coeffs, N, q);
      }  
      // b
      intel::hexl::EltwiseMultMod(tmp->coeffs, dec_in_a->coeffs, ks_key->s[i][j]->b->coeffs, N, q, 1);
      intel::hexl::EltwiseSubMod(out->b->coeffs, out->b->coeffs, tmp->coeffs, N, q);
    }
  }
  rlwe_from_DFT(out);
  intel::hexl::EltwiseAddMod(out->b->coeffs, out->b->coeffs, in->b->coeffs, N, q);

  free_polynomial(a_i);
  free_polynomial(tmp);
  free_polynomial(dec_in_a);
}

void rlwe_mod_switch(RLWE c, intel::hexl::NTT * ntt){
  const uint64_t p = c->ntt->GetModulus(), N = c->b->N;
  const uint64_t q = ntt->GetModulus();
  for (size_t i = 0; i < c->k; i++){
    array_mod_switch(c->a[i]->coeffs, c->a[i]->coeffs, p, q, N);
  }
  array_mod_switch(c->b->coeffs, c->b->coeffs, p, q, N);
  c->ntt = ntt;
}

// RLWE RNS functions

RNS_RLWE_Key rlwe_alloc_RNS_key(uint64_t Q, uint64_t l, uint64_t log_mod, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res;
  res = (RNS_RLWE_Key) safe_malloc(sizeof(*res));
  res->log_mod = log_mod;
  res->sigma = sigma;
  res->N = next_power_of_2(Q) << 1;
  res->Q = Q;
  res->l = l;
  res->s = polynomial_new_int_polynomial(Q);
  res->s_RNS = polynomial_new_RNS_polynomial(Q, l, ntt);
  return res;
}

void free_rlwe_RNS_key(RNS_RLWE_Key key){
  free_polynomial(key->s);
  free_RNS_polynomial(key->s_RNS);
  free(key);
}

RNS_RLWE_Key rlwe_new_RNS_key(uint64_t Q, uint64_t l, uint64_t log_mod, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res = rlwe_alloc_RNS_key(Q, l, log_mod, ntt, sigma);
  const uint64_t mod_mask = (1ULL << log_mod) - 1;
  generate_random_bytes(Q*sizeof(uint64_t), (uint8_t *) res->s->coeffs);
  for (size_t i = 0; i < Q; i++) res->s->coeffs[i] &= mod_mask; 
  polynomial_to_RNS(res->s_RNS, res->s);
  return res;
}

void rlwe_RNSc_extract_lwe_key(LWE_Key out, RNS_RLWE_Key in){
  const uint64_t Q = in->Q;
  assert(Q == out->n);
  out->q = in->s_RNS->ntt[0]->GetModulus();
  for (size_t i = 0; i < Q; i++){
    out->s[i] = in->s->coeffs[i];
  }
}

void rlwe_RNSc_extract_lwe(LWE out, RNSc_RLWE in, uint64_t idx){
  const uint64_t Q = in->a->Q;
  assert(in->a->l == 1);
  for (size_t i = 0; i <= idx; i++){
    out->a[i] = in->a->coeffs[0][idx - i];
  }
  for (size_t i = idx + 1; i < Q; i++){
    out->a[i] = in->a->coeffs[0][Q + idx - i];
  }
  out->b = in->b->coeffs[0][idx];
}

// negacyclic extract
void rlwe_RNSc_extract_lwe_nc(LWE out, RNSc_RLWE in, uint64_t idx){
  const uint64_t Q = in->a->Q;
  assert(in->a->l == 1);
  for (size_t i = 0; i <= idx; i++){
    out->a[i] = in->a->coeffs[0][idx - i];
  }
  for (size_t i = idx + 1; i < Q; i++){
    out->a[i] = -in->a->coeffs[0][Q + idx - i];
  }
  out->b = in->b->coeffs[0][idx];
}

RNS_RLWE_Key rlwe_new_RNS_gaussian_key(uint64_t Q, uint64_t l, double key_sigma, intel::hexl::NTT ** ntt, double sigma){
  RNS_RLWE_Key res = rlwe_alloc_RNS_key(Q, l, ceil(key_sigma), ntt, sigma);
  clwe_sample_gaussian(res->s->coeffs, Q, key_sigma);
  polynomial_to_RNS(res->s_RNS, res->s);
  return res;
}

RNS_RLWE rlwe_alloc_RNS_sample(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_RLWE res;
  res = (RNS_RLWE) safe_malloc(sizeof(*res));
  res->a = polynomial_new_RNS_polynomial(Q, l, ntt);
  res->b = polynomial_new_RNS_polynomial(Q, l, ntt);
  return res;
}

void free_RNS_rlwe_sample(RNS_RLWE c){
  free_RNS_polynomial(c->a);
  free_RNS_polynomial(c->b);
  free(c);
}

void rlwe_copy_RNS_sample(RNS_RLWE out, RNS_RLWE in){
  polynomial_copy_RNS_polynomial(out->a, in->a);
  polynomial_copy_RNS_polynomial(out->b, in->b);
}

void rlwe_copy_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in){
  rlwe_copy_RNS_sample((RNS_RLWE) out, (RNS_RLWE) in);
}

void free_rlwe_RNS_sample(void * p){
  const RNS_RLWE pp = (RNS_RLWE) p;
  free_RNS_polynomial(pp->a);
  free_RNS_polynomial(pp->b);
  free(pp);
}

void rlwe_RNSc_sample_of_zero(RNSc_RLWE out, RNS_RLWE_Key key){
  const uint64_t Q = key->Q, l = key->l;
  RNS_Polynomial tmp = polynomial_new_RNS_polynomial(Q, l, key->s_RNS->ntt);
  polynomial_gen_random_RNSc_polynomial(out->a);
  polynomial_RNSc_to_RNS(tmp, out->a);
  polynomial_mul_RNS_polynomial(out->b, key->s_RNS, tmp);
  polynomial_RNSc_add_noise(out->b, out->b, key->sigma);
  free_RNS_polynomial(tmp);
}


void rlwe_scale_RNSc_rlwe(RNSc_RLWE c, uint64_t scale){
  for (size_t j = 0; j < c->a->l; j++){
    const uint64_t mod = c->a->ntt[j]->GetModulus();
    for (size_t k = 0; k < c->a->Q; k++){
      c->a->coeffs[j][k] = intel::hexl::MultiplyMod(c->a->coeffs[j][k], scale, mod);
      c->b->coeffs[j][k] = intel::hexl::MultiplyMod(c->b->coeffs[j][k], scale, mod);
    } 
  }
}

void rlwe_RNSc_mod_switch(RNSc_RLWE c, uint64_t q){
  const uint64_t p = c->a->ntt[0]->GetModulus();
  array_mod_switch(c->a->coeffs[0], c->a->coeffs[0], p, q, c->a->Q);
  array_mod_switch(c->b->coeffs[0], c->b->coeffs[0], p, q, c->a->Q);
}

// Return Q_0 * X^m
RNSc_RLWE rlwe_new_RNSc_sample(RNS_RLWE_Key key, uint64_t m){
  RNSc_RLWE res = rlwe_new_RNSc_sample_of_zero(key);
  const uint64_t p = res->b->ntt[0]->GetModulus();
  res->b->coeffs[0][m] = intel::hexl::AddUIntMod(res->b->coeffs[0][m], 1, p);
  return res;
}

void rlwe_RNS_phase(RNS_Polynomial out, RNS_RLWE in, RNS_RLWE_Key key){
  polynomial_mul_RNS_polynomial((RNSc_Polynomial)out, in->a, key->s_RNS);
  polynomial_RNSc_to_RNS(out, (RNSc_Polynomial)out);
  polynomial_sub_RNS_polynomial(out, in->b, out);
}

RNSc_RLWE rlwe_new_RNSc_sample_of_zero(RNS_RLWE_Key key){
  RNSc_RLWE res = (RNSc_RLWE) rlwe_alloc_RNS_sample(key->Q, key->l, key->s_RNS->ntt);
  rlwe_RNSc_sample_of_zero(res, key);
  return res;
}

void rlwe_RNS_sample_of_zero(RNS_RLWE out, RNS_RLWE_Key key){
  rlwe_RNSc_sample_of_zero((RNSc_RLWE) out, key);
  rlwe_RNSc_to_RNS(out, (RNSc_RLWE) out);
}

RNS_RLWE rlwe_new_RNS_sample_of_zero(RNS_RLWE_Key key){
  RNS_RLWE res = rlwe_alloc_RNS_sample(key->Q, key->l, key->s_RNS->ntt);
  rlwe_RNS_sample_of_zero(res, key);
  return res;
}

RNS_RLWE rlwe_new_RNS_trivial_sample_of_zero(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_RLWE res = rlwe_alloc_RNS_sample(Q, l, ntt);
  for (size_t i = 0; i < res->a->l; i++){
    memset(res->a->coeffs[i], 0, sizeof(uint64_t)*res->a->N);
    memset(res->b->coeffs[i], 0, sizeof(uint64_t)*res->b->N);
  }
  return res;
}

void rlwe_RNS_trivial_sample_of_zero(RNS_RLWE out){
  for (size_t i = 0; i < out->a->l; i++){
    memset(out->a->coeffs[i], 0, sizeof(uint64_t)*out->a->N);
    memset(out->b->coeffs[i], 0, sizeof(uint64_t)*out->b->N);
  }
}

void rlwe_shrink_RNSc_sample(RNSc_RLWE c){
  polynomial_base_reduce_RNSc(c->a);
  polynomial_base_reduce_RNSc(c->b);
}

void rlwe_automorphism_RNSc(RNSc_RLWE out, RNSc_RLWE in, uint64_t gen, RNS_RLWE_KS_Key ksk){
  RNSc_RLWE tmp = (RNSc_RLWE) rlwe_alloc_RNS_sample(out->a->Q, out->a->l, out->a->ntt);
  polynomial_RNSc_permute(tmp->a, in->a, gen);  
  polynomial_RNSc_permute(tmp->b, in->b, gen);
  rlwe_RNSc_keyswitch(out, tmp, ksk);
  free_rlwe_RNS_sample(tmp);
}

void rlwe_addto_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in){
  polynomial_add_RNSc_polynomial(out->a, out->a, in->a);
  polynomial_add_RNSc_polynomial(out->b, out->b, in->b);
}

void rlwe_decomp_automorphism_RNSc(RNSc_RLWE out, RNSc_RLWE in, uint64_t gen, RNS_RLWE_KS_Key * ksk){
  const uint64_t Q = out->a->Q;
  RNSc_RLWE tmp = (RNSc_RLWE) rlwe_new_RNS_trivial_sample_of_zero(Q, out->a->l, out->a->ntt);
  assert(gen < Q);
  gen = discret_log_base_11_mod_12289[gen - 1];
  if(!(gen&1)) rlwe_copy_RNSc_sample(out, in);
  else{
    const uint64_t dec_gen = 11;
    polynomial_RNSc_permute(tmp->a, in->a, dec_gen);  
    polynomial_RNSc_permute(tmp->b, in->b, dec_gen);
    rlwe_RNSc_keyswitch(out, tmp, ksk[0]);
  }
  gen>>=1;
  for (size_t i = 1; gen; i++, gen>>=1){
    const uint64_t dec_gen = intel::hexl::PowMod(11, ((gen&1) << i), Q);
    if(dec_gen == 1) continue;
    polynomial_RNSc_permute(tmp->a, out->a, dec_gen);  
    polynomial_RNSc_permute(tmp->b, out->b, dec_gen);
    rlwe_RNSc_keyswitch(out, tmp, ksk[i]);
  }
  free_rlwe_RNS_sample(tmp);
}


void rlwe_RNSc_to_RNS(RNS_RLWE out, RNSc_RLWE in){
  polynomial_RNSc_to_RNS(out->a, in->a);
  polynomial_RNSc_to_RNS(out->b, in->b);
}

void rlwe_RNS_to_RNSc(RNSc_RLWE out, RNS_RLWE in){
  polynomial_RNS_to_RNSc(out->a, in->a);
  polynomial_RNS_to_RNSc(out->b, in->b);
}

RNS_RLWE_KS_Key rlwe_new_RNS_ks_key(RNS_RLWE_Key out_key, RNS_RLWE_Key in_key){
  const uint64_t l = in_key->l;
  RNS_RLWE_KS_Key res;
  res = (RNS_RLWE_KS_Key) safe_malloc(sizeof(*res));
  res->s = (RNS_RLWE*) safe_malloc(sizeof(RNS_RLWE)*l);
  res->l = l;
  for (size_t i = 0; i < l; i++){
    res->s[i] = (RNS_RLWE) rlwe_new_RNSc_sample_of_zero(out_key);
    const uint64_t p = res->s[i]->b->ntt[i]->GetModulus();
    for (size_t j = 0; j < in_key->Q; j++){
      res->s[i]->b->coeffs[i][j] = intel::hexl::AddUIntMod(res->s[i]->b->coeffs[i][j], (p+in_key->s->coeffs[j])%p, p);
    }
    rlwe_RNSc_to_RNS(res->s[i], (RNSc_RLWE) res->s[i]);
  }
  return res;
}

void naive_full_mul(uint64_t * out, uint64_t * in1, uint64_t * in2, size_t N){
  memset(out, 0, sizeof(uint64_t)*N*2);
  for (size_t i = 0; i < N; i++){
    for (size_t j = 0; j < N; j++){
      out[i + j] += in1[i]*in2[j];
    }
  }
}

void print_vec(const char * msg, uint64_t * v, size_t size);

RNS_RLWE_KS_Key * rlwe_new_RNS_priv_ks_key(RNS_RLWE_Key out_key, RNS_RLWE_Key in_key){
  RNS_RLWE_KS_Key * res = (RNS_RLWE_KS_Key *) safe_malloc(sizeof(RNS_RLWE_KS_Key)*2);
  // B 
  // print_vec("B key", tmp_key->s->coeffs, tmp_key->Q);
  const uint64_t l = in_key->l;
  const uint64_t sigma = out_key->sigma;
  out_key->sigma = 0;
  res[1] = (RNS_RLWE_KS_Key) safe_malloc(sizeof(*res[1]));
  res[1]->s = (RNS_RLWE*) safe_malloc(sizeof(RNS_RLWE)*l);
  res[1]->l = l;
  for (size_t i = 0; i < l; i++){
    res[1]->s[i] = (RNS_RLWE) rlwe_new_RNSc_sample_of_zero(out_key);
    const uint64_t p = res[1]->s[i]->a->ntt[i]->GetModulus();
    res[1]->s[i]->a->coeffs[i][0] = intel::hexl::AddUIntMod(res[1]->s[i]->a->coeffs[i][0], 1, p);
    rlwe_RNSc_to_RNS(res[1]->s[i], (RNSc_RLWE) res[1]->s[i]);
  }
  // A
  res[0] = (RNS_RLWE_KS_Key) safe_malloc(sizeof(*res[0]));
  res[0]->s = (RNS_RLWE*) safe_malloc(sizeof(RNS_RLWE)*l);
  res[0]->l = l;
  for (size_t i = 0; i < l; i++){
    res[0]->s[i] = (RNS_RLWE) rlwe_new_RNSc_sample_of_zero(out_key);
    polynomial_mul_RNS_polynomial((RNSc_Polynomial) res[0]->s[i]->a, res[1]->s[i]->a, in_key->s_RNS);
    polynomial_mul_RNS_polynomial((RNSc_Polynomial) res[0]->s[i]->b, res[1]->s[i]->b, in_key->s_RNS);
  }
  for (size_t i = 0; i < l; i++){
    rlwe_RNS_to_RNSc((RNSc_RLWE) res[1]->s[i], res[1]->s[i]);
    polynomial_RNSc_add_noise((RNSc_Polynomial) res[1]->s[i]->b, (RNSc_Polynomial) res[1]->s[i]->b, sigma);
    polynomial_RNSc_add_noise((RNSc_Polynomial) res[0]->s[i]->b, (RNSc_Polynomial) res[0]->s[i]->b, sigma);
    rlwe_RNSc_to_RNS(res[1]->s[i], (RNSc_RLWE) res[1]->s[i]);
    rlwe_RNSc_to_RNS(res[0]->s[i], (RNSc_RLWE) res[0]->s[i]);
  }
  out_key->sigma = sigma;
  return res;
}

void free_rlwe_RNS_ks_key(RNS_RLWE_KS_Key key){
  for (size_t i = 0; i < key->l; i++){
    free_rlwe_RNS_sample(key->s[i]);
  }
  free(key->s);
  free(key);
}

RNS_RLWE_KS_Key rlwe_new_RNS_automorphism_key(RNS_RLWE_Key key, uint64_t gen){
  RNS_RLWE_Key perm_key = rlwe_alloc_RNS_key(key->Q, key->l, key->log_mod, key->s_RNS->ntt, key->sigma);
  polynomial_int_permute_mod_Q(perm_key->s, key->s, gen);
  // It's not necessary to permute s_RNS
  RNS_RLWE_KS_Key res = rlwe_new_RNS_ks_key(key, perm_key);
  free_rlwe_RNS_key(perm_key);
  return res;
}

RNS_RLWE_KS_Key * rlwe_new_RNS_automorphism_keyset(RNS_RLWE_Key key){
  RNS_RLWE_Key perm_key = rlwe_alloc_RNS_key(key->Q, key->l, key->log_mod, key->s_RNS->ntt, key->sigma);
  const uint64_t Q = key->Q;
  RNS_RLWE_KS_Key * res = (RNS_RLWE_KS_Key *) safe_malloc(sizeof(RNS_RLWE_KS_Key)*Q);
  for (size_t i = 1; i < Q; i++){
    polynomial_int_permute_mod_Q(perm_key->s, key->s, i);
    res[i] = rlwe_new_RNS_ks_key(key, perm_key);
  }
  free_rlwe_RNS_key(perm_key);
  return res;
}

RNS_RLWE_KS_Key * rlwe_new_RNS_decomp_automorphism_keyset(RNS_RLWE_Key key, uint64_t unit_generator){
  RNS_RLWE_Key perm_key = rlwe_alloc_RNS_key(key->Q, key->l, key->log_mod, key->s_RNS->ntt, key->sigma);
  const uint64_t log_Q = (uint64_t) log2(key->Q);
  RNS_RLWE_KS_Key * res = (RNS_RLWE_KS_Key *) safe_malloc(sizeof(RNS_RLWE_KS_Key)*(log_Q +  1));
  for (size_t i = 0; i <= log_Q; i++){
    const uint64_t gen = intel::hexl::PowMod(unit_generator, 1ULL << i, key->Q);
    polynomial_int_permute_mod_Q(perm_key->s, key->s, gen);
    res[i] = rlwe_new_RNS_ks_key(key, perm_key);
  }
  free_rlwe_RNS_key(perm_key);
  return res;
}

RNS_RLWE_KS_Key * rlwe_new_RNS_automorphism_keyset_rou(RNS_RLWE_Key key, uint64_t root_of_unity, uint64_t N){
  RNS_RLWE_Key perm_key = rlwe_alloc_RNS_key(key->Q, key->l, key->log_mod, key->s_RNS->ntt, key->sigma);
  const uint64_t Q = key->Q;
  RNS_RLWE_KS_Key * res = (RNS_RLWE_KS_Key *) safe_malloc(sizeof(RNS_RLWE_KS_Key)*2*N);
  for (size_t i = 0; i < 2*N; i++){
    const uint64_t gen = intel::hexl::PowMod(root_of_unity, i, Q);
    polynomial_int_permute_mod_Q(perm_key->s, key->s, gen);
    res[i] = rlwe_new_RNS_ks_key(key, perm_key);
  }
  free_rlwe_RNS_key(perm_key);
  return res;
}

void rlwe_RNSc_mul_by_xai(RNSc_RLWE out, RNSc_RLWE in, uint64_t a){
  polynomial_RNSc_mul_by_xai(out->a, in->a, a);
  polynomial_RNSc_mul_by_xai(out->b, in->b, a);
}

void rlwe_RNSc_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key ksk){
  assert(in != out);
  rlwe_RNS_trivial_sample_of_zero((RNS_RLWE)out);
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->Q, in->a->l, in->a->ntt);
  RNSc_Polynomial tmp2 = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->Q, in->a->l, in->a->ntt);
  const uint64_t l = in->a->l;
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in->a->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], ksk->s[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], ksk->s[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  polynomial_RNS_to_RNSc(out->a, (RNS_Polynomial) out->a);
  polynomial_RNS_to_RNSc(out->b, (RNS_Polynomial) out->b);
  polynomial_RNSc_reduce_mod_XQ_minus_1(out->a);
  polynomial_RNSc_reduce_mod_XQ_minus_1(out->b);
  polynomial_add_RNSc_polynomial(out->b, out->b, in->b);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
}

void rlwe_RNSc_priv_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key kska, RNS_RLWE_KS_Key kskb){
  rlwe_RNS_trivial_sample_of_zero((RNS_RLWE)out);
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->Q, in->a->l, in->a->ntt);
  RNSc_Polynomial tmp2 = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in->a->Q, in->a->l, in->a->ntt);
  const uint64_t l = in->a->l;
  // a
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in->a->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kska->s[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kska->s[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseSubMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  // b
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in->b->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kskb->s[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], kskb->s[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  // reduce
  polynomial_RNS_to_RNSc(out->a, (RNS_Polynomial) out->a);
  polynomial_RNS_to_RNSc(out->b, (RNS_Polynomial) out->b);
  polynomial_RNSc_reduce_mod_XQ_minus_1(out->a);
  polynomial_RNSc_reduce_mod_XQ_minus_1(out->b);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
}
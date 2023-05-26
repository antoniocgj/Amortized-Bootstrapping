// RLWE
#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <hexl/hexl.hpp>
#include <cstdint>
#include <iostream>
#include <vector>
#ifndef PORTABLE_BUILD
#include <x86intrin.h>
#endif

#include "../src/pre_comp_data.h"

// Note: RNS_* and RNSc_* must have the same memory layout

/* Polynomials */
typedef struct _IntPolynomial {
  uint64_t * coeffs;
  uint64_t N;
} * IntPolynomial;

typedef struct _RNS_Polynomial {
  uint64_t ** coeffs;
  intel::hexl::NTT ** ntt = NULL; 
  uint64_t N, l, Q;
} * RNS_Polynomial;

/* RNS polynomial in coefficient representation*/
typedef struct _RNSc_Polynomial {
  uint64_t ** coeffs;
  intel::hexl::NTT ** ntt = NULL; 
  uint64_t N, l, Q;
} * RNSc_Polynomial;


/* LWE */

typedef struct _LWE_Key {
  uint64_t * s;
  uint64_t n, q;
  double sigma;
} * LWE_Key;

typedef struct _LWE {
  uint64_t * a, b;
  uint64_t n, q;
} * LWE;

typedef struct _LWE_KS_Key {
  LWE ***s;
  uint64_t t, base_bit;
} * LWE_KS_Key;


/* RLWE */

typedef struct _RLWE_Key {
  IntPolynomial * s, * s_DFT;
  uint64_t N, q, k;
  intel::hexl::NTT * ntt = NULL;
  double sigma;
} * RLWE_Key;

typedef struct _RLWE {
  IntPolynomial * a, b;
  uint64_t k;
  intel::hexl::NTT * ntt = NULL;
} * RLWE;

typedef struct _RLWE_KS_Key {
  RLWE **s;
  uint64_t t, base_bit;
} * RLWE_KS_Key;

/* RLWE RNS */

typedef struct _RNS_RLWE_Key {
  uint64_t log_mod;
  IntPolynomial s;
  RNS_Polynomial s_RNS;
  uint64_t N, Q, l;
  double sigma;
} * RNS_RLWE_Key;

typedef struct _RNS_RLWE {
  RNS_Polynomial a, b;
} * RNS_RLWE;

typedef struct _RNSc_RLWE {
  RNSc_Polynomial a, b;
} * RNSc_RLWE;

typedef struct _RNS_RLWE_KS_Key {
  RNS_RLWE * s;
  uint64_t l;
} * RNS_RLWE_KS_Key;


/* GSW */
typedef struct _GSW_Key {
  RNS_RLWE_Key rlwe_key;
  uint64_t l;
} * RNS_GSW_Key;

typedef struct _RNS_GSW {
  RNS_RLWE * samples;
  uint64_t l;
} * RNS_GSW;

typedef struct _RNSc_GSW {
  RNSc_RLWE * samples;
  uint64_t l;
} * RNSc_GSW;

// half GSW is a GSW ciphertext with only the \\ell last lines from a complete GSW (which has (k+1)\\ell lines), encrypting the decomposed message. 
typedef struct _RNS_half_GSW {
  RNS_RLWE * samples;
  uint64_t l;
} * RNS_half_GSW;

typedef struct _RNSc_half_GSW {
  RNSc_RLWE * samples;
  uint64_t l;
} * RNSc_half_GSW;

// hrd

typedef struct _RLWE_Bootstrap_KeySet {
  RNS_GSW * bk, * s_tmp;
  RNS_RLWE_KS_Key * priv_ksk, * aut_ksk;
  LWE_KS_Key lwe_ksk;
  RLWE_KS_Key packing_ksk;
  RNSc_RLWE * tmp_out;
  LWE * tmp_lwe_out, tmp_extracted_in;
  RLWE tmp_packed_in;
  uint64_t Q, Q2, l, rou_2nth;
  bool use_two_steps_dimred;
} * RLWE_Bootstrap_KeySet;

// clwe sampling
void clwe_sample_uniform(uint64_t * out, uint64_t size, uint64_t p);
void clwe_sample_gaussian(uint64_t * out, uint64_t size, double sigma);
void clwe_sample_gaussian_addto(uint64_t * out, uint64_t size, double sigma);


// polynomial
IntPolynomial polynomial_new_int_polynomial(uint64_t N);
RNS_Polynomial polynomial_new_RNS_polynomial(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt);
void polynomial_RNS_zero(RNS_Polynomial p);
RNS_Polynomial * polynomial_new_array_of_RNS_polynomials(uint64_t Q, uint64_t l, uint64_t size, intel::hexl::NTT ** ntt);
void polynomial_to_RNS(RNS_Polynomial out, IntPolynomial in);
void polynomial_gen_random_RNSc_polynomial(RNSc_Polynomial out);
void polynomial_mul_RNS_polynomial(RNSc_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2);
void polynomial_sub_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2);
void polynomial_sub_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2);
void polynomial_RNSc_to_RNS(RNS_Polynomial out, RNSc_Polynomial in);
void polynomial_RNS_to_RNSc(RNSc_Polynomial out, RNS_Polynomial in);
void polynomial_RNSc_add_noise(RNSc_Polynomial out, RNSc_Polynomial in, double sigma);
void polynomial_base_reduce_RNSc(RNSc_Polynomial out);
void polynomial_base_reduce_RNSc_wo_free(RNSc_Polynomial out);
void polynomial_base_extend_RNSc(RNSc_Polynomial out, uint64_t * in, uint64_t p);
void polynomial_RNSc_reduce_mod_XQ_minus_1(RNSc_Polynomial out);
void polynomial_RNSc_permute(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t gen);
void free_RNS_polynomial(void * p);
void polynomial_RNSc_negate(RNSc_Polynomial out, RNSc_Polynomial in);
void polynomial_add_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2);
void polynomial_RNSc_permute(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t gen);
void polynomial_int_permute_mod_Q(IntPolynomial out, IntPolynomial in, uint64_t gen);
void polynomial_copy_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in);
void polynomial_reduce_mod_Xd_minus_1_mod_q(uint64_t * p, int64_t d, uint64_t q, size_t size);
void polynomial_reduce_mod_Xd_minus_1(uint64_t * p, int64_t d, size_t size);
void polynomial_RNSc_mul_by_xai(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t a);
void polynomial_int_decompose_i(IntPolynomial out, IntPolynomial in, uint64_t Bg_bit, uint64_t l, uint64_t q, uint64_t bit_size, uint64_t i);
void free_polynomial(void * p);

// lwe
LWE_Key lwe_alloc_key(uint64_t n);
void free_lwe_sample(LWE c);
LWE lwe_alloc_sample(uint64_t n, uint64_t q);
LWE_Key lwe_new_key(uint64_t n, uint64_t q, double sec_sigma, double err_sigma);
LWE_Key lwe_new_sparse_ternary_key(uint64_t n, uint64_t q, uint64_t h, double err_sigma);
void lwe_sample(LWE c, uint64_t m, LWE_Key key);
LWE lwe_new_sample(uint64_t m, LWE_Key key);
LWE lwe_new_trivial_sample(uint64_t m, uint64_t n, uint64_t q);
uint64_t lwe_phase(LWE c, LWE_Key key);
LWE_KS_Key lwe_new_KS_key(LWE_Key out_key, LWE_Key in_key, uint64_t t, uint64_t base_bit);
void lwe_keyswitch(LWE out, LWE in, LWE_KS_Key ks_key);


// rlwe
RLWE_Key rlwe_alloc_key(uint64_t N, uint64_t k,  intel::hexl::NTT * ntt);
RLWE rlwe_alloc_sample(uint64_t N, uint64_t k, intel::hexl::NTT * ntt);
RLWE_Key rlwe_new_key(uint64_t N, uint64_t q, uint64_t k, double sec_sigma, double err_sigma, intel::hexl::NTT * ntt);
RLWE_Key rlwe_new_sparse_ternary_key(uint64_t N, uint64_t q, uint64_t k, uint64_t h, double err_sigma, intel::hexl::NTT * ntt);
void rlwe_sample(RLWE c, IntPolynomial m, RLWE_Key key);
RLWE rlwe_new_sample(IntPolynomial m, RLWE_Key key);
void rlwe_trivial_sample(RLWE c, IntPolynomial m);
RLWE rlwe_new_trivial_sample(IntPolynomial m, uint64_t N, uint64_t k, intel::hexl::NTT * ntt);
void rlwe_phase(IntPolynomial out, RLWE c, RLWE_Key key);
void rlwe_to_DFT(RLWE c);
void rlwe_from_DFT(RLWE c);
RLWE_KS_Key rlwe_new_full_packing_KS_key(RLWE_Key out_key, LWE_Key in_key, uint64_t t, uint64_t base_bit);
void rlwe_full_packing_keyswitch(RLWE out, LWE * in, uint64_t size, RLWE_KS_Key ks_key);
RLWE_KS_Key rlwe_new_KS_key(RLWE_Key out_key, RLWE_Key in_key, uint64_t t, uint64_t base_bit);
void rlwe_keyswitch(RLWE out, RLWE in, RLWE_KS_Key ks_key);
void rlwe_mod_switch(RLWE c, intel::hexl::NTT * ntt);
void rlwe_key_mod_switch(RLWE_Key key, intel::hexl::NTT * ntt);


// rlwe rns
RNS_RLWE_Key rlwe_alloc_RNS_key(uint64_t Q, uint64_t l, uint64_t log_mod, intel::hexl::NTT ** ntt, double sigma);
void free_RNS_rlwe_sample(RNS_RLWE c);
RNS_RLWE_Key rlwe_new_RNS_key(uint64_t Q, uint64_t l, uint64_t log_mod, intel::hexl::NTT ** ntt, double sigma);
RNS_RLWE_Key rlwe_new_RNS_gaussian_key(uint64_t Q, uint64_t l, double key_sigma, intel::hexl::NTT ** ntt, double sigma);
RNS_RLWE rlwe_alloc_RNS_sample(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt);
void rlwe_RNS_sample_of_zero(RNS_RLWE out, RNS_RLWE_Key key);
void rlwe_RNSc_sample_of_zero(RNSc_RLWE out, RNS_RLWE_Key key);
RNS_RLWE rlwe_new_RNS_sample_of_zero(RNS_RLWE_Key key);
RNSc_RLWE rlwe_new_RNSc_sample_of_zero(RNS_RLWE_Key key);
RNS_RLWE rlwe_new_RNS_trivial_sample_of_zero(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt);
void rlwe_RNS_phase(RNS_Polynomial out, RNS_RLWE in, RNS_RLWE_Key key);
void rlwe_RNSc_to_RNS(RNS_RLWE out, RNSc_RLWE in);
void rlwe_RNS_to_RNSc(RNSc_RLWE out, RNS_RLWE in);
RNSc_RLWE rlwe_new_RNSc_sample(RNS_RLWE_Key key, uint64_t m);
void rlwe_RNS_trivial_sample_of_zero(RNS_RLWE out);
void rlwe_copy_RNS_sample(RNS_RLWE out, RNS_RLWE in);
void rlwe_copy_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in);
void rlwe_shrink_RNSc_sample(RNSc_RLWE c);
void rlwe_RNSc_mul_by_xai(RNSc_RLWE out, RNSc_RLWE in, uint64_t a);
void rlwe_RNSc_extract_lwe_key(LWE_Key out, RNS_RLWE_Key in);
void rlwe_RNSc_extract_lwe(LWE out, RNSc_RLWE in, uint64_t idx);
void rlwe_RNSc_extract_lwe_nc(LWE out, RNSc_RLWE in, uint64_t idx);
void rlwe_scale_RNSc_rlwe(RNSc_RLWE c, uint64_t scale);
void rlwe_RNSc_mod_switch(RNSc_RLWE c, uint64_t q);
void rlwe_addto_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in);
void rlwe_decomp_automorphism_RNSc(RNSc_RLWE out, RNSc_RLWE in, uint64_t gen, RNS_RLWE_KS_Key * ksk);
RNS_RLWE_KS_Key * rlwe_new_RNS_decomp_automorphism_keyset(RNS_RLWE_Key key, uint64_t unit_generator);



// gsw
RNS_GSW_Key gsw_new_RNS_key(RNS_RLWE_Key rlwe_key);
RNS_GSW gsw_alloc_RNS_sample(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt);
RNSc_GSW gsw_new_RNSc_sample(RNS_GSW_Key key, uint64_t m);
void gsw_mul_RNSc_rlwe(RNSc_RLWE out, RNS_GSW in1, RNSc_RLWE in2);
void gsw_RNSc_to_RNS(RNS_GSW out, RNSc_GSW in);
void gsw_RNS_to_RNSc(RNSc_GSW out, RNS_GSW in);
void gsw_mul_RNSc_gsw(RNSc_GSW out, RNS_GSW in1, RNSc_GSW in2);
RNSc_GSW gsw_new_RNSc_ks_key(RNS_GSW_Key key);
void gsw_shrink_RNSc_gsw_sample(RNSc_GSW c);
void gsw_RNS_copy(RNS_GSW out, RNS_GSW in);
void gsw_scale_RNSc_gsw(RNSc_GSW c, uint64_t scale);


// half gsw
RNS_half_GSW gsw_alloc_half_RNS_sample(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt);
RNS_half_GSW * gsw_alloc_half_RNS_sample_array(uint64_t size, uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt);
void gsw_mul_RNSc_half_gsw(RNSc_half_GSW out, RNS_GSW in1, RNSc_half_GSW in2);
void gsw_automorphism_RNSc_half_gsw(RNSc_half_GSW out, RNSc_half_GSW in, uint64_t gen, RNS_RLWE_KS_Key ksk);
void gsw_RNS_to_half_RNS_gsw(RNS_half_GSW out, RNS_GSW in);
void gsw_RNS_to_half_RNSc_gsw(RNSc_half_GSW out, RNS_GSW in);
void gsw_half_RNSc_to_RNS_gsw(RNS_GSW out, RNSc_half_GSW in, RNS_RLWE_KS_Key * ksk);
void gsw_shrink_RNSc_half_gsw_sample(RNSc_GSW c);
void gsw_half_RNSc_to_RNSc_gsw(RNSc_GSW out, RNSc_half_GSW in, RNS_RLWE_KS_Key * ksk);
void gsw_RNSc_half_gsw_set_ell(RNSc_half_GSW c, uint64_t ell);
void gsw_scale_RNSc_half_gsw(RNSc_half_GSW c, uint64_t scale);
void gsw_decomp_automorphism_RNSc_half_gsw(RNSc_half_GSW out, RNSc_half_GSW in, uint64_t gen, RNS_RLWE_KS_Key * ksk);

// rlwe ks
RNS_RLWE_KS_Key rlwe_new_RNS_ks_key(RNS_RLWE_Key out_key, RNS_RLWE_Key in_key);
void rlwe_RNSc_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key ksk);
RNS_RLWE_KS_Key rlwe_new_RNS_automorphism_key(RNS_RLWE_Key key, uint64_t gen);
void rlwe_automorphism_RNSc(RNSc_RLWE out, RNSc_RLWE in, uint64_t gen, RNS_RLWE_KS_Key ksk);
void rlwe_RNSc_priv_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key kska, RNS_RLWE_KS_Key kskb);
RNS_RLWE_KS_Key * rlwe_new_RNS_priv_ks_key(RNS_RLWE_Key out_key, RNS_RLWE_Key in_key);
RNS_RLWE_KS_Key * rlwe_new_RNS_automorphism_keyset(RNS_RLWE_Key key);
RNS_RLWE_KS_Key * rlwe_new_RNS_automorphism_keyset_rou(RNS_RLWE_Key key, uint64_t root_of_unity, uint64_t N);

// hrd
void ntt_forward(uint64_t * out, uint64_t * in, uint64_t * ws, uint64_t Q, uint64_t n);
void nc_ntt_forward(uint64_t * out, uint64_t * in, uint64_t rou_2nth, uint64_t q, uint64_t n);
void calcRootsOfUnit(uint64_t * out, uint64_t min_RoU, uint64_t q, uint64_t n);
RNS_GSW * new_rlwe_bootstrap_key(RNS_GSW_Key out_key, RLWE_Key in_key, uint64_t rou_2nth);
void packed_bootstrap_wo_extract(RNSc_RLWE * out, RNSc_RLWE * tv, RLWE in, uint64_t rou_2nth, RNS_GSW * s, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk);
void rlwe_bootstrap_and_ks(RLWE out, RLWE in, RNSc_RLWE tv, RLWE_Bootstrap_KeySet bks);
void lwe_amortized_bootstrap(LWE * out, LWE * in, RNSc_RLWE tv, RLWE_Bootstrap_KeySet bks);
void test_vector_from_array(RNSc_RLWE tv, uint64_t * in, uint64_t p, uint64_t size);
RLWE_Bootstrap_KeySet new_rlwe_bootstrap_keyset(RLWE_Key in_key, RNS_GSW_Key gsw_key, LWE_Key lwe_key, uint64_t ks_l, uint64_t ks_base_bit);
void intt_reference(RNS_GSW * p, const uint64_t n, uint64_t root_of_unity, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk);


// LMK+
void blind_rotate_LMK(RNSc_RLWE tv, uint64_t * a, RNS_GSW * s, RNS_RLWE_KS_Key * ak, uint64_t size);
void bootstrap_LMK_wo_extract(RNSc_RLWE out, RNSc_RLWE tv, LWE in, RNS_GSW * s, RNS_RLWE_KS_Key * ak);
void bootstrap_LMK(LWE out, RNSc_RLWE tv, LWE in, RNS_GSW * s, RNS_RLWE_KS_Key * ak);
void bootstrap_LMK_and_ks(LWE out, RNSc_RLWE tv, LWE in, RNS_GSW * s, RNS_RLWE_KS_Key * ak, LWE_KS_Key ksk);
RNS_GSW * new_LMK_bootstrap_key(RNS_GSW_Key out_key, LWE_Key in_key);


// misc
void gen_sparse_ternary_array(uint64_t * out, uint64_t size, uint64_t h, uint64_t q);
uint64_t next_power_of_2(uint64_t x);
void array_reduce_mod_N(uint64_t * out,  uint64_t * in, uint64_t size, uint64_t p);
void array_mod_switch(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n);
void array_mod_switch_from_2k(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n);
uint64_t int_mod_switch(uint64_t in, uint64_t p, uint64_t q);
intel::hexl::NTT ** new_ntt_list(uint64_t * primes, uint64_t N, uint64_t l);
uint64_t double2int(double x);
uint64_t __debug_get_exp_message_from_noisy_RNS(RNS_Polynomial in, uint64_t * p);
void compute_RNS_Qhat_array(uint64_t * out, uint64_t * p, uint64_t l);
double RNS_compose_double(uint64_t * in, uint64_t * p, uint64_t * qi_hat, uint64_t l);
uint64_t __debug_get_exp_message_from_noisy_RNSc(RNSc_Polynomial in, uint64_t * p);
uint64_t __debug_get_exp_message_from_noisy_RNS_2(RNS_Polynomial in, uint64_t * p, uint64_t * Q_hat);
void array_additive_inverse_mod_switch(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n);
uint64_t mod_dist(uint64_t a, uint64_t b, uint64_t q);
void print_array(const char * msg, uint64_t * v, size_t size);

// Misc from third party
uint64_t RoundqQ(uint64_t v, uint64_t q, uint64_t Q);
void generate_random_bytes(uint64_t amount, uint8_t * pointer);
double generate_normal_random(double sigma);
void * safe_malloc(size_t size);
void * safe_aligned_malloc(size_t size);

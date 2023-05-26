#include "rns-gsw.h"

#define INTT_MAX_DEPTH 2
#define INTT_M 64
#define SHRINK_LEVELS 0

// Note: Decomposed automorphism will be added in a future version
#ifdef USE_DECOMP_AUTOMORPHISM
#define MACRO_rlwe_automorphism(OUT, IN, GEN) rlwe_decomp_automorphism_RNSc(OUT, IN, GEN, aut_ksk);
#define MACRO_hgsw_automorphism(OUT, IN, GEN) gsw_decomp_automorphism_RNSc_half_gsw(OUT, IN, GEN, aut_ksk);
#else
#define MACRO_rlwe_automorphism(OUT, IN, GEN) rlwe_automorphism_RNSc(OUT, IN, GEN, aut_ksk[GEN]);
#define MACRO_hgsw_automorphism(OUT, IN, GEN) gsw_automorphism_RNSc_half_gsw(OUT, IN, GEN, aut_ksk[GEN]);
#endif

static RNSc_half_GSW * acc, tmp;
// The intt algorithm does not run in place, so we need a buffer for temporary results. 
void setup_tmp_buffer(uint64_t n, uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  acc = (RNSc_half_GSW *) gsw_alloc_half_RNS_sample_array(n, Q, l, ntt);
  tmp = (RNSc_half_GSW) gsw_alloc_half_RNS_sample(Q, l, ntt);
}

// Trivial INTT. 
// In the paper: Lines 5 to 7 in Algorithm 13.
void intt_trivial(RNS_GSW * p, uint64_t * in_consts,  uint64_t * in_consts_inv, RNSc_half_GSW * acc, uint64_t * w, uint64_t * w_inv, uint64_t w_power, uint64_t n, int depth, uint64_t Ni, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  const uint64_t Q = p[0]->samples[0]->b->Q, mask_mod = Ni - 1;
  for (size_t i = 0; i < n; i++){
    gsw_RNS_to_half_RNSc_gsw(acc[i], p[0]);
    const uint64_t next_power = (w_power*i*1) & mask_mod;
    const uint64_t gen = intel::hexl::MultiplyMod(w[0], w_inv[next_power], Q);
    const uint64_t gen_c = intel::hexl::MultiplyMod(gen, in_consts[0], Q);
    const uint64_t gen_c2 = intel::hexl::MultiplyMod(gen_c, in_consts_inv[1], Q);
    MACRO_hgsw_automorphism(acc[i], acc[i], gen_c2);  
    for (size_t j = 1; j < n - 1; j++){
      const uint64_t power = (w_power*i*j) & mask_mod;
      const uint64_t next_power = (w_power*i*(j + 1)) & mask_mod;
      const uint64_t gen = intel::hexl::MultiplyMod(w[power], w_inv[next_power], Q);
      const uint64_t gen_c = intel::hexl::MultiplyMod(gen, in_consts[j], Q);
      const uint64_t gen_c2 = intel::hexl::MultiplyMod(gen_c, in_consts_inv[j+1], Q);
      gsw_mul_RNSc_half_gsw(tmp, p[j], acc[i]);
      MACRO_hgsw_automorphism(acc[i], tmp, gen_c2);
    }
    const uint64_t power = (w_power*i*(n - 1)) & mask_mod;
    const uint64_t gen_c3 = intel::hexl::MultiplyMod(w[power], in_consts[n-1], Q);
    gsw_mul_RNSc_half_gsw(tmp, p[n - 1], acc[i]);
    MACRO_hgsw_automorphism(acc[i], tmp, gen_c3);
  }
  for (size_t i = 0; i < n; i++){
    gsw_half_RNSc_to_RNS_gsw(p[i], acc[i], priv_ksk); 
  } 
}

// Recursive INTT algorithm.
// In the paper: Algorithm 13.
// Uses tv as output
void _intt_rlwe(RNSc_RLWE * tv, RNS_GSW * p, uint64_t * in_consts, uint64_t * in_consts_inv, uint64_t * out_consts, uint64_t * out_consts_inv, RNSc_half_GSW * acc, uint64_t * w, uint64_t * w_inv, uint64_t w_power, uint64_t n, int depth, uint64_t Ni, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  const uint64_t Q = p[0]->samples[0]->b->Q, mask_mod = Ni - 1;
  const uint64_t m = INTT_M, nom = n/m, l_ini = p[0]->l;
  RNS_GSW g[m][nom];
  uint64_t g_in_consts[nom], g_in_consts_inv[nom];
  for (size_t i = 0; i < m; i++){
    for (size_t j = 0; j < nom; j++){
      g[i][j] = p[m*j + i];
      g_in_consts[j] = in_consts[m*j + i];
      g_in_consts_inv[j] = in_consts_inv[m*j + i];
    }
    intt_trivial(g[i], g_in_consts, g_in_consts_inv, &acc[i*nom], w, w_inv, w_power*m, nom, depth + 1, Ni, aut_ksk, priv_ksk);
  }
  const uint64_t inv_alpha = l_ini != acc[0]->l ? acc[0]->samples[acc[0]->l]->b->ntt[acc[0]->l]->GetModulus() : 1;
  const uint64_t lm1 = 0, l2m1 = g[0][0]->l;
  printf("\nacc->l: %lu p->l %lu depth: %d inv_alpha: %lu\n", acc[0]->l, p[0]->l, depth, inv_alpha);
  for (size_t k1 = 0; k1 < nom; k1++){
    for (size_t k2 = 0; k2 < m; k2++){
      // first iteration
      const uint64_t gen = out_consts_inv[k1 + nom*k2];
      MACRO_rlwe_automorphism(tv[k1 + nom*k2], tv[k1 + nom*k2], gen);
      for (size_t i = 0; i < m - 1; i++){
        const uint64_t power = (w_power*i*k1 + w_power*nom*i*k2) & mask_mod;
        const uint64_t next_power = (w_power*(i + 1)*k1 + w_power*(nom)*(i+1)*k2) & mask_mod;
        const uint64_t gen = intel::hexl::MultiplyMod(w[power], w_inv[next_power], Q);
        gsw_mul_RNSc_rlwe(tmp->samples[lm1], g[i][k1], tv[k1 + nom*k2]);
        MACRO_rlwe_automorphism(tv[k1 + nom*k2], tmp->samples[lm1], gen);
      } 
      // last iteration
      const uint64_t power = (w_power*(m - 1)*k1 + w_power*(nom)*(m - 1)*k2) & mask_mod;
      gsw_mul_RNSc_rlwe(tmp->samples[lm1], g[m - 1][k1], tv[k1 + nom*k2]);
      const uint64_t gen_co = intel::hexl::MultiplyMod(w[power], out_consts[k1 + nom*k2], Q);
      MACRO_rlwe_automorphism(tv[k1 + nom*k2], tmp->samples[lm1], gen_co);
    }
  }
}

// Packed version of the blind rotate.
void packed_blind_rotate(RNSc_RLWE * tv, uint64_t * a_ntt, RNS_GSW * s, const uint64_t n, uint64_t root_of_unity_2nth, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  uint64_t rou[n], irou[n];
  const uint64_t Q = s[0]->samples[0]->a->Q, Ni = n, l = s[0]->l;
  intel::hexl::NTT ** ntt = s[0]->samples[0]->a->ntt;

  const uint64_t root_of_unity = intel::hexl::MultiplyMod(root_of_unity_2nth, root_of_unity_2nth, Q);
  calcRootsOfUnit(rou, root_of_unity, Q, Ni);
  calcRootsOfUnit(irou, intel::hexl::InverseMod(root_of_unity, Q), Q, Ni);

  const uint64_t inv_root_of_unity_2nth = intel::hexl::InverseMod(root_of_unity_2nth, Q);
  uint64_t in_cons_inv[n], out_cons[n], out_cons_inv[n];
  for (size_t i = 0; i < n; i++) in_cons_inv[i] = intel::hexl::InverseMod(a_ntt[i], Q);
  calcRootsOfUnit(out_cons, inv_root_of_unity_2nth, Q, n);
  const uint64_t inv_n = intel::hexl::InverseMod(n, Q);
  for (size_t i = 0; i < n; i++) out_cons[i] = intel::hexl::MultiplyMod(out_cons[i], inv_n, Q);
  for (size_t i = 0; i < n; i++) out_cons_inv[i] = intel::hexl::InverseMod(out_cons[i], Q);
  _intt_rlwe(tv, s, a_ntt, in_cons_inv, out_cons, out_cons_inv, acc, irou, rou, 1, n, 1, n, aut_ksk, priv_ksk);
}

RNS_GSW * new_rlwe_bootstrap_key(RNS_GSW_Key out_key, RLWE_Key in_key, uint64_t rou_2nth){
  assert(in_key->k == 1);
  const uint64_t N = in_key->N, Q = out_key->rlwe_key->Q, Q2 = in_key->q;
  uint64_t s_dft[N];
  array_additive_inverse_mod_switch(s_dft, in_key->s[0]->coeffs, Q2, Q, N);
  nc_ntt_forward(s_dft, s_dft, rou_2nth, Q, N);
  RNS_GSW * s = (RNS_GSW *) safe_malloc(sizeof(RNS_GSW)*N);
  for (size_t i = 0; i < N; i++){
    s[i] = (RNS_GSW) gsw_new_RNSc_sample(out_key, s_dft[i]);
    gsw_RNSc_to_RNS(s[i], (RNSc_GSW) s[i]);
  }
  return s;
} 

void packed_bootstrap_wo_extract(RNSc_RLWE * out, RNSc_RLWE * tv, RLWE in, uint64_t rou_2nth, RNS_GSW * s, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  assert(in->k == 1);
  const uint64_t n = in->a[0]->N, Q = out[0]->a->Q, Q2 = in->ntt->GetModulus(), l = s[0]->l;
  if(!acc) setup_tmp_buffer(n, Q, l, out[0]->a->ntt);
  array_mod_switch(in->b->coeffs, in->b->coeffs, Q2, Q, n);
  array_mod_switch(in->a[0]->coeffs, in->a[0]->coeffs, Q2, Q, n);
  for (size_t i = 0; i < n; i++){
    rlwe_RNSc_mul_by_xai(out[i], tv[i], Q - in->b->coeffs[i]);
  }
  uint64_t a_ntt[n];
  nc_ntt_forward(a_ntt, in->a[0]->coeffs, rou_2nth, Q, n);
  packed_blind_rotate(out, a_ntt, s, n, rou_2nth, aut_ksk, priv_ksk);
}

// generates the complete keyset for the rlwe bootstrap and ks
RLWE_Bootstrap_KeySet new_rlwe_bootstrap_keyset(RLWE_Key in_key, RNS_GSW_Key gsw_key, LWE_Key lwe_key, uint64_t ks_l, uint64_t ks_base_bit){
  // parameters
  // ls and bases
  // const uint64_t ks_l = 23, ks_base_bit = 1;
  // paramters extracted from keys
  RNS_RLWE_Key rlwe_key = gsw_key->rlwe_key;
  const uint64_t Q = rlwe_key->Q, Ni = in_key->N, l = gsw_key->l;
  const uint64_t Q2 = in_key->q;
  intel::hexl::NTT * h_ntt = in_key->ntt;
  intel::hexl::NTT ** ntt = rlwe_key->s_RNS->ntt;
  const uint64_t p0 = ntt[0]->GetModulus();
  uint64_t rou_2nth = intel::hexl::MinimalPrimitiveRoot(2*Ni, Q);
  assert(rou_2nth != 0);
  std::cout << "rou 2nth: " << rou_2nth << "\n";

  // generate RNS KS keys
  RNS_RLWE_KS_Key * priv_ksk = rlwe_new_RNS_priv_ks_key(rlwe_key, rlwe_key);
  RNS_RLWE_KS_Key * aut_ksk = rlwe_new_RNS_automorphism_keyset(rlwe_key);

  // generate RNS Bootstrapping key
  RNS_GSW * bk = new_rlwe_bootstrap_key(gsw_key, in_key, rou_2nth);

  // LWE ks key
  LWE_Key extracted_key = lwe_alloc_key(Q);
  rlwe_RNSc_extract_lwe_key(extracted_key, rlwe_key);
  extracted_key->q = Q2;
  LWE_KS_Key lwe_ksk = lwe_new_KS_key(lwe_key, extracted_key, ks_l, ks_base_bit);

  // packing key
  assert(in_key->q == lwe_key->q);
  RLWE_KS_Key packing_ksk = rlwe_new_full_packing_KS_key(in_key, lwe_key, ks_l, ks_base_bit);

  // create struct containing the keys
  RLWE_Bootstrap_KeySet res;
  res = (RLWE_Bootstrap_KeySet) safe_malloc(sizeof(*res));
  res->aut_ksk = aut_ksk;
  res->priv_ksk = priv_ksk;
  res->bk = bk;
  res->lwe_ksk = lwe_ksk;
  res->packing_ksk = packing_ksk;
  // allocate temporaries
  res->tmp_packed_in = rlwe_alloc_sample(Ni, 1, h_ntt);
  res->tmp_out = (RNSc_RLWE *) safe_malloc(sizeof(RNSc_RLWE)*Ni);
  for (size_t i = 0; i < Ni; i++) res->tmp_out[i] = (RNSc_RLWE) rlwe_new_RNS_trivial_sample_of_zero(Q, l, ntt);
  res->tmp_lwe_out = (LWE *) safe_malloc(sizeof(LWE)*Ni);
  for (size_t i = 0; i < Ni; i++) res->tmp_lwe_out[i] = lwe_alloc_sample(lwe_key->n, Q2);
  res->tmp_extracted_in = lwe_alloc_sample(Q, Q2); 
  res->s_tmp = (RNS_GSW *) safe_malloc(sizeof(RNS_GSW)*Ni);
  for (size_t i = 0; i < Ni; i++) res->s_tmp[i] = (RNS_GSW) gsw_new_RNSc_sample(gsw_key, 0);
  // save some parameters
  res->l = l;
  res->Q2 = Q2;
  res->Q = Q;
  res->rou_2nth = rou_2nth;
  return res;
}

void rlwe_bootstrap_and_ks(RLWE out, RLWE in, RNSc_RLWE tv, RLWE_Bootstrap_KeySet bks){
  const uint64_t Ni = in->b->N;
  RNSc_RLWE tv_vec[Ni];
  for (size_t i = 0; i < Ni; i++) tv_vec[i] = tv;

  // copy bk
  for (size_t i = 0; i < Ni; i++) gsw_RNS_copy(bks->s_tmp[i], bks->bk[i]);

  // run packed bootstrap
  packed_bootstrap_wo_extract(bks->tmp_out, tv_vec, in, bks->rou_2nth, bks->s_tmp, bks->aut_ksk, bks->priv_ksk);

  // mod switch to Q2, extract, and key switch to lwe_key
  for (size_t i = 0; i < Ni; i++){
    for (size_t j = 0; j < bks->l - 1 - SHRINK_LEVELS; j++){
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->a);
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->b);
    }
    rlwe_RNSc_mod_switch(bks->tmp_out[i], bks->Q2);
    rlwe_RNSc_extract_lwe(bks->tmp_extracted_in, bks->tmp_out[i], 0);
    lwe_keyswitch(bks->tmp_lwe_out[i], bks->tmp_extracted_in, bks->lwe_ksk);
    // reconstruct tmp_out buffer for the next boostrap
    bks->tmp_out[i]->a->l = bks->l;
    bks->tmp_out[i]->b->l = bks->l;
  }

  // packing
  rlwe_full_packing_keyswitch(out, bks->tmp_lwe_out, Ni, bks->packing_ksk);
}

void lwe_amortized_bootstrap(LWE * out, LWE * in, RNSc_RLWE tv, RLWE_Bootstrap_KeySet bks){
  const uint64_t Ni = bks->tmp_packed_in->b->N;

  // packing
  rlwe_full_packing_keyswitch(bks->tmp_packed_in, in, Ni, bks->packing_ksk);

  // same TV for all lwes
  RNSc_RLWE tv_vec[Ni];
  for (size_t i = 0; i < Ni; i++) tv_vec[i] = tv;

  // copy bk
  for (size_t i = 0; i < Ni; i++) gsw_RNS_copy(bks->s_tmp[i], bks->bk[i]);

  // run packed bootstrap
  packed_bootstrap_wo_extract(bks->tmp_out, tv_vec, bks->tmp_packed_in, bks->rou_2nth, bks->s_tmp, bks->aut_ksk, bks->priv_ksk);

  // mod switch to Q2, extract, and key switch to lwe_key
  for (size_t i = 0; i < Ni; i++){
    for (size_t j = 0; j < bks->l - 1 - SHRINK_LEVELS; j++){
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->a);
      polynomial_base_reduce_RNSc_wo_free(bks->tmp_out[i]->b);
    }
    rlwe_RNSc_mod_switch(bks->tmp_out[i], bks->Q2);
    rlwe_RNSc_extract_lwe(bks->tmp_extracted_in, bks->tmp_out[i], 0);
    lwe_keyswitch(out[i], bks->tmp_extracted_in, bks->lwe_ksk);
    // reconstruct tmp_out buffer for the next boostrap
    bks->tmp_out[i]->a->l = bks->l;
    bks->tmp_out[i]->b->l = bks->l;
  }
}

void test_vector_from_array(RNSc_RLWE tv, uint64_t * in, uint64_t p, uint64_t size){
  const uint64_t p0 = tv->b->ntt[0]->GetModulus();
  uint64_t pi_prod = 1;
  for (size_t i = 1; i < tv->b->l; i++) pi_prod = intel::hexl::MultiplyMod(pi_prod, tv->b->ntt[i]->GetModulus(), p0);

  const uint64_t slot_size = tv->b->Q/size;
  for (size_t i = 0; i < tv->b->Q; i++){
    const uint64_t in_ms = i/slot_size < size ? int_mod_switch(in[i/slot_size], p, p0) : 0;
    tv->b->coeffs[0][i] += intel::hexl::MultiplyMod(in_ms, pi_prod, p0);
  }   
}
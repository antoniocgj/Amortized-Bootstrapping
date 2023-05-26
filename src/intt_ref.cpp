#include "rns-gsw.h"

// rho
#define INTT_MAX_DEPTH 2
// M
#define INTT_M 32
// Use shrink
#define SHRINK

static RNSc_half_GSW * acc, tmp;
// The intt algorithm does not run in place, so we need a buffer for temporary results. 
static void setup_tmp_buffer(uint64_t n, uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  acc = (RNSc_half_GSW *) gsw_alloc_half_RNS_sample_array(n, Q, l, ntt);
  tmp = (RNSc_half_GSW) gsw_alloc_half_RNS_sample(Q, l, ntt);
}

// Trivial INTT. 
// In the paper: Lines 5 to 7 in Algorithm 13.
static void intt_trivial(RNS_GSW * p, RNSc_half_GSW * acc, uint64_t * w, uint64_t * w_inv, uint64_t w_power, uint64_t n, int depth, uint64_t Ni, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  const uint64_t Q = p[0]->samples[0]->b->Q, mask_mod = Ni - 1;
#ifdef SHRINK
  gsw_RNSc_half_gsw_set_ell(tmp, p[0]->l);
#endif
  for (size_t i = 0; i < n; i++){
    gsw_RNS_to_half_RNSc_gsw(acc[i], p[0]);
    const uint64_t next_power = (w_power*i*1) & mask_mod;
    const uint64_t gen = intel::hexl::MultiplyMod(w[0], w_inv[next_power], Q);
    gsw_automorphism_RNSc_half_gsw(acc[i], acc[i], gen, aut_ksk[next_power]);  
    for (size_t j = 1; j < n - 1; j++){
      const uint64_t power = (w_power*i*j) & mask_mod;
      const uint64_t next_power = (w_power*i*(j + 1)) & mask_mod;
      const uint64_t gen = intel::hexl::MultiplyMod(w[power], w_inv[next_power], Q);
      gsw_mul_RNSc_half_gsw(tmp, p[j], acc[i]);
      gsw_automorphism_RNSc_half_gsw(acc[i], tmp, gen, aut_ksk[(2*Ni + next_power - power)&mask_mod]);
    }
    const uint64_t power = (w_power*i*(n - 1)) & mask_mod;
    gsw_mul_RNSc_half_gsw(tmp, p[n - 1], acc[i]);
    gsw_automorphism_RNSc_half_gsw(acc[i], tmp, w[power], aut_ksk[(2*Ni - power)&mask_mod]);
  }
  for (size_t i = 0; i < n; i++){
    #ifdef SHRINK
    gsw_half_RNSc_to_RNSc_gsw((RNSc_GSW) p[i], acc[i], priv_ksk); 
    gsw_shrink_RNSc_gsw_sample((RNSc_GSW) p[i]);
    gsw_RNSc_half_gsw_set_ell(acc[i], acc[i]->l - 1);
    gsw_RNSc_to_RNS(p[i], (RNSc_GSW) p[i]);
    #else
    gsw_half_RNSc_to_RNS_gsw(p[i], acc[i], priv_ksk); 
    #endif
  } 
}

// Recursive INTT algorithm.
// In the paper: Algorithm 13.
static void _intt(RNS_GSW * p, RNSc_half_GSW * acc, uint64_t * w, uint64_t * w_inv, uint64_t w_power, uint64_t n, int depth, uint64_t Ni, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  const uint64_t Q = p[0]->samples[0]->b->Q, mask_mod = Ni - 1;
  const uint64_t m = INTT_M, nom = n/m, l_ini = p[0]->l;
  if(depth == INTT_MAX_DEPTH || (n%m)){
    intt_trivial(p, acc, w, w_inv, w_power, n, depth + 1, Ni, aut_ksk, priv_ksk);
    return;
  }
  RNS_GSW g[m][nom];
  for (size_t i = 0; i < m; i++){
    for (size_t j = 0; j < nom; j++){
      g[i][j] = p[m*j + i];
    }
    _intt(g[i], &acc[i*nom], w, w_inv, w_power*m, nom, depth + 1, Ni, aut_ksk, priv_ksk);
  }
#ifdef SHRINK
  gsw_RNSc_half_gsw_set_ell(tmp, g[0][0]->l);
#endif
  const uint64_t inv_alpha = l_ini != acc[0]->l ? acc[0]->samples[acc[0]->l]->b->ntt[acc[0]->l]->GetModulus() : 1;
  for (size_t k1 = 0; k1 < nom; k1++){
    for (size_t k2 = 0; k2 < m; k2++){
      // first iteration
      #ifdef SHRINK
      gsw_RNSc_half_gsw_set_ell(acc[k1 + nom*k2], g[0][0]->l);
      #endif
      gsw_RNS_to_half_RNSc_gsw(acc[k1 + nom*k2], g[0][k1]);
      const uint64_t next_power = (w_power*1*k1 + w_power*nom*1*k2) & mask_mod;
      const uint64_t gen = intel::hexl::MultiplyMod(w[0], w_inv[next_power], Q);
      gsw_automorphism_RNSc_half_gsw(acc[k1 + nom*k2], acc[k1 + nom*k2], gen, aut_ksk[next_power]);
      // iterations 1 to m - 2 
      for (size_t i = 1; i < m - 1; i++){
        const uint64_t power = (w_power*i*k1 + w_power*nom*i*k2) & mask_mod;
        const uint64_t next_power = (w_power*(i + 1)*k1 + w_power*(nom)*(i+1)*k2) & mask_mod;
        const uint64_t gen = intel::hexl::MultiplyMod(w[power], w_inv[next_power], Q);
        #ifdef SHRINK
        gsw_scale_RNSc_half_gsw(acc[k1 + nom*k2], inv_alpha);
        #endif
        gsw_mul_RNSc_half_gsw(tmp, g[i][k1], acc[k1 + nom*k2]);
        gsw_automorphism_RNSc_half_gsw(acc[k1 + nom*k2], tmp, gen, aut_ksk[(2*Ni + next_power - power)&mask_mod]);
      } 
      // last iteration
      const uint64_t power = (w_power*(m - 1)*k1 + w_power*(nom)*(m - 1)*k2) & mask_mod;
      #ifdef SHRINK
      gsw_scale_RNSc_half_gsw(acc[k1 + nom*k2], inv_alpha);
      #endif
      gsw_mul_RNSc_half_gsw(tmp, g[m - 1][k1], acc[k1 + nom*k2]);
      gsw_automorphism_RNSc_half_gsw(acc[k1 + nom*k2], tmp, w[power], aut_ksk[(2*Ni - power)&mask_mod]);
    }
  }
  printf("\nacc->l: %lu p->l %lu d: %d\n", acc[0]->l, p[0]->l, depth);

  // result always in (full) gsw
  for (size_t i = 0; i < n; i++){
    gsw_half_RNSc_to_RNSc_gsw((RNSc_GSW) p[i], acc[i], priv_ksk); 
    #ifdef SHRINK
    if(depth > 1){
      gsw_shrink_RNSc_gsw_sample((RNSc_GSW) p[i]);
      gsw_RNSc_half_gsw_set_ell(acc[i], acc[i]->l - 1);
    }
    #endif
    gsw_RNSc_to_RNS(p[i], (RNSc_GSW) p[i]);
  } 
}

// High level call for the reference INTT
// p: Array of GSW samples encrypting the input polynomial in point-value representation.
// n: size of p.
// root_of_unity: the n-th root of unity modulo Q.
// aut_ksk and priv_ksk: key switching keys.
void intt_reference(RNS_GSW * p, const uint64_t n, uint64_t root_of_unity, RNS_RLWE_KS_Key * aut_ksk, RNS_RLWE_KS_Key * priv_ksk){
  uint64_t rou[n], irou[n];
  const uint64_t Q = p[0]->samples[0]->a->Q, Ni = n;
  if(!acc) setup_tmp_buffer(n, Q, p[0]->l, p[0]->samples[0]->a->ntt);
  for (size_t i = 0; i < Ni; i++){
    rou[i] = intel::hexl::PowMod(root_of_unity, i, Q);
    irou[i] = intel::hexl::InverseMod(rou[i], Q);
  }
  _intt(p, acc, irou, rou, 1, n, 1, n, aut_ksk, priv_ksk);
}
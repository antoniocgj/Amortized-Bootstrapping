#include "rns-gsw.h"

void blind_rotate_LMK(RNSc_RLWE tv, uint64_t * a, RNS_GSW * s, RNS_RLWE_KS_Key * ak, uint64_t size){
  const uint64_t Q = tv->b->Q;
  const uint64_t inv_a0 = intel::hexl::InverseMod(a[0], Q);
  RNSc_RLWE rotated_tv = (RNSc_RLWE) rlwe_alloc_RNS_sample(Q, tv->a->l, tv->a->ntt);
  rlwe_automorphism_RNSc(tv, tv, inv_a0, ak[inv_a0]);
  for (size_t i = 0; i < size - 1; i++){
    const uint64_t inv_a_ip1 = intel::hexl::InverseMod(a[i+1], Q);
    const uint64_t gen = intel::hexl::MultiplyMod(a[i], inv_a_ip1, Q);
    gsw_mul_RNSc_rlwe(rotated_tv, s[i], tv);
    rlwe_automorphism_RNSc(tv, rotated_tv, gen, ak[gen]);
  }
  const uint64_t gen = a[size - 1];
  gsw_mul_RNSc_rlwe(rotated_tv, s[size - 1], tv);
  rlwe_automorphism_RNSc(tv, rotated_tv, gen, ak[gen]);
  free_RNS_rlwe_sample((RNS_RLWE)rotated_tv);
}

void bootstrap_LMK_wo_extract(RNSc_RLWE out, RNSc_RLWE tv, LWE in, RNS_GSW * s, RNS_RLWE_KS_Key * ak){
  const uint64_t p = in->q, q = out->a->Q;
  uint64_t a[in->n];
  array_mod_switch(a, in->a, p, q, in->n);
  rlwe_RNSc_mul_by_xai(out, tv, intel::hexl::SubUIntMod(out->a->Q, int_mod_switch(in->b, p, q), out->a->Q));
  blind_rotate_LMK(out, a, s, ak, in->n);
}

void bootstrap_LMK(LWE out, RNSc_RLWE tv, LWE in, RNS_GSW * s, RNS_RLWE_KS_Key * ak){
  RNSc_RLWE tmp = (RNSc_RLWE) rlwe_alloc_RNS_sample(tv->a->Q, tv->a->l, tv->a->ntt);
  bootstrap_LMK_wo_extract(tmp, tv, in, s, ak);
  for (size_t i = 0; i < tv->a->l - 1; i++) rlwe_shrink_RNSc_sample(tmp);
  rlwe_RNSc_extract_lwe(out, tmp, 0);
  free_RNS_rlwe_sample((RNS_RLWE)tmp);
}

void bootstrap_LMK_and_ks(LWE out, RNSc_RLWE tv, LWE in, RNS_GSW * s, RNS_RLWE_KS_Key * ak, LWE_KS_Key ksk){
  LWE tmp = lwe_alloc_sample(tv->a->Q, out->q);
  bootstrap_LMK(tmp, tv, in, s, ak);
  array_mod_switch(tmp->a, tmp->a, tv->a->ntt[0]->GetModulus(), out->q, tmp->n);
  tmp->b = int_mod_switch(tmp->b, tv->a->ntt[0]->GetModulus(), out->q);
  lwe_keyswitch(out, tmp, ksk);
  free_lwe_sample(tmp);
}


RNS_GSW * new_LMK_bootstrap_key(RNS_GSW_Key out_key, LWE_Key in_key){
  const uint64_t n = in_key->n;
  RNS_GSW * s = (RNS_GSW *) safe_malloc(sizeof(RNS_GSW)*n);
  uint64_t s_in[n];
  array_additive_inverse_mod_switch(s_in, in_key->s, in_key->q, out_key->rlwe_key->Q, n);
  for (size_t i = 0; i < n; i++){
    s[i] = (RNS_GSW) gsw_new_RNSc_sample(out_key, s_in[i]);
    gsw_RNSc_to_RNS(s[i], (RNSc_GSW) s[i]);
  }
  return s;
}
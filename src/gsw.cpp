#include "rns-gsw.h"

RNS_GSW_Key gsw_new_RNS_key(RNS_RLWE_Key rlwe_key){
  RNS_GSW_Key res;
  res = (RNS_GSW_Key) safe_malloc(sizeof(*res));
  res->rlwe_key = rlwe_key;
  res->l = rlwe_key->l;
  return res;
}

RNS_GSW gsw_alloc_RNS_sample(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_GSW res;
  res = (RNS_GSW) safe_malloc(sizeof(*res));
  res->samples = (RNS_RLWE *) safe_malloc(sizeof(RNS_RLWE)*2*l);
  for (size_t i = 0; i < 2*l; i++){
    res->samples[i] = rlwe_alloc_RNS_sample(Q, l, ntt);
  }
  res->l = l;
  return res;
}

RNS_half_GSW gsw_alloc_half_RNS_sample(uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_half_GSW res;
  res = (RNS_half_GSW) safe_malloc(sizeof(*res));
  res->samples = (RNS_RLWE *) safe_malloc(sizeof(RNS_RLWE)*2*l);
  for (size_t i = 0; i < l; i++){
    res->samples[i] = rlwe_alloc_RNS_sample(Q, l, ntt);
  }
  res->l = l;
  return res;
}

RNS_half_GSW * gsw_alloc_half_RNS_sample_array(uint64_t size, uint64_t Q, uint64_t l, intel::hexl::NTT ** ntt){
  RNS_half_GSW * res = (RNS_half_GSW *) safe_malloc(sizeof(RNS_half_GSW)*size);
  for (size_t i = 0; i < size; i++){
    res[i] = gsw_alloc_half_RNS_sample(Q, l, ntt);
  }
  return res;
}

void gsw_shrink_RNSc_gsw_sample(RNSc_GSW c){
  for (size_t i = 0; i < c->l - 1; i++){
    rlwe_shrink_RNSc_sample(c->samples[i]);
    rlwe_shrink_RNSc_sample(c->samples[i + c->l]);
  }
  free_RNS_rlwe_sample((RNS_RLWE) c->samples[c->l - 1]);
  free_RNS_rlwe_sample((RNS_RLWE) c->samples[2*c->l - 1]);
  for (size_t i = c->l - 1; i < 2*c->l - 2; i++){
    c->samples[i] = c->samples[i + 1];
  }
  c->l -= 1;
}

void gsw_shrink_RNSc_half_gsw_sample(RNSc_half_GSW c){
  for (size_t i = 0; i < c->l - 1; i++){
    rlwe_shrink_RNSc_sample(c->samples[i]);
  }
  free_RNS_rlwe_sample((RNS_RLWE) c->samples[c->l - 1]);
  c->l -= 1;
}

void gsw_RNSc_half_gsw_set_ell(RNSc_half_GSW c, uint64_t ell){
  for (size_t i = 0; i < ell; i++){
    c->samples[i]->a->l = ell;
    c->samples[i]->b->l = ell;
  }
  c->l = ell;
}

void gsw_RNSc_to_RNS(RNS_GSW out, RNSc_GSW in){
  const uint64_t l = in->samples[0]->a->l;
  for (size_t i = 0; i < 2*l; i++){
    rlwe_RNSc_to_RNS(out->samples[i], in->samples[i]);
  }
}

void gsw_RNS_to_RNSc(RNSc_GSW out, RNS_GSW in){
  const uint64_t l = in->samples[0]->a->l;
  for (size_t i = 0; i < 2*l; i++){
    rlwe_RNS_to_RNSc(out->samples[i], in->samples[i]);
  }
}

/* Returns X^m */
RNSc_GSW gsw_new_RNSc_sample(RNS_GSW_Key key, uint64_t m){
  const uint64_t l = key->rlwe_key->l;
  RNSc_GSW res = (RNSc_GSW) gsw_alloc_RNS_sample(key->rlwe_key->Q, key->rlwe_key->l, key->rlwe_key->s_RNS->ntt);
  for (size_t i = 0; i < l; i++){
    const uint64_t p = res->samples[i]->a->ntt[i]->GetModulus();
    // a
    rlwe_RNSc_sample_of_zero(res->samples[i], key->rlwe_key);
    res->samples[i]->a->coeffs[i][m] = intel::hexl::AddUIntMod(res->samples[i]->a->coeffs[i][m], 1, p);
    // b
    rlwe_RNSc_sample_of_zero(res->samples[i + l], key->rlwe_key);
    res->samples[i + l]->b->coeffs[i][m] = intel::hexl::AddUIntMod(res->samples[i + l]->b->coeffs[i][m], 1, p);
  }
  return res;
}

/* Returns GSW(-rlwe_key)*/
// RNSc_GSW gsw_new_RNSc_ks_key(RNS_GSW_Key key){
//   const uint64_t l = key->rlwe_key->l;
//   RNSc_GSW res = (RNSc_GSW) gsw_alloc_RNS_sample(key->rlwe_key->Q, key->rlwe_key->l, key->rlwe_key->s_RNS->ntt);
//   for (size_t i = 0; i < l; i++){
//     const uint64_t p = key->rlwe_key->s_RNS->ntt[i]->GetModulus();
//     // a
//     rlwe_RNSc_sample_of_zero(res->samples[i], key->rlwe_key);
//     for (size_t j = 0; j < key->rlwe_key->Q; j++){
//       res->samples[i]->a->coeffs[i][j] += intel::hexl::SubUIntMod(0, key->rlwe_key->s->coeffs[j], p);
//     }
    
//     // b
//     rlwe_RNSc_sample_of_zero(res->samples[i + l], key->rlwe_key);
//     for (size_t j = 0; j < key->rlwe_key->Q; j++){
//       res->samples[i + l]->b->coeffs[i][j] += intel::hexl::SubUIntMod(0, key->rlwe_key->s->coeffs[j], p);
//     }
//   }
//   return res;
// }

void gsw_mul_RNSc_rlwe(RNSc_RLWE out, RNS_GSW in1, RNSc_RLWE in2){
  assert(in1->l >= in2->a->l);
  assert(out != in2);
  const uint64_t l = in2->a->l;
  RNSc_Polynomial tmp = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in2->a->Q, l, in2->a->ntt);
  RNSc_Polynomial tmp2 = (RNSc_Polynomial) polynomial_new_RNS_polynomial(in2->a->Q, l, in2->a->ntt);
  // static uint64_t max_l = l;
  // assert(l <= max_l);
  // tmp->l = l;
  // tmp2->l = l;

  const uint64_t p = in2->a->ntt[0]->GetModulus();
  polynomial_base_extend_RNSc(tmp, in2->a->coeffs[0], p);
  polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
  for (size_t j = 0; j < l; j++){
    const uint64_t q = tmp->ntt[j]->GetModulus();
    intel::hexl::EltwiseMultMod(out->a->coeffs[j], in1->samples[0]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
    intel::hexl::EltwiseMultMod(out->b->coeffs[j], in1->samples[0]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
  }

  // a
  for (size_t i = 1; i < l; i++){
    const uint64_t p = in2->a->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in2->a->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();     
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], in1->samples[i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], in1->samples[i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }
  // b
  for (size_t i = 0; i < l; i++){
    const uint64_t p = in2->b->ntt[i]->GetModulus();
    polynomial_base_extend_RNSc(tmp, in2->b->coeffs[i], p);
    polynomial_RNSc_to_RNS((RNS_Polynomial) tmp, tmp);
    for (size_t j = 0; j < l; j++){
      const uint64_t q = tmp->ntt[j]->GetModulus();
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], in1->samples[tmp->l + i]->a->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->a->coeffs[j], out->a->coeffs[j], tmp2->coeffs[j], tmp->N, q);
      intel::hexl::EltwiseMultMod(tmp2->coeffs[j], in1->samples[tmp->l + i]->b->coeffs[j], tmp->coeffs[j], tmp->N, q, 1);
      intel::hexl::EltwiseAddMod(out->b->coeffs[j], out->b->coeffs[j], tmp2->coeffs[j], tmp->N, q);
    }
  }

  polynomial_RNS_to_RNSc(out->a, (RNS_Polynomial) out->a);
  polynomial_RNS_to_RNSc(out->b, (RNS_Polynomial) out->b);
  polynomial_RNSc_reduce_mod_XQ_minus_1(out->a);
  polynomial_RNSc_reduce_mod_XQ_minus_1(out->b);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
}


void gsw_mul_RNSc_gsw(RNSc_GSW out, RNS_GSW in1, RNSc_GSW in2){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < 2*l; i++){
    gsw_mul_RNSc_rlwe(out->samples[i], in1, in2->samples[i]);
  }
}

void gsw_scale_RNSc_gsw(RNSc_GSW c, uint64_t scale){
  for (size_t i = 0; i < 2*c->l; i++){
    for (size_t j = 0; j < c->l; j++){
      const uint64_t mod = c->samples[0]->a->ntt[j]->GetModulus();
      for (size_t k = 0; k < c->samples[i]->a->Q; k++){
        c->samples[i]->a->coeffs[j][k] = intel::hexl::MultiplyMod(c->samples[i]->a->coeffs[j][k], scale, mod);
        c->samples[i]->b->coeffs[j][k] = intel::hexl::MultiplyMod(c->samples[i]->b->coeffs[j][k], scale, mod);
      } 
    }
  }
}


void gsw_scale_RNSc_half_gsw(RNSc_half_GSW c, uint64_t scale){
  for (size_t i = 0; i < c->l; i++){
    for (size_t j = 0; j < c->l; j++){
      const uint64_t mod = c->samples[0]->a->ntt[j]->GetModulus();
      for (size_t k = 0; k < c->samples[i]->a->Q; k++){
        c->samples[i]->a->coeffs[j][k] = intel::hexl::MultiplyMod(c->samples[i]->a->coeffs[j][k], scale, mod);
        c->samples[i]->b->coeffs[j][k] = intel::hexl::MultiplyMod(c->samples[i]->b->coeffs[j][k], scale, mod);
      } 
    }
  }
}


void gsw_RNS_to_half_RNS_gsw(RNS_half_GSW out, RNS_GSW in){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    rlwe_copy_RNS_sample(out->samples[i], in->samples[i + l]);
  }
}

void gsw_RNS_copy(RNS_GSW out, RNS_GSW in){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < 2*l; i++){
    rlwe_copy_RNS_sample(out->samples[i], in->samples[i]);
  }
}

void gsw_RNS_to_half_RNSc_gsw(RNSc_half_GSW out, RNS_GSW in){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    rlwe_RNS_to_RNSc(out->samples[i], in->samples[i + l]);
  }
}

void gsw_half_RNSc_to_RNS_gsw(RNS_GSW out, RNSc_half_GSW in, RNS_RLWE_KS_Key * ksk){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    rlwe_RNSc_to_RNS(out->samples[i + l], in->samples[i]);
    rlwe_RNSc_priv_keyswitch((RNSc_RLWE) out->samples[i], in->samples[i], ksk[0], ksk[1]);
    rlwe_RNSc_to_RNS(out->samples[i], (RNSc_RLWE) out->samples[i]);
  }
}

void gsw_half_RNSc_to_RNSc_gsw(RNSc_GSW out, RNSc_half_GSW in, RNS_RLWE_KS_Key * ksk){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    rlwe_copy_RNS_sample((RNS_RLWE) out->samples[i + l], (RNS_RLWE) in->samples[i]);
    rlwe_RNSc_priv_keyswitch(out->samples[i], in->samples[i], ksk[0], ksk[1]);
  }
}


void gsw_mul_RNSc_half_gsw(RNSc_half_GSW out, RNS_GSW in1, RNSc_half_GSW in2){
  assert(out != in2);
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    gsw_mul_RNSc_rlwe(out->samples[i], in1, in2->samples[i]);
  }
}

void gsw_automorphism_RNSc_half_gsw(RNSc_half_GSW out, RNSc_half_GSW in, uint64_t gen, RNS_RLWE_KS_Key ksk){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    rlwe_automorphism_RNSc(out->samples[i], in->samples[i], gen, ksk);  
  }
}

void gsw_decomp_automorphism_RNSc_half_gsw(RNSc_half_GSW out, RNSc_half_GSW in, uint64_t gen, RNS_RLWE_KS_Key * ksk){
  const uint64_t l = out->samples[0]->a->l;
  for (size_t i = 0; i < l; i++){
    rlwe_decomp_automorphism_RNSc(out->samples[i], in->samples[i], gen, ksk);  
  }
}



#include "rns-gsw.h"
#include "benchmark_util.h"

void test_LMK_wo_extract(){
  // Parameters
  const uint64_t N = 896, l = 2, bit_size = 49, p = 12289, key_hw = 256; 
  const uint64_t N_prime = next_power_of_2(p)<<1;
 
  // Generate Primes
  const uint64_t p_star = 1965443; // a random 21-bit prime for the output modulus
  auto Q = intel::hexl::GeneratePrimes(l, bit_size, true, N_prime);

  
  // Rings
  intel::hexl::NTT ** ntt = new_ntt_list(Q.data(), N_prime, l);

  // generate keys
  std::cout << "Generating private keys...\n";
  LWE_Key in_key = lwe_new_sparse_ternary_key(N, p_star, key_hw, 1);
  RNS_RLWE_Key rlwe_key = rlwe_new_RNS_gaussian_key(p, l, 3.2, ntt, 1);
  RNS_GSW_Key gsw_key = gsw_new_RNS_key(rlwe_key);
  LWE_Key extracted_key = lwe_alloc_key(p);
  rlwe_RNSc_extract_lwe_key(extracted_key, rlwe_key);
  extracted_key->q = p_star;

  // Generate Input
  LWE in = lwe_new_sample(int_mod_switch(5000, p, p_star), in_key);
  printf("in: %ld\n", int_mod_switch(lwe_phase(in, in_key), p_star, p));
 
  // generate bootstrapping keys
  std::cout << "Generating bootstrapping keys...\n";
  RNS_GSW * bk = new_LMK_bootstrap_key(gsw_key, in_key);
  RNS_RLWE_KS_Key * aut_ksk = rlwe_new_RNS_automorphism_keyset(rlwe_key);
  LWE_KS_Key ksk = lwe_new_KS_key(in_key, extracted_key, 21, 1);

  // generate test_vector
  uint64_t lut[p];
  for (size_t i = 0; i < p; i++){
    lut[i] = i;
  }
  RNSc_RLWE tv = (RNSc_RLWE) rlwe_new_RNS_trivial_sample_of_zero(p, l, ntt);
  test_vector_from_array(tv, lut, p, p);

  // create output
  LWE out = lwe_alloc_sample(N, p_star);

  // Run LMK 10 times and gets the average
  std::cout << "Bootstrapping...\n";
  MEASURE_TIME("LMK", 10, "LMK", 
    bootstrap_LMK_and_ks(out, tv, in, bk, aut_ksk, ksk);
  );

  std::cout << "Decrypting. \n";
  // Decrypt
  std::cout << "Out: " << int_mod_switch(lwe_phase(out, in_key), p_star, p) << "\n";
}

int main(int argc, char const *argv[])
{
  test_LMK_wo_extract();
  return 0;
}

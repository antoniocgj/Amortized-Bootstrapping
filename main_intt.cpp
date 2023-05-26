#include "rns-gsw.h"
#include "benchmark_util.h"

void test_intt_reference(){
  // Parameters
  const uint64_t N = 512, l = 4, bit_size = 49, p = 12289;
  const uint64_t N_prime = next_power_of_2(p)<<1;
  std::cout << "p: " << p << "\n";
  std::cout << "N': " << N_prime << "\n";
  // Generate primes
  auto Q = intel::hexl::GeneratePrimes(l, bit_size, true, N_prime);
  std::cout << "p: ";
  for (size_t i = 0; i < l; i++){
    std::cout << Q.data()[i] << ", ";
  }
  std::cout << "\n";

  // Ring
  intel::hexl::NTT ** ntt = new_ntt_list(Q.data(), N_prime, l);

  // Generate private key
  RNS_RLWE_Key rlwe_key = rlwe_new_RNS_gaussian_key(p, l, 3.2, ntt, 2);
  RNS_GSW_Key gsw_key = gsw_new_RNS_key(rlwe_key);

  // Prepare Clear-text Forward NTT
  intel::hexl::NTT * h_ntt = new intel::hexl::NTT(N, p);
  uint64_t rou = h_ntt->GetMinimalRootOfUnity();
  assert(rou != 0);
  std::cout << "2N-th RoU: " << rou << "\n";
  rou = intel::hexl::MultiplyMod(rou, rou, p);
  std::cout << "N-th RoU: " << rou << "\n";

  const uint64_t inv_n = intel::hexl::InverseMod(N, p);

  // Generate evaluation keys
  RNS_RLWE_KS_Key * priv_ksk = rlwe_new_RNS_priv_ks_key(rlwe_key, rlwe_key);
  RNS_RLWE_KS_Key * aut_ksk = rlwe_new_RNS_automorphism_keyset_rou(rlwe_key, rou, N);
  
  // Generate Input
  uint64_t a[N], a_dft[N], ws[N];
  for (size_t i = 0; i < N; i++){
    a[i] = i;
  }
  calcRootsOfUnit(ws, rou, p, N);
  ntt_forward(a_dft, a, ws, p, N);

  // encrypt
  std::cout << "Encrypting input. \n";
  RNS_GSW a_enc[N];
  for (size_t i = 0; i < N; i++){
    a_enc[i] = (RNS_GSW) gsw_new_RNSc_sample(gsw_key, a_dft[i]);
    gsw_RNSc_to_RNS(a_enc[i], (RNSc_GSW) a_enc[i]);
  }

  // intt
  std::cout << "INTT. \n";
  MEASURE_TIME("Intt", 1, "INTT execution time", 
    intt_reference(a_enc, N, rou, aut_ksk, priv_ksk);
  );

  // decrypt
  std::cout << "Decrypting. \n";
  RNS_Polynomial m = polynomial_new_RNS_polynomial(p, a_enc[0]->l, ntt);
  uint64_t res, err = 0;
  printf("a: ");
  for (size_t i = 0; i < N; i++){
    rlwe_RNS_phase(m, a_enc[i]->samples[a_enc[0]->l], rlwe_key);
    res = __debug_get_exp_message_from_noisy_RNS(m, Q.data());
    res = intel::hexl::MultiplyMod(res, inv_n, p);
    printf("%lu ", res);
    if(i!=res) err++;
  }
  printf("\n");
  if(err) printf("Failed: %lu\n", err);
  else printf("Pass.\n");
}

int main(int argc, char const *argv[])
{
  test_intt_reference();
  return 0;
}

#include "rns-gsw.h"
#include "benchmark_util.h"

void test_amortized_bootstrap(){
  // Parameters
  const uint64_t N = 1024, l = 3, bit_size = 49, p = 12289, p_start_size = 24, key_hw = 256;
  const uint64_t N_prime = next_power_of_2(p)<<1;

  // Generate primes
  auto p_star_list = intel::hexl::GeneratePrimes(1, p_start_size, true, N);
  const uint64_t p_star = p_star_list[0]; // p* is the output modulus
  auto Q = intel::hexl::GeneratePrimes(l, bit_size, true, N_prime);
  
  // Rings
  intel::hexl::NTT ** ntt = new_ntt_list(Q.data(), N_prime, l);
  intel::hexl::NTT * h_ntt = new intel::hexl::NTT(N, p_star);

  // generate keys
  std::cout << "Generating private keys...\n";
  RLWE_Key in_key = rlwe_new_sparse_ternary_key(N, p_star, 1, key_hw, 1, h_ntt);
  RNS_RLWE_Key rlwe_key = rlwe_new_RNS_gaussian_key(p, l, 3.2, ntt, 1);
  RNS_GSW_Key gsw_key = gsw_new_RNS_key(rlwe_key);
  LWE_Key lwe_key = lwe_new_sparse_ternary_key(N, p_star, key_hw, 1);


  // input message precision
  const uint64_t input_prec = 8,
                 input_base = 1<<input_prec,
                 message_scale = p_star/input_base;

  // generate input
  std::cout << "Generating Input...\n";
  LWE in[N];
  // Message is multiplied by the scaling. 
  // We add scale/2 so that the bootstrap evaluates a rounding instead of a floor operation.
  for (size_t i = 0; i < N; i++){
    in[i] = lwe_new_sample((i%input_base)*message_scale + message_scale/2, lwe_key);
  }

  // generate test_vector
  uint64_t lut[input_base];
  for (size_t i = 0; i < input_base; i++){
    lut[i] = i*message_scale;
  }
  RNSc_RLWE tv = (RNSc_RLWE) rlwe_new_RNS_trivial_sample_of_zero(p, l, ntt);
  test_vector_from_array(tv, lut, p_star, input_base);

  // generate bootstrapping key
  std::cout << "Generating bootstrapping keys...\n";
  RLWE_Bootstrap_KeySet bks = new_rlwe_bootstrap_keyset(in_key, gsw_key, lwe_key, p_start_size, 1);

  // run rlwe bootstrap
  std::cout << "Bootstrapping...\n";
  MEASURE_TIME("", 1, "RLWE bootstrap", 
    lwe_amortized_bootstrap(in, in, tv, bks);
  )

  // decrypt
  std::cout << "Decrypting...\n";
  uint64_t m_out[N];

  // round
  for (size_t i = 0; i < N; i++){
    m_out[i] = ((lwe_phase(in[i], lwe_key) + message_scale/2) / message_scale)%input_base;
  } 
  // print result
  print_array("Out:", m_out, N);
  // calculate error standard deviation
  double err = 0;
  for (size_t i = 0; i < N; i++){
    err += pow(mod_dist(m_out[i], i%input_base, input_base), 2);
  }
  printf("Error: 2^%f \n", log2(sqrt(err/N)));

  // Add an offset = message_scale/2 to bootstrap again. We add scale/2 so that the bootstrap evaluates a rounding instead of a floor operation.
  for (size_t i = 0; i < N; i++){
    in[i]->b = intel::hexl::AddUIntMod(in[i]->b, message_scale/2, p_star);
  }

  // run rlwe bootstrap again
  std::cout << "Bootstrapping again...\n";
  MEASURE_TIME("", 1, "RLWE bootstrap", 
    lwe_amortized_bootstrap(in, in, tv, bks);
  )

  // decrypt
  std::cout << "Decrypting...\n";

  // round
  for (size_t i = 0; i < N; i++){
    m_out[i] = ((lwe_phase(in[i], lwe_key) + message_scale/2) / message_scale)%input_base;
  } 
  // print result
  print_array("Out:", m_out, N);
  // calculate error standard deviation
  err = 0;
  for (size_t i = 0; i < N; i++){
    err += pow(mod_dist(m_out[i], i%input_base, input_base), 2);
  }
  printf("Error: 2^%f \n", log2(sqrt(err/N)));
}

int main(int argc, char const *argv[])
{
  // In this test, we use the amortized bootstrapping to run an RLWE bootstrapping because we want measure error (and probability of failure) after the packing. See the file "src/amortized_boostrapping.cpp" for running the amortized boostrapping without the repacking to RLWE at the end. 
  test_amortized_bootstrap();
  return 0;
}

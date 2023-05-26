# Amortized Bootstrapping Revisited: Simpler, Asymptotically-faster, Implemented

This repository holds the source code for reproducing the results of the paper *[Amortized bootstrapping revisited: Simpler, asymptotically-faster, implemented](https://eprint.iacr.org/2023/014)*.

## Citation

```
@misc{cryptoeprint:2023/014,
      author = {Antonio Guimarães and Hilder V. L. Pereira and Barry van Leeuwen},
      title = {Amortized Bootstrapping Revisited: Simpler, Asymptotically-faster, Implemented},
      howpublished = {Cryptology ePrint Archive, Paper 2023/014},
      year = {2023},
      note = {\url{https://eprint.iacr.org/2023/014}},
      url = {https://eprint.iacr.org/2023/014}
}
```

## Examples

We provide three main examples:

- [main_non_amortized](./main_non_amortized.cpp): The non-amortized version of our bootstrapping. It's essentially a CLWE RNS variant of the [LMK+22](https://eprint.iacr.org/2022/198.pdf) bootstrapping method.
- [main_amortized](./main_amortized.cpp): Our amortized bootstrapping. This implementation specializes in rho = 2.
- [main_rlwe](./main_rlwe.cpp): An RLWE bootstrapping using our amortized bootstrapping. We use it to measure the noise and probability of failure for the amortized bootstrapping. 
- [main_intt](./main_intt.cpp): Executes a reference implementation of our homomorphic INTT. This implementation accepts arbitrary rho and shrinking. Check the parameters at the beginning of [intt_ref.cpp](./src/intt_ref.cpp).

## Compiling and Running

Notice that this implementation does not use any compression techniques. Thus, memory requirements for our examples are from **80GB up to 150 GB**, depending on the example and the parameter set.

1. Compile [Intel HEXL](https://github.com/intel/hexl):
```console
$ make hexl
```
2. Compile this code:
```console
$ make
```
3. Run our examples:
```console
$ ./main_non_amortized
$ ./main_amortized
$ ./main_rlwe
$ ./main_intt
```
**NOTE**: In case of Assertion Failure, please run the code again. This proof-of-concept implementation does not treat a few corner cases (e.g. when NTT(a) contains a zero). These corner cases do not affect performance and would be easy to treat in a production implementation. 

For reproducing the results of the paper, we recommend using an [m6i.metal](https://instances.vantage.sh/aws/ec2/m6i.metal) instance or similar.

## License

[Apache License Version 2.0](./LICENSE)

We include code from the following third-party libraries:

- [Intel HEXL](https://github.com/intel/hexl): [Apache License 2.0](https://github.com/intel/hexl/blob/development/LICENSE), Copyright 2020 Intel Corporation
- [FIPS202 from Kyber](https://github.com/pq-crystals/kyber/blob/master/ref/fips202.c): [Public Domain](https://creativecommons.org/share-your-work/public-domain/cc0/)

And small snippets of code from:
- [MOSFHET](https://github.com/antoniocgj/MOSFHET): [Apache License 2.0](https://github.com/antoniocgj/MOSFHET/blob/main/LICENSE), Copyright 2022 Antonio Guimarães et al.
- [OpenFHE](https://github.com/openfheorg/openfhe-development): [BSD-2-Clause license](https://github.com/openfheorg/openfhe-development/blob/main/LICENSE), Copyright (c) 2022, OpenFHE

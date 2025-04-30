# TwoFloat: Double-Word Arithmetics in C++

Double-word arithmetic is a technique used to represent a real number as an unevaluated sum of two floating-point numbers, commonly referred to as "double-double" arithmetic. This representation allows to represent numbers with twice the precision of the underlying floating-point type. 

This C++-library implements two recently proposed double-word arithmetics.
- `twofloat::doubleword`: Double-word arithmetic  proposed by Joldes et al. 2017 ([Tight and rigorous error bounds for basic building blocks of double-word arithmetic](https://doi.org/10.1145/3121432))
- `twofloat::pair`: Pair arithmetic proposed by Lange and Rump in 2020 ([Faithfully Rounded Floating-point Computations](https://doi.org/10.1145/3290955)). This arithmetic skips the normalization step after each operation, which makes it faster than the double-word arithmetic. However, faithful rounding is only guaranteed up to a certain number of operations on the same number.


## FMA support
Some operations that involve multiplications or divisions can benefit from fused multiply-add (FMA) instructions. **If this instruction is available on the target architecture, the FMA version of an operation should be preferred for better performance and accuracy.** If FMA instructions are unavailable, FMA implementations should not be used because emulating FMA instructions is considerably slower than using non-FMA implementations. 

FMA instructions have been available on Intel CPUs since Haswell (2013) and AMD CPUs since Piledriver (2012). However, many FPUs in embedded CPUs do not support FMA instructions.


## Usage
This library is a header-only library requiring at least C++17. To use it, include the header file of the desired arithmetic:
```cpp
// Pair arithmetic by Lange and Rump in 2020
#include <twofloat/arithmetics/pair-arithmetic.hpp>

// Double-word arithmetic by Joldes et al. 2017
#include <twofloat/arithmetics/double-word-arithmetic.hpp>
```

All arithmetics operate on the `two<T>` type. This type is templated with the type of the underlying floating point numbers. The operators of each arithmetic are defined in their respective namespace. The following code shows how to add two numbers using the pair arithmetic:

```cpp
using namespace twofloat;

two<float> a(1.0);   // a = { 1.0,  0.0 }
two<float> b(1e-9);  // b = { 1e-9, 0.0 }
auto c = pair::add(a, b); // c = { 1.0, 1e-9 }
```

For some operations, different algorithms are implemented. The algorithm is chosen by template parameters. For example, the `doubleword::mul` operation provides FMA and non-FMA implementations. When using a non-FMA implementation, the user must choose between a fast and an accurate algorithm. The documentation provides more information about the different algorithms that each operation implements. The following code shows how to select different algorithms:
    
```cpp
using namespace twofloat;

two<float> a(1.0);
two<float> b(1.0);

// Fast FMA implementation
two<float> c = doubleword::mul<doubleword::Mode::Fast, true>(a, b);

// Accurate FMA implementation
two<float> d = doubleword::mul<doubleword::Mode::Accurate, true>(a, b);

// Fast non-FMA implementation
two<float> e = doubleword::mul<doubleword::Mode::Fast, false>(a, b);
```

We do not provide operator overloads because different algorithms are implemented for some operations, and we do not want to choose a default algorithm for the user.

## Precision 
Generally speaking, the precision of the double-word arithmetic is twice the precision of the underlying floating-point type. However, the range stays the same. The following table shows the precision of the different arithmetics: 
| Type | Base | Precision Bits | Roundoff Error Unit | Decimal digits | Range |
| ---- | --------- | -------------- | ------------------- | -------------- | ----- | 
| `float` | β=2 | p=23+1 | u<sub>float</sub> = β<sup>1-p</sup> = 5.96e-8 | log<sub>10</sub>(2<sup>p-1</sup>) = 6.9 | 1.2e-38 to 3.4e38 |
| `double` | β=2 | p=52+1 | u<sub>double</sub> = β<sup>1-p</sup> = 1.11e-16 | log<sub>10</sub>(2<sup>p-1</sup>) = 15.7 | 2.2e-308 to 1.8e308 |
| `two<float>` | β=2 | p=2*24 | u<sub>float</sub><sup>2</sup> = 3.55e-15 |  log<sub>10</sub>(2<sup>p-1</sup>) = 14.1 | 1.2e-38 to 3.4e38 |
| `two<double>` | β=2 | p=2*53 | u<sub>double</sub><sup>2</sup> = 1.23e-31 | log<sub>10</sub>(2<sup>p-1</sup>) = 31.6 | 2.2e-308 to 1.8e308 |


This information is available in C++ using the specializations of the `std::numeric_limits` template provided in `libtwofloat/limits.hpp`:
```cpp
#include <twofloat/limits.hpp>

std::numeric_limits<two<float>>::digits10; // digits10 (two<float>) = 14
std::numeric_limits<two<double>>::digits10; // digits10 (two<double>) = 31
```

## Runtime and error bounds
### Double-word arithmetic (Joldes et al. 2017)
The double-word arithmetic by Joldes et al. provides error bounds for each operation. The error bounds are given in units u of the roundoff error of the underlying floating-point type (see table above). For example, when using `two<float>`, u is equal to u<sub>float</sub>. 

The number of floating-point operations (FP ops) required for each operation is different to Table 1 in Joldes et al. 2017 for several reasons:
- We take negations and comparisons into account
- The non-FMA algorithms use the non-FMA version of the `TwoProd` algorithm, which requires significantly more FP ops than the FMA version (called `FastTwoProd`). This is more realistic because if an FMA is available, the user will likely use the FMA version.

| Operation | FMA | Mode | Error bound formally proved | # of FP ops | Name in Joldes et al. 2017 | 
| --------- | --- | ---- | ----------- | --------------------------- | -------------------------- |
| DW **+** FP | No | | 2u<sup>2</sup> | 10 | Algorithm 4 (DWPlusFP) |
| DW **+** DW | No | Sloppy | N/A | 11 | Algorithm 5 (SloppyDWPlusDW) |
| DW **+** DW | No | Accurate | 3u<sup>2</sup>+13u<sup>3</sup> | 20 | Algorithm 6 (AccurateDWPlusDW) |
| DW **x** FP | No | Accurate | 1.5u<sup>2</sup>+4u<sup>3</sup>  | 29 | Algorithm 7 (DWTimesFP1) |
| DW **x** FP | No | Fast | 3u<sup>2</sup> | 23 | Algorithm 8 (DWTimesFP2) |
| DW **x** FP | Yes | Accurate | 2u<sup>2</sup> | 7 | Algorithm 9 (DWTimesFP3) |
| DW **x** DW | No | Fast | 7u<sup>2</sup> | 28 | Algorithm 10 (DWTimesDW1) |
| DW **x** DW | Yes | Fast | 6u<sup>2</sup> | 9 | Algorithm 11 (DWTimesDW2) |
| DW **x** DW | Yes | Accurate | 5u<sup>2</sup> | 10 | Algorithm 12 (DWTimesDW3) |
| DW **:** FP | No |  | 3u<sup>2</sup> | 29 | Algorithm 15 (DWDivFP3) |
| DW **:** FP | Yes |  | 3u<sup>2</sup> | 11 | Algorithm 15 (DWDivFP3) |
| DW **:** DW | No | Fast | 15u<sup>2</sup>+56u<sup>3</sup> | 36 | Algorithm 17 (DWDivDW2) |
| DW **:** DW | Yes | Fast | 15u<sup>2</sup>+56u<sup>3</sup> | 14 | Algorithm 17 (DWDivDW2) |
| DW **:** DW | Yes | Accurate | 9.8u<sup>2</sup> | 34 | Algorithm 18 (DWDivDW3) |

DW: Double-word (`two<T>`), FP: Floating-point (`T`)

### Pair arithmetic (Lange and Rump 2020)

| Operation | FMA | # of FP ops | Name in Lange and Rump 2020 |
| --------- | --- | ----------- | -------------------------- |
| DW **+** FP | No | 7 | derived from CPairSum |
| DW **+** DW | No | 8 | CPairSum |
| DW **x** FP | No | 23 | derived from CPairMul |
| DW **x** FP | Yes | 5 | derived from CPairMul |
| DW **x** DW | No | 25 | CPairMul |
| DW **x** DW | Yes | 7 | CPairMul |
| DW **:** DW | No | 8 | CPairDiv |

## Documentation
The documentation can be build using Doxygen. We are working on providing a hosted version of the documentation.

## Testing
The repository contains some basic tests for the different arithmetics. To run them, checkout the repository and use the following commands:

```bash
cmake -B build
cmake --build build
cmake --build build --target test
```

## Runtime of basic algorithms
| Algorithm | # of FP ops |
| --------- | ----------- |
| Fast2Sum | 3 |
| TwoSum | 6 |
| Split | 6 |
| TwoProd w/o FMA | 9 + 2*Split = 21 |
| TwoProd w/ FMA | 3 |

## Errata
`pair::sub` seems to be unprecise. Tested in the double precision matrix residual implementation.

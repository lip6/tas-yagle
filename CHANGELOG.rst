
===================
Tas/Yagle Changelog
===================


GCC 15 compaibility -- July 27, 2025
====================================

Fix compilation of old code under GCC 15 (for SUSE Tumbleweed).

Note: Use "feature_test_macros" to have the list of GLIBC/GCC
      configuration flags.

This code is old and have lots of contradictory constraints to compile.
* vfork() -> needs _POSIX_C_SOURCE < 200809L (ANSI C C99)
* signal/SA_RESTART -> non-POSIX.
* caddr_t -> deprecated in newer ANSI norms.
* inline -> needs C99.
* C++ comments (//) in C code -> needs C99.
* C23 is no longer backward compatible with non-prototype functions
  declarations.

To have the old code working, we need the compiler to:
* Use the -std=gnu99 old dialect.
* Do not choke on callbacks parameters: -Wno-incompatible-pointer-types.

Mofifications that were made to the code, while looking for a solution,
but which may not have been needed :

* Change back C++ style comment into C (// comment -> /* comment */).
* Write complete function prototypes, completed up to Xyagle (in order
  of compilation). Would be great to fully finish it.
* Remove of inline call. This I may need to revert if performance issue
  arise.


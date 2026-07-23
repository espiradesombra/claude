//
// prime_sieve.c
//
// Copyright (C) July 2002, Tomás Oliveira e Silva
//
// e-mail: tos@ua.pt
// www:    http://www.ieeta.pt/~tos
//
// Comparison of two simple (but serious) implementations of the segmented sieve of
// Eratosthenes. The second implementation can generate primes reasonably fast near
// 1e18 (using around 400Mbytes of memory).
//
// _implementation_=0 gives a classical segmented sieve
// _implementation_=1 gives a cache-friendly segmented sieve
//
// See timing results for the two implementations at the end.
//
// Main idea: use one linked list for each interval of the segmented sieve, putting in it
//   the primes than have an odd multiple in that interval (but not in a previous interval);
//   this allows a better utilization of the processor data caches, giving significant time
//   savings (up to a factor of 6) when working near 1e18. The amount of memory used is
//   approximately 8*pi(sqrt(N)), where N is the last number of the interval, and pi(x) is
//   the usual prime counting function.
//
// Assumptions: pointers have 4 bytes, gcc compiler
//
//
// Released under the GNU general public license (version 2 or any later version); see the
// gpl.txt file (or the page http://www.gnu.org/licenses/gpl.html) for details.
//
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


//
// configuration specification
//
// default parameters optimized for integers near 1e18; see tables at the end
//

#ifndef _implementation_
# define _implementation_     1
#endif
#ifndef _sieve_bits_log2_
# if _implementation_ == 0
#  define _sieve_bits_log2_  25
# else
#  define _sieve_bits_log2_  19
# endif
#endif
#ifndef _bucket_size_log2_
# if _implementation_ == 0
#  define _bucket_size_log2_ 20
# else
#  define _bucket_size_log2_ 10
# endif
#endif

#if _implementation_ == 1
# if _bucket_size_log2_ > 16
#  error "_bucket_size_log2_ is too large"
# endif
#endif

//
// basic type definitions
//

typedef unsigned char      u08;
typedef unsigned int       u32;
typedef unsigned long long u64;


//
// memory allocation
//

static void *get_memory(u32 size)
{
  size_t m;

  m = (size_t)malloc(size + 255); // this assumes that sizeof(void *) = sizeof(size_t)
  if((void *)m == NULL)
    exit(1);
  m = (m + (size_t)255) & ~(size_t)255;
  return (void *)m; // pointer aligned on a 256 byte boundary
}


//
// count the number of zeros
//

static u32 count_zero_bits(u08 *addr,u32 size)
{
  static u32 data[256];
  u32 i,j;

  if(data[1] == 0)
    for(i = 0;i < 256;i++)
      for(j = i ^ 255;j;j >>= 1)
        if(j & 1)
          data[i]++;
  j = 0;
  for(i = 0;i < size;i++)
    j += data[(u32)addr[i] & 255];
  return j;
}


//
// generation of the (small) primes used by the main sieve
//

#define number_of_small_primes 6541
static u32 small_primes[number_of_small_primes];
static u32 small_sieve[1024];
static u32 small_base;

static void update_small_sieve(void)
{
  u32 i,j;

  if(small_primes[0] == 0)
  { // initialize the small_primes array
    for(j = 0;j < 1024;j++)
      small_sieve[j] = 0;
    for(i = 3;i < 256;i += 2)
      if((small_sieve[i >> 6] & (1 << ((i >> 1) & 31))) == 0)
        for(j = (i * i) >> 1;j < 32768;j += i)
          small_sieve[j >> 5] |= 1 << (j & 31);
    j = 0;
    for(i = 3;i < 65536;i += 2)
      if((small_sieve[i >> 6] & (1 << ((i >> 1) & 31))) == 0)
        small_primes[j++] = i;
    if(j != number_of_small_primes)
      exit(2); // this should never happen
  }
  for(j = 0;j < 1024;j++)
    small_sieve[j] = 0;
  for(i = 0;i < number_of_small_primes;i++)
  {
    j = small_primes[i] * small_primes[i];
    if(j >= small_base + 65536)
      break;
    if(j < small_base)
    {
      j = small_base / small_primes[i];
      j *= small_primes[i];
      if(j < small_base)
        j += small_primes[i];
      if((j & 1) == 0)
        j += small_primes[i];
    }
    for(j = (j - small_base) >> 1;j < 32768;j += small_primes[i])
      small_sieve[j >> 5] |= 1 << (j & 31);
  }
}

#if 0

static void test_small_sieve(void)
{
#define _test_small_sieve_
  u32 pi,i,j;

  pi = 0;
  small_base = 0;
  i = 1000000;
  while(small_base < 4000000)
  {
    update_small_sieve();
    j = (i - small_base) >> 4;
    if(j > 4096)
      j = 4096;
    pi += count_zero_bits((u08 *)small_sieve,j);
    if(small_base + (j << 4) == i)
    {
      printf("%9u %u\n",i,pi);
      i += 1000000;
    }
    if(j < 4096)
      pi += count_zero_bits((u08 *)small_sieve + j,4096 - j);
    small_base += 65536;
  }
}

#endif


//
// main sieve
//
// the following structure is used to record the
// information required to sieve an interval
//
// when _implementation_ = 0 the value of _bucket_size_log2_ should
//   be rather large
// when _implementation_ = 1 the value of _bucket_size_log2_ should
//   be small (and a multiple of the L1 or L2 data cache line size)
//

#define primes_per_bucket ((1 << (_bucket_size_log2_ - 3)) - 1)

typedef struct bucket
{
  struct bucket *next; // pointer to next bucket
  u32 count;           // count of the number of primes in this bucket
  struct
  {
    u32 p; // prime
    u32 o; // the bit number of the first odd multiple (>= main_base) of the prime
  }
  data[primes_per_bucket];
}
bucket;

static u32 main_sieve[1 << (_sieve_bits_log2_ - 5)];
static u64 main_base,main_limit;
static u32 next_prime;

#if _implementation_ == 0
static bucket *main_list;
#else
static bucket **main_lists,*available_buckets;
static u32 list_size_log2,current_list;
#endif

#if _implementation_ == 0
# define new_bucket() do { bucket *b; b = get_memory(sizeof(bucket)); \
    b->next = main_list; b->count = 0; main_list = b; } while(0)
#else
# define more_buckets() do { u32 i,j; i = 1 << (20 - _bucket_size_log2_); \
    available_buckets = (bucket *)get_memory(i * sizeof(bucket)); for(j = 0;j < i;j++) \
    available_buckets[j].next = (j < i - 1) ? &available_buckets[j + 1] : NULL; } while(0)
# define new_bucket(k) do { bucket *b; if(available_buckets == NULL) more_buckets(); \
    b = available_buckets; available_buckets = available_buckets->next; \
    b->next = main_lists[k]; main_lists[k] = b; b->count = 0; } while(0)
#endif

static void init_main_sieve(void)
{
  u64 t,end;
  u32 i,j;
#if _implementation_ == 1
  u32 k;
#endif

  if(next_prime == 0)
  {
#if _implementation_ == 0
    main_list = NULL;
    new_bucket();
#else
    i = 1 + (u32)ceil(sqrt((double)main_limit));
    i = 2 + (i >> _sieve_bits_log2_);
    for(list_size_log2 = 2;(1 << list_size_log2) < i;list_size_log2++)
      ;
    current_list = 0;
    available_buckets = NULL;
    main_lists = (bucket **)get_memory((1 << list_size_log2) * sizeof(bucket *));
    for(k = 0;k < (1 << list_size_log2);k++)
    {
      main_lists[k] = NULL;
      new_bucket(k);
    }
#endif
    small_base = 0;
    update_small_sieve();
    next_prime = 3;
  }
  end = main_base + (u64)(2 << _sieve_bits_log2_);
  while((t = (u64)next_prime * (u64)next_prime) < end)
  {
    if(next_prime >= small_base + 65536)
    {
      small_base += 65536;
      update_small_sieve();
    }
    i = (next_prime - small_base) >> 1;
    if((small_sieve[i >> 5] & (1 << (i & 31))) == 0)
    {
      if(t < main_base)
      {
        t = main_base / (u64)next_prime;
        t *= (u64)next_prime;
        if(t < main_base)
          t += (u64)next_prime;
        if(((u32)t & 1) == 0)
          t += (u64)next_prime;
      }
      i = (u32)((t - main_base) >> 1); // bit number
#if _implementation_ == 0
      if(main_list->count == primes_per_bucket)
        new_bucket();
      j = main_list->count++;
      main_list->data[j].p = next_prime;
      main_list->data[j].o = i;
#else
      k = (current_list + (i >> _sieve_bits_log2_)) & ((1 << list_size_log2) - 1);
      if(main_lists[k]->count == primes_per_bucket)
        new_bucket(k);
      j = main_lists[k]->count++;
      main_lists[k]->data[j].p = next_prime;
      main_lists[k]->data[j].o = i & ((1 << _sieve_bits_log2_) - 1);
#endif
    }
    next_prime += 2;
  }
}

static void do_main_sieve(void)
{
  bucket *b;
#if _implementation_ == 1
  bucket *c;
  u32 j,k;
#endif
  u32 i,p,o;

  init_main_sieve();
  for(i = 0;i < (1 << (_sieve_bits_log2_ - 5));i++)
    main_sieve[i] = 0;
#if _implementation_ == 0
  for(b = main_list;b != NULL;b = b->next)
    for(i = 0;i < b->count;i++)
    {
      p = b->data[i].p;
      for(o = b->data[i].o;o < (1 << _sieve_bits_log2_);o += p)
        main_sieve[o >> 5] |= 1 << (o & 31);
      b->data[i].o = o - (1 << _sieve_bits_log2_);
    }
#else
  b = main_lists[current_list];
  while(b != NULL)
  {
    for(i = 0;i < b->count;i++)
    {
      p = b->data[i].p;
      for(o = b->data[i].o;o < (1 << _sieve_bits_log2_);o += p)
        main_sieve[o >> 5] |= 1 << (o & 31);
      k = (current_list + (o >> _sieve_bits_log2_)) & ((1 << list_size_log2) - 1);
      if(main_lists[k]->count == primes_per_bucket)
        new_bucket(k);
      j = main_lists[k]->count++;
      main_lists[k]->data[j].p = p;
      main_lists[k]->data[j].o = o & ((1 << _sieve_bits_log2_) - 1);
    }
    c = b;
    b = b->next;
    c->next = available_buckets;
    available_buckets = c;
  }
  main_lists[current_list] = NULL;
  new_bucket(current_list);
  current_list = (current_list + 1) & ((1 << list_size_log2) - 1);
#endif
}

#if 0

static void test_main_sieve(void)
{
#define _test_main_sieve_
#define _first_   0ull
#define _last_    4ull
#define _step_    1000000000ull
#define _init_pi_ 0ull
  u64 pi,i;
  u32 j,k;

  if(sizeof(bucket) != (1 << _bucket_size_log2_))
    exit(3);
  k = 1 << (_sieve_bits_log2_ - 3);
  if((u32)_step_ & 15 || _step_ < ((u64)k << 4))
    exit(4);
  pi = _init_pi_;
  main_base = _first_ * _step_;
  main_limit = (_last_ + 1ull) * _step_;
  i = (_first_ + 1ull) * _step_;
  next_prime = 0;
  while(main_base < _last_ * _step_)
  {
    do_main_sieve();
    j = (u32)(i - main_base) >> 4;
    if(j > k)
      j = k;
    pi += (u64)count_zero_bits((u08 *)main_sieve,j);
    if(main_base + (u64)(j << 4) == i)
    {
      printf("%10llu %llu\n",i,pi);
      fflush(stdout);
      i += _step_;
    }
    if(j < k)
      pi += (u64)count_zero_bits((u08 *)main_sieve + j,k - j);
    main_base += (u64)k << 4;
  }
  exit(0);
#undef _first_
#undef _last_
#undef _step_
#undef _init_pi_
}

#endif


//
// main program
//

int main(int argc,char **argv)
{
  double t;
  u32 i,j;
  u64 pi;

#ifdef _test_small_sieve_
  test_small_sieve();
#endif
#ifdef _test_main_sieve_
  test_main_sieve();
#endif

  if(argc == 1)
    i = 15;
  else
    i = atoi(argv[1]);
  if(i < 6)
    i = 6;
  if(i > 18)
    i = 18;

  printf("%u %2u %2u",_implementation_,_sieve_bits_log2_,_bucket_size_log2_);
  main_base = 1ull;
  for(j = 0;j < i;j++)
    main_base *= 10ull;
  main_limit = main_base + 2000000000ull;
  next_prime = 0;
  printf(" %2d",i);
  t = (double)clock();
  init_main_sieve();
  t = ((double)clock() - t) / (double)CLOCKS_PER_SEC;
  printf(" %6.2f",t);
  j = 1 << (_sieve_bits_log2_ - 3);
  pi = 0ull;
  main_limit = main_base + 1000000000ull;
  if(((u32)main_base | (u32)main_limit) & 63)
  {
    fprintf(stderr,"Warning: prime number counts may be incorrect\n");
    fprintf(stderr,"         main_base and main_limit should be multiples of 64\n");
  }
  t = (double)clock();
  for(;;)
  {
    do_main_sieve();
    i = (u32)(main_limit - main_base) >> 4;
    if(i <= j)
      break;
    pi += (u64)count_zero_bits((u08 *)main_sieve,j);
#if 0
    //
    // example code to print the prime numbers between
    // main_base and main_base+2*j
    //
    {
      u32 k;

      for(k = 0;k < (j << 3);k++)
        if((main_sieve[k >> 5] & (1 << (k & 31))) == 0)
          printf("%llu\n",main_base + (u32)(2 * k + 1));
    }
#endif
    main_base += (u64)j << 4;
  }
  pi += (u64)count_zero_bits((u08 *)main_sieve,i);
  t = ((double)clock() - t) / (double)CLOCKS_PER_SEC;
  printf(" %7.2f %8llu\n",t,pi);
  return 0;
}


//
// speed measurements on an 900MHz Athlon (64k two-way L1 data cache)
//
// I  = _implementation_
// SB = _sieve_bits_log2_
// BS = _bucket_size_log2_
//
// Best times (in seconds) to sieve an interval of 1e9 integers starting at
// 1e12, 1e14, 1e16 and 1e18; the last column gives, for comparison, the time
// used by the author's fastest implementation of a segmented sieve (in the
// same processor), using a mod 30 wheel and assembly language
//
//   n I SB BS   time      n I SB BS   time     fastest
//  -- - -- -- ------     -- - -- -- ------     -------
//  12 0 19 20  13.31     12 1 19 10  14.08        2.63
//  14 0 21 20  25.79     14 1 19 10  20.62        4.05
//  16 0 21 20  97.39     16 1 19 10  26.26        5.58
//  18 0 25 20 208.13     18 1 20 10  32.36        9.61
//
//
// _implementation_=0,  _bucket_size_log2_=20
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   25.07  141.98 1166.44       ?      n SB   time
//   18   17.25   76.04  584.33       ?     -- -- ------
//   19   13.31   43.36  304.62       ?     12 19  13.31
//   20   14.36   30.33  163.52 1282.99     14 21  25.79
//   21   15.47   25.79   97.39  672.26     16 21  97.39
//   22   61.43   73.58  116.73  408.58     18 25 208.13
//   23   87.25  103.27  130.35  278.88
//   24  103.92  123.04  145.30  225.99
//   25  117.13  138.90  159.45  208.13
//   26  130.87  155.76  176.47  211.21
// ---- ------- ------- ------- -------
// init    0.04    0.37   3.73    39.14
// ---- ------- ------- ------- -------
//
// _implementation_=1,  _bucket_size_log2_=6
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   19.52   28.44   37.21   54.70      n SB   time
//   18   17.30   25.41   32.80   43.15     -- -- ------
//   19   14.83   22.70   29.66   37.35     12 19  14.83
//   20   15.28   22.97   29.90   36.91     14 19  22.70
//   21   16.02   24.77   34.28   43.45     16 19  29.66
//   22   61.54   75.82   88.95  102.21     18 20  36.99
// ---- ------- ------- ------- -------
// init    0.05    0.42    4.17   42.98
// ---- ------- ------- ------- -------
//
// _implementation_=1,  _bucket_size_log2_=7
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   18.44   26.69   34.78   53.05      n SB   time
//   18   16.56   24.07   30.65   40.64     -- -- ------
//   19   14.42   21.60   27.73   34.65     12 19  14.42
//   20   15.04   22.05   28.13   34.29     14 19  21.60
//   21   15.90   23.67   31.59   39.30     16 19  27.73
//   22   61.37   74.94   86.86   98.15     18 20  34.29
// ---- ------- ------- ------- -------
// init    0.04    0.41    4.04   41.83
// ---- ------- ------- ------- -------
//
// _implementation_=1,  _bucket_size_log2_=8
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   17.83   25.68   33.57   52.16      n SB   time
//   18   16.14   23.32   29.49   39.40     -- -- ------
//   19   14.20   21.05   26.79   33.30     12 19  14.20
//   20   14.95   21.66   27.33   32.98     14 19  21.05
//   21   15.87   23.12   30.32   37.31     16 19  26.79
//   22   61.68   74.56   85.69   96.41     18 20  32.98
// ---- ------- ------- ------- -------
// init    0.04    0.42    4.07   42.29
// ---- ------- ------- ------- -------
//
// _implementation_=1,  _bucket_size_log2_=9
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   17.64   25.36   33.25   52.22      n SB   time
//   18   15.97   22.97   28.98   38.99     -- -- ------
//   19   14.12   20.76   26.38   32.82     12 19  14.12
//   20   14.96   21.41   26.99   32.41     14 19  20.76
//   21   15.91   22.83   29.78   36.36     16 19  26.38
//   22   61.07   73.91   84.65   93.54     18 20  32.41
// ---- ------- ------- ------- -------
// init    0.04    0.41    4.04   42.19
// ---- ------- ------- ------- -------
//
// _implementation_=1,  _bucket_size_log2_=10
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   17.59   25.42   33.35   52.50      n SB   time
//   18   15.89   22.85   28.92   38.89     -- -- ------
//   19   14.08   20.62   26.26   32.75     12 19  14.08
//   20   14.91   21.26   26.86   32.36     14 19  20.62
//   21   15.84   22.78   29.51   35.93     16 19  26.26
//   22   61.19   74.00   84.71   95.03     18 20  32.36
// ---- ------- ------- ------- -------
// init    0.04    0.41    4.03   42.32
// ---- ------- ------- ------- -------
//
// _implementation_=1,  _bucket_size_log2_=11
//
//   SB   t(12)   t(14)   t(16)   t(18)        best data
// ---- ------- ------- ------- -------
//   17   17.75   25.85   33.61   53.11      n SB   time
//   18   15.91   22.90   29.12   39.09     -- -- ------
//   19   14.07   20.59   26.30   33.01     12 19  14.07
//   20   14.91   21.27   26.73   32.37     14 19  20.59
//   21   15.96   22.78   29.58   36.13     16 19  26.30
//   22   60.37   74.06   85.05   94.91     18 20  32.37
// ---- ------- ------- ------- -------
// init    0.04    0.41    4.01   43.19
// ---- ------- ------- ------- -------
//

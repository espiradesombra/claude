//
// fast_sieve.c
//
// Copyright (C) August 2010, Tomás Oliveira e Silva
//
// e-mail: tos@ua.pt
// www:    http://www.ieeta.pt/~tos
//
// simple single-threaded cache-friendly implementation of the segmented sieve of Eratosthenes
// capable of generating primes up to 2^64
//
// main idea: use one linked list for each interval of the segmented sieve, putting in it the
// primes than have an odd multiple in that interval (but not in a previous interval); this uses
// the processor data caches in a better way, giving significant time savings when working near
// 1e18. The amount of memory used is approximately 8*pi(sqrt(N)), where N is the last number of
// the interval, and pi(x) is the usual prime counting function.
//
// assumptions: Linux, GNU libc, C99-capable compiler, linear address space
//
// see some speed measurements at the end of the source code
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


//
// sieve configuration
//

// number of bytes of a pointer (either 4 or 8)
#ifndef _pointer_size_
# if defined(__i386__)
#  define _pointer_size_      4u
# elif defined(__x86_64__)
#  define _pointer_size_      8u
# else
#  define _pointer_size_      4u
# endif
#endif

// base-2 logarithm of the number of bits of the sieve array
// should be adjusted for better performance (see some speed measurements at the end)
// this default gives very good performance of an Intel Core2 Quad Q9550 processor
#ifndef _sieve_bits_log2_
# define _sieve_bits_log2_   23u
#endif

// base-2 logarithm of the size (in bytes) of the bucket structure
#ifndef _bucket_size_log2_
# define _bucket_size_log2_  12u
#endif

// threshold used to separate the small from the large primes (can be set to 0xFFFFFFFFu)
// primes larger than _prime_threshold_ are guaranteed to have at most one multiple in the
// sieve array (use 0xFFFFFFFFu for the classical segmented Eratosthenes sieve)
#ifndef _prime_threshold_
# define _prime_threshold_  (1u << _sieve_bits_log2_)
#endif
#if _prime_threshold_ < (1u << _sieve_bits_log2_)
# error "_prime_threshold_ cannot be smaller than (1u << _sieve_bits_log2_)"
#endif
#if _prime_threshold_ == 0xFFFFFFFFu
# undef _prime_threshold_
# define _prime_threshold_ 0xFFFFFFFEu
#endif


//
// useful definitions and data types
//

typedef unsigned      char u08;  //  8-bit unsigned integer
typedef unsigned       int u32;  // 32-bit unsigned integer
typedef unsigned long long u64;  // 64-bit unsigned integer
typedef             double r64;  // 64-bit floating point number

#define _sieve_span_  (2ull << _sieve_bits_log2_)
#if _pointer_size_ == 4u
  typedef u32 uXX; // for a 32-bit architecture
# define _sieve_words_  (1u << (_sieve_bits_log2_ - 5u))
#elif _pointer_size_ == 8u
  typedef u64 uXX; // for a 64-bit architecture
# define _sieve_words_  (1u << (_sieve_bits_log2_ - 6u))
#else
# error "bad _pointer_size_ value"
#endif


//
// support stuff
//

typedef struct memory_block
{
  struct memory_block *next_block;
  void *address;
}
memory_block;

static memory_block *all_memory_blocks = NULL;

static void *get_memory(u32 size,u32 alignment)
{
  memory_block *mb;
  void *a;
  uXX m;

  if((u32)sizeof(void *) != _pointer_size_)
  { // sanity check
    fprintf(stderr,"get_memory: bad pointer size\n");
    exit(1);
  }
  if(alignment != (alignment & -alignment) || alignment < 64)
  { // the requested alignment must be a power of 2 (we also want it to be at least 64)
    fprintf(stderr,"get_memory: bad alignment requested\n");
    exit(1);
  }
  if(size < alignment)
    size = alignment;
  if((size & 63u) != 0u)
  { // the requested size must be a multiple of 64
    fprintf(stderr,"get_memory: bad size requested\n");
    exit(1);
  }
  a = malloc((size_t)(size + --alignment) + sizeof(memory_block));
  if(a == NULL)
  {
    fprintf(stderr,"get_memory: out of memory\n");
    exit(1);
  }
  m = ((uXX)a + (uXX)alignment) & ~(uXX)alignment;
  mb = (memory_block *)(m + size);
  mb->next_block = all_memory_blocks;
  mb->address = a;
  all_memory_blocks = mb;
  return (void *)m;
}

static void free_all_memory(void)
{
  memory_block *mb1,*mb2;

  for(mb1 = all_memory_blocks;mb1 != NULL;mb1 = mb2)
  {
    mb2 = mb1->next_block;
    free(mb1->address);
  }
  all_memory_blocks = NULL;
}

static r64 cpu_time(void)
{
  struct rusage r;
  r64 t;

  getrusage(RUSAGE_SELF,&r);
  t = (r64)(r.ru_utime.tv_sec + r.ru_stime.tv_sec);
  t += 1.0e-6 * (r64)(r.ru_utime.tv_usec + r.ru_stime.tv_usec);
  return t;
}

static u32 count_zero_bits(u08 *addr,u32 n_bytes)
{ //                        00           01           10           11
#define EXPAND1(x)          x+2,         x+1 ,        x+1 ,        x+0
#define EXPAND2(x)  EXPAND1(x+2),EXPAND1(x+1),EXPAND1(x+1),EXPAND1(x+0)
#define EXPAND3(x)  EXPAND2(x+2),EXPAND2(x+1),EXPAND2(x+1),EXPAND2(x+0)
#define EXPAND4(x)  EXPAND3(x+2),EXPAND3(x+1),EXPAND3(x+1),EXPAND3(x+0)
  static const u32 data[256] = { EXPAND4(0) };
#undef EXPAND1
#undef EXPAND2
#undef EXPAND3
#undef EXPAND4
  u32 i,j;

  for(i = j = 0u;i < n_bytes;i++)
    j += data[(u32)addr[i]];
  return j;
}

#if _pointer_size_ == 4u
static const u32 mark_mask[32u] =
{
  0x00000001u,0x00000002u,0x00000004u,0x00000008u,
  0x00000010u,0x00000020u,0x00000040u,0x00000080u,
  0x00000100u,0x00000200u,0x00000400u,0x00000800u,
  0x00001000u,0x00002000u,0x00004000u,0x00008000u,
  0x00010000u,0x00020000u,0x00040000u,0x00080000u,
  0x00100000u,0x00200000u,0x00400000u,0x00800000u,
  0x01000000u,0x02000000u,0x04000000u,0x08000000u,
  0x10000000u,0x20000000u,0x40000000u,0x80000000u
};
# define mark_1(s,o)  (s)[(o) >> 5u] |= (1u << ((o) & 31u))
# define test_1(s,o)  (((s)[(o) >> 5u] & (1u << ((o) & 31u))) == 0u)
# define mark_2(s,o)  (s)[(o) >> 5u] |= mark_mask[(o) & 31u]
# define test_2(s,o)  (((s)[(o) >> 5u] & mark_mask[(o) & 31u]) == 0u)
# define unmark(s,o)  (s)[(o) >> 5u] &= ~mark_mask[(o) & 31u]
#elif _pointer_size_ == 8u
static const u64 mark_mask[64u] =
{
  0x0000000000000001ull,0x0000000000000002ull,0x0000000000000004ull,0x0000000000000008ull,
  0x0000000000000010ull,0x0000000000000020ull,0x0000000000000040ull,0x0000000000000080ull,
  0x0000000000000100ull,0x0000000000000200ull,0x0000000000000400ull,0x0000000000000800ull,
  0x0000000000001000ull,0x0000000000002000ull,0x0000000000004000ull,0x0000000000008000ull,
  0x0000000000010000ull,0x0000000000020000ull,0x0000000000040000ull,0x0000000000080000ull,
  0x0000000000100000ull,0x0000000000200000ull,0x0000000000400000ull,0x0000000000800000ull,
  0x0000000001000000ull,0x0000000002000000ull,0x0000000004000000ull,0x0000000008000000ull,
  0x0000000010000000ull,0x0000000020000000ull,0x0000000040000000ull,0x0000000080000000ull,
  0x0000000100000000ull,0x0000000200000000ull,0x0000000400000000ull,0x0000000800000000ull,
  0x0000001000000000ull,0x0000002000000000ull,0x0000004000000000ull,0x0000008000000000ull,
  0x0000010000000000ull,0x0000020000000000ull,0x0000040000000000ull,0x0000080000000000ull,
  0x0000100000000000ull,0x0000200000000000ull,0x0000400000000000ull,0x0000800000000000ull,
  0x0001000000000000ull,0x0002000000000000ull,0x0004000000000000ull,0x0008000000000000ull,
  0x0010000000000000ull,0x0020000000000000ull,0x0040000000000000ull,0x0080000000000000ull,
  0x0100000000000000ull,0x0200000000000000ull,0x0400000000000000ull,0x0800000000000000ull,
  0x1000000000000000ull,0x2000000000000000ull,0x4000000000000000ull,0x8000000000000000ull
};
# define mark_1(s,o)  (s)[(o) >> 6u] |= (1ull << ((o) & 63u))
# define test_1(s,o)  (((s)[(o) >> 6u] & (1ull << ((o) & 63u))) == 0ull)
# define mark_2(s,o)  (s)[(o) >> 6u] |= mark_mask[(o) & 63u]
# define test_2(s,o)  (((s)[(o) >> 6u] & mark_mask[(o) & 63u]) == 0ull)
# define unmark(s,o)  (s)[(o) >> 6u] &= ~mark_mask[(o) & 63u]
#endif

static uXX pattern[3u * 5u * 7u * 11u * 13u];  // sieve initialization pattern

static void init_pattern(void)
{
  u32 i;

  for(i = 0u;i < 3u * 5u * 7u * 11u * 13u;i++)
    pattern[i] = (uXX)0u;
  for(i = (3u >> 1u);i < 3u * 5u * 7u * 11u * 13u * 8u * _pointer_size_;i += 3u)
    mark_2(pattern,i);
  for(i = (5u >> 1u);i < 3u * 5u * 7u * 11u * 13u * 8u * _pointer_size_;i += 5u)
    mark_2(pattern,i);
  for(i = (7u >> 1u);i < 3u * 5u * 7u * 11u * 13u * 8u * _pointer_size_;i += 7u)
    mark_2(pattern,i);
  for(i = (11u >> 1u);i < 3u * 5u * 7u * 11u * 13u * 8u * _pointer_size_;i += 11u)
    mark_2(pattern,i);
  for(i = (13u >> 1u);i < 3u * 5u * 7u * 11u * 13u * 8u * _pointer_size_;i += 13u)
    mark_2(pattern,i);
}

static u64 convert(u32 m,u32 n)
{
  u64 c;

  for(c = (u64)m;n > 0u;n--)
  {
    if(c > (0xFFFFFFFFFFFFFFFFull / 10ull))
    {
      fprintf(stderr,"convert: arithmetic overflow\n");
      exit(0);
    }
    c *= 10ull;
  }
  return c;
}


//
// the bucket data structure
//
// each bucket must be properly aligned; its starting address MUST be a multiple of the size of
// the bucket (which itself must be a power of two)
//
// to avoid touching non-consecutive memory, when a bucket is partially full, only the pointer to
// the next free data[] entry is kept; the memory alignment of each bucket makes possible a very
// simple and efficient way to test when a bucket is full; there is no need to keep track of the
// number of stored data[] entries in each bucket
//

#define _primes_per_bucket_  ((1u << (_bucket_size_log2_ - 3u)) - 1u)
#define _bucket_addr_mask_   ((uXX)((1u << _bucket_size_log2_) - 1u))

typedef struct
{
  u32 p;  // prime number
  u32 o;  // base+2*o+1 is the first odd multiple of p larger than base
}
po;

typedef struct bucket
{
  struct bucket *next;           // pointer to the end of the next bucket
#if _pointer_size_ == 4u
  u32 dummy;                     // to ensure that everything is properly aligned
#endif
  po data[_primes_per_bucket_];  // prime and bit number data for this bucket
}
bucket;

static bucket *available_buckets = NULL; // list of free buckets

static void more_buckets(void)
{
  u32 i,j;

  i = (u32)sizeof(bucket);
  if(i != (1u << _bucket_size_log2_))
  { // sanity check
    fprintf(stderr,"more_buckets: bad bucket size\n");
    exit(1);
  }
  if(available_buckets != NULL)
  { // another sanity check
    fprintf(stderr,"more_buckets: no need to call this function now\n");
    exit(1);
  }
  if((uXX)&available_buckets->data[0] != (uXX)8u)
  { // final sanity check
    fprintf(stderr,"more_buckets: bad data offset\n");
    exit(1);
  }
  j = 1u << (24u - _bucket_size_log2_); // number of buckets to allocate (16MiB)
  available_buckets = (bucket *)get_memory(j * i,i);
  for(i = 0u;i < j - 1u;i++)
    available_buckets[i].next = &available_buckets[i + 1u];
  available_buckets[i].next = NULL;
}


//
// auxiliary segmented Eratosthenes sieve used to generate odd primes up to 2^32
//

#define _number_of_aux_primes_  6536u
#define _aux_sieve_span_        (1u << 19u)
#define _aux_sieve_words_       (_aux_sieve_span_ / (2u * 8u * _pointer_size_))
#if _aux_sieve_span_ < 65536u
# error "_aux_sieve_span_ is too small"
#endif

static u32 aux_primes[_number_of_aux_primes_];  // the odd primes used in the auxiliary sieve
static uXX aux_sieve[_aux_sieve_words_];        // the sieve data (one bit per odd integer)
static u32 aux_base;                            // must be a multiple of _aux_sieve_span_

static void update_aux_sieve(void)
{
  u32 i,j,k,o;

  if((aux_base & (_aux_sieve_span_ - 1u)) != 0u)
  { // sanity check
    fprintf(stderr,"update_aux_sieve: illegal aux_base\n");
    exit(1);
  }
  // if necessary, initialize the aux_primes array
  if(aux_primes[0u] == 0u)
  {
    for(i = 0u;i < 65536u / (2u * 8u * _pointer_size_);i++)
      aux_sieve[i] = (uXX)0u;
    for(i = 3u;i < 256u;i += 2u)
    {
      j = i >> 1u;
      if(test_2(aux_sieve,j))
        for(j = (i * i) >> 1u;j < 32768u;j += i)
          mark_2(aux_sieve,j);
    }
    j = 0u;
    for(i = (17u >> 1u);i < 32768u;i++) // start at 17
      if(test_2(aux_sieve,i))
        aux_primes[j++] = 2u * i + 1u;
    if(j != _number_of_aux_primes_ || aux_primes[_number_of_aux_primes_ - 1u] != 65521u)
    {
      fprintf(stderr,"update_aux_sieve: initialization of the aux_primes array failed\n");
      exit(1);
    }
    init_pattern();
  }
  // initialize the aux_sieve array
  o = aux_base % (3u * 5u * 7u * 11u * 13u);
# if _pointer_size_ == 4u
  o = (o + ((o * 41u) & 63u) * 3u * 5u * 7u * 11u * 13u) >> 6u;   // 41 = -1/15015 mod 64
# elif _pointer_size_ == 8u
  o = (o + ((o * 105u) & 127u) * 3u * 5u * 7u * 11u * 13u) >> 7u; // 105 = -1/15015 mod 128
# endif
  for(i = 0u;i < _aux_sieve_words_;i += k)
  {
    j = _aux_sieve_words_ - i;        // remaining sieve words
    k = 3u * 5u * 7u * 11u * 13u - o; // remaining pattern words
    if(j < k)
      k = j;
    for(j = 0;j < k;j++)
      aux_sieve[i + j] = pattern[o + j];
    o = 0u;
  }
  if(aux_base == 0u)
  {
    aux_sieve[0] |= mark_mask[0]; // mark 1 as not prime, and mark 3, 5, 7, 11, and 13 as prime
    aux_sieve[0] &= ~(mark_mask[1] | mark_mask[2] | mark_mask[3] | mark_mask[5] | mark_mask[6]);
  }
  // mark the multiples of the primes larger than 13
  for(i = 0u;i < _number_of_aux_primes_;i++)
  {
    j = aux_primes[i] * aux_primes[i];
    if(j > aux_base + (_aux_sieve_span_ - 1u))
      break;
    if(j > aux_base)
      j = (j - aux_base) >> 1u;
    else
    {
      j = aux_primes[i] - aux_base % aux_primes[i];
      if((j & 1u) == 0u)
        j += aux_primes[i];
      j >>= 1u;
    }
    while(j < (_aux_sieve_span_ / 2u))
    {
      mark_2(aux_sieve,j);
      j += aux_primes[i];
    }
  }
}

static void test_aux_sieve(void)
{
  static const u32 pi[32u] =
  { // pi(2^k)
          1u,      2u,      4u,       6u,      11u,      18u,       31u,       54u,
         97u,    172u,    309u,     564u,    1028u,    1900u,     3512u,     6542u,
      12251u,  23000u,  43390u,   82025u,  155611u,  295947u,   564163u,  1077871u,
    2063689u,3957809u,7603553u,14630843u,28192750u,54400028u,105097565u,203280221u
  };
  u32 i,j,c,x,l0,l1;
  r64 dt;

  dt = cpu_time();
  x = 1u;  // current value of x
  c = 1u;  // value of pi(x+1)
  l0 = 0u; // base-2 logarithm of the next threshold minus 1
  l1 = 1u; // next threshold; equal to 2^(l0+1)-1
  fprintf(stderr,"test_aux_sieve: ");
  for(aux_base = 0u;l0 < 32u;aux_base += _aux_sieve_span_)
  {
    update_aux_sieve();
    for(i = 0u;i < _aux_sieve_words_;)
    {
      j = (l1 - x) / (16u * _pointer_size_); // number of words that can be processed quickly
      if(j > 0u)
      {
        if(j > _aux_sieve_words_ - i)
          j = _aux_sieve_words_ - i;
        c += count_zero_bits((u08 *)&aux_sieve[i],j * _pointer_size_);
        i += j;
        x += j * (16u * _pointer_size_);
      }
      else
      {
        for(j = 0u;j < 8u * _pointer_size_;j++)
        {
          if((aux_sieve[i] & mark_mask[j]) == (uXX)0u)
            c++;
          if(x == l1)
          {
            if(c != pi[l0])
            {
              fprintf(stderr,"bad   pi(2^%u)=%u != %u\n",l0 + 1u,pi[l0],c);
              exit(1);
            }
            fprintf(stderr,"%2u\b\b",++l0);
            l1 = 2u * l1 + 1u;
          }
          x += 2u;
        }
        i++;
      }
    }
  }
  dt = cpu_time() - dt;
  fprintf(stderr,"good (%.3f seconds)\n",dt);
}


//
// main segmented Eratosthenes sieve used to generate odd primes up to 2^64
//
// the aux_sieve variables are used by the update_sieve function and must not be modified
//

static uXX sieve[_sieve_words_]; // the sieve data (one bit per odd integer)
static u64 sieve_base;           // must be even and not larger than 2^64-_sieve_span_
static u64 expected_sieve_base;  // expected value of sieve_base
static u64 sieve_limit;          // maximum value of the prime numbers to be generated
static po *small_list;           // first free po entry in the small primes list
static po **large_lists;         // first free po entry in the head of each list of buckets
static u32 next_prime;           // first "prime" not in the buckets
static u32 list_size_log2;       // base-2 logarithm of the number of bucket lists
static u32 list_mask;            // equal to (1u<<list_size_log2)-1u
static u32 current_list;         // current list number

static void init_sieve(u64 limit)
{
  bucket *b;
  u32 i,j;

  free_all_memory();
  if((limit & 1ull) == 0ull)
    limit = limit - 1ull; // make it odd, to avoid arithmetic overflow latter on
  if(limit < (_sieve_span_ - 1ull))
  {
    fprintf(stderr,"init_sieve: bad sieve limit\n");
    exit(1);
  }
  if(_prime_threshold_ > 0xFFFFFFFFu - (1u << _sieve_bits_log2_))
    if(limit >= (u64)_prime_threshold_ * (u64)_prime_threshold_)
    { // almost certain arithmetic overflow when sieve_base is close to sieve_limit
      fprintf(stderr,"init_sieve: limit incompatible with _prime_threshold_\n");
      exit(1);
    }
  sieve_limit = limit;
  expected_sieve_base = 1ull; // illegal value
  available_buckets = NULL;
  // prepare the large lists
  i = 0u;
  for(j = 1u << 31u;j != 0u;j >>= 1u)
    if((u64)(i + j) * (u64)(i + j) <= sieve_limit)
      i += j; // at the end, i will be the integer part of the square root of sieve_limit
  if(i <= _prime_threshold_)
  { // all primes not larger than _prime_threshold_ will be placed on small_list
    large_lists = NULL;
  }
  else
  {
    i = 2u + (i >> _sieve_bits_log2_); // minimum number of lists
    for(list_size_log2 = 2u;(1u << list_size_log2) < i;list_size_log2++)
      ;
    list_mask = (1u << list_size_log2) - 1u;
    large_lists = (po **)get_memory((1u << list_size_log2) * (u32)sizeof(po *),64u);
    for(i = 0u;i < (1u << list_size_log2);i++)
    {
      if(available_buckets == NULL)
        more_buckets();
      b = available_buckets;
      available_buckets = available_buckets->next;
      large_lists[i] = &b->data[0];
      b->next = NULL;
    }
    current_list = 0u;
  }
  // prepare the small list
  if(available_buckets == NULL)
    more_buckets();
  b = available_buckets;
  available_buckets = available_buckets->next;
  small_list = &b->data[0];
  b->next = NULL;
  // prepare the auxiliary sieve
  aux_base = 0u;
  update_aux_sieve();
  next_prime = 17u;
}

#define find_next_prime()                               \
  do                                                    \
  {                                                     \
    u32 ii;                                             \
                                                        \
    next_prime += 2u;                                   \
    if(next_prime == 0xFFFFFFFFu)                       \
      break;                                            \
    if(next_prime > aux_base + (_aux_sieve_span_ - 1u)) \
    {                                                   \
      aux_base += _aux_sieve_span_;                     \
      update_aux_sieve();                               \
    }                                                   \
    ii = (next_prime - aux_base) >> 1u;                 \
    if(test_2(aux_sieve,ii))                            \
      break;                                            \
  }                                                     \
  while(1)

#define insert_po(pp,oo,ooo)                                               \
  do                                                                       \
  {                                                                        \
    bucket *bb1,*bb2;                                                      \
    u32 ll;                                                                \
    po *zz;                                                                \
                                                                           \
    ll = (current_list + ((oo) >> _sieve_bits_log2_) + (ooo)) & list_mask; \
    zz = large_lists[ll];                                                  \
    if(((uXX)zz & _bucket_addr_mask_) == (uXX)0u)                          \
    {                                                                      \
      if(available_buckets == NULL)                                        \
        more_buckets();                                                    \
      bb1 = (bucket *)zz; /* end of the full bucket */                     \
      bb2 = available_buckets; /* start of the new empty bucket */         \
      available_buckets = available_buckets->next;                         \
      bb2->next = &bb1[-1];                                                \
      zz = &bb2->data[0];                                                  \
    }                                                                      \
    zz->p = (pp);                                                          \
    zz->o = (oo) & ((1u << _sieve_bits_log2_) - 1u);                       \
    large_lists[ll] = &zz[1];                                              \
  }                                                                        \
  while(0)

static void update_sieve(void)
{
  bucket *b1,*b2;
  u32 i,j,k,p,o;
  u64 e,t;
  po *z;

  if(expected_sieve_base == 1ull)
    expected_sieve_base = sieve_base & ~1ull;
  if(sieve_base != expected_sieve_base || sieve_base > sieve_limit - (_sieve_span_ - 1ull))
  {
    fprintf(stderr,"update_sieve: illegal sieve_base\n");
    exit(1);
  }
  e = sieve_base + (_sieve_span_ - 1ull);
  // if needed, add primes to small_list or to large_lists
  while(next_prime < 0xFFFFFFFFu)
  {
    t = (u64)next_prime * (u64)next_prime;
    if(t > e)
      break;
    if(t > sieve_base)
      o = (u32)(t - sieve_base) >> 1u;
    else
    {
      o = next_prime - (u32)(sieve_base % (u64)next_prime);
      if((o & 1u) == 0u)
        o = (o >> 1u) + (next_prime >> 1u); // arithmetic overflow not possible
      else
        o = (o >> 1u);
    }
    if(next_prime <= _prime_threshold_)
    { // insert it in small_list
      if(((uXX)small_list & _bucket_addr_mask_) == (uXX)0u)
      { // the current bucket is full
        if(available_buckets == NULL)
          more_buckets();
        b1 = (bucket *)small_list; // end of the full bucket
        b2 = available_buckets;    // start of the new empty bucket
        available_buckets = available_buckets->next;
        b2->next = &b1[-1];
        small_list = &b2->data[0];
      }
      small_list->p = next_prime;
      small_list->o = o;
      small_list++;
    }
    else
    { // insert it in large_lists
      insert_po(next_prime,o,0);
    }
    find_next_prime();
  }
  // initialize the sieve array
  o = (u32)(sieve_base % (u64)(3u * 5u * 7u * 11u * 13u));
# if _pointer_size_ == 4u
  o = (o + ((o * 41u) & 63u) * 3u * 5u * 7u * 11u * 13u) >> 6u;   // 41 = -1/15015 mod 64
# elif _pointer_size_ == 8u
  o = (o + ((o * 105u) & 127u) * 3u * 5u * 7u * 11u * 13u) >> 7u; // 105 = -1/15015 mod 128
# endif
  for(i = 0u;i < _sieve_words_;i += k)
  {
    j = _sieve_words_ - i;            // remaining sieve words
    k = 3u * 5u * 7u * 11u * 13u - o; // remaining pattern words
    if(j < k)
      k = j;
    for(j = 0u;j < k;j++)
      sieve[i + j] = pattern[o + j];
    o = 0u;
  }
  if(sieve_base == 0ull)
    sieve[0] |= mark_mask[0]; // mark 1 as not prime
  if(sieve_base < 4ull)
    unmark(sieve,(3u - (u32)sieve_base) >> 1u); // mark 3 as prime
  if(sieve_base < 6ull)
    unmark(sieve,(5u - (u32)sieve_base) >> 1u); // mark 5 as prime
  if(sieve_base < 8ull)
    unmark(sieve,(7u - (u32)sieve_base) >> 1u); // mark 7 as prime
  if(sieve_base < 12ull)
    unmark(sieve,(11u - (u32)sieve_base) >> 1u); // mark 11 as prime
  if(sieve_base < 14ull)
    unmark(sieve,(13u - (u32)sieve_base) >> 1u); // mark 13 as prime
  // sieve with the small primes
  for(z = &small_list[-1];((uXX)z & _bucket_addr_mask_) != (uXX)0u;z--)
  { // deal with the (possibly partially full) last bucket
    p = z->p;
    for(o = z->o;o < (1u << _sieve_bits_log2_);o += p)
      mark_2(sieve,o);
    z->o = o - (1u << _sieve_bits_log2_);
  }
  b1 = (bucket *)z;
  while((b1 = b1->next) != NULL)
    for(i = 0;i < _primes_per_bucket_;i++)
    { // deal with the rest (full buckets)
      p = b1->data[i].p;
      for(o = b1->data[i].o;o < (1u << _sieve_bits_log2_);o += p)
        mark_2(sieve,o);
      b1->data[i].o = o - (1u << _sieve_bits_log2_);
    }
  // sieve with the large primes
  if(large_lists != NULL)
  {
    for(z = &large_lists[current_list][-1];((uXX)z & _bucket_addr_mask_) != (uXX)0u;z--)
    { // deal with the (possibly partially full) last bucket
      p = z->p;
      o = z->o;
      mark_2(sieve,o);
      o += p;
      if(o < (1u << _sieve_bits_log2_))
        insert_po(p,o,(1u << (32u - _sieve_bits_log2_))); // overflow
      else
        insert_po(p,o,0u); // no overflow
    }
    b1 = (bucket *)z;
    b2 = b1->next;
    b1->next = available_buckets;
    available_buckets = b1;
    while(b2 != NULL)
    {
      for(i = 0;i < _primes_per_bucket_;i++)
      { // deal with the rest (full buckets)
        p = b2->data[i].p;
        o = b2->data[i].o;
        mark_2(sieve,o);
        o += p;
        if(o < (1u << _sieve_bits_log2_))
          insert_po(p,o,(1u << (32u - _sieve_bits_log2_))); // overflow in the new o value
        else
          insert_po(p,o,0u); // no overflow in the new o value
      }
      b1 = b2;
      b2 = b2->next;
      b1->next = available_buckets;
      available_buckets = b1;
    }
    b1 = available_buckets;
    available_buckets = available_buckets->next;
    large_lists[current_list] = &b1->data[0];
    b1->next = NULL;
  }
  //
  // final updates
  //
  current_list = (current_list + 1u) & list_mask;
  expected_sieve_base = sieve_base + _sieve_span_;
}

#undef find_next_prime
#undef insert_po

static void test_sieve(int test_mask)
{
  static const u32 pi[36u] =
  { // pi(2^k)
            1ull,        2ull,         4ull,         6ull,
           11ull,       18ull,        31ull,        54ull,
           97ull,      172ull,       309ull,       564ull,
         1028ull,     1900ull,      3512ull,      6542ull,
        12251ull,    23000ull,     43390ull,     82025ull,
       155611ull,   295947ull,    564163ull,   1077871ull,
      2063689ull,  3957809ull,   7603553ull,  14630843ull,
     28192750ull, 54400028ull, 105097565ull, 203280221ull,
    393615806ull,762939111ull,1480206279ull,2874398515ull
  };
  static const u32 d_pi[5u] =
  {
    271431096u,271441873u,271428657u,271436608u,271437821u
  };
  u32 i,j,c,l0;
  u64 x,l1;
  FILE *fp;
  r64 dt;

  if((test_mask & 0x01u) != 0u)
  { // first test
    init_sieve((1ull << 36u) + 2ull * _sieve_span_);
    dt = cpu_time();
    x = 1ull;  // current value of x
    c = 1u;    // value of pi(x+1)
    l0 = 0u;   // base-2 logarithm of the next threshold minus 1
    l1 = 1ull; // next threshold; equal to 2^(l0+1)-1
    fprintf(stderr,"test_sieve (1):  0\b\b");
    for(sieve_base = 0ull;l0 < 36u;sieve_base = expected_sieve_base)
    {
      update_sieve();
      for(i = 0u;i < _sieve_words_;)
      {
        j = (u32)((l1 - x) / (u64)(16u * _pointer_size_));
        if(j > 0u)
        {
          if(j > _sieve_words_ - i)
            j = _sieve_words_ - i;
          c += count_zero_bits((u08 *)&sieve[i],j * _pointer_size_);
          i += j;
          x += (u64)(j * (16u * _pointer_size_));
        }
        else
        {
          for(j = 0u;j < 8u * _pointer_size_;j++)
          {
            if((sieve[i] & mark_mask[j]) == (uXX)0u)
              c++;
            if(x == l1)
            {
              if(c != pi[l0])
              {
                fprintf(stderr,"bad   pi(2^%u)=%u != %u\n",l0 + 1u,pi[l0],c);
                exit(1);
              }
              fprintf(stderr,"%2u\b\b",++l0);
              l1 = 2ull * l1 + 1ull;
            }
            x += 2ull;
          }
          i++;
        }
      }
    }
    dt = cpu_time() - dt;
    fprintf(stderr,"good (%.3f seconds)\n",dt);
  }
  if((test_mask & 0x02u) != 0u)
  { // second test
    init_sieve(1000000ull * 10000000000ull + 2ull * _sieve_span_);
    dt = cpu_time();
    l0 = 999995u;
    x = (u64)l0 * 10000000000ull + 1ull;
    c = 0u;
    l1 = (u64)(l0 + 1u) * 10000000000ull + 1ull;
    fprintf(stderr,"test_sieve (2): %7u\b\b\b\b\b\b\b",l0);
    for(sieve_base = (u64)l0 * 10000000000ull;l0 < 1000000u;sieve_base = expected_sieve_base)
    {
      update_sieve();
      for(i = 0u;i < _sieve_words_;)
      {
        j = (u32)((l1 - x) / (u64)(16u * _pointer_size_));
        if(j > 0u)
        {
          if(j > _sieve_words_ - i)
            j = _sieve_words_ - i;
          c += count_zero_bits((u08 *)&sieve[i],j * _pointer_size_);
          i += j;
          x += (u64)(j * (16u * _pointer_size_));
        }
        if(x == l1)
        {
          if(c != d_pi[l0 - 999995u])
          {
            fprintf(stderr,"bad   d_pi(%u)=%u != %u\n",l0,d_pi[l0 - 999995u],c);
            exit(1);
          }
          fprintf(stderr,"%7u\b\b\b\b\b\b\b",++l0);
          c = 0u;
          l1 = (u64)(l0 + 1u) * 10000000000ull + 1ull;
        }
      }
    }
    dt = cpu_time() - dt;
    fprintf(stderr,"good (%.3f seconds)\n",dt);
  }
  if((test_mask & 0x04u) != 0u)
  { // third test
    fprintf(stderr,"test_sieve (3)\n");
    init_sieve(0xFFFFFFFFFFFFFFFFull);
    sieve_base = sieve_limit - (1u << 30u) - _sieve_span_ + 1ull;
    dt = cpu_time();
    update_sieve(); // the first time is slow
    dt = cpu_time() - dt;
    fprintf(stderr,"  init .... %.3f seconds\n",dt);
    dt = cpu_time();
    for(i = 0u;i < (1u << (29u - _sieve_bits_log2_));i++)
    {
      sieve_base = expected_sieve_base;
      update_sieve();
    }
    dt = cpu_time() - dt;
    fprintf(stderr,"  sieve ... %.3f seconds\n",dt);
    if(expected_sieve_base != 0ull)
    {
      fprintf(stderr,"test_sieve (3): unexpected expected_sieve_base\n");
      exit(1);
    }
    fp = fopen("data.gp","wb");
    if(fp == NULL)
    {
      fprintf(stderr,"test_sieve (3): unable to create pari-gp script\n");
      exit(1);
    }
    fprintf(fp,"p=nextprime(%llu);v=1;\n",sieve_base);
    fprintf(fp,"f(q)=if(v==1,if(p==q,p=nextprime(p+1);,v=0;print(p,\" != \",q);););\n");
    for(i = 0u;i < (1u << _sieve_bits_log2_);i++)
      if(test_2(sieve,i))
        fprintf(fp,"f(%llu);\n",sieve_base + (u64)(2u * i + 1u));
    fprintf(fp,"if(p<2^64,v=0;);if(v==1,v=\"good\";,v=\"bad\";);\n");
    fprintf(fp,"print(\"test_sieve (3): \",v);\n");
    fprintf(fp,"quit();\n");
    fclose(fp);
    fprintf(stderr,"  use \"gp -q data.gp\" to finish this test\n");
  }
}


//
// measurements
//

static void measure_sieve_speed(u32 m,u32 n)
{
  r64 dt1,dt2;
  u64 c;
  u32 i;

  c = convert(m,n);
  init_sieve(c + (1ull << 30u));
  sieve_base = c - (1ull << 29u) - _sieve_span_;
  dt1 = cpu_time();
  update_sieve();
  dt1 = cpu_time() - dt1;
  dt2 = cpu_time();
  for(i = 0u;i < (1u << (29u - _sieve_bits_log2_));i++)
  {
    sieve_base = expected_sieve_base;
    update_sieve();
  }
  dt2 = cpu_time() - dt2;
  printf("%u %u %.3f %.3f\n",_sieve_bits_log2_,_bucket_size_log2_,dt1,dt2);
}


//
// find the best _sieve_bits_log2_ parameter
//

static r64 try(u32 sb,u32 bs,u32 m,u32 n,char *cflags,r64 t_best)
{
  char cmd[1024];
  r64 t1,t2;
  u32 a1,a2;
  FILE *fp;

  sprintf(cmd,"cc -D_sieve_bits_log2_=%u -D_bucket_size_log2_=%u %s fast_sieve.c",sb,bs,cflags);
  if(system(cmd) != 0)
  {
    fprintf(stderr,"compile: \"%s\" failed\n",cmd);
    exit(0);
  }
  sprintf(cmd,"./a.out -m %ue%u",m,n);
  fp = popen(cmd,"r");
  if(fp == NULL)
  {
    remove("a.out");
    fprintf(stderr,"compile: \"%s\" failed\n",cmd);
    exit(0);
  }
  if(fscanf(fp,"%u %u %lf %lf",&a1,&a2,&t1,&t2) != 4 || a1 != sb || a2 != bs)
  {
    pclose(fp);
    remove("a.out");
    fprintf(stderr,"compile: \"%s\" failed\n",cmd);
    exit(0);
  }
  pclose(fp);
  remove("a.out");
  printf("  %ue%u %2u %2u %8.3f (%7.3f)%s\n",m,n,sb,bs,t2,t1,(t2 < t_best) ? " *" : "");
  fflush(stdout);
  return t2;
}

static void find_best_sieve_size_log2(u32 bs,char *cflags)
{
  u32 m,n,sb,sb_best,sb_prev_best;
  r64 t,t_best;

  printf("Compilation flags: %s\n",cflags);
  sb_prev_best = 19u;
  for(n = 10u;n <= 19u;n++)
    for(m = 1u;m < 10u;m = 2u * m + (m >> 1u))
      if(n < 19u || m == 1u)
      {
        printf("Finding the best _sieve_size_log2_ for %ue%u ...\n",m,n);
        fflush(stdout);
        sb_best = sb_prev_best;
        t_best = try(sb_best,bs,m,n,cflags,1.0e100);
        for(sb = sb_prev_best - 1u;sb >= 10u;sb--)
        { // try smaller values
          t = try(sb,bs,m,n,cflags,t_best);
          if(t < t_best)
          {
            sb_best = sb;
            t_best = t;
          }
          else if(t > 1.1 * t_best)
            break;
        }
        for(sb = sb_prev_best + 1u;sb <= 28u;sb++)
        { // try larger values
          t = try(sb,bs,m,n,cflags,t_best);
          if(t < t_best)
          {
            sb_best = sb;
            t_best = t;
          }
          else if(t > 1.1 * t_best)
            break;
        }
        printf("  %ue%u %2u %2u %8.3f best\n",m,n,sb_prev_best = sb_best,bs,t_best);
        fflush(stdout);
      }
}


//
// example application: find prime gaps larger than or equal to 1000
//

static void kilo_prime_gaps(u32 m,u32 n)
{
#define _g_  1000u
  u32 i,li,ci,b1,b2;
  uXX lw,cw,t;
  u64 c,p;

  c = convert(m,n);
  c &= ~1ull;
  if(c < 1000000000000000ull)
    c = 1000000000000000ull;
  init_sieve(c + 1000000000000000ull);
  li = _sieve_words_;
  lw = (uXX)1u;
  sieve_base = c - 2048ull;
  while(sieve_base < sieve_limit - _sieve_span_)
  {
    update_sieve();
    for(i = 0u;i < _sieve_words_;i++)
    {
      cw = ~sieve[i];
      if(cw != (uXX)0u)
      {
        ci = i + _sieve_words_;
        if(ci - li >= _g_ / (16u * _pointer_size_))
        {
          b1 = 8u * _pointer_size_ * li - 1u;
          for(t = lw;t != (uXX)0u;t >>= 1u)
            b1++;
          b2 = 8u * _pointer_size_ * ci;
          for(t = cw;(t & (uXX)1u) == (uXX)0u;t >>= 1u)
            b2++;
          if(b2 - b1 >= _g_ / 2u)
          {
            p = sieve_base + 2ull * (u64)b2 + 1ull - _sieve_span_;
            if(p > c)
            {
              printf("%4u %llu\n",2u * (b2 - b1),p);
              fflush(stdout);
            }
          }
        }
        lw = cw;
        li = ci;
      }
    }
    li -= _sieve_words_;
    sieve_base = expected_sieve_base;
  }
#undef _g_
}


//
// main program
//

int main(int argc,char **argv)
{
  u32 m,n;

  if((argc == 2 || argc == 3) && strcmp(argv[1],"-b") == 0)
  {
    find_best_sieve_size_log2(12u,(argc == 2) ? "-Wall -O2" : argv[2]);
    return 0;
  }
  if(argc == 3 && strcmp(argv[1],"-g") == 0 && sscanf(argv[2],"%ue%u",&m,&n) == 2)
  {
    kilo_prime_gaps(m,n);
    return 0;
  }
  if(argc == 3 && strcmp(argv[1],"-m") == 0 && sscanf(argv[2],"%ue%u",&m,&n) == 2)
  {
    measure_sieve_speed(m,n);
    return 0;
  }
  if(argc == 2 && strcmp(argv[1],"-ta") == 0)
  {
    test_aux_sieve();
    return 0;
  }
  if(argc == 3 && strcmp(argv[1],"-ts") == 0)
  {
    test_sieve(atoi(argv[2]));
    return 0;
  }
  if(argc == 2 && strcmp(argv[1],"-v") == 0)
  {
    printf("_pointer_size_ ........ %u\n",_pointer_size_);
    printf("_bucket_size_log2_ .... %u\n",_bucket_size_log2_);
    printf("_sieve_bits_log2_ ..... %u\n",_sieve_bits_log2_);
    printf("_prime_threshold_ ..... %u\n",_prime_threshold_);
    return 0;
  }
  fprintf(stderr,"usage: %s -b [cflags]  (find the best _sieve_size_log2_)\n",argv[0]);
  fprintf(stderr,"       %s -g MeN       (example: kilo prime gaps after M*10^N)\n",argv[0]);
  fprintf(stderr,"       %s -m MeN       (measure the sieve speed near M*10^N)\n",argv[0]);
  fprintf(stderr,"       %s -ta          (test the aux_sieve functions)\n",argv[0]);
  fprintf(stderr,"       %s -ts mask     (test the sieve functions)\n",argv[0]);
  fprintf(stderr,"         0x01 ... 2^36 test\n");
  fprintf(stderr,"         0x02 ... 10^16 test\n");
  fprintf(stderr,"         0x04 ... 2^64 test (requires about 1.6GiB of memory)\n");
  fprintf(stderr,"       %s -v           (dump compilation-time parameters)\n",argv[0]);
  return 1;
}


//
// all time measurements for _bucket_size_log2_ equal to 12
//
// measurement interval: 2^30
//
// brand    processor name     freq   L1    L2                     main memory
// -----  ----------------  -------  ---  ----  ------------------------------
//   AMD          Athlon64  2.20GHz  64k  512k    2G 333MHz DDR (dual channel)
// Intel  Core2 Quad Q9550  2.83GHz  32k    6M  8G 1100MHz DDR2 (dual channel)
//
//
// time measurements for a 2.20GHz Athlon 64
//
//  MeN     16     17     18     19     20     21     22
// ---- ------ ------ ------ ------ ------ ------ ------
// 1e10  1.968  1.766  1.656  1.615  2.282
// 2e10         1.889  1.727  1.667  2.347
// 5e10         2.062  1.858  1.752  2.428
// 1e11         2.189  1.971  1.825  2.517
// 2e11         2.320  2.099  1.928  2.600
// 5e11         2.508  2.283  2.120  2.763
// 1e12         2.665  2.442  2.274  2.887
// 2e12         2.809  2.587  2.422  3.051
// 5e12         2.994  2.770  2.600  3.233
// 1e13         3.134  2.906  2.735  3.377
// 2e13         3.278  3.049  2.867  3.510
// 5e13         3.484  3.224  3.038  3.684
// 1e14         3.650  3.368  3.164  3.817
// 2e14         3.846  3.512  3.288  3.946
// 5e14         4.130  3.711  3.461  4.115
// 1e15         4.381  3.865  3.587  4.237
// 2e15         4.660  4.031  3.725  4.364
// 5e15         5.058  4.276  3.902  4.529
// 1e16                4.484  4.044  4.656
// 2e16                4.726  4.198  4.778
// 5e16                5.075  4.414  4.950
// 1e17                5.369  4.606  5.088
// 2e17                5.701  4.817  5.232  5.579
// 5e17                6.209  5.135  5.434  5.768
// 1e18                6.662  5.399  5.604  5.910  8.083
// 2e18                7.193  5.697  5.794  6.061  8.359
// 5e18                8.027  6.148  6.095  6.272  8.699
// 1e19                8.782  6.554  6.339  6.453  8.970
// ---- ------ ------ ------ ------ ------ ------ ------
//  MeN     16     17     18     19     20     21     22
//
// time measurements for a 2.20GHz Athlon 64 with _prime_threshold_ = 0xFFFFFFFF
//
//  MeN     16     17     18     19     20     21     22     23     24     25     26     27     28
// ---- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
// 1e10  2.092  1.931  1.831  1.800  2.445
// 2e10         2.042  1.904  1.854  2.513
// 5e10         2.260  2.033  1.939  2.597
// 1e11         2.505  2.127  2.011  2.690
// 2e11         2.820  2.295  2.114  2.768
// 5e11                2.638  2.285  2.929
// 1e12                3.077  2.518  3.064
// 2e12                3.545  2.783  3.192
// 5e12                4.409  3.271  3.475  3.626
// 1e13                5.298  3.754  3.768  3.789  4.178
// 2e13                6.529  4.409  4.137  4.035  4.418 19.834
// 5e13                       5.626  4.813  4.500  4.875 20.680
// 1e14                              5.542  4.991  5.352 21.242
// 2e14                              6.568  5.658  5.983 21.631
// 5e14                              8.555  6.918  7.005 22.811
// 1e15                             10.739  8.266  7.995 23.595
// 2e15                             13.732 10.029  9.217 24.573
// 5e15                                    13.317 11.290 26.144
// 1e16                                    16.668 13.349 27.395
// 2e16                                    21.203 16.028 29.165
// 5e16                                    29.792 20.802 31.914
// 1e17                                    38.998 25.944 34.803
// 2e17                                    51.065 32.758 38.575
// 5e17                                    73.346 45.172 45.488 45.994 46.965 47.387 47.422 47.401
// 1e18                                    96.879 58.065 52.767 49.975 49.318 49.020 48.724 48.510
// 2e18                                                         55.315 52.286 50.911 50.124 49.685
// 5e18                                                                57.701 54.027 52.301 51.411
// 1e19                                                                63.358 57.191 54.268 52.882
// ---- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
//  MeN     16     17     18     19     20     21     22     23     24     25     26     27     28
//
// time measurements for a 2.83GHz Q9550
//
//  MeN     15     16     17     18     19     20     21     22     23     24
// ---- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
// 1e10         1.417  1.303  1.251  1.396  1.429
// 2e10         1.489  1.370  1.305  1.437
// 5e10         1.582  1.463  1.403  1.506  1.525  1.531  1.526  1.551
// 1e11         1.650  1.532  1.467  1.560  1.571  1.572  1.566  1.591  1.846
// 2e11         1.716  1.597  1.534  1.638  1.622  1.612  1.605  1.625  1.895
// 5e11         1.801  1.682  1.620  1.722  1.700  1.671  1.657  1.680  1.961
// 1e12         1.864  1.745  1.684  1.787  1.770  1.723  1.698  1.721  2.008
// 2e12         1.926  1.805  1.748  1.852  1.833  1.781  1.740  1.758  2.059
// 5e12         2.005  1.884  1.828  1.933  1.916  1.867  1.804  1.812  2.125
// 1e13         2.066  1.942  1.888  1.995  1.977  1.928  1.869  1.857  2.176
// 2e13         2.124  2.000  1.946  2.054  2.038  1.990  1.928  1.908  2.227
// 5e13         2.203  2.073  2.020  2.131  2.116  2.067  2.007  1.997  2.306
// 1e14         2.264  2.131  2.078  2.189  2.176  2.129  2.070  2.053  2.368
// 2e14  2.520  2.331  2.192  2.139  2.253  2.235  2.195  2.135  2.122  2.453
// 5e14         2.423  2.274  2.217  2.331  2.320  2.276  2.219  2.202  2.536
// 1e15  2.725  2.488  2.324  2.274  2.390  2.378  2.336  2.278  2.263  2.600
// 2e15         2.557  2.395  2.331  2.446  2.436  2.396  2.336  2.320  2.662
// 5e15         2.654  2.474  2.402  2.512  2.503  2.473  2.413  2.376  2.742
// 1e16         2.732  2.538  2.467  2.576  2.569  2.529  2.467  2.450  2.807
// 2e16         2.821  2.602  2.516  2.637  2.627  2.586  2.521  2.502  2.870
// 5e16         2.968  2.698  2.606  2.714  2.701  2.659  2.593  2.578  2.953
// 1e17         3.094  2.776  2.668  2.775  2.763  2.715  2.647  2.611  3.017
// 2e17         3.240  2.869  2.733  2.836  2.822  2.771  2.701  2.682  3.081
// 5e17         3.467  3.015  2.829  2.922  2.907  2.849  2.771  2.755  3.164
// 1e18                3.160  2.919  2.994  2.973  2.909  2.823  2.807  3.231
// 2e18                3.346  3.035  3.078  3.047  2.969  2.877  2.860  3.292
// 5e18                              3.233  3.158  3.052  2.949  2.929  3.376
// 1e19                              3.393  3.264  3.119  3.005  2.975  3.437
// ---- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
//  MeN     15     16     17     18     19     20     21     22     23     24
//
// time measurements for a 2.83GHz Q9550 with _prime_threshold_ = 0xFFFFFFFF
//
//  MeN     16     17     18     19     20     21     22     23     24     25     26
// ---- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
// 1e10  1.435  1.305  1.256  1.400  1.432
// 2e10  1.591  1.368  1.310  1.443
// 5e10  1.856  1.505  1.409  1.514  1.529  1.532  1.526  1.550
// 1e11         1.667  1.471  1.568  1.577  1.574  1.565  1.581  1.848
// 2e11         1.865  1.586  1.649  1.629  1.615  1.606  1.627  1.898
// 5e11         2.163  1.808  1.741  1.708  1.675  1.657  1.679  1.964
// 1e12                       1.873  1.780  1.729  1.698  1.722  2.011
// 2e12                       2.047  1.851  1.787  1.742  1.761  2.063
// 5e12                              2.032  1.874  1.807  1.810  2.132
// 1e13                              2.211  1.960  1.875  1.857  2.182
// 2e13                                     2.100  1.933  1.892  2.233
// 5e13                                     2.334  2.055  2.002  2.314
// 1e14                                     2.544  2.212  2.064  2.375
// 2e14                                            2.401  2.177  2.463
// 5e14                                            2.665  2.385  2.560  2.748
// 1e15                                            2.901  2.574  2.705  2.893
// 2e15                                            3.204  2.784  2.882  2.978  9.169
// 5e15                                            3.728  3.114  3.183  3.254  9.465
// 1e16                                            4.297  3.399  3.461  3.538  9.725
// 2e16                                            5.093  3.890  3.807  3.886 10.159
// 5e16                                            6.656  4.753  4.404  4.447 10.806
// 1e17                                                   5.693  5.024  4.968 11.306
// 2e17                                                   6.985  5.823  5.590 11.939
// 5e17                                                   9.438  7.257  6.606 12.634
// 1e18                                                          8.713  7.578 13.308
// 2e18                                                         10.639  8.769 14.200
// 5e18                                                         14.199 10.853 15.650
// 1e19                                                         18.130 13.000 17.018
// ---- ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
//  MeN     16     17     18     19     20     21     22     23     24     25     26
//
// Single measurement at 1e18 with SB=23 on a Q5550 ............. 2.809
// Four concurrent measurements at 1e18 with SB=23 on a Q9550 ... 3.821 3.828 4.095 4.096
//
//
// time measurements summary (measurement interval: 2^30); measurements with _prime_threshold_ set
// to 0xFFFFFFFFu are measurements of the classical segmented Eratosthenes sieve
//
//
//        2.20GHz Athlon64   2.20GHz Athlon64     2.83GHz Q9550      2.83GHz Q9550
//                            with threshold                        with threshold
// ----   -- ------ ------   -- ------ ------   -- ------ ------   -- ------ ------
//  MeN   SB  sieve   init   SB  sieve   init   SB  sieve   init   SB  sieve   init
// ----   -- ------ ------   -- ------ ------   -- ------ ------   -- ------ ------
// 1e10   19  1.615  0.002   19  1.800  0.003   18  1.251  0.001   18  1.256  0.001
// 2e10   19  1.667  0.003   19  1.854  0.003   18  1.305  0.001   18  1.310  0.002
// 5e10   19  1.752  0.003   19  1.939  0.004   18  1.403  0.002   18  1.409  0.001
// 1e11   19  1.825  0.005   19  2.011  0.004   18  1.467  0.002   18  1.471  0.002
// 2e11   19  1.928  0.006   19  2.114  0.005   18  1.534  0.002   18  1.586  0.002
// 5e11   19  2.120  0.008   19  2.285  0.009   18  1.620  0.004   22  1.657  0.017
// 1e12   19  2.274  0.011   19  2.518  0.011   18  1.684  0.005   22  1.698  0.018
// 2e12   19  2.422  0.014   19  2.783  0.015   22  1.740  0.019   22  1.742  0.020
// 5e12   19  2.600  0.021   19  3.271  0.022   22  1.804  0.023   22  1.807  0.024
// 1e13   19  2.735  0.030   19  3.754  0.030   23  1.857  0.043   23  1.857  0.044
// 2e13   19  2.867  0.040   21  4.035  0.053   23  1.908  0.048   23  1.892  0.049
// 5e13   19  3.038  0.062   21  4.500  0.076   23  1.997  0.060   23  2.002  0.062
// 1e14   19  3.164  0.086   21  4.991  0.101   23  2.053  0.073   23  2.064  0.075
// 2e14   19  3.288  0.121   21  5.658  0.135   23  2.122  0.090   23  2.177  0.094
// 5e14   19  3.461  0.190   21  6.918  0.205   23  2.202  0.124   23  2.385  0.132
// 1e15   19  3.587  0.267   22  7.995  0.314   23  2.263  0.163   23  2.574  0.173
// 2e15   19  3.725  0.392   22  9.217  0.442   23  2.320  0.222   23  2.784  0.237
// 5e15   19  3.902  0.625   22 11.290  0.659   23  2.376  0.327   23  3.114  0.350
// 1e16   19  4.044  0.867   22 13.349  0.919   23  2.450  0.451   23  3.399  0.482
// 2e16   19  4.198  1.225   22 16.028  1.276   23  2.502  0.622   24  3.807  0.726
// 5e16   19  4.414  1.933   22 20.802  1.979   23  2.578  0.961   24  4.404  1.096
// 1e17   19  4.606  2.747   22 25.944  2.783   23  2.611  1.349   25  4.968  1.671
// 2e17   19  4.817  3.884   22 32.758  3.884   23  2.682  1.889   25  5.590  2.267
// 5e17   19  5.135  6.176   22 45.172  6.071   23  2.755  2.960   25  6.606  3.437
// 1e18   19  5.399  8.778   28 48.510 32.372   23  2.807  4.170   25  7.578  4.748
// 2e18   19  5.697 12.490   28 49.685 36.265   23  2.860  5.877   25  8.769  6.589
// 5e18   20  6.095 19.577   28 51.411 43.765   23  2.929  9.290   25 10.853 10.241
// 1e19   20  6.339 27.903   28 52.882 51.957   23  2.975 13.138   25 13.000 14.347
// ----   -- ------ ------   -- ------ ------   -- ------ ------   -- ------ ------
//  MeN   SB  sieve   init   SB  sieve   init   SB  sieve   init   SB  sieve   init
//

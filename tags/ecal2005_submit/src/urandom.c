/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <errno.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "transsys.h"


static int ulong_random_seed;
static unsigned long ulong_rnd_state[64];



/*** adapted from the random() sources by Jan T. Kim ***/
/*
 * modifications:
 *    - setstate / initstate interface modified to operate with
 *      unsigned long * rather than with char *, as the original
 *      char * variant causes 32 / 64 bit incompatibility.
 *    - ulong_random always uses type 4 R.N.G. (state array of 64
 *      unsigned longs), shorter R.N.G.s are not supported.
 *      Therefore, the length argument has been removed from
 *      the initstate interface.
 */

/*
 * future plan: use struct to hold both state and seed to allow
 *     transparent save & restore operations in which seed is
 *     documented.
 */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *	@(#)random.c	5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 */


/* An improved random number generation package.  In addition to the standard
   rand()/srand() like interface, this package also has a special state info
   interface.  The initstate() routine is called with a seed, an array of
   bytes, and a count of how many bytes are being passed in; this array is
   then initialized to contain information for random number generation with
   that much state information.  Good sizes for the amount of state
   information are 32, 64, 128, and 256 bytes.  The state can be switched by
   calling the setstate() function with the same array as was initiallized
   with initstate().  By default, the package runs with 128 bytes of state
   information and generates far better random numbers than a linear
   congruential generator.  If the amount of state information is less than
   32 bytes, a simple linear congruential R.N.G. is used.  Internally, the
   state information is treated as an array of longs; the zeroeth element of
   the array is the type of R.N.G. being used (small integer); the remainder
   of the array is the state information for the R.N.G.  Thus, 32 bytes of
   state information will give 7 longs worth of state information, which will
   allow a degree seven polynomial.  (Note: The zeroeth word of state
   information also has some other information stored in it; see setstate
   for details).  The random number generation technique is a linear feedback
   shift register approach, employing trinomials (since there are fewer terms
   to sum up that way).  In this approach, the least significant bit of all
   the numbers in the state table will act as a linear feedback shift register,
   and will have period 2^deg - 1 (where deg is the degree of the polynomial
   being used, assuming that the polynomial is irreducible and primitive).
   The higher order bits will have longer periods, since their values are
   also influenced by pseudo-random carries out of the lower bits.  The
   total period of the generator is approximately deg*(2**deg - 1); thus
   doubling the amount of state information has a vast influence on the
   period of the generator.  Note: The deg*(2**deg - 1) is an approximation
   only good for large deg, when the period of the shift register is the
   dominant factor.  With deg equal to seven, the period is actually much
   longer than the 7*(2**7 - 1) predicted by this formula.  */



/* For each of the currently supported random number generators, we have a
   break value on the amount of state information (you need at least thi
   bytes of state info to support this random number generator), a degree for
   the polynomial (actually a trinomial) that the R.N.G. is based on, and
   separation between the two lower order coefficients of the trinomial.  */

/* x**63 + x + 1.  */
#define	TYPE_4		4
#define	BREAK_4		256
#define	DEG_4		63
#define	SEP_4		1
#define MAX_TYPES	5

long int ulong_random(void);


/* Array versions of the above information to make code run faster.
   Relies on fact that TYPE_i == i.  */


static unsigned long int randtbl[DEG_4 + 1] = {
	0x00000004, 0x07a6c921, 0x75247e98, 0x958685e2,
	0xeeb14548, 0xd436aed1, 0xeeadcc28, 0xd5a6d6b0, 
	0x882e20e0, 0xc1165f36, 0x2fed1390, 0x44ca0be4,
	0x86da3f50, 0x8d76cf46, 0xcd46b0b8, 0x2d717398, 
	0x422d0ad8, 0xde7747c7, 0x5097ce88, 0xded84046,
	0xdde7cd78, 0xe47c24af, 0x194a39e8, 0xf4d833a8, 
	0x5ffe8570, 0x02787a04, 0x070b6cc0, 0x09f6a498,
	0x9513abc0, 0x2ea921bc, 0xe21aa4f8, 0x30df9bf0, 
	0x7dbc80a8, 0x33786e9d, 0x121d5b78, 0x8893446a,
	0xad28f228, 0x10db401d, 0x80b8ef28, 0x7c226980, 
	0x4f7bb080, 0xcb0bb342, 0x2a73d7f0, 0x2db9414c,
	0x1e9c84b0, 0x94232c82, 0x55f835b8, 0x44543d68, 
	0x7dd87cf8, 0xd78cf6a3, 0xffdd0f68, 0xb5011ece,
	0xcc390358, 0xe52c329b, 0xb47b15e8, 0xb2ef08b8, 
	0x52dc7210, 0xeeab5470, 0xf7a53f20, 0x7d7ace80,
	0x1e483a20, 0x53d0e918, 0x1dc50cf8, 0xaf493880, 
};

/* Fvoid * and Rvoid * are two pointers into the state info, a front and a rear
   pointer.  These two pointers are always rand_sep places aparts, as they
   cycle through the state information.  (Yes, this does mean we could get
   away with just one pointer, but the code for random is more efficient
   this way).  The pointers are left positioned as they would be from the call:
	initstate(1, randtbl, 128);
   (The position of the rear pointer, rptr, is really 0 (as explained above
   in the initialization of randtbl) because the state table pointer is set
   to point to randtbl[1] (as explained below).)  */

static unsigned long int *fptr = &randtbl[SEP_4 + 1];
static unsigned long int *rptr = &randtbl[1];



/* The following things are the pointer to the state information table,
   the type of the current generator, the degree of the current polynomial
   being used, and the separation between the two pointers.
   Note that for efficiency of random, we remember the first location of
   the state information, not the zeroeth.  Hence it is valid to access
   state[-1], which is used to store the type of the R.N.G.
   Also, we remember the last location, since this is more efficient than
   indexing every time to find the address of the last element to see if
   the front and rear pointers have wrapped.  */

static unsigned long int *state = randtbl + 1;

const static int rand_type = TYPE_4;
const static int rand_deg = DEG_4;
const static int rand_sep = SEP_4;

static unsigned long int *end_ptr = &randtbl[sizeof(randtbl) / sizeof(randtbl[0])];

/* Initialize the random number generator based on the given seed.  If the
   type is the trivial no-state-information type, just remember the seed.
   Otherwise, initializes state[] based on the given "seed" via a linear
   congruential generator.  Then, the pointers are set to known locations
   that are exactly rand_sep places apart.  Lastly, it cycles the state
   information a given number of times to get rid of any initial dependencies
   introduced by the L.C.R.N.G.  Note that the initialization of randtbl[]
   for default usage relies on values produced by this routine.  */

void ulong_srandom( unsigned int x)
{
  register long int i;
  state[0] = x;
  /* fprintf(stderr, "srandom: state[%3d] = %12ld (%10lx)\n", 0, state[0], state[0]); */
  for (i = 1; i < rand_deg; ++i)
  {
#ifdef LINUX_RANDOM
    state[i] = 1103515145L * state[i - 1] + 12345;
#else
    state[i] = 1103515245L * state[i - 1] + 12345;
#endif
    /* fprintf(stderr, "srandom: state[%3ld] = %12lu (%10lx)\n", i, state[i], state[i]); */
  }
  fptr = &state[rand_sep];
  rptr = &state[0];
  for (i = 0; i < 10 * rand_deg; ++i)
    (void) ulong_random();
}



/* Initialize the state information in the given array of N bytes for
   future random number generation.  Based on the number of bytes we
   are given, and the break values for the different R.N.G.'s, we choose
   the best (largest) one we can and set things up for it.  srandom is
   then called to initialize the state information.  Note that on return
   from srandom, we set state[-1] to be the type multiplexed with the current
   value of the rear pointer; this is so successive calls to initstate won't
   lose this information and will be able to restart with setstate.
   Note: The first thing we do is save the current state, if any, just like
   setstate so that it doesn't matter when initstate is called.
   Returns a pointer to the old state.  */

unsigned long *ulong_initstate(unsigned int seed, unsigned long *arg_state)
{
  unsigned long *ostate = &state[-1];

  state[-1] = (MAX_TYPES * (rptr - state)) + rand_type;
  state = arg_state + 1;	/* First location.  */
  /* Must set END_void * before srandom.  */
  end_ptr = state + rand_deg;
  ulong_srandom(seed);
  state[-1] = (MAX_TYPES * (rptr - state)) + rand_type;
  return ostate;
}



/* Restore the state from the given state array.
   Note: It is important that we also remember the locations of the pointers
   in the current state information, and restore the locations of the pointers
   from the old state information.  This is done by multiplexing the pointer
   location into the zeroeth word of the state information. Note that due
   to the order in which things are done, it is OK to call setstate with the
   same state as the current state
   Returns a pointer to the old state information.  */

unsigned long *ulong_setstate(unsigned long *arg_state)
{
  register unsigned long int *new_state = arg_state;
  register int type = new_state[0] % MAX_TYPES;
  register int rear = new_state[0] / MAX_TYPES;
  void *ostate = (void *) &state[-1];

  state[-1] = (MAX_TYPES * (rptr - state)) + rand_type;
  switch (type)
  {
  case TYPE_4:
    break;
  default:
    /* State info munged.  */
    fprintf(stderr, "ulong_setstate(): state info (type %d) munged -- not changed\n", type);
    errno = EINVAL;
    return NULL;
  }
  state = new_state + 1;
  rptr = state + rear;
  fptr = state + (rear + rand_sep) % rand_deg;
  /* Set end_ptr too.  */
  end_ptr = state + rand_deg;
  return ostate;
}



/* If we are using the trivial TYPE_0 R.N.G., just do the old linear
   congruential bit.  Otherwise, we do our fancy trinomial stuff, which is the
   same in all ther other cases due to all the global variables that have been
   set up.  The basic operation is to add the number at the rear pointer into
   the one at the front pointer.  Then both pointers are advanced to the next
   location cyclically in the table.  The value returned is the sum generated,
   reduced to 31 bits by throwing away the "least random" low bit.
   Note: The code takes advantage of the fact that both the front and
   rear pointers can't wrap on the same call by not testing the rear
   pointer if the front one has wrapped.  Returns a 31-bit random number.  */

long int ulong_random(void)
{
  long int i;

  *fptr += *rptr;
  *fptr &= 0xffffffff;
  /* Chucking least random bit.  */
  i = (*fptr >> 1) & 0x7fffffffUL;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return i;
}

/*** end of code taken from random() sources ***/


unsigned long urandom_long(unsigned long range)
{
  unsigned long max;
  unsigned long r;
  static unsigned long num_urandom_long_calls = 0;

  if (range == 0)
    return (0);
  max = 0x80000000UL - 0x80000000UL % range;
  do
  {
    r = ulong_random();
/*
    if (max < r)
      fprintf(stderr, "urandom_long: repeat at call %lu (range %lu, max %lu, r %lu)\n", num_urandom_long_calls, range, max, r);
*/
  }
  while (max < r);
  num_urandom_long_calls++;
  return (r % range);
}


double urandom_double(void)
{
  return ((double) (ulong_random() & 0x7fffff) / 0x800000);
}


/*
 * adapted from "Numerical Recipes in C", p. 289
 * uses ulong_random() rather than ran1() from Numerical Recipes
 *
 * Hack/deficiency: random_gauss does not respond to setting state
 * appropriately. Ought to manage iset and gset in state struct
 * of some sort...  see "future plan".
 */

double urandom_gauss(void)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if (iset == 0)
  {
    do
    {
      v1 = 2.0 * ulong_random() / 0x7fffffff - 1.0;
      v2 = 2.0 * ulong_random() / 0x7fffffff - 1.0;
      rsq = v1 * v1 + v2 * v2;
    }
    while ((rsq >= 1.0) || (rsq == 0.0));
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    return (v2 * fac);
  }
  else
  {
    iset = 0;
    return (gset);
  }
}


int write_urandom_state(FILE *f)
{
  int i;
  unsigned long junk[64], *rndstate;

  /* fprintf(stderr, "write_urandom_state: starting\n"); */
  /*** Note: it is necessary to call ulong_initstate() before saving ulong_rnd_state ***/
  /*** This call modifies ulong_rnd_state, without this modification resuming        ***/
  /*** does not work!                                                                ***/
  rndstate = ulong_initstate(4711, junk);
  fprintf(f, "ulong_rndstate 64\n");
  for (i = 0; i < 64; i++)
    fprintf(f, "%lu\n", rndstate[i]);
  ulong_setstate(rndstate);
  return (0);
}


int read_urandom_state(FILE *f)
{
  int i;
  char buf[514];
  unsigned long junk[64];

  /* fprintf(stderr, "read_urandom_state: starting\n"); */
  ulong_initstate(4711, junk);
  if (fgets(buf, 514, f) != buf)
  {
    fprintf(stderr, "read_urandom_state: error reading state information\n");
    return (-1);
  }
  if (strcmp(buf, "ulong_rndstate 64\n"))
  {
    fprintf(stderr, "read_urandom_state: no \"ulong_rndstate 64\" label\n");
    fprintf(stderr, "\"%s\"", buf);
    return (-1);
  }
  ulong_random_seed = strtol(buf, (char **) NULL, 10);
  buf[0] = '\0';
  for (i = 0; i < 64; i++)
  {
    if (fgets(buf, 514, f) != buf)
    {
      fprintf(stderr, "read_urandom_state: error reading state information (state array)\n");
      return (-1);
    }
    ulong_rnd_state[i] = strtoul(buf, NULL, 10);
  }
  ulong_setstate(ulong_rnd_state);
  return (0);
}


/*
 * Randomly shuffles the num elements of array a
 */

void *urandom_shuffle(size_t num, size_t s, void *a)
{
  long i, r;
  char *buf;
  char *a1 = a;

  buf = malloc(s);
  if (buf == NULL)
  {
    fprintf(stderr, "random_shuffle: malloc failed\n");
    return (NULL);
  }
  for (i = 0; i < num; i++)
  {
    r = urandom_long(num);
    memcpy(buf, a1 + r * s, s);
    memcpy(a1 + r * s, a1 + i * s, s);
    memcpy(a1 + i * s, buf, s);
  }
  free(buf);
  return (a);
}


#ifdef TEST

int main(int argc, char *argv[])
{
  char rstate[256];
  long i, lr, r, l1[10], l2[10];
  unsigned long virginstate[64], junk[64];
  FILE *f;

  
  ulong_initstate(1, virginstate);
  ulong_initstate(1, junk);
  f = fopen("vstate.dat", "w");
  for (i = 0; i < 64; i += 8)
  {
    for (r = i; r < i + 8; r++)
      fprintf(f, "0x%08lx, ", virginstate[r]);
    fprintf(f, "\n");
  }
  fclose(f);
  initstate(1, rstate, 256);
  ulong_srandom(1);
  printf("ulong_random random\n");
  for (i = 0; i < 10000000; i++)
  {
    lr = ulong_random();
    r = random();
    if (lr != r)
      printf("%3ld: %12ld %12ld\n", i, lr, r);
  }
  printf("\n");
  ulong_srandom(12345);
  printf("urandom_long values after seeding with 12335\n");
  for (i = 0; i < 10; i++)
  {
    lr = urandom_long(100);
    printf("%3ld: %10ld\n", i, lr);
  }
  printf("\n");
  f = fopen("l.rnd", "w");
  write_urandom_state(f);
  fclose(f);
  printf("state written\n");
  for (i = 0; i < 10; i++)
  {
    l1[i] = urandom_long(100);
    printf("%3ld: %10ld\n", i, l1[i]);
  }
  printf("\n");
  for (i = 0; i < 10; i++)
  {
    lr = urandom_long(100);
    printf("%3ld: %10ld\n", i, lr);
  }
  printf("\n");
  f = fopen("l.rnd", "r");
  read_urandom_state(f);
  fclose(f);
  printf("state read\n");
  f = fopen("lx.rnd", "w");
  write_urandom_state(f);
  fclose(f);
  for (i = 0; i < 10; i++)
  {
    l2[i] = urandom_long(100);
    printf("%3ld: %10ld\n", i, l2[i]);
  }
  printf("\n");
  for (i = 0; i < 10; i++)
  {
    lr = urandom_long(100);
    printf("%3ld: %10ld\n", i, lr);
  }
  printf("\n");
  return (EXIT_SUCCESS);
}


#endif


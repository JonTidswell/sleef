//   Copyright Naoki Shibata and contributors 2010 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


// Horner's Scheme
//   * good for throughput and accuracy
//   * can be slow for high degree polynomials (high latency)
// Estrin's Method
//   * faster than Horner's Scheme (improved latency)
//   * sacrifices some of both throughput and accuracy

#define USE_ESTRIN_NOT_HORNER 0

#if !defined(USE_ESTRIN_NOT_HORNER)
#define USE_ESTRIN_NOT_HORNER 1
#endif

// simple default name 
#if (USE_ESTRIN_NOT_HORNER)
#define POLY2 EstrinPoly2
#define POLY3 EstrinPoly3
#define POLY4 EstrinPoly4
#define POLY5 EstrinPoly5
#define POLY6 EstrinPoly6
#define POLY7 EstrinPoly7
#define POLY8 EstrinPoly8
#define POLY9 EstrinPoly9
#define POLY10 EstrinPoly10
#define POLY11 EstrinPoly11
#define POLY12 EstrinPoly12
#define POLY13 EstrinPoly13
#define POLY14 EstrinPoly14
#define POLY15 EstrinPoly15
#define POLY16 EstrinPoly16
#define POLY17 EstrinPoly17
#define POLY18 EstrinPoly18
#define POLY19 EstrinPoly19
#define POLY20 EstrinPoly20
#define POLY21 EstrinPoly21
#define POLY22 EstrinPoly22
#define POLY23 EstrinPoly23
#define POLY24 EstrinPoly24
#else
#define POLY2 EXHornerPoly2
#define POLY3 EXHornerPoly3
#define POLY4 EXHornerPoly4
#define POLY5 EXHornerPoly5
#define POLY6 EXHornerPoly6
#define POLY7 EXHornerPoly7
#define POLY8 EXHornerPoly8
#define POLY9 EXHornerPoly9
#define POLY10 EXHornerPoly10
#define POLY11 EXHornerPoly11
#define POLY12 EXHornerPoly12
#define POLY13 EXHornerPoly13
#define POLY14 EXHornerPoly14
#define POLY15 EXHornerPoly15
#define POLY16 EXHornerPoly16
#define POLY17 EXHornerPoly17
#define POLY18 EXHornerPoly18
#define POLY19 EXHornerPoly19
#define POLY20 EXHornerPoly20
#define POLY21 EXHornerPoly21
#define POLY22 EXHornerPoly22
#define POLY23 EXHornerPoly23
#define POLY24 EXHornerPoly24
#endif


// Assumptions:
//       MLA(x,y,z) == x * y + z
//       C2V(c) casts c from scalar to vector (if appropriate) 


// These are macros for evaluating polynomials using Estrin's method

#define EstrinPoly2(x, c1, c0) MLA(x, C2V(c1), C2V(c0))
#define EstrinPoly3(x, x2, c2, c1, c0) MLA(x2, C2V(c2), MLA(x, C2V(c1), C2V(c0)))
#define EstrinPoly4(x, x2, c3, c2, c1, c0) MLA(x2, MLA(x, C2V(c3), C2V(c2)), MLA(x, C2V(c1), C2V(c0)))
#define EstrinPoly5(x, x2, x4, c4, c3, c2, c1, c0) MLA(x4, C2V(c4), EstrinPoly4(x, x2, c3, c2, c1, c0))
#define EstrinPoly6(x, x2, x4, c5, c4, c3, c2, c1, c0) MLA(x4, EstrinPoly2(x, c5, c4), EstrinPoly4(x, x2, c3, c2, c1, c0))
#define EstrinPoly7(x, x2, x4, c6, c5, c4, c3, c2, c1, c0) MLA(x4, EstrinPoly3(x, x2, c6, c5, c4), EstrinPoly4(x, x2, c3, c2, c1, c0))
#define EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0) MLA(x4, EstrinPoly4(x, x2, c7, c6, c5, c4), EstrinPoly4(x, x2, c3, c2, c1, c0))
#define EstrinPoly9(x, x2, x4, x8, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, C2V(c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly10(x, x2, x4, x8, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly2(x, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly11(x, x2, x4, x8, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly3(x, x2, ca, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly12(x, x2, x4, x8, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly4(x, x2, cb, ca, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly13(x, x2, x4, x8, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly5(x, x2, x4, cc, cb, ca, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly14(x, x2, x4, x8, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly6(x, x2, x4, cd, cc, cb, ca, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly15(x, x2, x4, x8, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly7(x, x2, x4, ce, cd, cc, cb, ca, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x8, EstrinPoly8(x, x2, x4, cf, ce, cd, cc, cb, ca, c9, c8), EstrinPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly17(x, x2, x4, x8, x16, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, C2V(d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly18(x, x2, x4, x8, x16, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly2(x, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly19(x, x2, x4, x8, x16, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly3(x, x2, d2, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly20(x, x2, x4, x8, x16, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly4(x, x2, d3, d2, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly21(x, x2, x4, x8, x16, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly5(x, x2, x4, d4, d3, d2, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly22(x, x2, x4, x8, x16, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly6(x, x2, x4, d5, d4, d3, d2, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly23(x, x2, x4, x8, x16, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly7(x, x2, x4, d6, d5, d4, d3, d2, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))
#define EstrinPoly24(x, x2, x4, x8, x16, d7, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0)\
  MLA(x16, EstrinPoly8(x, x2, x4, d7, d6, d5, d4, d3, d2, d1, d0), EstrinPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0))



// These are macros for evaluating polynomials using Horner's method

// EXtended parameter list to match EstrinPoly macro parameter list
#define EXHornerPoly2(x, c1, c0) \
            HornerPoly2(x, c1, c0)
#define EXHornerPoly3(x, x2, c2, c1, c0) \
            ( ((void)x2), \
            HornerPoly3(x, c2, c1, c0) ) 
#define EXHornerPoly4(x, x2, c3, c2, c1, c0) \
            ( ((void)x2), \
            HornerPoly4(x, c3, c2, c1, c0) )
#define EXHornerPoly5(x, x2, x4, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), \
            HornerPoly5(x, c4, c3, c2, c1, c0) )
#define EXHornerPoly6(x, x2, x4, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), \
            HornerPoly6(x, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly7(x, x2, x4, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), \
            HornerPoly7(x, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly8(x, x2, x4, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), \
            HornerPoly8(x, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly9(x, x2, x4, x8, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly9(x, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly10(x, x2, x4, x8, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly10(x, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly11(x, x2, x4, x8, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly11(x, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly12(x, x2, x4, x8, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly12(x, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly13(x, x2, x4, x8, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly13(x, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly14(x, x2, x4, x8, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly14(x, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly15(x, x2, x4, x8, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly15(x, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly16(x, x2, x4, x8, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), \
            HornerPoly16(x, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly17(x, x2, x4, x8, x16, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly17(x, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly18(x, x2, x4, x8, x16, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly18(x, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly19(x, x2, x4, x8, x16, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly19(x, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly20(x, x2, x4, x8, x16, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly20(x, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly21(x, x2, x4, x8, x16, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly21(x, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly22(x, x2, x4, x8, x16, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly22(x, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly23(x, x2, x4, x8, x16, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly23(x, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )
#define EXHornerPoly24(x, x2, x4, x8, x16, d7, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            ( ((void)x2), ((void)x4), ((void)x8), ((void)x16), \
            HornerPoly24(x, d7, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) )

  
// Standard Parameters versions to match EstrinPoly macros
#define HornerPoly2(x, c1, c0) \
            MLA(x, C2V(c1), C2V(c0))
#define HornerPoly3(x, c2, c1, c0) \
            MLA(x, HornerPoly2(x, c2, c1), C2V(c0))
#define HornerPoly4(x, c3, c2, c1, c0) \
            MLA(x, HornerPoly3(x, c3, c2, c1), C2V(c0))
#define HornerPoly5(x, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly4(x, c4, c3, c2, c1), C2V(c0))
#define HornerPoly6(x, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly5(x, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly7(x, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly6(x, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly8(x, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly7(x, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly9(x, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly8(x, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly10(x, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly9(x, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly11(x, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly10(x, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly12(x, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly11(x, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly13(x, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly12(x, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly14(x, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly13(x, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly15(x, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly14(x, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly16(x, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly15(x, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly17(x, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly16(x, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly18(x, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly17(x, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly19(x, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly18(x, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly20(x, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly19(x, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly21(x, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly20(x, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly22(x, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly21(x, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly23(x, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly22(x, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))
#define HornerPoly24(x, d7, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) \
            MLA(x, HornerPoly23(x, d7, d6, d5, d4, d3, d2, d1, d0, cf, ce, cd, cc, cb, ca, c9, c8, c7, c6, c5, c4, c3, c2, c1), C2V(c0))


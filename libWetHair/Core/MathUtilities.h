//
// This file is part of the libWetHair open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// The Clear BSD License
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef __MATH_UTILITIES_H__
#define __MATH_UTILITIES_H__

#include <Eigen/Core>
#include <iostream>
#include <algorithm>
#include "ThreadUtils.h"
#include "MathDefs.h"

template< int DIM >
class TwoDScene;

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989


#ifndef M_PI
const double M_PI = PI;
#endif

#ifdef WIN32
#undef min
#undef max
#endif

using std::min;
using std::max;
using std::swap;


namespace mathutils
{

  template<class T>
  inline T sqr(const T& x)
  { return x*x; }
  
  template<class T>
  inline T cube(const T& x)
  { return x*x*x; }
  
  template<class T>
  inline T min(T a1, T a2, T a3)
  { return min(a1, min(a2, a3)); }
  
  template<class T>
  inline T min(T a1, T a2, T a3, T a4)
  { return min(min(a1, a2), min(a3, a4)); }
  
  template<class T>
  inline T min(T a1, T a2, T a3, T a4, T a5)
  { return min(min(a1, a2), min(a3, a4), a5); }
  
  template<class T>
  inline T min(T a1, T a2, T a3, T a4, T a5, T a6)
  { return min(min(a1, a2), min(a3, a4), min(a5, a6)); }
  
  template<class T>
  inline T max(T a1, T a2, T a3)
  { return max(a1, max(a2, a3)); }
  
  template<class T>
  inline T max(T a1, T a2, T a3, T a4)
  { return max(max(a1, a2), max(a3, a4)); }
  
  template<class T>
  inline T max(T a1, T a2, T a3, T a4, T a5)
  { return max(max(a1, a2), max(a3, a4),  a5); }
  
  template<class T>
  inline T max(T a1, T a2, T a3, T a4, T a5, T a6)
  { return max(max(a1, a2), max(a3, a4),  max(a5, a6)); }
  
  template<class T>
  inline void minmax(T a1, T a2, T& amin, T& amax)
  {
    if(a1<a2){
      amin=a1;
      amax=a2;
    }else{
      amin=a2;
      amax=a1;
    }
  }
  
  template<class T>
  inline void minmax(T a1, T a2, T a3, T& amin, T& amax)
  {
    if(a1<a2){
      if(a1<a3){
        amin=a1;
        if(a2<a3) amax=a3;
        else amax=a2;
      }else{
        amin=a3;
        if(a1<a2) amax=a2;
        else amax=a1;
      }
    }else{
      if(a2<a3){
        amin=a2;
        if(a1<a3) amax=a3;
        else amax=a1;
      }else{
        amin=a3;
        amax=a1;
      }
    }
  }
  
  template<class T>
  inline void minmax(T a1, T a2, T a3, T a4, T& amin, T& amax)
  {
    if(a1<a2){
      if(a3<a4){
        amin=min(a1,a3);
        amax=max(a2,a4);
      }else{
        amin=min(a1,a4);
        amax=max(a2,a3);
      }
    }else{
      if(a3<a4){
        amin=min(a2,a3);
        amax=max(a1,a4);
      }else{
        amin=min(a2,a4);
        amax=max(a1,a3);
      }
    }
  }
  
  template<class T>
  inline void minmax(T a1, T a2, T a3, T a4, T a5, T& amin, T& amax)
  {
    //@@@ the logic could be shortcircuited a lot!
    amin=min(a1,a2,a3,a4,a5);
    amax=max(a1,a2,a3,a4,a5);
  }
  
  template<class T>
  inline void minmax(T a1, T a2, T a3, T a4, T a5, T a6, T& amin, T& amax)
  {
    //@@@ the logic could be shortcircuited a lot!
    amin=min(a1,a2,a3,a4,a5,a6);
    amax=max(a1,a2,a3,a4,a5,a6);
  }
  
  template<class T>
  inline void update_minmax(T a1, T& amin, T& amax)
  {
    if(a1<amin) amin=a1;
    else if(a1>amax) amax=a1;
  }
  
  // swap so that a<b
  template<class T>
  inline void sort(T &a, T &b)
  {
    if(a>b) swap(a,b);
  }
  
  // swap so that a<b<c
  template<class T>
  inline void sort(T &a, T &b, T &c)
  {
    if(a>b) swap(a,b);
    if(a>c) swap(a,c);
    if(b>c) swap(b,c);
  }
  
  // swap so that a<b<c<d
  template<class T>
  inline void sort(T &a, T &b, T &c, T &d)
  {
    if(a>b) swap(a,b);
    if(c>d) swap(c,d);
    if(a>c) swap(a,c);
    if(b>d) swap(b,d);
    if(b>c) swap(b,c);
  }
  
  template<class T>
  inline T clamp(T a, T lower, T upper)
  {
    if(a<lower) return lower;
    else if(a>upper) return upper;
    else return a;
  }
  
  // only makes sense with T=float or double
  template<class T>
  inline T smooth_step(T r)
  {
    if(r<0) return 0;
    else if(r>1) return 1;
    return r*r*r*(10+r*(-15+r*6));
  }
  
  template<class T>
  inline int bipart_closest(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v, T val)
  {
    std::vector<T> buff(v.size());
    
    memcpy(&buff[0], v.data(), v.size() * sizeof(T));
    
    std::sort(buff.begin(), buff.end());
    
    auto const it = std::lower_bound(buff.begin(), buff.end(), val);
    if (it == buff.end()) { return v.size() - 1; }
    
    return std::distance(buff.begin(), it);
  }
  
  // only makes sense with T=float or double
  template<class T>
  inline T smooth_step(T r, T r_lower, T r_upper, T value_lower, T value_upper)
  { return value_lower + smooth_step((r-r_lower)/(r_upper-r_lower)) * (value_upper-value_lower); }
  
  // only makes sense with T=float or double
  template<class T>
  inline T ramp(T r)
  { return smooth_step((r+1)/2)*2-1; }
  
#ifdef WIN32
  inline int lround(double x)
  {
    if(x>0)
      return (x-floor(x)<0.5) ? (int)floor(x) : (int)ceil(x);
    else
      return (x-floor(x)<=0.5) ? (int)floor(x) : (int)ceil(x);
  }
  
  inline double remainder(double x, double y)
  {
    return x-std::floor(x/y+0.5)*y;
  }
#endif
  
  inline unsigned int round_up_to_power_of_two(unsigned int n)
  {
    int exponent=0;
    --n;
    while(n){
      ++exponent;
      n>>=1;
    }
    return 1<<exponent;
  }
  
  inline unsigned int round_down_to_power_of_two(unsigned int n)
  {
    int exponent=0;
    while(n>1){
      ++exponent;
      n>>=1;
    }
    return 1<<exponent;
  }
  
  // Transforms even the sequence 0,1,2,3,... into reasonably good random numbers
  // Challenge: improve on this in speed and "randomness"!
  // This seems to pass several statistical tests, and is a bijective map (of 32-bit unsigned ints)
  inline unsigned int randhash(unsigned int seed)
  {
    unsigned int i=(seed^0xA3C59AC3u)*2654435769u;
    i^=(i>>16);
    i*=2654435769u;
    i^=(i>>16);
    i*=2654435769u;
    return i;
  }
  
  // the inverse of randhash
  inline unsigned int unhash(unsigned int h)
  {
    h*=340573321u;
    h^=(h>>16);
    h*=340573321u;
    h^=(h>>16);
    h*=340573321u;
    h^=0xA3C59AC3u;
    return h;
  }
  
  // returns repeatable stateless pseudo-random number in [0,1]
  inline double randhashd(unsigned int seed)
  { return randhash(seed)/(double)UINT_MAX; }
  inline float randhashf(unsigned int seed)
  { return randhash(seed)/(float)UINT_MAX; }
  
  // returns repeatable stateless pseudo-random number in [a,b]
  inline double randhashd(unsigned int seed, double a, double b)
  { return (b-a)*randhash(seed)/(double)UINT_MAX + a; }
  inline float randhashf(unsigned int seed, float a, float b)
  { return ( (b-a)*randhash(seed)/(float)UINT_MAX + a); }
  
  inline int intlog2(int x)
  {
    int exp=-1;
    while(x){
      x>>=1;
      ++exp;
    }
    return exp;
  }
  
  template<class T>
  inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
  {
    T s=std::floor(x);
    i=(int)s;
    if(i<i_low){
      i=i_low;
      f=0;
    }else if(i>i_high-2){
      i=i_high-2;
      f=1;
    }else
      f=(T)(x-s);
  }
  
  template<class S, class T>
  inline S lerp(const S& value0, const S& value1, T f)
  { return (1-f)*value0 + f*value1; }
  
  template<typename S, typename T>
  inline S lerp_weno( const S value[], T f )
  {
    S p1 = value[0] + (value[1] - value[0]) * (f + 2.0) + (value[2] - 2.0 * value[1] + value[0]) * (f + 2.0) * (f + 1.0) * 0.5
    + (value[3] - 3.0 * value[2] + 3.0 * value[1] - value[0]) * (f + 2.0) * (f + 1.0) * f / 6.0;
    S p2 = value[1] + (value[2] - value[1]) * (f + 1.0) + (value[3] - 2.0 * value[2] + value[1]) * (f + 1.0) * f * 0.5
    + (value[4] - 3.0 * value[3] + 3.0 * value[2] - value[1]) * (f + 1.0) * f * (f - 1.0) / 6.0;
    S p3 = value[2] + (value[3] - value[2]) * f + (value[4] - 2.0 * value[3] + value[2]) * f * (f - 1.0) * 0.5
    + (value[5] - 3.0 * value[4] + 3.0 * value[3] - value[2]) * f * (f - 1.0) * (f - 2.0) / 6.0;
    
    T C1 = (2 - f) * (3 - f) / 20.0;
    T C2 = (3 - f) * (f + 2) / 10.0;
    T C3 = (f + 2) * (f + 1) / 20.0;
    
    T IS1 = (814.0 * value[3] * value[3] + 4326 * value[2] * value[2] + 2976 * value[1] * value[1] + 244 * value[0] * value[0] - 3579 * value[2] * value[3] - 6927 * value[2] * value[1]
             + 1854 * value[2] * value[0] + 2634 * value[3] * value[1] - 683 * value[3] * value[0] - 1659 * value[1] * value[0])   / 180.0;
    T IS2 = (1986 * value[3] * value[3] + 1986 * value[2] * value[2] + 244 * value[1] * value[1] + 244 * value[4] * value[4] + 1074 * value[2] * value[4] - 3777 * value[2] * value[3]
             - 1269 * value[2] * value[1] + 1074 * value[3] * value[1] - 1269 * value[4] * value[3] - 293 * value[4] * value[1])  / 180.0;
    T IS3 = (814 * value[2] * value[2] + 4326 * value[3] * value[3] + 2976 * value[4] * value[4] + 244 * value[5] * value[5] - 683 * value[2] * value[5] + 2634 * value[2] * value[4]
             - 3579 * value[2] * value[3] - 6927 * value[3] * value[4] + 1854 * value[3] * value[5] - 1659 * value[4] * value[5]) / 180.0;
    
    const T epsilon = 1e-6;
    T alpha1 = C1 / ((IS1 + epsilon) * (IS1 + epsilon));
    T alpha2 = C2 / ((IS2 + epsilon) * (IS2 + epsilon));
    T alpha3 = C3 / ((IS3 + epsilon) * (IS3 + epsilon));
    
    T sumalpha = alpha1 + alpha2 + alpha3;
    T w1 = alpha1 / sumalpha;
    T w2 = alpha2 / sumalpha;
    T w3 = alpha3 / sumalpha;
    
    return p1 * w1 + p2 * w2 + p3 * w3;
  }
  
  template<class S, class T>
  inline S bilerp(const S& v00, const S& v10,
                  const S& v01, const S& v11,
                  T fx, T fy)
  {
    return lerp(lerp(v00, v10, fx),
                lerp(v01, v11, fx),
                fy);
  }
  
  template<class T>
  inline Eigen::Matrix<T, 2, 1> grad_bilerp(const T& v00, const T& v10,
                                            const T& v01, const T& v11,
                                            T fx, T fy)
  {
    return Eigen::Matrix<T, 2, 1>(fy - 1.0, fx - 1.0) * v00 +
    Eigen::Matrix<T, 2, 1>(1.0 - fy, -fx) * v10 +
    Eigen::Matrix<T, 2, 1>(-fy, 1.0 - fx) * v01 +
    Eigen::Matrix<T, 2, 1>(fy, fx) * v11;
  }
  
  template<class T>
  inline Eigen::Matrix<T, 3, 1> grad_trilerp(const T& v000, const T& v100,
                                             const T& v010, const T& v110,
                                             const T& v001, const T& v101,
                                             const T& v011, const T& v111,
                                            T fx, T fy, T fz)
  {
    return
    Eigen::Matrix<T, 3, 1>(-(fy - 1.) * (fz - 1.), -(fx - 1.) * (fz - 1.), -(fx - 1.) * (fy - 1.)) * v000 +
    Eigen::Matrix<T, 3, 1>((fy - 1.) * (fz - 1.), fx * (fz - 1.), fx * (fy - 1.)) * v100 +
    Eigen::Matrix<T, 3, 1>(fy * (fz - 1.), (fx - 1.) * (fz - 1.), fy * (fx - 1.) ) * v010 +
    Eigen::Matrix<T, 3, 1>(-fy * (fz - 1.), -fx * (fz - 1.), -fx * fy ) * v110 +
    Eigen::Matrix<T, 3, 1>(fz * (fy - 1.), fz * (fx - 1.), (fx - 1.) * (fy - 1.)) * v001 +
    Eigen::Matrix<T, 3, 1>(-fz * (fy - 1.), -fx * fz, -fx * (fy - 1.)) * v101 +
    Eigen::Matrix<T, 3, 1>(-fy * fz, -fz * (fx - 1.), -fy * (fx - 1.) ) * v011 +
    Eigen::Matrix<T, 3, 1>(fy * fz, fx * fz, fx * fy ) * v111;
  }
  
  
  template<class S, class T>
  inline S trilerp(const S& v000, const S& v100,
                   const S& v010, const S& v110,
                   const S& v001, const S& v101,
                   const S& v011, const S& v111,
                   T fx, T fy, T fz)
  {
    return lerp(bilerp(v000, v100, v010, v110, fx, fy),
                bilerp(v001, v101, v011, v111, fx, fy),
                fz);
  }
  
  template<class S, class T>
  inline S quadlerp(const S& v0000, const S& v1000,
                    const S& v0100, const S& v1100,
                    const S& v0010, const S& v1010,
                    const S& v0110, const S& v1110,
                    const S& v0001, const S& v1001,
                    const S& v0101, const S& v1101,
                    const S& v0011, const S& v1011,
                    const S& v0111, const S& v1111,
                    T fx, T fy, T fz, T ft)
  {
    return lerp(trilerp(v0000, v1000, v0100, v1100, v0010, v1010, v0110, v1110, fx, fy, fz),
                trilerp(v0001, v1001, v0101, v1101, v0011, v1011, v0111, v1111, fx, fy, fz),
                ft);
  }
  
  // f should be between 0 and 1, with f=0.5 corresponding to balanced weighting between w0 and w2
  template<class T>
  inline void quadratic_bspline_weights(T f, T& w0, T& w1, T& w2)
  {
    w0=T(0.5)*sqr(f-1);
    w1=T(0.75)-sqr(f-T(0.5));;
    w2=T(0.5)*sqr(f);
  }
  
  // f should be between 0 and 1
  template<class T>
  inline void cubic_interp_weights(T f, T& wneg1, T& w0, T& w1, T& w2)
  {
    T f2(f*f), f3(f2*f);
    wneg1=-T(1./3)*f+T(1./2)*f2-T(1./6)*f3;
    w0=1-f2+T(1./2)*(f3-f);
    w1=f+T(1./2)*(f2-f3);
    w2=T(1./6)*(f3-f);
  }
  
  template<class S, class T>
  inline S cubic_interp(const S& value_neg1, const S& value0, const S& value1, const S& value2, T f)
  {
    T wneg1, w0, w1, w2;
    cubic_interp_weights(f, wneg1, w0, w1, w2);
    return wneg1*value_neg1 + w0*value0 + w1*value1 + w2*value2;
  }
  
  template<class T>
  void zero(std::vector<T>& v)
  { for(int i=(int)v.size()-1; i>=0; --i) v[i]=0; }
  
  template<class T>
  T abs_max(const std::vector<T>& v)
  {
    T m=0;
    for(int i=(int)v.size()-1; i>=0; --i){
      if(std::fabs(v[i])>m)
        m=std::fabs(v[i]);
    }
    return m;
  }
  
  template<class T>
  bool contains(const std::vector<T>& a, T e)
  {
    for(unsigned int i=0; i<a.size(); ++i)
      if(a[i]==e) return true;
    return false;
  }
  
  template<class T>
  void add_unique(std::vector<T>& a, T e)
  {
    for(unsigned int i=0; i<a.size(); ++i)
      if(a[i]==e) return;
    a.push_back(e);
  }
  
  template<class T>
  void insert(std::vector<T>& a, unsigned int index, T e)
  {
    a.push_back(a.back());
    for(unsigned int i=(unsigned int)a.size()-1; i>index; --i)
      a[i]=a[i-1];
    a[index]=e;
  }
  
  template<class T>
  void erase(std::vector<T>& a, unsigned int index)
  {
    for(unsigned int i=index; i<a.size()-1; ++i)
      a[i]=a[i+1];
    a.pop_back();
  }
  
  template<class T>
  void erase_swap(std::vector<T>& a, unsigned int index)
  {
    for(unsigned int i=index; i<a.size()-1; ++i)
      swap(a[i], a[i+1]);
    a.pop_back();
  }
  
  template<class T>
  void erase_unordered(std::vector<T>& a, unsigned int index)
  {
    a[index]=a.back();
    a.pop_back();
  }
  
  template<class T>
  void erase_unordered_swap(std::vector<T>& a, unsigned int index)
  {
    swap(a[index], a.back());
    a.pop_back();
  }
  
  template<class T>
  void find_and_erase_unordered(std::vector<T>& a, const T& doomed_element)
  {
    for(unsigned int i=0; i<a.size(); ++i)
      if(a[i]==doomed_element){
        erase_unordered(a, i);
        return;
      }
  }
  
  template<class T>
  void replace_once(std::vector<T>& a, const T& old_element, const T& new_element)
  {
    for(unsigned int i=0; i<a.size(); ++i)
      if(a[i]==old_element){
        a[i]=new_element;
        return;
      }
  }
  
  template<class T>
  void write_matlab(std::ostream& output, const std::vector<T>& a, const char *variable_name, bool column_vector=true, int significant_digits=18)
  {
    output<<variable_name<<"=[";
    std::streamsize old_precision=output.precision();
    output.precision(significant_digits);
    for(unsigned int i=0; i<a.size(); ++i){
      output<<a[i]<<" ";
    }
    output<<"]";
    if(column_vector)
      output<<"'";
    output<<";"<<std::endl;
    output.precision(old_precision);
  }
  
  template<class T>
  void write_matlab(std::ostream& output, const Eigen::SparseMatrix<T>& mat, const char *variable_name, int significant_digits=18)
  {
    output<<variable_name<<"=spconvert([";
    std::streamsize old_precision=output.precision();
    output.precision(significant_digits);
    for (int k=0; k<mat.outerSize(); ++k) {
      for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat,k); it; ++it)
      {
        if(it.value() != 0.0) output << (it.row()+1) << ", " << (it.col()+1) << ", " << it.value() << std::endl;
      }
    }
    output<<"]);"<<std::endl;
    output.precision(old_precision);
  }
  
  const static scalar bdf_weights[6][13] =
  {
    {1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {9.0/4.0, -6.0, 11.0/2.0, -2.0, 1.0/4.0, 0, 0, 0, 0, 0, 0, 0, 0},
    {121.0/36.0, -11.0, 29.0/2.0, -92.0/9.0, 51.0/12.0, -1.0, 1.0/9.0, 0, 0, 0, 0, 0, 0},
    {625.0/144.0, -50.0/3.0, 57.0/2.0, -266.0/9.0, 497.0/24.0, -10.0, 59.0/18.0, -2.0/3.0, 1.0/16.0, 0, 0, 0, 0},
    {18769.0/3600.0, -137.0/6.0, 287.0/6.0, -587.0/9.0, 1537.0/24.0, -3506.0/75.0, 461.0/18.0, -31.0/3.0, 139.0/48.0, -1.0/2.0, 1.0/25.0, 0, 0},
    {2401.0/400.0, -147.0/5.0, 291.0/4.0, -368.0/3.0, 1237.0/8.0, -3772.0/25.0, 5216.0/45.0, -70.0, 521.0/16.0, -101.0/9.0, 269.0/100.0, -2.0/5.0, 1.0/36.0}
  };
  
  const static scalar bdf_velocity_weights[6][7] =
  {
    {1.0, -1.0, 0, 0, 0, 0, 0},
    {3.0/2.0, -2.0, 1.0/2.0, 0, 0, 0, 0},
    {11.0/6.0, -3.0, 3.0/2.0, -1.0/3.0, 0, 0, 0},
    {25.0/12.0, -4.0, 3.0, -4.0/3.0, 1.0/4.0, 0, 0},
    {137.0/60.0, -5.0, 5.0, -10.0/3.0, 5.0/4.0, -1.0/5.0, 0},
    {49.0/20.0, -6.0, 15.0/2.0, -20.0/3.0, 15.0/4.0, -6.0/5.0, 1.0/6.0}
  };
  
  const static scalar gauss_legendre_point[5][5] =
  {
    {0.5, 0, 0, 0, 0},
    {0.2113248654051871344705659794271923601627349853515625, 0.7886751345948128655294340205728076398372650146484375, 0, 0, 0},
    {0.112701665379258297861042592558078467845916748046875, 0.5, 0.887298334620741702138957407441921532154083251953125, 0, 0},
    {0.069431844202973713731097404888714663684368133544921875, 0.33000947820757187134432797392946667969226837158203125, 0.66999052179242812865567202607053332030773162841796875, 0.9305681557970262307577513638534583151340484619140625, 0},
    {0.046910077030668018149839326724759303033351898193359375, 0.2307653449471585016539165735593996942043304443359375, 0.5, 0.7692346550528414983460834264406003057956695556640625, 0.953089922969331926339009442017413675785064697265625}
  };
  
  const static scalar gauss_legendre_weight[5][5] =
  {
    {1, 0, 0, 0, 0},
    {0.5, 0.5, 0, 0, 0},
    {5.0/18.0, 4.0/9.0, 5.0/18.0, 0, 0},
    {0.1739274225687269248563637802362791262567043304443359375, 0.326072577431273102899211835392634384334087371826171875, 0.326072577431273102899211835392634384334087371826171875, 0.1739274225687269248563637802362791262567043304443359375},
    {0.118463442528094542449679238416138105094432830810546875, 0.2393143352496832354514566532088792882859706878662109375, 64./255., 0.2393143352496832354514566532088792882859706878662109375, 0.118463442528094542449679238416138105094432830810546875}
  };
  
  inline int find_gauss_legendre_index(const scalar& alpha, int N_GAUSS)
  {
    int i = N_GAUSS - 1;
    for(; i >= 0; --i)
    {
      if(alpha >= gauss_legendre_point[N_GAUSS - 1][i]) return i;
    }
    
    return i;
  }
  
  inline int mod(int x, int m) {
    return (x%m + m)%m;
  }
  
  bool approxSymmetric( const MatrixXs& A, const scalar& eps );
  
  inline Vector2s rotate(const Vector2s& v, const scalar& a)
  {
    Matrix2s A;
    A << cos(a), -sin(a), sin(a), cos(a);
    
    return A * v;
  }
  
/*  needs update to work with twist vars
  template<typename T, int DIM>
  inline void swapVector(std::vector<T>& x)
  {
    int np = x.size();
    int npp = x.size() / DIM;
    for(int i = 0; i < npp; ++i)
    {
      std::swap(x[i], x[np - i - 1]);
    }
  }

  template<int DIM>
  inline void swapVectorNd(VectorXs& x)
  {
    int np = x.size() / DIM;
    int npp = np / DIM;
    Vectors<DIM> v;
    for(int i = 0; i < npp; ++i)
    {
      v = x.segment<DIM>(i * DIM);
      x.segment<DIM>(i * DIM) = x.segment<DIM>((np - i - 1) * DIM);
      x.segment<DIM>((np - i - 1) * DIM) = v;
    }
  }
*/
  template<typename T>
  inline T linear_kernel(const Eigen::Matrix<T, 2, 1>& d, const T& h) {
    return std::max( ((T) 1.0 - (T) fabs(d(0) / h)) * ((T) 1.0 - (T) fabs(d(1) / h)), (T) 0.0 );
  }
  
  template<typename T>
  inline T linear_kernel(const Eigen::Matrix<T, 3, 1>& d, const T& h) {
    return std::max( ((T) 1.0 - (T) fabs(d(0) / h)) * ((T) 1.0 - (T) fabs(d(1) / h)) * ((T) 1.0 - (T) fabs(d(2) / h)), (T) 0.0 );
  }
  
  template<typename T>
  inline T smooth_kernel( const T& r2, const T& h ) {
    return std::max( (T) pow((T) 1.0 - r2 / (h * h), (T) 3.0), (T) 0.0 );
  }
  
  template<typename T>
  inline T smooth_kernel_laplacian( const T& r2, const T& h ) {
    T x2 = (T) sqrt( r2 / (h*h));
    return x2 > (T) 1.0 ? (T) 0.0 : ((T) 1.0 - x2);
  }
  
  template<typename T>
  inline T sharp_kernel( const T& r2, const T& h ) {
    return std::max( (T) (h * h / std::max(r2, (T) 1.0e-5) - (T) 1.0), (T) 0.0 );
  }
  
  // Equ. (2) in [Akinci et al. 2013]
  template<typename T>
  inline T akinci_cohesion_kernel( const T& r, const T& h ) {
    T coeff = (T) (32. / M_PI) / pow(h, (T) 9.);
    
    if((T) 2. * r > h && r <= h)
    {
      return coeff * (h - r) * (h - r) * (h - r) * r * r * r;
    } else if(r > (T) 0. && (T) 2. * r <= h)
    {
      return coeff * ((T) 2. * (h - r) * (h - r) * (h - r) * r * r * r - pow(h, (T) 6.) / (T) 64.);
    } else {
      return (T) 0.0;
    }
  }
  
  template<typename T>
  inline T grad_akinci_cohesion_kernel( const T& r, const T& h ) {
    T coeff = (T) (32. / M_PI) / pow(h, (T) 9.);
    
    if((T) 2. * r > h && r <= h)
    {
      return coeff * (r*r) * (h-r) * (h-r) * (h-r*2.0) *3.0;
    } else if(r > (T) 0. && (T) 2. * r <= h)
    {
      return coeff * (r*r) * (h-r) * (h-r) * (h-r*2.0) *6.0;
    } else {
      return (T) 0.0;
    }
  }

  template<typename T>
  inline T poly6_kernel_2d ( const T& r2, const T& h ) {
    if(r2 <= h*h) return (T) 4.0 / (M_PI * pow(h, (T) 8.0)) * pow(h*h - r2, (T) 3.0);
    else return 0.0;
  }
  
  template<typename T>
  inline Eigen::Matrix<T, 2, 1> grad_poly6_kernel_2d ( const Eigen::Matrix<T, 2, 1>& d, const T& h ) {
    T r2 = d.squaredNorm();
    if(r2 <= h*h) return (T) -24.0 / (M_PI * pow(h, (T) 8.0)) * pow(h*h - r2, (T) 2.0) * d;
    else return Eigen::Matrix<T, 2, 1>::Zero();
  }
  
  template<typename T>
  inline void orientationMatrix(const Eigen::Matrix<T, 3, 1>& t, const T& r, Eigen::Matrix<T, 3, 3>& A)
  {
    const T t2 = t(0)*t(0);
    const T t3 = t(1)*t(1);
    const T t4 = t2+t3;
    if(t4 == 0.0) {
      const T t6 = t(2)*t(2);
      const T t7 = t2+t3+t6;
      if(t7 == 0.0) {
        return;
      }
      A(0, 2) = -r*1.0/sqrt(t6)*t(2);
      A(1, 1) = r;
      A(2, 0) = t(2);
      A(0, 0) = A(0, 1) = A(1, 0) = A(1, 2) = A(2, 1) = A(2, 2) = 0.0;
    } else {
      const T t5 = 1.0/sqrt(t4);
      const T t6 = t(2)*t(2);
      const T t7 = t2+t3+t6;
      const T t8 = 1.0/sqrt(t7);
      A(0, 0) = t(0);
      A(0, 1) = -r*t5*t(1);
      A(0, 2) = -r*t5*t8*t(0)*t(2);
      A(1, 0) = t(1);
      A(1, 1) = r*t5*t(0);
      A(1, 2) = -r*t5*t8*t(1)*t(2);
      A(2, 0) = t(2);
      A(2, 1) = 0.0;
      A(2, 2) = r*sqrt(t4)*t8;
    }
  }
  
  template<typename T>
  inline T softsqrt(const T& x, const T& h)
  {
    if(x < h) {
      return x * (15.0 * h * h - 10.0 * h * x + 3.0 * x * x) / (8.0 * sqrt(h * h * h * h * h));
    } else {
      return sqrt(x);
    }
  }
  
  template<typename T>
  inline T gradsoftsqrt(const T& x, const T& h)
  {
    if(x < h) {
      return (15.0 * h * h - 20.0 * h * x + 9.0 * x * x) / (8.0 * sqrt(h * h * h * h * h));
    } else {
      return 0.5 / sqrt(x);
    }
  }
  
  template<typename T>
  inline T hesssoftsqrt(const T& x, const T& h)
  {
    if(x < h) {
      return -(10.0 * h - 9.0 * x) / (4.0 * sqrt(h * h * h * h * h));
    } else {
      return -0.25 / sqrt(x * x * x);
    }
  }
  
  template<typename T>
  inline T softmin(const T& x, const T& xb, const T& h)
  {
    if (x < xb - h) {
      return x;
    } else if (x < xb + h) {
      return (x * ((h + xb) * (T) 2.0 - x) - (h - xb) * (h - xb)) / ((T) 4.0 * h);
    } else {
      return xb;
    }
  }
  
  template<typename T>
  inline T gradsoftmin(const T& x, const T& xb, const T& h)
  {
    if (x < xb - h) {
      return (T) 1.0;
    } else if (x < xb + h) {
      return (h - x + xb) / ((T) 2.0 * h);
    } else {
      return (T) 0;
    }
  }
  
  template<typename T>
  inline T softmax(const T& x, const T& xa, const T& h)
  {
    if(x < xa - h) {
      return xa;
    } else if (x < xa + h) {
      return (h * h + (T) 2.0 * h * x + (T) 2.0 * h * xa + x * x - (T) 2.0 * x * xa + xa * xa) / ((T) 4.0 * h);
    } else {
      return x;
    }
  }
  
  template<typename T>
  inline T gradsoftmax(const T& x, const T& xa, const T& h)
  {
    if(x < xa - h) {
      return (T) 0;
    } else if (x < xa + h) {
      return (h + x - xa) / ((T) 2.0 * h);
    } else {
      return (T) 1.0;
    }
  }
  
  template<typename T>
  inline T softclamp(const T& x, const T& xa, const T& xb, const T& h)
  {
    if(x < xa - h) {
      return xa;
    } else if (x < xa + h) {
      const T t2 = x*x;
      const T t3 = 1.0/(h*h*h);
      const T t4 = h*h;
      const T t5 = xa*xa;
      return t3*(t4*t5*6.0+(t4*t4)*3.0-t5*t5+h*t4*xa*8.0)*(1.0/1.6E1)-(t2*t2)*t3*(1.0/1.6E1)+t2*t3*(t4*6.0-t5*6.0)*(1.0/1.6E1)+t3*x*(h*t4*8.0-t4*xa*1.2E1+t5*xa*4.0)*(1.0/1.6E1)+t2*t3*x*xa*(1.0/4.0);
    } else if (x < xb - h) {
      return x;
    } else if (x < xb + h) {
      const T t2 = x*x;
      const T t3 = 1.0/(h*h*h);
      const T t4 = h*h;
      const T t5 = xb*xb;
      return t3*(t4*t5*6.0+(t4*t4)*3.0-t5*t5-h*t4*xb*8.0)*(-1.0/1.6E1)+(t2*t2)*t3*(1.0/1.6E1)-t2*t3*(t4*3.0-t5*3.0)*(1.0/8.0)+t3*x*(h*t4*2.0+t4*xb*3.0-t5*xb)*(1.0/4.0)-t2*t3*x*xb*(1.0/4.0);
    } else {
      return xb;
    }
  }
  
  template<typename T>
  inline T gradsoftclamp(const T& x, const T& xa, const T& xb, const T& h)
  {
    if(x < xa - h) {
      return (T) 0;
    } else if (x < xa + h) {
      return (h + x - xa) * (h + x - xa) * (2.0 * h - x + xa) / ((T) 4.0 * h * h * h);
    } else if (x < xb - h) {
      return (T) 1.0;
    } else if (x < xb + h) {
      return (h - x + xb) * (h - x + xb) * (2.0 * h + x - xb) / ((T) 4.0 * h * h * h);
    } else {
      return (T) 0;
    }
  }
  
  template<typename T>
  inline T hesssoftclamp(const T& x, const T& xa, const T& xb, const T& h)
  {
    if(x < xa - h) {
      return (T) 0;
    } else if (x < xa + h) {
      return 3.0 * (h * h - x * x + 2 * x * xa - xa * xa) / (4.0 * h * h * h);
    } else if (x < xb - h) {
      return (T) 0;
    } else if (x < xb + h) {
      return -3.0 * (h * h - x * x + 2 * x * xb - xb * xb) / (4.0 * h * h * h);
    } else {
      return (T) 0;
    }
  }
  
  template<typename T>
  inline T hardclamp(const T& x, const T& a, const T& b)
  {
    return (x < a) ? a : (x > b ? b : x);
  }
  
  template<typename T>
  inline T cross2(const Eigen::Matrix<T, 2, 1>& a, const Eigen::Matrix<T, 2, 1>& b)
  {
    return a(0) * b(1) - a(1) * b(0);
  }

  template<typename T>
  inline T unicross(const Eigen::Matrix<T, 2, 1>& a, const Eigen::Matrix<T, 2, 1>& b)
  {
    return fabs(cross2(a, b));
  }
  
  template<typename T>
  inline bool isnan(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v)
  {
    return isnan(v.sum());
  }
  
  template<typename T>
  inline void check_isnan(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, 1>& v)
  {
#ifndef NDEBUG
    if(std::isnan((double) v.sum()))
    {
      std::cerr << "NAN in " << name << std::endl;
      std::cerr << v << std::endl;
      exit(EXIT_FAILURE);
    }
#endif
  }
  
  
  template<typename T>
  inline T unicross(const Eigen::Matrix<T, 3, 1>& a, const Eigen::Matrix<T, 3, 1>& b)
  {
    return a.cross(b).norm();
  }
  
  template<typename T>
  inline Eigen::Matrix<T, 3, 1> unicross_vec(const Eigen::Matrix<T, 2, 1>& a, const Eigen::Matrix<T, 2, 1>& b)
  {
    return Eigen::Matrix<T, 3, 1>((T)0.0, (T)0.0, cross2(a, b));
  }
  
  template<typename T>
  inline Eigen::Matrix<T, 3, 1> unicross_vec(const Eigen::Matrix<T, 3, 1>& a, const Eigen::Matrix<T, 3, 1>& b)
  {
    return a.cross(b);
  } 
  
  template<typename T, int DIM>
  inline T pointedgedist(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3)
  {
    T alpha = (x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2);
    alpha = std::min((T) 1.0, std::max((T) 0.0, alpha));
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha*(x3-x2);
    Eigen::Matrix<T, DIM, 1> n = closest-x1;
    return n.norm();
  }
  
  template<typename T, int DIM>
  inline void gradsoftcapsule_phi(Eigen::Matrix<T, DIM * 3, 1>& grad, const Eigen::Matrix<T, DIM, 1>& p, const Eigen::Matrix<T, DIM, 1>& a, const Eigen::Matrix<T, DIM, 1>& b, const T& radius, const T& h) {
    Eigen::Matrix<T, DIM, 1> pa = p - a;
    Eigen::Matrix<T, DIM, 1> ba = b - a;
    T sigma = ba.dot(ba);
    T d = pa.dot(ba) / sigma;
    T alpha = softclamp(d, 0.0, 1.0, h);
    Eigen::Matrix<T, DIM, 1> u = pa - ba * alpha;
    Eigen::Matrix<T, DIM, DIM> dudp = Eigen::Matrix<T, DIM, DIM>::Identity() - gradsoftclamp(d, 0.0, 1.0, h) * (ba * ba.transpose()) / sigma;
    Eigen::Matrix<T, DIM, 1> ddda = (ba * (2.0 * d - 1.0) - pa) / sigma;
    Eigen::Matrix<T, DIM, 1> dddb = (pa - 2.0 * d * ba) / sigma;
    
    scalar dscdd = gradsoftclamp(d, 0.0, 1.0, h);
    
    Eigen::Matrix<T, DIM, DIM> duda = (alpha - 1.0) * Eigen::Matrix<T, DIM, DIM>::Identity() - ba * (dscdd * ddda.transpose());
    Eigen::Matrix<T, DIM, DIM> dudb = -alpha * Eigen::Matrix<T, DIM, DIM>::Identity() - ba * (dscdd * dddb.transpose());
    
    grad.segment(0,DIM) = gradsoftsqrt(u.squaredNorm(), h) * 2.0 * dudp.transpose() * u;
    grad.segment(DIM,DIM) = gradsoftsqrt(u.squaredNorm(), h) * 2.0 * duda.transpose() * u;
    grad.segment(DIM * 2,DIM) = gradsoftsqrt(u.squaredNorm(), h) * 2.0 * dudb.transpose() * u;
  }
  
  template<typename T>
  inline void hesssoftcapsule_phi(Eigen::Matrix<T, 6, 6>& hess, const Eigen::Matrix<T, 2, 1>& p, const Eigen::Matrix<T, 2, 1>& a, const Eigen::Matrix<T, 2, 1>& b, const T& radius, const T& h) {
    const T qx0 = p(0);
    const T qx1 = p(1);
    const T qa0 = a(0);
    const T qa1 = a(1);
    const T qb0 = b(0);
    const T qb1 = b(1);
    
    const T t2 = qa0-qb0;
    const T t3 = qa1-qb1;
    const T t5 = t2*t2;
    const T t6 = t3*t3;
    const T t7 = std::max(t5+t6, 1e-20);
    const T t8 = 1.0/t7;
    const T t9 = qa0-qx0;
    const T t10 = t2*t9;
    const T t11 = qa1-qx1;
    const T t12 = t3*t11;
    const T t13 = t10+t12;
    const T t14 = t8*t13;
    const T t15 = softclamp(t14, 0.0, 1.0, h);
    const T t18 = t2*t15;
    const T t4 = -qa0+qx0+t18;
    const T t22 = t3*t15;
    const T t16 = -qa1+qx1+t22;
    const T t19 = gradsoftclamp(t14, 0.0, 1.0, h);
    const T t26 = t5*t8*t19;
    const T t17 = t26-1.0;
    const T t20 = 1.0/(t7*t7);
    const T t21 = hesssoftclamp(t14, 0.0, 1.0, h);
    const T t23 = t4*t4;
    const T t24 = t16*t16;
    const T t25 = t23+t24;
    const T t30 = t4*t17*2.0;
    const T t31 = t2*t3*t8*t16*t19*2.0;
    const T t27 = t30+t31;
    const T t28 = gradsoftsqrt(t25, h);
    const T t29 = hesssoftsqrt(t25, h);
    const T t32 = t6*t8*t19;
    const T t33 = t32-1.0;
    const T t34 = qa0*2.0;
    const T t35 = qb0*2.0;
    const T t36 = t34-t35;
    const T t37 = t13*t20*t36;
    const T t38 = t19*t19;
    const T t39 = qb0+qx0-t34;
    const T t40 = t8*t39;
    const T t41 = t37+t40;
    const T t42 = qa1*2.0;
    const T t43 = qb1*2.0;
    const T t44 = t42-t43;
    const T t45 = qb1+qx1-t42;
    const T t46 = t8*t45;
    const T t47 = t13*t20*t44;
    const T t48 = t46+t47;
    const T t49 = t3*t19*t48;
    const T t50 = -t15+t49+1.0;
    const T t51 = t5*t19*t20*t36;
    const T t53 = t8*t9;
    const T t52 = t37-t53;
    const T t54 = t2*t3*t16*t19*t20*t36*2.0;
    const T t82 = t2*t19*t52;
    const T t55 = t15-t82;
    const T t56 = t5*t19*t20*t44;
    const T t58 = t8*t11;
    const T t57 = t47-t58;
    const T t59 = t2*t3*t16*t19*t20*t44*2.0;
    const T t87 = t3*t19*t57;
    const T t60 = t15-t87;
    const T t61 = t2*t3*t8*t17*t19*2.0;
    const T t62 = t2*t3*t8*t19*t33*2.0;
    const T t63 = t3*t4*t5*t20*t21*2.0;
    const T t64 = t2*t6*t16*t20*t21*2.0;
    const T t65 = t61+t62+t63+t64;
    const T t66 = t28*t65;
    const T t67 = t16*t33*2.0;
    const T t68 = t2*t3*t4*t8*t19*2.0;
    const T t69 = t67+t68;
    const T t70 = t27*t29*t69;
    const T t71 = t66+t70;
    const T t72 = t5*t6*t20*t38*2.0;
    const T t73 = t2*t19*t41;
    const T t74 = -t15+t73+1.0;
    const T t75 = t4*t74*2.0;
    const T t76 = t3*t16*t19*t41*2.0;
    const T t77 = t75+t76;
    const T t78 = t16*t50*2.0;
    const T t79 = t2*t4*t19*t48*2.0;
    const T t80 = t78+t79;
    const T t81 = t6*t19*t20*t36;
    const T t83 = t2*t3*t4*t19*t20*t36*2.0;
    const T t84 = t4*t55*2.0;
    const T t113 = t3*t16*t19*t52*2.0;
    const T t85 = t84-t113;
    const T t86 = t6*t19*t20*t44;
    const T t88 = t2*t3*t4*t19*t20*t44*2.0;
    const T t89 = t16*t60*2.0;
    const T t115 = t2*t4*t19*t57*2.0;
    const T t90 = t89-t115;
    const T t91 = t5*t8*t21*t41;
    const T t139 = t2*t20*t36;
    const T t92 = t8-t139;
    const T t93 = t2*t6*t8*t38*t41*2.0;
    const T t94 = t2*t3*t8*t16*t21*t41*2.0;
    const T t95 = t27*t29*t77;
    const T t96 = t3*t19*t33*t41*2.0;
    const T t97 = t2*t3*t8*t19*t74*2.0;
    const T t98 = t29*t69*t77;
    const T t99 = t41*t41;
    const T t100 = t8*2.0;
    const T t101 = t36*t36;
    const T t102 = 1.0/(t7*t7*t7);
    const T t103 = t13*t101*t102*2.0;
    const T t104 = t20*t36*t39*2.0;
    const T t110 = t13*t20*2.0;
    const T t105 = t100+t103+t104-t110;
    const T t106 = t20*t39*t44;
    const T t107 = t20*t36*t45;
    const T t108 = t13*t36*t44*t102*2.0;
    const T t109 = t106+t107+t108;
    const T t111 = t20*t36*t39;
    const T t150 = t9*t20*t36;
    const T t112 = t8+t103-t110+t111-t150;
    const T t168 = t11*t20*t36;
    const T t114 = t106+t108-t168;
    const T t116 = t2*t8*t19;
    const T t117 = t2*t17*t19*t48*2.0;
    const T t118 = t2*t3*t8*t19*t50*2.0;
    const T t119 = t27*t29*t80;
    const T t120 = t6*t8*t21*t48;
    const T t121 = t33*t50*2.0;
    const T t174 = t3*t20*t44;
    const T t122 = t8-t174;
    const T t123 = t3*t5*t8*t38*t48*2.0;
    const T t124 = t2*t3*t4*t8*t21*t48*2.0;
    const T t125 = t29*t69*t80;
    const T t126 = t19*t41;
    const T t127 = t2*t19*t48*t74*2.0;
    const T t128 = t3*t19*t41*t50*2.0;
    const T t129 = t29*t77*t80;
    const T t130 = t48*t48;
    const T t131 = t44*t44;
    const T t132 = t13*t102*t131*2.0;
    const T t133 = t20*t44*t45*2.0;
    const T t134 = t100-t110+t132+t133;
    const T t135 = t19*t52;
    const T t159 = t9*t20*t44;
    const T t136 = t107+t108-t159;
    const T t137 = t20*t44*t45;
    const T t182 = t11*t20*t44;
    const T t138 = t8-t110+t132+t137-t182;
    const T t140 = t2*t19*t92;
    const T t141 = t5*t8*t21*t52;
    const T t142 = t2*t6*t8*t38*t52*2.0;
    const T t143 = t2*t3*t8*t16*t21*t52*2.0;
    const T t144 = t27*t29*t85;
    const T t145 = t3*t8*t19;
    const T t146 = t2*t3*t19*t20*t36;
    const T t147 = t3*t19*t33*t52*2.0;
    const T t148 = t6*t16*t19*t20*t36*2.0;
    const T t149 = t29*t69*t85;
    const T t151 = t3*t16*t19*t112*2.0;
    const T t152 = t6*t38*t41*t52*2.0;
    const T t153 = t3*t16*t21*t41*t52*2.0;
    const T t154 = t126+t135-t2*t19*t112-t2*t21*t41*t52;
    const T t155 = t151+t152+t153-t55*t74*2.0-t4*t154*2.0;
    const T t156 = t29*t77*t85;
    const T t157 = t156-t28*t155;
    const T t158 = t19*t48;
    const T t160 = t3*t19*t50*t52*2.0;
    const T t161 = t29*t80*t85;
    const T t162 = t37-t53;
    const T t163 = t37-t53;
    const T t164 = t9*t20*t36*2.0;
    const T t165 = -t103+t110+t164;
    const T t166 = t37-t53;
    const T t167 = t19*t57;
    const T t169 = -t108+t159+t168;
    const T t170 = t2*t3*t19*t20*t44;
    const T t171 = t2*t17*t19*t57*2.0;
    const T t172 = t4*t5*t19*t20*t44*2.0;
    const T t173 = t27*t29*t90;
    const T t175 = t3*t19*t122;
    const T t176 = t6*t8*t21*t57;
    const T t177 = t3*t5*t8*t38*t57*2.0;
    const T t178 = t2*t3*t4*t8*t21*t57*2.0;
    const T t179 = t29*t69*t90;
    const T t180 = t2*t19*t57*t74*2.0;
    const T t181 = t29*t77*t90;
    const T t183 = t2*t4*t19*t138*2.0;
    const T t184 = t5*t38*t48*t57*2.0;
    const T t185 = t2*t4*t21*t48*t57*2.0;
    const T t186 = t158+t167-t3*t19*t138-t3*t21*t48*t57;
    const T t187 = t183+t184+t185-t50*t60*2.0-t16*t186*2.0;
    const T t188 = t29*t80*t90;
    const T t189 = t188-t28*t187;
    const T t190 = t2*t19*t55*t57*2.0;
    const T t191 = t3*t19*t52*t60*2.0;
    const T t192 = t29*t85*t90;
    const T t193 = t47-t58;
    const T t194 = t47-t58;
    const T t195 = t11*t20*t44*2.0;
    const T t196 = t110-t132+t195;
    const T t197 = t47-t58;
    hess(0, 0) = t28*(t72+(t17*t17)*2.0+t2*t4*t5*t20*t21*2.0+t3*t5*t16*t20*t21*2.0)+(t27*t27)*t29;
    hess(0, 1) = t71;
    hess(0, 2) = t95+t28*(t54+t93+t94+t17*(-t15+t2*t19*(t37+t8*(qa0*-2.0+qb0+qx0))+1.0)*2.0+t4*(t51+t91-t8*t19*t36)*2.0-t3*t8*t16*t19*2.0);
    hess(0, 3) = t119+t28*(t59+t117+t118+t4*(t56+t5*t8*t21*t48)*2.0-t2*t8*t16*t19*2.0+t2*t3*t8*t16*t21*t48*2.0);
    hess(0, 4) = t144-t28*(t54+t142+t143-t17*t55*2.0+t4*(t51+t141-t8*t19*t36)*2.0-t3*t8*t16*t19*2.0);
    hess(0, 5) = t173-t28*(t59+t171+t4*(t56+t5*t8*t21*t57)*2.0-t2*t8*t16*t19*2.0-t2*t3*t8*t19*t60*2.0+t2*t3*t8*t16*t21*t57*2.0);
    hess(1, 0) = t71;
    hess(1, 1) = t28*(t72+(t33*t33)*2.0+t2*t4*t6*t20*t21*2.0+t3*t6*t16*t20*t21*2.0)+t29*(t69*t69);
    hess(1, 2) = t98+t28*(t83+t96+t97+t16*(t81+t6*t8*t21*t41)*2.0-t3*t4*t8*t19*2.0+t2*t3*t4*t8*t21*t41*2.0);
    hess(1, 3) = t125+t28*(t88+t121+t123+t124+t16*(t86+t120-t8*t19*t44)*2.0-t2*t4*t8*t19*2.0);
    hess(1, 4) = t149-t28*(t83+t147+t16*(t81+t6*t8*t21*t52)*2.0-t3*t4*t8*t19*2.0-t2*t3*t8*t19*t55*2.0+t2*t3*t4*t8*t21*t52*2.0);
    hess(1, 5) = t179-t28*(t88+t177+t178-t33*t60*2.0+t16*(t86+t176-t8*t19*t44)*2.0-t2*t4*t8*t19*2.0);
    hess(2, 0) = t95+t28*(t93+t94+t17*t74*2.0-t4*(-t91+t116+t140)*2.0-t3*t16*t19*t92*2.0);
    hess(2, 1) = t98+t28*(t96+t97+t148+t4*(t146-t3*t8*t19+t2*t3*t8*t21*t41)*2.0+t6*t8*t16*t21*t41*2.0);
    hess(2, 2) = t28*((t74*t74)*2.0+t4*(t19*t41*-2.0+t2*t21*t99+t2*t19*t105)*2.0+t6*t38*t99*2.0+t3*t16*t21*t99*2.0+t3*t16*t19*t105*2.0)+t29*(t77*t77);
    hess(2, 3) = t129+t28*(t127+t128+t4*(-t19*t48+t2*t19*t109+t2*t21*t41*t48)*2.0-t16*t19*t41*2.0+t3*t16*t19*t109*2.0+t3*t16*t21*t41*t48*2.0);
    hess(2, 4) = t157;
    hess(2, 5) = t181-t28*(t180+t4*(-t19*t57+t2*t19*t114+t2*t21*t41*t57)*2.0-t16*t19*t41*2.0-t3*t19*t41*t60*2.0+t3*t16*t19*t114*2.0+t3*t16*t21*t41*t57*2.0);
    hess(3, 0) = t119+t28*(t117+t118+t172+t16*(-t116+t170+t2*t3*t8*t21*t48)*2.0+t4*t5*t8*t21*t48*2.0);
    hess(3, 1) = t125+t28*(t121+t123+t124-t16*(-t120+t145+t175)*2.0-t2*t4*t19*t122*2.0);
    hess(3, 2) = t129+t28*(t127+t128+t16*(-t126+t3*t19*t109+t3*t21*t41*t48)*2.0-t4*t19*t48*2.0+t2*t4*t19*t109*2.0+t2*t4*t21*t41*t48*2.0);
    hess(3, 3) = t28*((t50*t50)*2.0+t16*(t19*t48*-2.0+t3*t21*t130+t3*t19*t134)*2.0+t5*t38*t130*2.0+t2*t4*t21*t130*2.0+t2*t4*t19*t134*2.0)+t29*(t80*t80);
    hess(3, 4) = t161-t28*(t160+t16*(-t135+t3*t19*t136+t3*t21*t48*t52)*2.0-t4*t19*t48*2.0-t2*t19*t48*t55*2.0+t2*t4*t19*t136*2.0+t2*t4*t21*t48*t52*2.0);
    hess(3, 5) = t189;
    hess(4, 0) = t144+t28*(-t142-t143+t17*t55*2.0+t4*(t116+t140-t141)*2.0+t3*t16*t19*t92*2.0);
    hess(4, 1) = t149-t28*(t147+t148+t4*(-t145+t146+t2*t3*t8*t21*t52)*2.0-t2*t3*t8*t19*t55*2.0+t6*t8*t16*t21*t52*2.0);
    hess(4, 2) = t157;
    hess(4, 3) = t161-t28*(t160+t4*(-t158+t2*t19*t136+t2*t21*t48*t52)*2.0-t16*t19*t52*2.0-t2*t19*t48*t55*2.0+t3*t16*t19*t136*2.0+t3*t16*t21*t48*t52*2.0);
    hess(4, 4) = t28*(t4*(t19*t52*2.0-t2*t21*(t162*t162)+t2*t19*t165)*-2.0+(t55*t55)*2.0+t6*t38*(t163*t163)*2.0-t3*t16*t19*t165*2.0+t3*t16*t21*(t166*t166)*2.0)+t29*(t85*t85);
    hess(4, 5) = t192-t28*(t190+t191+t4*(t167+t2*t19*t169-t2*t21*t52*t57)*2.0+t16*t19*t52*2.0+t3*t16*t19*t169*2.0-t3*t16*t21*t52*t57*2.0);
    hess(5, 0) = t173-t28*(t171+t172+t16*(-t116+t170+t2*t3*t8*t21*t57)*2.0-t2*t3*t8*t19*t60*2.0+t4*t5*t8*t21*t57*2.0);
    hess(5, 1) = t179+t28*(-t177-t178+t33*t60*2.0+t16*(t145+t175-t176)*2.0+t2*t4*t19*t122*2.0);
    hess(5, 2) = t181-t28*(t180+t16*(-t126+t3*t19*t114+t3*t21*t41*t57)*2.0-t4*t19*t57*2.0-t3*t19*t41*t60*2.0+t2*t4*t19*t114*2.0+t2*t4*t21*t41*t57*2.0);
    hess(5, 3) = t189;
    hess(5, 4) = t192-t28*(t190+t191+t16*(t135+t3*t19*t169-t3*t21*t52*t57)*2.0+t4*t19*t57*2.0+t2*t4*t19*t169*2.0-t2*t4*t21*t52*t57*2.0);
    hess(5, 5) = t28*(t16*(t19*t57*2.0-t3*t21*(t193*t193)+t3*t19*t196)*-2.0+(t60*t60)*2.0+t5*t38*(t194*t194)*2.0-t2*t4*t19*t196*2.0+t2*t4*t21*(t197*t197)*2.0)+t29*(t90*t90);
  }
  
  template<typename T>
  inline void hesssoftcapsule_phi(Eigen::Matrix<T, 9, 9>& hess, const Eigen::Matrix<T, 3, 1>& p, const Eigen::Matrix<T, 3, 1>& a, const Eigen::Matrix<T, 3, 1>& b, const T& radius, const T& h)
  {
    const scalar t2 = a(0)-b(0);
    const scalar t3 = a(1)-b(1);
    const scalar t4 = a(2)-b(2);
    const scalar t6 = t2*t2;
    const scalar t7 = t3*t3;
    const scalar t8 = t4*t4;
    const scalar t9 = t6+t7+t8;
    const scalar t10 = 1.0/t9;
    const scalar t11 = a(0)-p(0);
    const scalar t12 = t2*t11;
    const scalar t13 = a(1)-p(1);
    const scalar t14 = t3*t13;
    const scalar t15 = a(2)-p(2);
    const scalar t16 = t4*t15;
    const scalar t17 = t12+t14+t16;
    const scalar t18 = t10*t17;
    const scalar t19 = softclamp(t18, 0.0, 1.0, h);
    const scalar t22 = t2*t19;
    const scalar t5 = -a(0)+p(0)+t22;
    const scalar t24 = t3*t19;
    const scalar t20 = -a(1)+p(1)+t24;
    const scalar t26 = t4*t19;
    const scalar t21 = -a(2)+p(2)+t26;
    const scalar t23 = t5*t5;
    const scalar t25 = t20*t20;
    const scalar t27 = t21*t21;
    const scalar t28 = t23+t25+t27;
    const scalar t29 = gradsoftclamp(t18, 0.0, 1.0, h);
    const scalar t30 = t6*t10*t29;
    const scalar t31 = t30-1.0;
    const scalar t32 = t29;
    const scalar t33 = t6*t10*t32;
    const scalar t34 = t33-1.0;
    const scalar t35 = gradsoftsqrt(t28, h);
    const scalar t36 = 1.0/(t9*t9);
    const scalar t37 = hesssoftclamp(t18, 0.0, 1.0, h);
    const scalar t38 = hesssoftsqrt(t28, h);
    const scalar t39 = t5*t34*2.0;
    const scalar t40 = t2*t3*t10*t20*t32*2.0;
    const scalar t41 = t2*t4*t10*t21*t32*2.0;
    const scalar t42 = t39+t40+t41;
    const scalar t43 = t7*t10*t32;
    const scalar t44 = t43-1.0;
    const scalar t45 = t20*t44*2.0;
    const scalar t46 = t2*t3*t5*t10*t32*2.0;
    const scalar t47 = t3*t4*t10*t21*t32*2.0;
    const scalar t48 = t45+t46+t47;
    const scalar t49 = t8*t10*t32;
    const scalar t50 = t49-1.0;
    const scalar t51 = t21*t50*2.0;
    const scalar t52 = t2*t4*t5*t10*t32*2.0;
    const scalar t53 = t3*t4*t10*t20*t32*2.0;
    const scalar t54 = t51+t52+t53;
    const scalar t55 = a(0)*2.0;
    const scalar t56 = b(0)*2.0;
    const scalar t57 = t55-t56;
    const scalar t58 = b(0)+p(0)-t55;
    const scalar t59 = t10*t58;
    const scalar t60 = t17*t36*t57;
    const scalar t61 = t59+t60;
    const scalar t62 = t2*t32*t61;
    const scalar t63 = -t19+t62+1.0;
    const scalar t64 = t5*t63*2.0;
    const scalar t65 = t3*t20*t32*t61*2.0;
    const scalar t66 = t4*t21*t32*t61*2.0;
    const scalar t67 = t64+t65+t66;
    const scalar t68 = a(1)*2.0;
    const scalar t69 = b(1)*2.0;
    const scalar t70 = t68-t69;
    const scalar t71 = t17*t36*t70;
    const scalar t72 = b(1)+p(1)-t68;
    const scalar t73 = t10*t72;
    const scalar t74 = t71+t73;
    const scalar t75 = t3*t32*t74;
    const scalar t76 = -t19+t75+1.0;
    const scalar t77 = t20*t76*2.0;
    const scalar t78 = t2*t5*t32*t74*2.0;
    const scalar t79 = t4*t21*t32*t74*2.0;
    const scalar t80 = t77+t78+t79;
    const scalar t81 = a(2)*2.0;
    const scalar t82 = b(2)*2.0;
    const scalar t83 = t81-t82;
    const scalar t84 = t17*t36*t83;
    const scalar t85 = b(2)+p(2)-t81;
    const scalar t86 = t10*t85;
    const scalar t87 = t84+t86;
    const scalar t88 = t4*t32*t87;
    const scalar t89 = -t19+t88+1.0;
    const scalar t90 = t21*t89*2.0;
    const scalar t91 = t2*t5*t32*t87*2.0;
    const scalar t92 = t3*t20*t32*t87*2.0;
    const scalar t93 = t90+t91+t92;
    const scalar t94 = t6*t29*t36*t57;
    const scalar t96 = t10*t11;
    const scalar t95 = t60-t96;
    const scalar t98 = t2*t32*t95;
    const scalar t97 = t19-t98;
    const scalar t99 = t3*t20*t32*t95*2.0;
    const scalar t100 = t4*t21*t32*t95*2.0;
    const scalar t102 = t5*t97*2.0;
    const scalar t101 = t99+t100-t102;
    const scalar t103 = t2*t3*t20*t29*t35*t36*t57*2.0;
    const scalar t104 = t2*t4*t21*t29*t35*t36*t57*2.0;
    const scalar t105 = t6*t29*t36*t70;
    const scalar t107 = t10*t13;
    const scalar t106 = t71-t107;
    const scalar t109 = t3*t32*t106;
    const scalar t108 = t19-t109;
    const scalar t110 = t2*t5*t32*t106*2.0;
    const scalar t111 = t4*t21*t32*t106*2.0;
    const scalar t113 = t20*t108*2.0;
    const scalar t112 = t110+t111-t113;
    const scalar t114 = t2*t3*t20*t29*t35*t36*t70*2.0;
    const scalar t115 = t2*t4*t21*t29*t35*t36*t70*2.0;
    const scalar t116 = t6*t29*t36*t83;
    const scalar t118 = t10*t15;
    const scalar t117 = t84-t118;
    const scalar t120 = t4*t32*t117;
    const scalar t119 = t19-t120;
    const scalar t121 = t2*t5*t32*t117*2.0;
    const scalar t122 = t3*t20*t32*t117*2.0;
    const scalar t124 = t21*t119*2.0;
    const scalar t123 = t121+t122-t124;
    const scalar t125 = t2*t3*t20*t29*t35*t36*t83*2.0;
    const scalar t126 = t2*t4*t21*t29*t35*t36*t83*2.0;
    const scalar t127 = t3*t5*t6*t35*t36*t37*2.0;
    const scalar t128 = t2*t7*t20*t35*t36*t37*2.0;
    const scalar t129 = t7*t10*t29;
    const scalar t130 = t129-1.0;
    const scalar t131 = t2*t3*t8*t29*t32*t35*t36*2.0;
    const scalar t132 = t2*t3*t4*t21*t35*t36*t37*2.0;
    const scalar t133 = t6*t7*t29*t32*t35*t36*2.0;
    const scalar t134 = t4*t10*t21*t29*t35*2.0;
    const scalar t135 = t7*t29*t36*t57;
    const scalar t136 = t2*t3*t5*t29*t35*t36*t57*2.0;
    const scalar t137 = t3*t4*t21*t29*t35*t36*t57*2.0;
    const scalar t138 = t7*t29*t36*t70;
    const scalar t139 = t2*t3*t5*t29*t35*t36*t70*2.0;
    const scalar t140 = t3*t4*t21*t29*t35*t36*t70*2.0;
    const scalar t141 = t7*t29*t36*t83;
    const scalar t142 = t2*t3*t5*t29*t35*t36*t83*2.0;
    const scalar t143 = t3*t4*t21*t29*t35*t36*t83*2.0;
    const scalar t144 = t4*t5*t6*t35*t36*t37*2.0;
    const scalar t145 = t2*t8*t21*t35*t36*t37*2.0;
    const scalar t146 = t8*t10*t29;
    const scalar t147 = t146-1.0;
    const scalar t148 = t2*t4*t7*t29*t32*t35*t36*2.0;
    const scalar t149 = t2*t3*t4*t20*t35*t36*t37*2.0;
    const scalar t150 = t4*t7*t20*t35*t36*t37*2.0;
    const scalar t151 = t3*t8*t21*t35*t36*t37*2.0;
    const scalar t152 = t3*t4*t6*t29*t32*t35*t36*2.0;
    const scalar t153 = t2*t3*t4*t5*t35*t36*t37*2.0;
    const scalar t154 = t6*t8*t29*t32*t35*t36*2.0;
    const scalar t155 = t7*t8*t29*t32*t35*t36*2.0;
    const scalar t156 = t2*t5*t10*t29*t35*2.0;
    const scalar t157 = t3*t10*t20*t29*t35*2.0;
    const scalar t158 = t8*t29*t36*t57;
    const scalar t159 = t2*t4*t5*t29*t35*t36*t57*2.0;
    const scalar t160 = t3*t4*t20*t29*t35*t36*t57*2.0;
    const scalar t161 = t8*t29*t36*t70;
    const scalar t162 = t2*t4*t5*t29*t35*t36*t70*2.0;
    const scalar t163 = t3*t4*t20*t29*t35*t36*t70*2.0;
    const scalar t164 = t8*t29*t36*t83;
    const scalar t165 = t2*t4*t5*t29*t35*t36*t83*2.0;
    const scalar t166 = t3*t4*t20*t29*t35*t36*t83*2.0;
    const scalar t167 = t10*t17*2.0;
    const scalar t168 = t167-1.0;
    const scalar t169 = t2*t168;
    const scalar t170 = -a(0)+p(0)+t169;
    const scalar t171 = t2*t10*t29*t170;
    const scalar t172 = -t19+t171+1.0;
    const scalar t173 = t6*t10*2.0;
    const scalar t174 = t173-1.0;
    const scalar t175 = t10*t58*2.0;
    const scalar t176 = t17*t36*t57*2.0;
    const scalar t177 = t175+t176;
    const scalar t178 = t2*t177;
    const scalar t179 = -t167+t178+2.0;
    const scalar t180 = t10*t72*2.0;
    const scalar t181 = t17*t36*t70*2.0;
    const scalar t182 = t180+t181;
    const scalar t183 = t10*t85*2.0;
    const scalar t184 = t17*t36*t83*2.0;
    const scalar t185 = t183+t184;
    const scalar t186 = t2*t29*t36*t57*t170;
    const scalar t190 = t10*t11*2.0;
    const scalar t187 = t176-t190;
    const scalar t188 = t2*t187;
    const scalar t189 = -t167+t188+1.0;
    const scalar t191 = t3*t20*t29*t35*t36*t57*t170*2.0;
    const scalar t192 = t4*t21*t29*t35*t36*t57*t170*2.0;
    const scalar t193 = t2*t29*t36*t70*t170;
    const scalar t194 = t3*t20*t29*t35*t36*t70*t170*2.0;
    const scalar t195 = t4*t21*t29*t35*t36*t70*t170*2.0;
    const scalar t197 = t10*t13*2.0;
    const scalar t196 = t181-t197;
    const scalar t198 = t2*t29*t36*t83*t170;
    const scalar t199 = t3*t20*t29*t35*t36*t83*t170*2.0;
    const scalar t200 = t4*t21*t29*t35*t36*t83*t170*2.0;
    const scalar t202 = t10*t15*2.0;
    const scalar t201 = t184-t202;
    const scalar t203 = t3*t168;
    const scalar t204 = -a(1)+p(1)+t203;
    const scalar t205 = t3*t10*t29*t204;
    const scalar t206 = -t19+t205+1.0;
    const scalar t207 = t2*t3*t4*t21*t29*t35*t36*4.0;
    const scalar t208 = t7*t10*2.0;
    const scalar t209 = t208-1.0;
    const scalar t210 = t3*t182;
    const scalar t211 = -t167+t210+2.0;
    const scalar t212 = t3*t29*t36*t57*t204;
    const scalar t213 = t2*t5*t29*t35*t36*t57*t204*2.0;
    const scalar t214 = t4*t21*t29*t35*t36*t57*t204*2.0;
    const scalar t215 = t3*t29*t36*t70*t204;
    const scalar t216 = t3*t196;
    const scalar t217 = -t167+t216+1.0;
    const scalar t218 = t2*t5*t29*t35*t36*t70*t204*2.0;
    const scalar t219 = t4*t21*t29*t35*t36*t70*t204*2.0;
    const scalar t220 = t3*t29*t36*t83*t204;
    const scalar t221 = t2*t5*t29*t35*t36*t83*t204*2.0;
    const scalar t222 = t4*t21*t29*t35*t36*t83*t204*2.0;
    const scalar t223 = t4*t168;
    const scalar t224 = -a(2)+p(2)+t223;
    const scalar t225 = t4*t10*t29*t224;
    const scalar t226 = -t19+t225+1.0;
    const scalar t227 = t2*t3*t4*t20*t29*t35*t36*4.0;
    const scalar t228 = t2*t3*t4*t5*t29*t35*t36*4.0;
    const scalar t229 = t8*t10*2.0;
    const scalar t230 = t229-1.0;
    const scalar t231 = t4*t185;
    const scalar t232 = -t167+t231+2.0;
    const scalar t233 = t4*t29*t36*t57*t224;
    const scalar t234 = t2*t5*t29*t35*t36*t57*t224*2.0;
    const scalar t235 = t3*t20*t29*t35*t36*t57*t224*2.0;
    const scalar t236 = t4*t29*t36*t70*t224;
    const scalar t237 = t2*t5*t29*t35*t36*t70*t224*2.0;
    const scalar t238 = t3*t20*t29*t35*t36*t70*t224*2.0;
    const scalar t239 = t4*t29*t36*t83*t224;
    const scalar t240 = t4*t201;
    const scalar t241 = -t167+t240+1.0;
    const scalar t242 = t2*t5*t29*t35*t36*t83*t224*2.0;
    const scalar t243 = t3*t20*t29*t35*t36*t83*t224*2.0;
    const scalar t244 = t2*t10*t29*t174;
    const scalar t245 = t2*t10*t17*2.0;
    const scalar t246 = -a(0)+p(0)+t245;
    const scalar t251 = t2*t10*t29*t246;
    const scalar t247 = t19-t251;
    const scalar t248 = t3*t10*t20*t29*t35*t174*2.0;
    const scalar t249 = t4*t10*t21*t29*t35*t174*2.0;
    const scalar t250 = t3*t6*t29*t36*2.0;
    const scalar t252 = t2*t7*t20*t29*t35*t36*4.0;
    const scalar t253 = t4*t6*t29*t36*2.0;
    const scalar t254 = t2*t8*t21*t29*t35*t36*4.0;
    const scalar t255 = t2*t10*t58*2.0;
    const scalar t256 = t2*t17*t36*t57*2.0;
    const scalar t257 = -t167+t255+t256+1.0;
    const scalar t258 = t2*t10*t72*2.0;
    const scalar t259 = t2*t17*t36*t70*2.0;
    const scalar t260 = t258+t259;
    const scalar t261 = t2*t10*t85*2.0;
    const scalar t262 = t2*t17*t36*t83*2.0;
    const scalar t263 = t261+t262;
    const scalar t264 = t2*t29*t36*t57*t246;
    const scalar t265 = t2*t10*t11*2.0;
    const scalar t266 = t167-t256+t265;
    const scalar t267 = t2*t29*t36*t70*t246;
    const scalar t268 = t10*t20*t29*t35*t246*2.0;
    const scalar t269 = t2*t29*t36*t83*t246;
    const scalar t270 = t10*t21*t29*t35*t246*2.0;
    const scalar t271 = t2*t7*t29*t36*2.0;
    const scalar t272 = t3*t10*t17*2.0;
    const scalar t273 = -a(1)+p(1)+t272;
    const scalar t274 = t3*t5*t6*t29*t35*t36*4.0;
    const scalar t276 = t3*t10*t29*t273;
    const scalar t275 = t19-t276;
    const scalar t277 = t3*t10*t29*t209;
    const scalar t278 = t2*t5*t10*t29*t35*t209*2.0;
    const scalar t279 = t4*t10*t21*t29*t35*t209*2.0;
    const scalar t280 = t4*t7*t29*t36*2.0;
    const scalar t281 = t3*t8*t21*t29*t35*t36*4.0;
    const scalar t282 = t3*t10*t58*2.0;
    const scalar t283 = t3*t17*t36*t57*2.0;
    const scalar t284 = t282+t283;
    const scalar t285 = t3*t10*t72*2.0;
    const scalar t286 = t3*t17*t36*t70*2.0;
    const scalar t287 = -t167+t285+t286+1.0;
    const scalar t288 = t3*t10*t85*2.0;
    const scalar t289 = t3*t17*t36*t83*2.0;
    const scalar t290 = t288+t289;
    const scalar t291 = t32*t95;
    const scalar t292 = t3*t29*t36*t57*t273;
    const scalar t293 = t5*t10*t29*t35*t273*2.0;
    const scalar t294 = t3*t29*t36*t70*t273;
    const scalar t295 = t3*t10*t13*2.0;
    const scalar t296 = t167-t286+t295;
    const scalar t297 = t3*t29*t36*t83*t273;
    const scalar t298 = t10*t21*t29*t35*t273*2.0;
    const scalar t299 = t2*t8*t29*t36*2.0;
    const scalar t300 = t4*t10*t17*2.0;
    const scalar t301 = -a(2)+p(2)+t300;
    const scalar t302 = t4*t5*t6*t29*t35*t36*4.0;
    const scalar t305 = t4*t10*t29*t301;
    const scalar t303 = t19-t305;
    const scalar t304 = t3*t8*t29*t36*2.0;
    const scalar t306 = t4*t7*t20*t29*t35*t36*4.0;
    const scalar t307 = t4*t10*t29*t230;
    const scalar t308 = t2*t5*t10*t29*t35*t230*2.0;
    const scalar t309 = t3*t10*t20*t29*t35*t230*2.0;
    const scalar t310 = t4*t10*t58*2.0;
    const scalar t311 = t4*t17*t36*t57*2.0;
    const scalar t312 = t310+t311;
    const scalar t313 = t4*t10*t72*2.0;
    const scalar t314 = t4*t17*t36*t70*2.0;
    const scalar t315 = t313+t314;
    const scalar t316 = t4*t10*t85*2.0;
    const scalar t317 = t4*t17*t36*t83*2.0;
    const scalar t318 = -t167+t316+t317+1.0;
    const scalar t319 = t4*t29*t36*t57*t301;
    const scalar t320 = t5*t10*t29*t35*t301*2.0;
    const scalar t321 = t32*t106;
    const scalar t322 = t4*t29*t36*t70*t301;
    const scalar t323 = t10*t20*t29*t35*t301*2.0;
    const scalar t324 = t4*t29*t36*t83*t301;
    const scalar t325 = t4*t10*t15*2.0;
    const scalar t326 = t167-t317+t325;
    hess(0, 0) = t133+t154+t31*t34*t35*2.0+t5*t31*t38*t42*2.0+t2*t5*t6*t35*t36*t37*2.0+t3*t6*t20*t35*t36*t37*2.0+t4*t6*t21*t35*t36*t37*2.0+t2*t3*t10*t20*t29*t38*t42*2.0+t2*t4*t10*t21*t29*t38*t42*2.0;
    hess(0, 1) = t127+t128+t131+t132+t5*t31*t38*t48*2.0+t2*t3*t10*t31*t32*t35*2.0+t2*t3*t10*t29*t35*t44*2.0+t2*t3*t10*t20*t29*t38*t48*2.0+t2*t4*t10*t21*t29*t38*t48*2.0;
    hess(0, 2) = t144+t145+t148+t149+t5*t31*t38*t54*2.0+t2*t4*t10*t31*t32*t35*2.0+t2*t4*t10*t29*t35*t50*2.0+t2*t3*t10*t20*t29*t38*t54*2.0+t2*t4*t10*t21*t29*t38*t54*2.0;
    hess(0, 3) = t103+t104+t5*t35*(t94-t10*t29*t57+t6*t10*t37*t61)*2.0+t31*t35*t63*2.0+t5*t31*t38*t67*2.0-t3*t10*t20*t29*t35*2.0-t4*t10*t21*t29*t35*2.0+t2*t3*t10*t20*t35*t37*t61*2.0+t2*t3*t10*t20*t29*t38*t67*2.0+t2*t4*t10*t21*t35*t37*t61*2.0+t2*t4*t10*t21*t29*t38*t67*2.0+t2*t7*t10*t29*t32*t35*t61*2.0+t2*t8*t10*t29*t32*t35*t61*2.0;
    hess(0, 4) = t114+t115+t5*t35*(t105+t6*t10*t37*(t71+t10*(a(1)*-2.0+b(1)+p(1))))*2.0+t5*t31*t38*t80*2.0-t2*t10*t20*t29*t35*2.0+t2*t31*t32*t35*t74*2.0+t2*t3*t10*t29*t35*t76*2.0+t2*t3*t10*t20*t35*t37*t74*2.0+t2*t3*t10*t20*t29*t38*t80*2.0+t2*t4*t10*t21*t35*t37*t74*2.0+t2*t4*t10*t21*t29*t38*t80*2.0+t2*t8*t10*t29*t32*t35*t74*2.0;
    hess(0, 5) = t125+t126+t5*t35*(t116+t6*t10*t37*(t84+t10*(a(2)*-2.0+b(2)+p(2))))*2.0+t5*t31*t38*t93*2.0-t2*t10*t21*t29*t35*2.0+t2*t31*t32*t35*t87*2.0+t2*t4*t10*t29*t35*t89*2.0+t2*t3*t10*t20*t35*t37*t87*2.0+t2*t3*t10*t20*t29*t38*t93*2.0+t2*t4*t10*t21*t35*t37*t87*2.0+t2*t4*t10*t21*t29*t38*t93*2.0+t2*t7*t10*t29*t32*t35*t87*2.0;
    hess(0, 6) = -t103-t104+t134+t157-t5*t35*(t94-t10*t29*t57+t6*t10*t37*t95)*2.0+t31*t35*t97*2.0-t5*t31*t38*t101*2.0-t2*t3*t10*t20*t35*t37*t95*2.0-t2*t3*t10*t20*t29*t38*t101*2.0-t2*t4*t10*t21*t35*t37*t95*2.0-t2*t4*t10*t21*t29*t38*t101*2.0-t2*t7*t10*t29*t32*t35*t95*2.0-t2*t8*t10*t29*t32*t35*t95*2.0;
    hess(0, 7) = -t114-t115-t5*t35*(t105+t6*t10*t37*t106)*2.0-t5*t31*t38*t112*2.0+t2*t10*t20*t29*t35*2.0-t2*t31*t32*t35*t106*2.0+t2*t3*t10*t29*t35*t108*2.0-t2*t3*t10*t20*t35*t37*t106*2.0-t2*t3*t10*t20*t29*t38*t112*2.0-t2*t4*t10*t21*t35*t37*t106*2.0-t2*t4*t10*t21*t29*t38*t112*2.0-t2*t8*t10*t29*t32*t35*t106*2.0;
    hess(0, 8) = -t125-t126-t5*t35*(t116+t6*t10*t37*t117)*2.0-t5*t31*t38*t123*2.0+t2*t10*t21*t29*t35*2.0-t2*t31*t32*t35*t117*2.0+t2*t4*t10*t29*t35*t119*2.0-t2*t3*t10*t20*t35*t37*t117*2.0-t2*t3*t10*t20*t29*t38*t123*2.0-t2*t4*t10*t21*t35*t37*t117*2.0-t2*t4*t10*t21*t29*t38*t123*2.0-t2*t7*t10*t29*t32*t35*t117*2.0;
    hess(1, 0) = t127+t128+t131+t132+t20*t38*t42*t130*2.0+t2*t3*t10*t29*t34*t35*2.0+t2*t3*t10*t32*t35*t130*2.0+t2*t3*t5*t10*t29*t38*t42*2.0+t3*t4*t10*t21*t29*t38*t42*2.0;
    hess(1, 1) = t133+t155+t35*t44*t130*2.0+t20*t38*t48*t130*2.0+t2*t5*t7*t35*t36*t37*2.0+t3*t7*t20*t35*t36*t37*2.0+t4*t7*t21*t35*t36*t37*2.0+t2*t3*t5*t10*t29*t38*t48*2.0+t3*t4*t10*t21*t29*t38*t48*2.0;
    hess(1, 2) = t150+t151+t152+t153+t20*t38*t54*t130*2.0+t3*t4*t10*t29*t35*t50*2.0+t3*t4*t10*t32*t35*t130*2.0+t2*t3*t5*t10*t29*t38*t54*2.0+t3*t4*t10*t21*t29*t38*t54*2.0;
    hess(1, 3) = t136+t137+t20*t35*(t135+t7*t10*t37*t61)*2.0+t20*t38*t67*t130*2.0-t3*t5*t10*t29*t35*2.0+t3*t32*t35*t61*t130*2.0+t2*t3*t10*t29*t35*t63*2.0+t2*t3*t5*t10*t35*t37*t61*2.0+t2*t3*t5*t10*t29*t38*t67*2.0+t3*t4*t10*t21*t35*t37*t61*2.0+t3*t4*t10*t21*t29*t38*t67*2.0+t3*t8*t10*t29*t32*t35*t61*2.0;
    hess(1, 4) = -t134+t139+t140+t20*t35*(t138-t10*t29*t70+t7*t10*t37*t74)*2.0+t35*t76*t130*2.0+t20*t38*t80*t130*2.0-t2*t5*t10*t29*t35*2.0+t2*t3*t5*t10*t35*t37*t74*2.0+t2*t3*t5*t10*t29*t38*t80*2.0+t3*t4*t10*t21*t35*t37*t74*2.0+t3*t4*t10*t21*t29*t38*t80*2.0+t3*t6*t10*t29*t32*t35*t74*2.0+t3*t8*t10*t29*t32*t35*t74*2.0;
    hess(1, 5) = t142+t143+t20*t35*(t141+t7*t10*t37*t87)*2.0+t20*t38*t93*t130*2.0-t3*t10*t21*t29*t35*2.0+t3*t32*t35*t87*t130*2.0+t3*t4*t10*t29*t35*t89*2.0+t2*t3*t5*t10*t35*t37*t87*2.0+t2*t3*t5*t10*t29*t38*t93*2.0+t3*t4*t10*t21*t35*t37*t87*2.0+t3*t4*t10*t21*t29*t38*t93*2.0+t3*t6*t10*t29*t32*t35*t87*2.0;
    hess(1, 6) = -t136-t137-t20*t35*(t135+t7*t10*t37*t95)*2.0-t20*t38*t101*t130*2.0+t3*t5*t10*t29*t35*2.0-t3*t32*t35*t95*t130*2.0+t2*t3*t10*t29*t35*t97*2.0-t2*t3*t5*t10*t35*t37*t95*2.0-t2*t3*t5*t10*t29*t38*t101*2.0-t3*t4*t10*t21*t35*t37*t95*2.0-t3*t4*t10*t21*t29*t38*t101*2.0-t3*t8*t10*t29*t32*t35*t95*2.0;
    hess(1, 7) = t134-t139-t140+t156-t20*t35*(t138-t10*t29*t70+t7*t10*t37*t106)*2.0+t35*t108*t130*2.0-t20*t38*t112*t130*2.0-t2*t3*t5*t10*t35*t37*t106*2.0-t2*t3*t5*t10*t29*t38*t112*2.0-t3*t4*t10*t21*t35*t37*t106*2.0-t3*t4*t10*t21*t29*t38*t112*2.0-t3*t6*t10*t29*t32*t35*t106*2.0-t3*t8*t10*t29*t32*t35*t106*2.0;
    hess(1, 8) = -t142-t143-t20*t35*(t141+t7*t10*t37*t117)*2.0-t20*t38*t123*t130*2.0+t3*t10*t21*t29*t35*2.0-t3*t32*t35*t117*t130*2.0+t3*t4*t10*t29*t35*t119*2.0-t2*t3*t5*t10*t35*t37*t117*2.0-t2*t3*t5*t10*t29*t38*t123*2.0-t3*t4*t10*t21*t35*t37*t117*2.0-t3*t4*t10*t21*t29*t38*t123*2.0-t3*t6*t10*t29*t32*t35*t117*2.0;
    hess(2, 0) = t144+t145+t148+t149+t21*t38*t42*t147*2.0+t2*t4*t10*t29*t34*t35*2.0+t2*t4*t10*t32*t35*t147*2.0+t2*t4*t5*t10*t29*t38*t42*2.0+t3*t4*t10*t20*t29*t38*t42*2.0;
    hess(2, 1) = t150+t151+t152+t153+t21*t38*t48*t147*2.0+t3*t4*t10*t29*t35*t44*2.0+t3*t4*t10*t32*t35*t147*2.0+t2*t4*t5*t10*t29*t38*t48*2.0+t3*t4*t10*t20*t29*t38*t48*2.0;
    hess(2, 2) = t154+t155+t35*t50*t147*2.0+t21*t38*t54*t147*2.0+t2*t5*t8*t35*t36*t37*2.0+t3*t8*t20*t35*t36*t37*2.0+t4*t8*t21*t35*t36*t37*2.0+t2*t4*t5*t10*t29*t38*t54*2.0+t3*t4*t10*t20*t29*t38*t54*2.0;
    hess(2, 3) = t159+t160+t21*t35*(t158+t8*t10*t37*t61)*2.0+t21*t38*t67*t147*2.0-t4*t5*t10*t29*t35*2.0+t4*t32*t35*t61*t147*2.0+t2*t4*t10*t29*t35*t63*2.0+t2*t4*t5*t10*t35*t37*t61*2.0+t2*t4*t5*t10*t29*t38*t67*2.0+t3*t4*t10*t20*t35*t37*t61*2.0+t3*t4*t10*t20*t29*t38*t67*2.0+t4*t7*t10*t29*t32*t35*t61*2.0;
    hess(2, 4) = t162+t163+t21*t35*(t161+t8*t10*t37*t74)*2.0+t21*t38*t80*t147*2.0-t4*t10*t20*t29*t35*2.0+t4*t32*t35*t74*t147*2.0+t3*t4*t10*t29*t35*t76*2.0+t2*t4*t5*t10*t35*t37*t74*2.0+t2*t4*t5*t10*t29*t38*t80*2.0+t3*t4*t10*t20*t35*t37*t74*2.0+t3*t4*t10*t20*t29*t38*t80*2.0+t4*t6*t10*t29*t32*t35*t74*2.0;
    hess(2, 5) = -t156-t157+t165+t166+t21*t35*(t164-t10*t29*t83+t8*t10*t37*t87)*2.0+t35*t89*t147*2.0+t21*t38*t93*t147*2.0+t2*t4*t5*t10*t35*t37*t87*2.0+t2*t4*t5*t10*t29*t38*t93*2.0+t3*t4*t10*t20*t35*t37*t87*2.0+t3*t4*t10*t20*t29*t38*t93*2.0+t4*t6*t10*t29*t32*t35*t87*2.0+t4*t7*t10*t29*t32*t35*t87*2.0;
    hess(2, 6) = -t159-t160-t21*t35*(t158+t8*t10*t37*t95)*2.0-t21*t38*t101*t147*2.0+t4*t5*t10*t29*t35*2.0-t4*t32*t35*t95*t147*2.0+t2*t4*t10*t29*t35*t97*2.0-t2*t4*t5*t10*t35*t37*t95*2.0-t2*t4*t5*t10*t29*t38*t101*2.0-t3*t4*t10*t20*t35*t37*t95*2.0-t3*t4*t10*t20*t29*t38*t101*2.0-t4*t7*t10*t29*t32*t35*t95*2.0;
    hess(2, 7) = -t162-t163-t21*t35*(t161+t8*t10*t37*t106)*2.0-t21*t38*t112*t147*2.0+t4*t10*t20*t29*t35*2.0-t4*t32*t35*t106*t147*2.0+t3*t4*t10*t29*t35*t108*2.0-t2*t4*t5*t10*t35*t37*t106*2.0-t2*t4*t5*t10*t29*t38*t112*2.0-t3*t4*t10*t20*t35*t37*t106*2.0-t3*t4*t10*t20*t29*t38*t112*2.0-t4*t6*t10*t29*t32*t35*t106*2.0;
    hess(2, 8) = t156+t157-t165-t166-t21*t35*(t164-t10*t29*t83+t8*t10*t37*t117)*2.0+t35*t119*t147*2.0-t21*t38*t123*t147*2.0-t2*t4*t5*t10*t35*t37*t117*2.0-t2*t4*t5*t10*t29*t38*t123*2.0-t3*t4*t10*t20*t35*t37*t117*2.0-t3*t4*t10*t20*t29*t38*t123*2.0-t4*t6*t10*t29*t32*t35*t117*2.0-t4*t7*t10*t29*t32*t35*t117*2.0;
    hess(3, 0) = t248+t249+t5*t35*(t244-t2*t10*t32+t6*t36*t37*t170)*2.0+t34*t35*t172*2.0+t5*t38*t42*t172*2.0+t2*t3*t20*t35*t36*t37*t170*2.0+t2*t4*t21*t35*t36*t37*t170*2.0+t2*t7*t29*t32*t35*t36*t170*2.0+t2*t8*t29*t32*t35*t36*t170*2.0+t3*t10*t20*t29*t38*t42*t170*2.0+t4*t10*t21*t29*t38*t42*t170*2.0;
    hess(3, 1) = t207+t252+t5*t35*(t250-t3*t10*t32+t2*t3*t36*t37*t170)*2.0+t5*t38*t48*t172*2.0+t2*t3*t10*t32*t35*t172*2.0+t3*t10*t29*t35*t44*t170*2.0+t7*t20*t35*t36*t37*t170*2.0+t3*t4*t21*t35*t36*t37*t170*2.0+t3*t8*t29*t32*t35*t36*t170*2.0+t3*t10*t20*t29*t38*t48*t170*2.0+t4*t10*t21*t29*t38*t48*t170*2.0;
    hess(3, 2) = t227+t254+t5*t35*(t253-t4*t10*t32+t2*t4*t36*t37*t170)*2.0+t5*t38*t54*t172*2.0+t2*t4*t10*t32*t35*t172*2.0+t4*t10*t29*t35*t50*t170*2.0+t8*t21*t35*t36*t37*t170*2.0+t3*t4*t20*t35*t36*t37*t170*2.0+t4*t7*t29*t32*t35*t36*t170*2.0+t3*t10*t20*t29*t38*t54*t170*2.0+t4*t10*t21*t29*t38*t54*t170*2.0;
    hess(3, 3) = t191+t192+t35*t63*t172*2.0+t5*t35*(t186-t32*t61-t10*t29*t170+t2*t10*t29*t179+t2*t10*t37*t61*t170)*2.0+t5*t38*t67*t172*2.0+t3*t10*t20*t29*t35*t179*2.0+t4*t10*t21*t29*t35*t179*2.0+t3*t10*t20*t35*t37*t61*t170*2.0+t3*t10*t20*t29*t38*t67*t170*2.0+t4*t10*t21*t35*t37*t61*t170*2.0+t4*t10*t21*t29*t38*t67*t170*2.0+t7*t10*t29*t32*t35*t61*t170*2.0+t8*t10*t29*t32*t35*t61*t170*2.0;
    hess(3, 4) = t194+t195+t5*t35*(t193-t32*t74+t6*t10*t29*t182+t2*t10*t37*t74*t170)*2.0+t5*t38*t80*t172*2.0-t10*t20*t29*t35*t170*2.0+t2*t32*t35*t74*t172*2.0+t3*t10*t29*t35*t76*t170*2.0+t2*t3*t10*t20*t29*t35*t182*2.0+t2*t4*t10*t21*t29*t35*t182*2.0+t3*t10*t20*t35*t37*t74*t170*2.0+t3*t10*t20*t29*t38*t80*t170*2.0+t4*t10*t21*t35*t37*t74*t170*2.0+t4*t10*t21*t29*t38*t80*t170*2.0+t8*t10*t29*t32*t35*t74*t170*2.0;
    hess(3, 5) = t199+t200+t5*t35*(t198-t32*t87+t6*t10*t29*t185+t2*t10*t37*t87*t170)*2.0+t5*t38*t93*t172*2.0-t10*t21*t29*t35*t170*2.0+t2*t32*t35*t87*t172*2.0+t4*t10*t29*t35*t89*t170*2.0+t2*t3*t10*t20*t29*t35*t185*2.0+t2*t4*t10*t21*t29*t35*t185*2.0+t3*t10*t20*t35*t37*t87*t170*2.0+t3*t10*t20*t29*t38*t93*t170*2.0+t4*t10*t21*t35*t37*t87*t170*2.0+t4*t10*t21*t29*t38*t93*t170*2.0+t7*t10*t29*t32*t35*t87*t170*2.0;
    hess(3, 6) = -t191-t192+t35*t97*t172*2.0-t5*t35*(t186-t32*t95-t10*t29*t170+t2*t10*t29*t189+t2*t10*t37*t95*t170)*2.0-t5*t38*t101*t172*2.0-t3*t10*t20*t29*t35*t189*2.0-t4*t10*t21*t29*t35*t189*2.0-t3*t10*t20*t35*t37*t95*t170*2.0-t3*t10*t20*t29*t38*t101*t170*2.0-t4*t10*t21*t35*t37*t95*t170*2.0-t4*t10*t21*t29*t38*t101*t170*2.0-t7*t10*t29*t32*t35*t95*t170*2.0-t8*t10*t29*t32*t35*t95*t170*2.0;
    hess(3, 7) = -t194-t195-t5*t35*(t193-t32*t106+t6*t10*t29*t196+t2*t10*t37*t106*t170)*2.0-t5*t38*t112*t172*2.0+t10*t20*t29*t35*t170*2.0-t2*t32*t35*t106*t172*2.0+t3*t10*t29*t35*t108*t170*2.0-t2*t3*t10*t20*t29*t35*t196*2.0-t2*t4*t10*t21*t29*t35*t196*2.0-t3*t10*t20*t35*t37*t106*t170*2.0-t3*t10*t20*t29*t38*t112*t170*2.0-t4*t10*t21*t35*t37*t106*t170*2.0-t4*t10*t21*t29*t38*t112*t170*2.0-t8*t10*t29*t32*t35*t106*t170*2.0;
    hess(3, 8) = -t199-t200-t5*t35*(t198-t32*t117+t6*t10*t29*t201+t2*t10*t37*t117*t170)*2.0-t5*t38*t123*t172*2.0+t10*t21*t29*t35*t170*2.0-t2*t32*t35*t117*t172*2.0+t4*t10*t29*t35*t119*t170*2.0-t2*t3*t10*t20*t29*t35*t201*2.0-t2*t4*t10*t21*t29*t35*t201*2.0-t3*t10*t20*t35*t37*t117*t170*2.0-t3*t10*t20*t29*t38*t123*t170*2.0-t4*t10*t21*t35*t37*t117*t170*2.0-t4*t10*t21*t29*t38*t123*t170*2.0-t7*t10*t29*t32*t35*t117*t170*2.0;
    hess(4, 0) = t207+t274+t20*t35*(t271-t2*t10*t32+t2*t3*t36*t37*t204)*2.0+t20*t38*t42*t206*2.0+t2*t3*t10*t32*t35*t206*2.0+t2*t10*t29*t34*t35*t204*2.0+t5*t6*t35*t36*t37*t204*2.0+t2*t5*t10*t29*t38*t42*t204*2.0+t2*t4*t21*t35*t36*t37*t204*2.0+t2*t8*t29*t32*t35*t36*t204*2.0+t4*t10*t21*t29*t38*t42*t204*2.0;
    hess(4, 1) = t278+t279+t20*t35*(t277-t3*t10*t32+t7*t36*t37*t204)*2.0+t35*t44*t206*2.0+t20*t38*t48*t206*2.0+t2*t3*t5*t35*t36*t37*t204*2.0+t2*t5*t10*t29*t38*t48*t204*2.0+t3*t4*t21*t35*t36*t37*t204*2.0+t3*t6*t29*t32*t35*t36*t204*2.0+t3*t8*t29*t32*t35*t36*t204*2.0+t4*t10*t21*t29*t38*t48*t204*2.0;
    hess(4, 2) = t228+t281+t20*t35*(t280-t4*t10*t32+t3*t4*t36*t37*t204)*2.0+t20*t38*t54*t206*2.0+t3*t4*t10*t32*t35*t206*2.0+t4*t10*t29*t35*t50*t204*2.0+t8*t21*t35*t36*t37*t204*2.0+t2*t4*t5*t35*t36*t37*t204*2.0+t2*t5*t10*t29*t38*t54*t204*2.0+t4*t6*t29*t32*t35*t36*t204*2.0+t4*t10*t21*t29*t38*t54*t204*2.0;
    hess(4, 3) = t213+t214+t20*t35*(t212-t32*t61+t7*t10*t29*t177+t3*t10*t37*t61*t204)*2.0+t20*t38*t67*t206*2.0-t5*t10*t29*t35*t204*2.0+t3*t32*t35*t61*t206*2.0+t2*t10*t29*t35*t63*t204*2.0+t2*t3*t5*t10*t29*t35*t177*2.0+t3*t4*t10*t21*t29*t35*t177*2.0+t2*t5*t10*t35*t37*t61*t204*2.0+t2*t5*t10*t29*t38*t67*t204*2.0+t4*t10*t21*t35*t37*t61*t204*2.0+t4*t10*t21*t29*t38*t67*t204*2.0+t8*t10*t29*t32*t35*t61*t204*2.0;
    hess(4, 4) = t218+t219+t35*t76*t206*2.0+t20*t35*(t215-t32*t74-t10*t29*t204+t3*t10*t29*t211+t3*t10*t37*t74*t204)*2.0+t20*t38*t80*t206*2.0+t2*t5*t10*t29*t35*t211*2.0+t4*t10*t21*t29*t35*t211*2.0+t2*t5*t10*t35*t37*t74*t204*2.0+t2*t5*t10*t29*t38*t80*t204*2.0+t4*t10*t21*t35*t37*t74*t204*2.0+t4*t10*t21*t29*t38*t80*t204*2.0+t6*t10*t29*t32*t35*t74*t204*2.0+t8*t10*t29*t32*t35*t74*t204*2.0;
    hess(4, 5) = t221+t222+t20*t35*(t220-t32*t87+t7*t10*t29*t185+t3*t10*t37*t87*t204)*2.0+t20*t38*t93*t206*2.0-t10*t21*t29*t35*t204*2.0+t3*t32*t35*t87*t206*2.0+t4*t10*t29*t35*t89*t204*2.0+t2*t3*t5*t10*t29*t35*t185*2.0+t3*t4*t10*t21*t29*t35*t185*2.0+t2*t5*t10*t35*t37*t87*t204*2.0+t2*t5*t10*t29*t38*t93*t204*2.0+t4*t10*t21*t35*t37*t87*t204*2.0+t4*t10*t21*t29*t38*t93*t204*2.0+t6*t10*t29*t32*t35*t87*t204*2.0;
    hess(4, 6) = -t213-t214-t20*t35*(t212-t32*t95+t7*t10*t29*t187+t3*t10*t37*t95*t204)*2.0-t20*t38*t101*t206*2.0+t5*t10*t29*t35*t204*2.0-t3*t32*t35*t95*t206*2.0+t2*t10*t29*t35*t97*t204*2.0-t2*t3*t5*t10*t29*t35*t187*2.0-t3*t4*t10*t21*t29*t35*t187*2.0-t2*t5*t10*t35*t37*t95*t204*2.0-t2*t5*t10*t29*t38*t101*t204*2.0-t4*t10*t21*t35*t37*t95*t204*2.0-t4*t10*t21*t29*t38*t101*t204*2.0-t8*t10*t29*t32*t35*t95*t204*2.0;
    hess(4, 7) = -t218-t219+t35*t108*t206*2.0-t20*t35*(t215-t32*t106-t10*t29*t204+t3*t10*t29*t217+t3*t10*t37*t106*t204)*2.0-t20*t38*t112*t206*2.0-t2*t5*t10*t29*t35*t217*2.0-t4*t10*t21*t29*t35*t217*2.0-t2*t5*t10*t35*t37*t106*t204*2.0-t2*t5*t10*t29*t38*t112*t204*2.0-t4*t10*t21*t35*t37*t106*t204*2.0-t4*t10*t21*t29*t38*t112*t204*2.0-t6*t10*t29*t32*t35*t106*t204*2.0-t8*t10*t29*t32*t35*t106*t204*2.0;
    hess(4, 8) = -t221-t222-t20*t35*(t220-t32*t117+t7*t10*t29*t201+t3*t10*t37*t117*t204)*2.0-t20*t38*t123*t206*2.0+t10*t21*t29*t35*t204*2.0-t3*t32*t35*t117*t206*2.0+t4*t10*t29*t35*t119*t204*2.0-t2*t3*t5*t10*t29*t35*t201*2.0-t3*t4*t10*t21*t29*t35*t201*2.0-t2*t5*t10*t35*t37*t117*t204*2.0-t2*t5*t10*t29*t38*t123*t204*2.0-t4*t10*t21*t35*t37*t117*t204*2.0-t4*t10*t21*t29*t38*t123*t204*2.0-t6*t10*t29*t32*t35*t117*t204*2.0;
    hess(5, 0) = t227+t302+t21*t35*(t299-t2*t10*t32+t2*t4*t36*t37*t224)*2.0+t21*t38*t42*t226*2.0+t2*t4*t10*t32*t35*t226*2.0+t2*t10*t29*t34*t35*t224*2.0+t5*t6*t35*t36*t37*t224*2.0+t2*t5*t10*t29*t38*t42*t224*2.0+t2*t3*t20*t35*t36*t37*t224*2.0+t2*t7*t29*t32*t35*t36*t224*2.0+t3*t10*t20*t29*t38*t42*t224*2.0;
    hess(5, 1) = t228+t306+t21*t35*(t304-t3*t10*t32+t3*t4*t36*t37*t224)*2.0+t21*t38*t48*t226*2.0+t3*t4*t10*t32*t35*t226*2.0+t3*t10*t29*t35*t44*t224*2.0+t7*t20*t35*t36*t37*t224*2.0+t2*t3*t5*t35*t36*t37*t224*2.0+t2*t5*t10*t29*t38*t48*t224*2.0+t3*t6*t29*t32*t35*t36*t224*2.0+t3*t10*t20*t29*t38*t48*t224*2.0;
    hess(5, 2) = t308+t309+t21*t35*(t307-t4*t10*t32+t8*t36*t37*t224)*2.0+t35*t50*t226*2.0+t21*t38*t54*t226*2.0+t2*t4*t5*t35*t36*t37*t224*2.0+t3*t4*t20*t35*t36*t37*t224*2.0+t2*t5*t10*t29*t38*t54*t224*2.0+t4*t6*t29*t32*t35*t36*t224*2.0+t4*t7*t29*t32*t35*t36*t224*2.0+t3*t10*t20*t29*t38*t54*t224*2.0;
    hess(5, 3) = t234+t235+t21*t35*(t233-t32*t61+t8*t10*t29*t177+t4*t10*t37*t61*t224)*2.0+t21*t38*t67*t226*2.0-t5*t10*t29*t35*t224*2.0+t4*t32*t35*t61*t226*2.0+t2*t10*t29*t35*t63*t224*2.0+t2*t4*t5*t10*t29*t35*t177*2.0+t3*t4*t10*t20*t29*t35*t177*2.0+t2*t5*t10*t35*t37*t61*t224*2.0+t2*t5*t10*t29*t38*t67*t224*2.0+t3*t10*t20*t35*t37*t61*t224*2.0+t3*t10*t20*t29*t38*t67*t224*2.0+t7*t10*t29*t32*t35*t61*t224*2.0;
    hess(5, 4) = t237+t238+t21*t35*(t236-t32*t74+t8*t10*t29*t182+t4*t10*t37*t74*t224)*2.0+t21*t38*t80*t226*2.0-t10*t20*t29*t35*t224*2.0+t4*t32*t35*t74*t226*2.0+t3*t10*t29*t35*t76*t224*2.0+t2*t4*t5*t10*t29*t35*t182*2.0+t3*t4*t10*t20*t29*t35*t182*2.0+t2*t5*t10*t35*t37*t74*t224*2.0+t2*t5*t10*t29*t38*t80*t224*2.0+t3*t10*t20*t35*t37*t74*t224*2.0+t3*t10*t20*t29*t38*t80*t224*2.0+t6*t10*t29*t32*t35*t74*t224*2.0;
    hess(5, 5) = t242+t243+t35*t89*t226*2.0+t21*t35*(t239-t32*t87-t10*t29*t224+t4*t10*t29*t232+t4*t10*t37*t87*t224)*2.0+t21*t38*t93*t226*2.0+t2*t5*t10*t29*t35*t232*2.0+t3*t10*t20*t29*t35*t232*2.0+t2*t5*t10*t35*t37*t87*t224*2.0+t2*t5*t10*t29*t38*t93*t224*2.0+t3*t10*t20*t35*t37*t87*t224*2.0+t3*t10*t20*t29*t38*t93*t224*2.0+t6*t10*t29*t32*t35*t87*t224*2.0+t7*t10*t29*t32*t35*t87*t224*2.0;
    hess(5, 6) = -t234-t235-t21*t35*(t233-t32*t95+t8*t10*t29*t187+t4*t10*t37*t95*t224)*2.0-t21*t38*t101*t226*2.0+t5*t10*t29*t35*t224*2.0-t4*t32*t35*t95*t226*2.0+t2*t10*t29*t35*t97*t224*2.0-t2*t4*t5*t10*t29*t35*t187*2.0-t3*t4*t10*t20*t29*t35*t187*2.0-t2*t5*t10*t35*t37*t95*t224*2.0-t2*t5*t10*t29*t38*t101*t224*2.0-t3*t10*t20*t35*t37*t95*t224*2.0-t3*t10*t20*t29*t38*t101*t224*2.0-t7*t10*t29*t32*t35*t95*t224*2.0;
    hess(5, 7) = -t237-t238-t21*t35*(t236-t32*t106+t8*t10*t29*t196+t4*t10*t37*t106*t224)*2.0-t21*t38*t112*t226*2.0+t10*t20*t29*t35*t224*2.0-t4*t32*t35*t106*t226*2.0+t3*t10*t29*t35*t108*t224*2.0-t2*t4*t5*t10*t29*t35*t196*2.0-t3*t4*t10*t20*t29*t35*t196*2.0-t2*t5*t10*t35*t37*t106*t224*2.0-t2*t5*t10*t29*t38*t112*t224*2.0-t3*t10*t20*t35*t37*t106*t224*2.0-t3*t10*t20*t29*t38*t112*t224*2.0-t6*t10*t29*t32*t35*t106*t224*2.0;
    hess(5, 8) = -t242-t243+t35*t119*t226*2.0-t21*t35*(t239-t32*t117-t10*t29*t224+t4*t10*t29*t241+t4*t10*t37*t117*t224)*2.0-t21*t38*t123*t226*2.0-t2*t5*t10*t29*t35*t241*2.0-t3*t10*t20*t29*t35*t241*2.0-t2*t5*t10*t35*t37*t117*t224*2.0-t2*t5*t10*t29*t38*t123*t224*2.0-t3*t10*t20*t35*t37*t117*t224*2.0-t3*t10*t20*t29*t38*t123*t224*2.0-t6*t10*t29*t32*t35*t117*t224*2.0-t7*t10*t29*t32*t35*t117*t224*2.0;
    hess(6, 0) = -t248-t249-t5*t35*(t244-t2*t10*t32+t6*t36*t37*t246)*2.0+t34*t35*t247*2.0+t5*t38*t42*t247*2.0-t2*t3*t20*t35*t36*t37*t246*2.0-t2*t4*t21*t35*t36*t37*t246*2.0-t2*t7*t29*t32*t35*t36*t246*2.0-t2*t8*t29*t32*t35*t36*t246*2.0-t3*t10*t20*t29*t38*t42*t246*2.0-t4*t10*t21*t29*t38*t42*t246*2.0;
    hess(6, 1) = -t207-t252-t5*t35*(t250-t3*t10*t32+t2*t3*t36*t37*t246)*2.0+t5*t38*t48*t247*2.0+t2*t3*t10*t32*t35*t247*2.0-t3*t10*t29*t35*t44*t246*2.0-t7*t20*t35*t36*t37*t246*2.0-t3*t4*t21*t35*t36*t37*t246*2.0-t3*t8*t29*t32*t35*t36*t246*2.0-t3*t10*t20*t29*t38*t48*t246*2.0-t4*t10*t21*t29*t38*t48*t246*2.0;
    hess(6, 2) = -t227-t254-t5*t35*(t253-t4*t10*t32+t2*t4*t36*t37*t246)*2.0+t5*t38*t54*t247*2.0+t2*t4*t10*t32*t35*t247*2.0-t4*t10*t29*t35*t50*t246*2.0-t8*t21*t35*t36*t37*t246*2.0-t3*t4*t20*t35*t36*t37*t246*2.0-t4*t7*t29*t32*t35*t36*t246*2.0-t3*t10*t20*t29*t38*t54*t246*2.0-t4*t10*t21*t29*t38*t54*t246*2.0;
    hess(6, 3) = t35*t63*t247*2.0-t5*t35*(t264-t32*t61-t10*t29*t246+t2*t10*t29*t257+t2*t10*t37*t61*t246)*2.0+t5*t38*t67*t247*2.0-t3*t10*t20*t29*t35*t257*2.0-t4*t10*t21*t29*t35*t257*2.0-t3*t10*t20*t35*t37*t61*t246*2.0-t3*t10*t20*t29*t38*t67*t246*2.0-t4*t10*t21*t35*t37*t61*t246*2.0-t4*t10*t21*t29*t38*t67*t246*2.0-t7*t10*t29*t32*t35*t61*t246*2.0-t8*t10*t29*t32*t35*t61*t246*2.0-t3*t20*t29*t35*t36*t57*t246*2.0-t4*t21*t29*t35*t36*t57*t246*2.0;
    hess(6, 4) = t268-t5*t35*(t267-t32*t74+t2*t10*t29*t260+t2*t10*t37*t74*t246)*2.0+t5*t38*t80*t247*2.0+t2*t32*t35*t74*t247*2.0-t3*t10*t20*t29*t35*t260*2.0-t4*t10*t21*t29*t35*t260*2.0-t3*t10*t29*t35*t76*t246*2.0-t3*t10*t20*t35*t37*t74*t246*2.0-t3*t10*t20*t29*t38*t80*t246*2.0-t4*t10*t21*t35*t37*t74*t246*2.0-t4*t10*t21*t29*t38*t80*t246*2.0-t8*t10*t29*t32*t35*t74*t246*2.0-t3*t20*t29*t35*t36*t70*t246*2.0-t4*t21*t29*t35*t36*t70*t246*2.0;
    hess(6, 5) = t270-t5*t35*(t269-t32*t87+t2*t10*t29*t263+t2*t10*t37*t87*t246)*2.0+t5*t38*t93*t247*2.0+t2*t32*t35*t87*t247*2.0-t3*t10*t20*t29*t35*t263*2.0-t4*t10*t21*t29*t35*t263*2.0-t4*t10*t29*t35*t89*t246*2.0-t3*t10*t20*t35*t37*t87*t246*2.0-t3*t10*t20*t29*t38*t93*t246*2.0-t4*t10*t21*t35*t37*t87*t246*2.0-t4*t10*t21*t29*t38*t93*t246*2.0-t7*t10*t29*t32*t35*t87*t246*2.0-t3*t20*t29*t35*t36*t83*t246*2.0-t4*t21*t29*t35*t36*t83*t246*2.0;
    hess(6, 6) = t5*t35*(-t264+t291+t10*t29*t246+t2*t10*t29*t266-t2*t10*t37*t95*t246)*-2.0+t35*t97*t247*2.0-t5*t38*t101*t247*2.0-t3*t10*t20*t29*t35*t266*2.0-t4*t10*t21*t29*t35*t266*2.0+t3*t10*t20*t29*t38*t246*(t99+t100-t102)*2.0+t4*t10*t21*t29*t38*t246*(t99+t100-t102)*2.0+t3*t20*t29*t35*t36*t57*t246*2.0+t4*t21*t29*t35*t36*t57*t246*2.0+t3*t10*t20*t35*t37*t246*(t60-t96)*2.0+t4*t10*t21*t35*t37*t246*(t60-t96)*2.0+t7*t10*t29*t32*t35*t246*(t60-t96)*2.0+t8*t10*t29*t32*t35*t246*(t60-t96)*2.0;
    hess(6, 7) = -t268+t5*t35*(t267-t32*t106+t2*t10*t29*(t259-t2*t10*t13*2.0)+t2*t10*t37*t246*(t71-t107))*2.0-t5*t38*t112*t247*2.0-t2*t32*t35*t106*t247*2.0-t3*t10*t29*t35*t108*t246*2.0+t3*t10*t20*t29*t35*(t259-t2*t10*t13*2.0)*2.0+t4*t10*t21*t29*t35*(t259-t2*t10*t13*2.0)*2.0+t3*t10*t20*t29*t38*t246*(t110+t111-t113)*2.0+t4*t10*t21*t29*t38*t246*(t110+t111-t113)*2.0+t3*t20*t29*t35*t36*t70*t246*2.0+t4*t21*t29*t35*t36*t70*t246*2.0+t3*t10*t20*t35*t37*t246*(t71-t107)*2.0+t4*t10*t21*t35*t37*t246*(t71-t107)*2.0+t8*t10*t29*t32*t35*t246*(t71-t107)*2.0;
    hess(6, 8) = -t270+t5*t35*(t269-t32*t117+t2*t10*t29*(t262-t2*t10*t15*2.0)+t2*t10*t37*t246*(t84-t118))*2.0-t5*t38*t123*t247*2.0-t2*t32*t35*t117*t247*2.0-t4*t10*t29*t35*t119*t246*2.0+t3*t10*t20*t29*t35*(t262-t2*t10*t15*2.0)*2.0+t4*t10*t21*t29*t35*(t262-t2*t10*t15*2.0)*2.0+t3*t10*t20*t29*t38*t246*(t121+t122-t124)*2.0+t4*t10*t21*t29*t38*t246*(t121+t122-t124)*2.0+t3*t20*t29*t35*t36*t83*t246*2.0+t4*t21*t29*t35*t36*t83*t246*2.0+t3*t10*t20*t35*t37*t246*(t84-t118)*2.0+t4*t10*t21*t35*t37*t246*(t84-t118)*2.0+t7*t10*t29*t32*t35*t246*(t84-t118)*2.0;
    hess(7, 0) = -t207-t274-t20*t35*(t271-t2*t10*t32+t2*t3*t36*t37*t273)*2.0+t20*t38*t42*t275*2.0+t2*t3*t10*t32*t35*t275*2.0-t2*t10*t29*t34*t35*t273*2.0-t5*t6*t35*t36*t37*t273*2.0-t2*t5*t10*t29*t38*t42*t273*2.0-t2*t4*t21*t35*t36*t37*t273*2.0-t2*t8*t29*t32*t35*t36*t273*2.0-t4*t10*t21*t29*t38*t42*t273*2.0;
    hess(7, 1) = -t278-t279-t20*t35*(t277-t3*t10*t32+t7*t36*t37*t273)*2.0+t35*t44*t275*2.0+t20*t38*t48*t275*2.0-t2*t3*t5*t35*t36*t37*t273*2.0-t2*t5*t10*t29*t38*t48*t273*2.0-t3*t4*t21*t35*t36*t37*t273*2.0-t3*t6*t29*t32*t35*t36*t273*2.0-t3*t8*t29*t32*t35*t36*t273*2.0-t4*t10*t21*t29*t38*t48*t273*2.0;
    hess(7, 2) = -t228-t281-t20*t35*(t280-t4*t10*t32+t3*t4*t36*t37*t273)*2.0+t20*t38*t54*t275*2.0+t3*t4*t10*t32*t35*t275*2.0-t4*t10*t29*t35*t50*t273*2.0-t8*t21*t35*t36*t37*t273*2.0-t2*t4*t5*t35*t36*t37*t273*2.0-t2*t5*t10*t29*t38*t54*t273*2.0-t4*t6*t29*t32*t35*t36*t273*2.0-t4*t10*t21*t29*t38*t54*t273*2.0;
    hess(7, 3) = t293-t20*t35*(t292-t32*t61+t3*t10*t29*t284+t3*t10*t37*t61*t273)*2.0+t20*t38*t67*t275*2.0+t3*t32*t35*t61*t275*2.0-t2*t5*t10*t29*t35*t284*2.0-t4*t10*t21*t29*t35*t284*2.0-t2*t10*t29*t35*t63*t273*2.0-t2*t5*t10*t35*t37*t61*t273*2.0-t2*t5*t10*t29*t38*t67*t273*2.0-t2*t5*t29*t35*t36*t57*t273*2.0-t4*t10*t21*t35*t37*t61*t273*2.0-t4*t10*t21*t29*t38*t67*t273*2.0-t8*t10*t29*t32*t35*t61*t273*2.0-t4*t21*t29*t35*t36*t57*t273*2.0;
    hess(7, 4) = t35*t76*t275*2.0-t20*t35*(t294-t32*t74-t10*t29*t273+t3*t10*t29*t287+t3*t10*t37*t74*t273)*2.0+t20*t38*t80*t275*2.0-t2*t5*t10*t29*t35*t287*2.0-t4*t10*t21*t29*t35*t287*2.0-t2*t5*t10*t35*t37*t74*t273*2.0-t2*t5*t10*t29*t38*t80*t273*2.0-t2*t5*t29*t35*t36*t70*t273*2.0-t4*t10*t21*t35*t37*t74*t273*2.0-t4*t10*t21*t29*t38*t80*t273*2.0-t6*t10*t29*t32*t35*t74*t273*2.0-t8*t10*t29*t32*t35*t74*t273*2.0-t4*t21*t29*t35*t36*t70*t273*2.0;
    hess(7, 5) = t298-t20*t35*(t297-t32*t87+t3*t10*t29*t290+t3*t10*t37*t87*t273)*2.0+t20*t38*t93*t275*2.0+t3*t32*t35*t87*t275*2.0-t2*t5*t10*t29*t35*t290*2.0-t4*t10*t21*t29*t35*t290*2.0-t4*t10*t29*t35*t89*t273*2.0-t2*t5*t10*t35*t37*t87*t273*2.0-t2*t5*t10*t29*t38*t93*t273*2.0-t2*t5*t29*t35*t36*t83*t273*2.0-t4*t10*t21*t35*t37*t87*t273*2.0-t4*t10*t21*t29*t38*t93*t273*2.0-t6*t10*t29*t32*t35*t87*t273*2.0-t4*t21*t29*t35*t36*t83*t273*2.0;
    hess(7, 6) = -t293+t20*t35*(-t291+t292+t3*t10*t29*(t283-t3*t10*t11*2.0)+t3*t10*t37*t273*(t60-t96))*2.0-t20*t38*t101*t275*2.0-t3*t32*t35*t95*t275*2.0-t2*t10*t29*t35*t97*t273*2.0+t2*t5*t10*t29*t35*(t283-t3*t10*t11*2.0)*2.0+t4*t10*t21*t29*t35*(t283-t3*t10*t11*2.0)*2.0+t2*t5*t10*t29*t38*t273*(t99+t100-t102)*2.0+t4*t10*t21*t29*t38*t273*(t99+t100-t102)*2.0+t2*t5*t29*t35*t36*t57*t273*2.0+t4*t21*t29*t35*t36*t57*t273*2.0+t2*t5*t10*t35*t37*t273*(t60-t96)*2.0+t4*t10*t21*t35*t37*t273*(t60-t96)*2.0+t8*t10*t29*t32*t35*t273*(t60-t96)*2.0;
    hess(7, 7) = t20*t35*(-t294+t321+t10*t29*t273+t3*t10*t29*t296-t3*t10*t37*t106*t273)*-2.0+t35*t108*t275*2.0-t20*t38*t112*t275*2.0-t2*t5*t10*t29*t35*t296*2.0-t4*t10*t21*t29*t35*t296*2.0+t2*t5*t10*t29*t38*t273*(t110+t111-t113)*2.0+t4*t10*t21*t29*t38*t273*(t110+t111-t113)*2.0+t2*t5*t29*t35*t36*t70*t273*2.0+t4*t21*t29*t35*t36*t70*t273*2.0+t2*t5*t10*t35*t37*t273*(t71-t107)*2.0+t4*t10*t21*t35*t37*t273*(t71-t107)*2.0+t6*t10*t29*t32*t35*t273*(t71-t107)*2.0+t8*t10*t29*t32*t35*t273*(t71-t107)*2.0;
    hess(7, 8) = -t298+t20*t35*(t297-t32*t117+t3*t10*t29*(t289-t3*t10*t15*2.0)+t3*t10*t37*t273*(t84-t118))*2.0-t20*t38*t123*t275*2.0-t3*t32*t35*t117*t275*2.0-t4*t10*t29*t35*t119*t273*2.0+t2*t5*t10*t29*t35*(t289-t3*t10*t15*2.0)*2.0+t4*t10*t21*t29*t35*(t289-t3*t10*t15*2.0)*2.0+t2*t5*t10*t29*t38*t273*(t121+t122-t124)*2.0+t4*t10*t21*t29*t38*t273*(t121+t122-t124)*2.0+t2*t5*t29*t35*t36*t83*t273*2.0+t4*t21*t29*t35*t36*t83*t273*2.0+t2*t5*t10*t35*t37*t273*(t84-t118)*2.0+t4*t10*t21*t35*t37*t273*(t84-t118)*2.0+t6*t10*t29*t32*t35*t273*(t84-t118)*2.0;
    hess(8, 0) = -t227-t302-t21*t35*(t299-t2*t10*t32+t2*t4*t36*t37*t301)*2.0+t21*t38*t42*t303*2.0+t2*t4*t10*t32*t35*t303*2.0-t2*t10*t29*t34*t35*t301*2.0-t5*t6*t35*t36*t37*t301*2.0-t2*t5*t10*t29*t38*t42*t301*2.0-t2*t3*t20*t35*t36*t37*t301*2.0-t2*t7*t29*t32*t35*t36*t301*2.0-t3*t10*t20*t29*t38*t42*t301*2.0;
    hess(8, 1) = -t228-t306-t21*t35*(t304-t3*t10*t32+t3*t4*t36*t37*t301)*2.0+t21*t38*t48*t303*2.0+t3*t4*t10*t32*t35*t303*2.0-t3*t10*t29*t35*t44*t301*2.0-t7*t20*t35*t36*t37*t301*2.0-t2*t3*t5*t35*t36*t37*t301*2.0-t2*t5*t10*t29*t38*t48*t301*2.0-t3*t6*t29*t32*t35*t36*t301*2.0-t3*t10*t20*t29*t38*t48*t301*2.0;
    hess(8, 2) = -t308-t309-t21*t35*(t307-t4*t10*t32+t8*t36*t37*t301)*2.0+t35*t50*t303*2.0+t21*t38*t54*t303*2.0-t2*t4*t5*t35*t36*t37*t301*2.0-t3*t4*t20*t35*t36*t37*t301*2.0-t2*t5*t10*t29*t38*t54*t301*2.0-t4*t6*t29*t32*t35*t36*t301*2.0-t4*t7*t29*t32*t35*t36*t301*2.0-t3*t10*t20*t29*t38*t54*t301*2.0;
    hess(8, 3) = t320-t21*t35*(t319-t32*t61+t4*t10*t29*t312+t4*t10*t37*t61*t301)*2.0+t21*t38*t67*t303*2.0+t4*t32*t35*t61*t303*2.0-t2*t5*t10*t29*t35*t312*2.0-t3*t10*t20*t29*t35*t312*2.0-t2*t10*t29*t35*t63*t301*2.0-t2*t5*t10*t35*t37*t61*t301*2.0-t2*t5*t10*t29*t38*t67*t301*2.0-t2*t5*t29*t35*t36*t57*t301*2.0-t3*t10*t20*t35*t37*t61*t301*2.0-t3*t10*t20*t29*t38*t67*t301*2.0-t7*t10*t29*t32*t35*t61*t301*2.0-t3*t20*t29*t35*t36*t57*t301*2.0;
    hess(8, 4) = t323-t21*t35*(t322-t32*t74+t4*t10*t29*t315+t4*t10*t37*t74*t301)*2.0+t21*t38*t80*t303*2.0+t4*t32*t35*t74*t303*2.0-t2*t5*t10*t29*t35*t315*2.0-t3*t10*t20*t29*t35*t315*2.0-t3*t10*t29*t35*t76*t301*2.0-t2*t5*t10*t35*t37*t74*t301*2.0-t2*t5*t10*t29*t38*t80*t301*2.0-t2*t5*t29*t35*t36*t70*t301*2.0-t3*t10*t20*t35*t37*t74*t301*2.0-t3*t10*t20*t29*t38*t80*t301*2.0-t6*t10*t29*t32*t35*t74*t301*2.0-t3*t20*t29*t35*t36*t70*t301*2.0;
    hess(8, 5) = t35*t89*t303*2.0-t21*t35*(t324-t32*t87-t10*t29*t301+t4*t10*t29*t318+t4*t10*t37*t87*t301)*2.0+t21*t38*t93*t303*2.0-t2*t5*t10*t29*t35*t318*2.0-t3*t10*t20*t29*t35*t318*2.0-t2*t5*t10*t35*t37*t87*t301*2.0-t2*t5*t10*t29*t38*t93*t301*2.0-t2*t5*t29*t35*t36*t83*t301*2.0-t3*t10*t20*t35*t37*t87*t301*2.0-t3*t10*t20*t29*t38*t93*t301*2.0-t6*t10*t29*t32*t35*t87*t301*2.0-t7*t10*t29*t32*t35*t87*t301*2.0-t3*t20*t29*t35*t36*t83*t301*2.0;
    hess(8, 6) = -t320+t21*t35*(-t291+t319+t4*t10*t29*(t311-t4*t10*t11*2.0)+t4*t10*t37*t301*(t60-t96))*2.0-t21*t38*t101*t303*2.0-t4*t32*t35*t95*t303*2.0-t2*t10*t29*t35*t97*t301*2.0+t2*t5*t10*t29*t35*(t311-t4*t10*t11*2.0)*2.0+t3*t10*t20*t29*t35*(t311-t4*t10*t11*2.0)*2.0+t2*t5*t10*t29*t38*t301*(t99+t100-t102)*2.0+t3*t10*t20*t29*t38*t301*(t99+t100-t102)*2.0+t2*t5*t29*t35*t36*t57*t301*2.0+t3*t20*t29*t35*t36*t57*t301*2.0+t2*t5*t10*t35*t37*t301*(t60-t96)*2.0+t3*t10*t20*t35*t37*t301*(t60-t96)*2.0+t7*t10*t29*t32*t35*t301*(t60-t96)*2.0;
    hess(8, 7) = -t323+t21*t35*(-t321+t322+t4*t10*t29*(t314-t4*t10*t13*2.0)+t4*t10*t37*t301*(t71-t107))*2.0-t21*t38*t112*t303*2.0-t4*t32*t35*t106*t303*2.0-t3*t10*t29*t35*t108*t301*2.0+t2*t5*t10*t29*t35*(t314-t4*t10*t13*2.0)*2.0+t3*t10*t20*t29*t35*(t314-t4*t10*t13*2.0)*2.0+t2*t5*t10*t29*t38*t301*(t110+t111-t113)*2.0+t3*t10*t20*t29*t38*t301*(t110+t111-t113)*2.0+t2*t5*t29*t35*t36*t70*t301*2.0+t3*t20*t29*t35*t36*t70*t301*2.0+t2*t5*t10*t35*t37*t301*(t71-t107)*2.0+t3*t10*t20*t35*t37*t301*(t71-t107)*2.0+t6*t10*t29*t32*t35*t301*(t71-t107)*2.0;
    hess(8, 8) = t21*t35*(-t324+t32*t117+t10*t29*t301+t4*t10*t29*t326-t4*t10*t37*t117*t301)*-2.0+t35*t119*t303*2.0-t21*t38*t123*t303*2.0-t2*t5*t10*t29*t35*t326*2.0-t3*t10*t20*t29*t35*t326*2.0+t2*t5*t10*t29*t38*t301*(t121+t122-t124)*2.0+t3*t10*t20*t29*t38*t301*(t121+t122-t124)*2.0+t2*t5*t29*t35*t36*t83*t301*2.0+t3*t20*t29*t35*t36*t83*t301*2.0+t2*t5*t10*t35*t37*t301*(t84-t118)*2.0+t3*t10*t20*t35*t37*t301*(t84-t118)*2.0+t6*t10*t29*t32*t35*t301*(t84-t118)*2.0+t7*t10*t29*t32*t35*t301*(t84-t118)*2.0;
  }
  
  template<typename T, int DIM>
  inline T softcapsule_phi(const Eigen::Matrix<T, DIM, 1>& p, const Eigen::Matrix<T, DIM, 1>& a, const Eigen::Matrix<T, DIM, 1>& b, const T& radius, const T& h) {
    Eigen::Matrix<T, DIM, 1> pa = p - a;
    Eigen::Matrix<T, DIM, 1> ba = b - a;
    T alpha = softclamp(pa.dot(ba) / ba.dot(ba), 0.0, 1.0, h);
    return softsqrt((pa - ba * alpha).squaredNorm(), h) - radius;
  }
  
  template<typename T, int DIM>
  inline T softcircle_phi(const Eigen::Matrix<T, DIM, 1>& position, const Eigen::Matrix<T, DIM, 1>& centre, const T& radius, const T& h) {
    return softsqrt((position - centre).squaredNorm(), h) - radius;
  }
  
  template<typename T, int DIM>
  inline void gradsoftcircle_phi(Eigen::Matrix<T, DIM * 2, 1>& grad, const Eigen::Matrix<T, DIM, 1>& position, const Eigen::Matrix<T, DIM, 1>& centre, const T& radius, const T& h) {
    grad.segment(0, DIM) = (position - centre) * 2.0 * gradsoftsqrt((position - centre).squaredNorm(), h);
    grad.segment(DIM, DIM) = -grad.segment(0, DIM);
  }
  
  template<typename T, int DIM>
  inline void hesssoftcircle_phi(Eigen::Matrix<T, DIM * 2, DIM * 2>& hess, const Eigen::Matrix<T, DIM, 1>& position, const Eigen::Matrix<T, DIM, 1>& centre, const T& radius, const T& h) {
    Eigen::Matrix<T, DIM, 1> dp = (position - centre);
    T sndp = dp.squaredNorm();
    Eigen::Matrix<T, DIM, DIM> M = 2.0 * Eigen::Matrix<T, DIM, DIM>::Identity() * gradsoftsqrt(sndp, h) + dp * dp.transpose() * hesssoftsqrt(sndp, h) * 4.0;
    hess.block(0, 0, DIM, DIM) = M;
    hess.block(0, DIM, DIM, DIM) = -M;
    hess.block(DIM, 0, DIM, DIM) = -M;
    hess.block(DIM, DIM, DIM, DIM) = M;
  }
  
  template<typename T, int DIM>
  inline T pointedgedist(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, T& alpha)
  {
    alpha = std::min((T) 1.0, std::max((T) 0.0, (x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2)));
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha*(x3-x2);
    Eigen::Matrix<T, DIM, 1> n = closest-x1;
    return n.norm();
  }
  
  template<typename T, int DIM>
  inline T softpointedgedist(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, T& alpha, const T& h)
  {
    alpha = softclamp((x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2), (T) 0.0, (T) 1.0, h);
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha*(x3-x2);
    Eigen::Matrix<T, DIM, 1> n = closest-x1;
    return n.norm();
  }
  
  template<typename T>
  inline void gradsoftpointedgedist(
                                 const Eigen::Matrix<T, 2, 1>& x0,
                                 const Eigen::Matrix<T, 2, 1>& x1,
                                 const Eigen::Matrix<T, 2, 1>& x2,
                                 Vector6s& gradE,
                                 const T& h)
  {
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t4 = t2*t2;
    const scalar t5 = t3*t3;
    const scalar t6 = t4+t5;
    const scalar t7 = 1.0/t6;
    const scalar t8 = x0(0)-x1(0);
    const scalar t9 = t2*t8;
    const scalar t10 = x0(1)-x1(1);
    const scalar t11 = t3*t10;
    const scalar t12 = t9+t11;
    const scalar t14 = t7*t12;
    const scalar t13 = gradsoftclamp(-t14, 0.0, 1.0, h);
    const scalar t15 = softclamp(-t14, 0.0, 1.0, h);
    const scalar t19 = t2*t15;
    const scalar t16 = t19+x0(0)-x1(0);
    const scalar t17 = t3*t15;
    const scalar t18 = t17+x0(1)-x1(1);
    const scalar t20 = t16*t16;
    const scalar t21 = t18*t18;
    const scalar t22 = t20+t21;
    const scalar t23 = 1.0/sqrt(t22);
    const scalar t24 = x1(0)*2.0;
    const scalar t35 = x2(0)*2.0;
    const scalar t25 = t24-t35;
    const scalar t26 = 1.0/(t6*t6);
    const scalar t27 = -t24+x0(0)+x2(0);
    const scalar t28 = t7*t27;
    const scalar t36 = t12*t25*t26;
    const scalar t29 = t28-t36;
    const scalar t30 = x1(1)*2.0;
    const scalar t38 = x2(1)*2.0;
    const scalar t31 = t30-t38;
    const scalar t32 = -t30+x0(1)+x2(1);
    const scalar t33 = t7*t32;
    const scalar t39 = t12*t26*t31;
    const scalar t34 = t33-t39;
    const scalar t37 = t36-t7*t8;
    const scalar t40 = t39-t7*t10;
    gradE(0) += t23*(t16*(t4*t7*t13-1.0)*2.0+t2*t3*t7*t13*t18*2.0)*(-1.0/2.0);
    gradE(1) += t23*(t18*(t5*t7*t13-1.0)*2.0+t2*t3*t7*t13*t16*2.0)*(-1.0/2.0);
    gradE(2) += t23*(t16*(-t15+t2*t13*t29+1.0)*2.0+t3*t13*t18*t29*2.0)*(-1.0/2.0);
    gradE(3) += t23*(t18*(-t15+t3*t13*t34+1.0)*2.0+t2*t13*t16*t34*2.0)*(-1.0/2.0);
    gradE(4) += t23*(t16*(t15+t2*t13*t37)*2.0+t3*t13*t18*t37*2.0)*(-1.0/2.0);
    gradE(5) += t23*(t18*(t15+t3*t13*t40)*2.0+t2*t13*t16*t40*2.0)*(-1.0/2.0);
  }
  
  template<typename T>
  inline void gradsoftpointedgedist(
                                    const Eigen::Matrix<T, 3, 1>& x0,
                                    const Eigen::Matrix<T, 3, 1>& x1,
                                    const Eigen::Matrix<T, 3, 1>& x2,
                                    Vector9s& grad,
                                    const T& h)
  {
    
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t4 = x1(2)-x2(2);
    const scalar t6 = t2*t2;
    const scalar t7 = t3*t3;
    const scalar t8 = t4*t4;
    const scalar t9 = t6+t7+t8;
    const scalar t10 = 1.0/t9;
    const scalar t11 = x0(0)-x1(0);
    const scalar t12 = t2*t11;
    const scalar t13 = x0(1)-x1(1);
    const scalar t14 = t3*t13;
    const scalar t15 = x0(2)-x1(2);
    const scalar t16 = t4*t15;
    const scalar t17 = t12+t14+t16;
    const scalar t18 = t10*t17;
    const scalar t19 = softclamp(-t18, 0.0, 1.0, h);
    const scalar t22 = t2*t19;
    const scalar t5 = t22+x0(0)-x1(0);
    const scalar t24 = t3*t19;
    const scalar t20 = t24+x0(1)-x1(1);
    const scalar t25 = t4*t19;
    const scalar t21 = t25+x0(2)-x1(2);
    const scalar t23 = gradsoftclamp(-t18, 0.0, 1.0, h);
    const scalar t26 = t5*t5;
    const scalar t27 = t20*t20;
    const scalar t28 = t21*t21;
    const scalar t29 = t26+t27+t28;
    const scalar t30 = 1.0/sqrt(t29);
    const scalar t31 = x1(0)*2.0;
    const scalar t36 = x2(0)*2.0;
    const scalar t32 = t31-t36;
    const scalar t33 = 1.0/(t9*t9);
    const scalar t34 = -t31+x0(0)+x2(0);
    const scalar t35 = t10*t34;
    const scalar t38 = t17*t32*t33;
    const scalar t37 = t35-t38;
    const scalar t39 = x1(1)*2.0;
    const scalar t43 = x2(1)*2.0;
    const scalar t40 = t39-t43;
    const scalar t41 = -t39+x0(1)+x2(1);
    const scalar t42 = t10*t41;
    const scalar t45 = t17*t33*t40;
    const scalar t44 = t42-t45;
    const scalar t46 = x1(2)*2.0;
    const scalar t50 = x2(2)*2.0;
    const scalar t47 = t46-t50;
    const scalar t48 = -t46+x0(2)+x2(2);
    const scalar t49 = t10*t48;
    const scalar t52 = t17*t33*t47;
    const scalar t51 = t49-t52;
    const scalar t54 = t10*t11;
    const scalar t53 = t38-t54;
    const scalar t56 = t10*t13;
    const scalar t55 = t45-t56;
    const scalar t58 = t10*t15;
    const scalar t57 = t52-t58;
    grad(0) = t30*(t5*(t6*t10*t23-1.0)*2.0+t2*t3*t10*t20*t23*2.0+t2*t4*t10*t21*t23*2.0)*(-1.0/2.0);
    grad(1) = t30*(t20*(t7*t10*t23-1.0)*2.0+t2*t3*t5*t10*t23*2.0+t3*t4*t10*t21*t23*2.0)*(-1.0/2.0);
    grad(2) = t30*(t21*(t8*t10*t23-1.0)*2.0+t2*t4*t5*t10*t23*2.0+t3*t4*t10*t20*t23*2.0)*(-1.0/2.0);
    grad(3) = t30*(t5*(-t19+t2*t23*t37+1.0)*2.0+t3*t20*t23*t37*2.0+t4*t21*t23*t37*2.0)*(-1.0/2.0);
    grad(4) = t30*(t20*(-t19+t3*t23*t44+1.0)*2.0+t2*t5*t23*t44*2.0+t4*t21*t23*t44*2.0)*(-1.0/2.0);
    grad(5) = t30*(t21*(-t19+t4*t23*t51+1.0)*2.0+t2*t5*t23*t51*2.0+t3*t20*t23*t51*2.0)*(-1.0/2.0);
    grad(6) = t30*(t5*(t19+t2*t23*t53)*2.0+t3*t20*t23*t53*2.0+t4*t21*t23*t53*2.0)*(-1.0/2.0);
    grad(7) = t30*(t20*(t19+t3*t23*t55)*2.0+t2*t5*t23*t55*2.0+t4*t21*t23*t55*2.0)*(-1.0/2.0);
    grad(8) = t30*(t21*(t19+t4*t23*t57)*2.0+t2*t5*t23*t57*2.0+t3*t20*t23*t57*2.0)*(-1.0/2.0);
  }
  
  template<typename T>
  inline void hesssoftpointedgedist(
                                    const Eigen::Matrix<T, 2, 1>& x0,
                                    const Eigen::Matrix<T, 2, 1>& x1,
                                    const Eigen::Matrix<T, 2, 1>& x2,
                                    Matrix6s& hess,
                                    const T& h)
  {
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t5 = t2*t2;
    const scalar t6 = t3*t3;
    const scalar t7 = t5+t6;
    const scalar t8 = 1.0/t7;
    const scalar t9 = x0(0)-x1(0);
    const scalar t10 = t2*t9;
    const scalar t11 = x0(1)-x1(1);
    const scalar t12 = t3*t11;
    const scalar t13 = t10+t12;
    const scalar t14 = t8*t13;
    const scalar t15 = softclamp(-t14, 0.0, 1.0, h);
    const scalar t20 = t2*t15;
    const scalar t4 = t20+x0(0)-x1(0);
    const scalar t22 = t3*t15;
    const scalar t16 = t22+x0(1)-x1(1);
    const scalar t18 = gradsoftclamp(-t14, 0.0, 1.0, h);
    const scalar t23 = t5*t8*t18;
    const scalar t17 = t23-1.0;
    const scalar t19 = 1.0/(t7*t7);
    const scalar t21 = hesssoftclamp(-t14, 0.0, 1.0, h);
    const scalar t29 = t4*t17*2.0;
    const scalar t30 = t2*t3*t8*t16*t18*2.0;
    const scalar t24 = t29+t30;
    const scalar t25 = t4*t4;
    const scalar t26 = t16*t16;
    const scalar t27 = t25+t26;
    const scalar t28 = 1.0/sqrt(t27);
    const scalar t31 = t6*t8*t18;
    const scalar t32 = t31-1.0;
    const scalar t33 = 1.0/pow(t27,3.0/2.0);
    const scalar t34 = x1(0)*2.0;
    const scalar t36 = x2(0)*2.0;
    const scalar t35 = t34-t36;
    const scalar t37 = -t34+x0(0)+x2(0);
    const scalar t38 = t8*t37;
    const scalar t41 = t13*t19*t35;
    const scalar t39 = t38-t41;
    const scalar t40 = t18*t18;
    const scalar t42 = t2*t18*t39;
    const scalar t43 = -t15+t42+1.0;
    const scalar t44 = x1(1)*2.0;
    const scalar t48 = x2(1)*2.0;
    const scalar t45 = t44-t48;
    const scalar t46 = -t44+x0(1)+x2(1);
    const scalar t47 = t8*t46;
    const scalar t50 = t13*t19*t45;
    const scalar t49 = t47-t50;
    const scalar t51 = t3*t18*t49;
    const scalar t52 = -t15+t51+1.0;
    const scalar t53 = t5*t18*t19*t35;
    const scalar t55 = t8*t9;
    const scalar t54 = t41-t55;
    const scalar t56 = t2*t3*t16*t18*t19*t35*2.0;
    const scalar t57 = t2*t18*t54;
    const scalar t58 = t15+t57;
    const scalar t59 = t5*t18*t19*t45;
    const scalar t61 = t8*t11;
    const scalar t60 = t50-t61;
    const scalar t62 = t2*t3*t16*t18*t19*t45*2.0;
    const scalar t63 = t3*t18*t60;
    const scalar t64 = t15+t63;
    const scalar t65 = t3*t4*t5*t19*t21*2.0;
    const scalar t66 = t2*t6*t16*t19*t21*2.0;
    const scalar t67 = t2*t3*t8*t17*t18*2.0;
    const scalar t68 = t2*t3*t8*t18*t32*2.0;
    const scalar t69 = t65+t66+t67+t68;
    const scalar t70 = t28*t69*(1.0/2.0);
    const scalar t71 = t16*t32*2.0;
    const scalar t72 = t2*t3*t4*t8*t18*2.0;
    const scalar t73 = t71+t72;
    const scalar t74 = t70-t24*t33*t73*(1.0/4.0);
    const scalar t75 = t5*t6*t19*t40*2.0;
    const scalar t76 = t4*t43*2.0;
    const scalar t77 = t3*t16*t18*t39*2.0;
    const scalar t78 = t76+t77;
    const scalar t79 = t16*t52*2.0;
    const scalar t80 = t2*t4*t18*t49*2.0;
    const scalar t81 = t79+t80;
    const scalar t82 = t6*t18*t19*t35;
    const scalar t83 = t2*t3*t4*t18*t19*t35*2.0;
    const scalar t84 = t4*t58*2.0;
    const scalar t85 = t3*t16*t18*t54*2.0;
    const scalar t86 = t84+t85;
    const scalar t87 = t6*t18*t19*t45;
    const scalar t88 = t2*t3*t4*t18*t19*t45*2.0;
    const scalar t89 = t16*t64*2.0;
    const scalar t90 = t2*t4*t18*t60*2.0;
    const scalar t91 = t89+t90;
    const scalar t92 = t5*t8*t21*t39;
    const scalar t93 = t17*t43*2.0;
    const scalar t137 = t2*t19*t35;
    const scalar t94 = t8-t137;
    const scalar t95 = t2*t6*t8*t39*t40*2.0;
    const scalar t96 = t2*t3*t8*t16*t21*t39*2.0;
    const scalar t97 = t3*t18*t32*t39*2.0;
    const scalar t98 = t2*t3*t8*t18*t43*2.0;
    const scalar t99 = t39*t39;
    const scalar t100 = t13*t19*2.0;
    const scalar t101 = t8*2.0;
    const scalar t102 = t19*t35*t37*2.0;
    const scalar t103 = t35*t35;
    const scalar t104 = 1.0/(t7*t7*t7);
    const scalar t109 = t13*t103*t104*2.0;
    const scalar t105 = t100+t101+t102-t109;
    const scalar t106 = t19*t37*t45;
    const scalar t107 = t19*t35*t46;
    const scalar t113 = t13*t35*t45*t104*2.0;
    const scalar t108 = t106+t107-t113;
    const scalar t110 = t19*t35*t37;
    const scalar t111 = t9*t19*t35;
    const scalar t112 = t8+t100-t109+t110+t111;
    const scalar t114 = t11*t19*t35;
    const scalar t115 = t106-t113+t114;
    const scalar t116 = t2*t8*t18;
    const scalar t117 = t2*t17*t18*t49*2.0;
    const scalar t118 = t2*t3*t8*t18*t52*2.0;
    const scalar t119 = t6*t8*t21*t49;
    const scalar t120 = t32*t52*2.0;
    const scalar t170 = t3*t19*t45;
    const scalar t121 = t8-t170;
    const scalar t122 = t3*t5*t8*t40*t49*2.0;
    const scalar t123 = t2*t3*t4*t8*t21*t49*2.0;
    const scalar t124 = t18*t39;
    const scalar t125 = t2*t18*t43*t49*2.0;
    const scalar t126 = t3*t18*t39*t52*2.0;
    const scalar t127 = t49*t49;
    const scalar t128 = t19*t45*t46*2.0;
    const scalar t129 = t45*t45;
    const scalar t133 = t13*t104*t129*2.0;
    const scalar t130 = t100+t101+t128-t133;
    const scalar t131 = t9*t19*t45;
    const scalar t132 = t107-t113+t131;
    const scalar t134 = t19*t45*t46;
    const scalar t135 = t11*t19*t45;
    const scalar t136 = t8+t100-t133+t134+t135;
    const scalar t138 = t2*t18*t94;
    const scalar t139 = t2*t18*(t41-t55);
    const scalar t140 = t15+t139;
    const scalar t141 = t3*t8*t18;
    const scalar t142 = t2*t3*t18*t19*t35;
    const scalar t143 = t3*t18*t32*t54*2.0;
    const scalar t144 = t6*t16*t18*t19*t35*2.0;
    const scalar t145 = t4*t140*2.0;
    const scalar t146 = t85+t145;
    const scalar t147 = t2*t21*t39*t54;
    const scalar t148 = t124+t147-t18*t54-t2*t18*t112;
    const scalar t149 = t4*t148*2.0;
    const scalar t150 = t6*t39*t40*t54*2.0;
    const scalar t151 = t3*t16*t21*t39*t54*2.0;
    const scalar t152 = t18*t49;
    const scalar t153 = t3*t18*t52*t54*2.0;
    const scalar t162 = t3*t16*t18*(t41-t55)*2.0;
    const scalar t154 = t145+t162;
    const scalar t155 = t41-t55;
    const scalar t156 = t41-t55;
    const scalar t157 = t9*t19*t35*2.0;
    const scalar t158 = t100-t109+t157;
    const scalar t159 = t41-t55;
    const scalar t160 = t18*(t50-t61);
    const scalar t161 = -t113+t114+t131;
    const scalar t163 = t3*t18*(t50-t61);
    const scalar t164 = t15+t163;
    const scalar t165 = t2*t3*t18*t19*t45;
    const scalar t166 = t2*t17*t18*t60*2.0;
    const scalar t167 = t4*t5*t18*t19*t45*2.0;
    const scalar t168 = t16*t164*2.0;
    const scalar t169 = t90+t168;
    const scalar t171 = t3*t18*t121;
    const scalar t172 = t2*t18*t43*t60*2.0;
    const scalar t173 = t3*t21*t49*t60;
    const scalar t174 = t152+t173-t18*t60-t3*t18*t136;
    const scalar t175 = t16*t174*2.0;
    const scalar t176 = t5*t40*t49*t60*2.0;
    const scalar t177 = t2*t4*t21*t49*t60*2.0;
    const scalar t178 = t18*(t41-t55);
    const scalar t179 = t2*t18*t140*(t50-t61)*2.0;
    const scalar t180 = t3*t18*t164*(t41-t55)*2.0;
    const scalar t181 = t168+t2*t4*t18*(t50-t61)*2.0;
    const scalar t182 = t50-t61;
    const scalar t183 = t50-t61;
    const scalar t184 = t11*t19*t45*2.0;
    const scalar t185 = t100-t133+t184;
    const scalar t186 = t50-t61;
    hess(0, 0) = t28*(t75+(t17*t17)*2.0+t2*t4*t5*t19*t21*2.0+t3*t5*t16*t19*t21*2.0)*(1.0/2.0)-(t24*t24)*t33*(1.0/4.0);
    hess(0, 1) = t74;
    hess(0, 2) = t28*(t56+t93+t95+t96+t4*(t53+t92-t8*t18*t35)*2.0-t3*t8*t16*t18*2.0)*(1.0/2.0)-t24*t33*t78*(1.0/4.0);
    hess(0, 3) = t28*(t62+t117+t118+t4*(t59+t5*t8*t21*t49)*2.0-t2*t8*t16*t18*2.0+t2*t3*t8*t16*t21*t49*2.0)*(1.0/2.0)-t24*t33*t81*(1.0/4.0);
    hess(0, 4) = t28*(-t56+t17*t58*2.0+t4*(-t53+t8*t18*t35+t5*t8*t21*t54)*2.0+t3*t8*t16*t18*2.0+t2*t6*t8*t40*t54*2.0+t2*t3*t8*t16*t21*t54*2.0)*(1.0/2.0)-t24*t33*t86*(1.0/4.0);
    hess(0, 5) = t28*(-t62+t166-t4*(t59-t5*t8*t21*t60)*2.0+t2*t8*t16*t18*2.0+t2*t3*t8*t18*t64*2.0+t2*t3*t8*t16*t21*t60*2.0)*(1.0/2.0)-t24*t33*t91*(1.0/4.0);
    hess(1, 0) = t74;
    hess(1, 1) = t28*(t75+(t32*t32)*2.0+t2*t4*t6*t19*t21*2.0+t3*t6*t16*t19*t21*2.0)*(1.0/2.0)-t33*(t73*t73)*(1.0/4.0);
    hess(1, 2) = t28*(t83+t97+t98+t16*(t82+t6*t8*t21*t39)*2.0-t3*t4*t8*t18*2.0+t2*t3*t4*t8*t21*t39*2.0)*(1.0/2.0)-t33*t73*t78*(1.0/4.0);
    hess(1, 3) = t28*(t88+t120+t122+t123+t16*(t87+t119-t8*t18*t45)*2.0-t2*t4*t8*t18*2.0)*(1.0/2.0)-t33*t73*t81*(1.0/4.0);
    hess(1, 4) = t28*(-t83+t143-t16*(t82-t6*t8*t21*t54)*2.0+t3*t4*t8*t18*2.0+t2*t3*t8*t18*t58*2.0+t2*t3*t4*t8*t21*t54*2.0)*(1.0/2.0)-t33*t73*t86*(1.0/4.0);
    hess(1, 5) = t28*(-t88+t32*t64*2.0+t16*(-t87+t8*t18*t45+t6*t8*t21*t60)*2.0+t2*t4*t8*t18*2.0+t3*t5*t8*t40*t60*2.0+t2*t3*t4*t8*t21*t60*2.0)*(1.0/2.0)-t33*t73*t91*(1.0/4.0);
    hess(2, 0) = t28*(t93+t95+t96-t4*(-t92+t116+t138)*2.0-t3*t16*t18*t94*2.0)*(1.0/2.0)-t24*t33*t78*(1.0/4.0);
    hess(2, 1) = t28*(t97+t98+t144+t4*(t142-t3*t8*t18+t2*t3*t8*t21*t39)*2.0+t6*t8*t16*t21*t39*2.0)*(1.0/2.0)-t33*t73*t78*(1.0/4.0);
    hess(2, 2) = t28*((t43*t43)*2.0+t4*(t18*t39*-2.0+t2*t21*t99+t2*t18*t105)*2.0+t6*t40*t99*2.0+t3*t16*t21*t99*2.0+t3*t16*t18*t105*2.0)*(1.0/2.0)-t33*(t78*t78)*(1.0/4.0);
    hess(2, 3) = t28*(t125+t126+t4*(-t18*t49+t2*t18*t108+t2*t21*t39*t49)*2.0-t16*t18*t39*2.0+t3*t16*t18*t108*2.0+t3*t16*t21*t39*t49*2.0)*(1.0/2.0)-t33*t78*t81*(1.0/4.0);
    hess(2, 4) = t28*(t149+t150+t151+t43*t58*2.0-t3*t16*t18*t112*2.0)*(1.0/2.0)-t33*t78*t86*(1.0/4.0);
    hess(2, 5) = t28*(t172-t4*(t160+t2*t18*t115-t2*t21*t39*t60)*2.0+t16*t18*t39*2.0+t3*t18*t39*t64*2.0-t3*t16*t18*t115*2.0+t3*t16*t21*t39*t60*2.0)*(1.0/2.0)-t33*t78*t91*(1.0/4.0);
    hess(3, 0) = t28*(t117+t118+t167+t16*(-t116+t165+t2*t3*t8*t21*t49)*2.0+t4*t5*t8*t21*t49*2.0)*(1.0/2.0)-t24*t33*t81*(1.0/4.0);
    hess(3, 1) = t28*(t120+t122+t123-t16*(-t119+t141+t171)*2.0-t2*t4*t18*t121*2.0)*(1.0/2.0)-t33*t73*t81*(1.0/4.0);
    hess(3, 2) = t28*(t125+t126+t16*(-t124+t3*t18*t108+t3*t21*t39*t49)*2.0-t4*t18*t49*2.0+t2*t4*t18*t108*2.0+t2*t4*t21*t39*t49*2.0)*(1.0/2.0)-t33*t78*t81*(1.0/4.0);
    hess(3, 3) = t28*((t52*t52)*2.0+t16*(t18*t49*-2.0+t3*t18*t130+t3*t21*t127)*2.0+t5*t40*t127*2.0+t2*t4*t18*t130*2.0+t2*t4*t21*t127*2.0)*(1.0/2.0)-t33*(t81*t81)*(1.0/4.0);
    hess(3, 4) = t28*(t153-t16*(t178+t3*t18*t132-t3*t21*t49*t54)*2.0+t4*t18*t49*2.0+t2*t18*t49*t58*2.0-t2*t4*t18*t132*2.0+t2*t4*t21*t49*t54*2.0)*(1.0/2.0)-t33*t81*t86*(1.0/4.0);
    hess(3, 5) = t28*(t175+t176+t177+t52*t64*2.0-t2*t4*t18*t136*2.0)*(1.0/2.0)-t33*t81*t91*(1.0/4.0);
    hess(4, 0) = t28*(t17*t140*2.0+t4*(t116+t138+t5*t8*t21*(t41-t55))*2.0+t3*t16*t18*t94*2.0+t2*t6*t8*t40*(t41-t55)*2.0+t2*t3*t8*t16*t21*(t41-t55)*2.0)*(1.0/2.0)-t24*t33*t146*(1.0/4.0);
    hess(4, 1) = t28*(t143-t144+t4*(t141-t142+t2*t3*t8*t21*t54)*2.0+t6*t8*t16*t21*t54*2.0+t2*t3*t8*t18*t140*2.0)*(1.0/2.0)-t33*t73*t146*(1.0/4.0);
    hess(4, 2) = t28*(t149+t150+t151+t43*t140*2.0-t3*t16*t18*t112*2.0)*(1.0/2.0)-t33*t78*t146*(1.0/4.0);
    hess(4, 3) = t28*(t153+t4*(t152-t2*t18*t132+t2*t21*t49*t54)*2.0-t16*t18*(t41-t55)*2.0-t3*t16*t18*t132*2.0+t2*t18*t49*t140*2.0+t3*t16*t21*t49*t54*2.0)*(1.0/2.0)-t33*t81*t146*(1.0/4.0);
    hess(4, 4) = t33*(t154*t154)*(-1.0/4.0)+t28*(t4*(t18*(t41-t55)*2.0+t2*t21*(t155*t155)+t2*t18*t158)*2.0+(t140*t140)*2.0+t6*t40*(t156*t156)*2.0+t3*t16*t18*t158*2.0+t3*t16*t21*(t159*t159)*2.0)*(1.0/2.0);
    hess(4, 5) = t28*(t179+t180+t4*(t160+t2*t18*t161+t2*t21*(t41-t55)*(t50-t61))*2.0+t16*t18*(t41-t55)*2.0+t3*t16*t18*t161*2.0+t3*t16*t21*(t41-t55)*(t50-t61)*2.0)*(1.0/2.0)-t33*t154*t169*(1.0/4.0);
    hess(5, 0) = t28*(t166-t167+t16*(t116-t165+t2*t3*t8*t21*t60)*2.0+t4*t5*t8*t21*t60*2.0+t2*t3*t8*t18*t164*2.0)*(1.0/2.0)-t24*t33*t169*(1.0/4.0);
    hess(5, 1) = t28*(t32*t164*2.0+t16*(t141+t171+t6*t8*t21*(t50-t61))*2.0+t2*t4*t18*t121*2.0+t3*t5*t8*t40*(t50-t61)*2.0+t2*t3*t4*t8*t21*(t50-t61)*2.0)*(1.0/2.0)-t33*t73*t169*(1.0/4.0);
    hess(5, 2) = t28*(t172+t16*(t124-t3*t18*t115+t3*t21*t39*t60)*2.0-t4*t18*(t50-t61)*2.0-t2*t4*t18*t115*2.0+t3*t18*t39*t164*2.0+t2*t4*t21*t39*t60*2.0)*(1.0/2.0)-t33*t78*t169*(1.0/4.0);
    hess(5, 3) = t28*(t175+t176+t177+t52*t164*2.0-t2*t4*t18*t136*2.0)*(1.0/2.0)-t33*t81*t169*(1.0/4.0);
    hess(5, 4) = t28*(t179+t180+t16*(t178+t3*t18*t161+t3*t21*(t41-t55)*(t50-t61))*2.0+t4*t18*(t50-t61)*2.0+t2*t4*t18*t161*2.0+t2*t4*t21*(t41-t55)*(t50-t61)*2.0)*(1.0/2.0)-t33*t154*t169*(1.0/4.0);
    hess(5, 5) = t33*(t181*t181)*(-1.0/4.0)+t28*(t16*(t18*(t50-t61)*2.0+t3*t21*(t182*t182)+t3*t18*t185)*2.0+(t164*t164)*2.0+t5*t40*(t183*t183)*2.0+t2*t4*t18*t185*2.0+t2*t4*t21*(t186*t186)*2.0)*(1.0/2.0);

  }
  
  template<typename T>
  inline void hesssoftpointedgedist(
                                    const Eigen::Matrix<T, 3, 1>& x0,
                                    const Eigen::Matrix<T, 3, 1>& x1,
                                    const Eigen::Matrix<T, 3, 1>& x2,
                                    Matrix9s& hess,
                                    const T& h)
  {
    
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t4 = x1(2)-x2(2);
    const scalar t6 = t2*t2;
    const scalar t7 = t3*t3;
    const scalar t8 = t4*t4;
    const scalar t9 = t6+t7+t8;
    const scalar t10 = 1.0/t9;
    const scalar t11 = x0(0)-x1(0);
    const scalar t12 = t2*t11;
    const scalar t13 = x0(1)-x1(1);
    const scalar t14 = t3*t13;
    const scalar t15 = x0(2)-x1(2);
    const scalar t16 = t4*t15;
    const scalar t17 = t12+t14+t16;
    const scalar t18 = t10*t17;
    const scalar t19 = softclamp(-t18, 0.0, 1.0, h);
    const scalar t22 = t2*t19;
    const scalar t5 = t22+x0(0)-x1(0);
    const scalar t24 = t3*t19;
    const scalar t20 = t24+x0(1)-x1(1);
    const scalar t25 = t4*t19;
    const scalar t21 = t25+x0(2)-x1(2);
    const scalar t23 = gradsoftclamp(-t18, 0.0, 1.0, h);
    const scalar t31 = t6*t10*t23;
    const scalar t32 = t31-1.0;
    const scalar t38 = t5*t32*2.0;
    const scalar t39 = t2*t3*t10*t20*t23*2.0;
    const scalar t40 = t2*t4*t10*t21*t23*2.0;
    const scalar t26 = t38+t39+t40;
    const scalar t27 = t5*t5;
    const scalar t28 = t20*t20;
    const scalar t29 = t21*t21;
    const scalar t30 = t27+t28+t29;
    const scalar t33 = 1.0/(t9*t9);
    const scalar t34 = t23*t23;
    const scalar t35 = hesssoftclamp(-t18, 0.0, 1.0, h);
    const scalar t36 = 1.0/sqrt(t30);
    const scalar t37 = 1.0/pow(t30,3.0/2.0);
    const scalar t41 = t7*t10*t23;
    const scalar t42 = t41-1.0;
    const scalar t43 = t8*t10*t23;
    const scalar t44 = t43-1.0;
    const scalar t45 = x1(0)*2.0;
    const scalar t47 = x2(0)*2.0;
    const scalar t46 = t45-t47;
    const scalar t48 = -t45+x0(0)+x2(0);
    const scalar t49 = t10*t48;
    const scalar t51 = t17*t33*t46;
    const scalar t50 = t49-t51;
    const scalar t52 = x1(1)*2.0;
    const scalar t54 = x2(1)*2.0;
    const scalar t53 = t52-t54;
    const scalar t55 = -t52+x0(1)+x2(1);
    const scalar t56 = t10*t55;
    const scalar t58 = t17*t33*t53;
    const scalar t57 = t56-t58;
    const scalar t59 = t3*t23*t57;
    const scalar t60 = -t19+t59+1.0;
    const scalar t61 = x1(2)*2.0;
    const scalar t63 = x2(2)*2.0;
    const scalar t62 = t61-t63;
    const scalar t64 = -t61+x0(2)+x2(2);
    const scalar t65 = t10*t64;
    const scalar t67 = t17*t33*t62;
    const scalar t66 = t65-t67;
    const scalar t68 = t4*t23*t66;
    const scalar t69 = -t19+t68+1.0;
    const scalar t70 = t6*t23*t33*t46;
    const scalar t72 = t10*t11;
    const scalar t71 = t51-t72;
    const scalar t73 = t2*t3*t20*t23*t33*t46*2.0;
    const scalar t74 = t2*t4*t21*t23*t33*t46*2.0;
    const scalar t75 = t2*t23*t71;
    const scalar t76 = t19+t75;
    const scalar t77 = t6*t23*t33*t53;
    const scalar t79 = t10*t13;
    const scalar t78 = t58-t79;
    const scalar t80 = t2*t3*t20*t23*t33*t53*2.0;
    const scalar t81 = t2*t4*t21*t23*t33*t53*2.0;
    const scalar t82 = t3*t23*t78;
    const scalar t83 = t19+t82;
    const scalar t84 = t6*t23*t33*t62;
    const scalar t86 = t10*t15;
    const scalar t85 = t67-t86;
    const scalar t87 = t2*t3*t20*t23*t33*t62*2.0;
    const scalar t88 = t2*t4*t21*t23*t33*t62*2.0;
    const scalar t89 = t4*t23*t85;
    const scalar t90 = t19+t89;
    const scalar t91 = t2*t3*t10*t23*t32*2.0;
    const scalar t92 = t2*t3*t10*t23*t42*2.0;
    const scalar t93 = t2*t3*t8*t33*t34*2.0;
    const scalar t94 = t3*t5*t6*t33*t35*2.0;
    const scalar t95 = t2*t7*t20*t33*t35*2.0;
    const scalar t96 = t2*t3*t4*t21*t33*t35*2.0;
    const scalar t97 = t91+t92+t93+t94+t95+t96;
    const scalar t98 = t36*t97*(1.0/2.0);
    const scalar t99 = t20*t42*2.0;
    const scalar t100 = t2*t3*t5*t10*t23*2.0;
    const scalar t101 = t3*t4*t10*t21*t23*2.0;
    const scalar t102 = t99+t100+t101;
    const scalar t103 = t98-t26*t37*t102*(1.0/4.0);
    const scalar t104 = t6*t7*t33*t34*2.0;
    const scalar t105 = t21*t44*2.0;
    const scalar t106 = t2*t4*t5*t10*t23*2.0;
    const scalar t107 = t3*t4*t10*t20*t23*2.0;
    const scalar t108 = t105+t106+t107;
    const scalar t109 = t2*t23*t50;
    const scalar t110 = -t19+t109+1.0;
    const scalar t111 = t5*t110*2.0;
    const scalar t112 = t3*t20*t23*t50*2.0;
    const scalar t113 = t4*t21*t23*t50*2.0;
    const scalar t114 = t111+t112+t113;
    const scalar t115 = t4*t10*t21*t23*2.0;
    const scalar t116 = t20*t60*2.0;
    const scalar t117 = t2*t5*t23*t57*2.0;
    const scalar t118 = t4*t21*t23*t57*2.0;
    const scalar t119 = t116+t117+t118;
    const scalar t120 = t21*t69*2.0;
    const scalar t121 = t2*t5*t23*t66*2.0;
    const scalar t122 = t3*t20*t23*t66*2.0;
    const scalar t123 = t120+t121+t122;
    const scalar t124 = t7*t23*t33*t46;
    const scalar t125 = t2*t3*t5*t23*t33*t46*2.0;
    const scalar t126 = t3*t4*t21*t23*t33*t46*2.0;
    const scalar t127 = t5*t76*2.0;
    const scalar t128 = t3*t20*t23*t71*2.0;
    const scalar t129 = t4*t21*t23*t71*2.0;
    const scalar t130 = t127+t128+t129;
    const scalar t131 = t7*t23*t33*t53;
    const scalar t132 = t2*t3*t5*t23*t33*t53*2.0;
    const scalar t133 = t3*t4*t21*t23*t33*t53*2.0;
    const scalar t134 = t20*t83*2.0;
    const scalar t135 = t2*t5*t23*t78*2.0;
    const scalar t136 = t4*t21*t23*t78*2.0;
    const scalar t137 = t134+t135+t136;
    const scalar t138 = t7*t23*t33*t62;
    const scalar t139 = t2*t3*t5*t23*t33*t62*2.0;
    const scalar t140 = t3*t4*t21*t23*t33*t62*2.0;
    const scalar t141 = t21*t90*2.0;
    const scalar t142 = t2*t5*t23*t85*2.0;
    const scalar t143 = t3*t20*t23*t85*2.0;
    const scalar t144 = t141+t142+t143;
    const scalar t145 = t2*t4*t10*t23*t32*2.0;
    const scalar t146 = t2*t4*t10*t23*t44*2.0;
    const scalar t147 = t2*t4*t7*t33*t34*2.0;
    const scalar t148 = t4*t5*t6*t33*t35*2.0;
    const scalar t149 = t2*t8*t21*t33*t35*2.0;
    const scalar t150 = t2*t3*t4*t20*t33*t35*2.0;
    const scalar t151 = t145+t146+t147+t148+t149+t150;
    const scalar t152 = t36*t151*(1.0/2.0);
    const scalar t153 = t152-t26*t37*t108*(1.0/4.0);
    const scalar t154 = t3*t4*t10*t23*t42*2.0;
    const scalar t155 = t3*t4*t10*t23*t44*2.0;
    const scalar t156 = t3*t4*t6*t33*t34*2.0;
    const scalar t157 = t4*t7*t20*t33*t35*2.0;
    const scalar t158 = t3*t8*t21*t33*t35*2.0;
    const scalar t159 = t2*t3*t4*t5*t33*t35*2.0;
    const scalar t160 = t154+t155+t156+t157+t158+t159;
    const scalar t161 = t36*t160*(1.0/2.0);
    const scalar t162 = t161-t37*t102*t108*(1.0/4.0);
    const scalar t163 = t6*t8*t33*t34*2.0;
    const scalar t164 = t7*t8*t33*t34*2.0;
    const scalar t165 = t2*t5*t10*t23*2.0;
    const scalar t166 = t3*t10*t20*t23*2.0;
    const scalar t167 = t8*t23*t33*t46;
    const scalar t168 = t2*t4*t5*t23*t33*t46*2.0;
    const scalar t169 = t3*t4*t20*t23*t33*t46*2.0;
    const scalar t170 = t8*t23*t33*t53;
    const scalar t171 = t2*t4*t5*t23*t33*t53*2.0;
    const scalar t172 = t3*t4*t20*t23*t33*t53*2.0;
    const scalar t173 = t8*t23*t33*t62;
    const scalar t174 = t2*t4*t5*t23*t33*t62*2.0;
    const scalar t175 = t3*t4*t20*t23*t33*t62*2.0;
    const scalar t176 = t6*t10*t35*t50;
    const scalar t178 = t2*t33*t46;
    const scalar t177 = t10-t178;
    const scalar t179 = t2*t7*t10*t34*t50*2.0;
    const scalar t180 = t2*t8*t10*t34*t50*2.0;
    const scalar t181 = t2*t3*t10*t20*t35*t50*2.0;
    const scalar t182 = t2*t4*t10*t21*t35*t50*2.0;
    const scalar t183 = t3*t23*t42*t50*2.0;
    const scalar t184 = t2*t3*t10*t23*t110*2.0;
    const scalar t185 = t3*t8*t10*t34*t50*2.0;
    const scalar t186 = t3*t4*t10*t21*t35*t50*2.0;
    const scalar t187 = t4*t23*t44*t50*2.0;
    const scalar t188 = t2*t4*t10*t23*t110*2.0;
    const scalar t189 = t4*t7*t10*t34*t50*2.0;
    const scalar t190 = t3*t4*t10*t20*t35*t50*2.0;
    const scalar t191 = t50*t50;
    const scalar t192 = t17*t33*2.0;
    const scalar t193 = t10*2.0;
    const scalar t194 = t46*t46;
    const scalar t195 = 1.0/(t9*t9*t9);
    const scalar t196 = t33*t46*t48*2.0;
    const scalar t198 = t17*t194*t195*2.0;
    const scalar t197 = t192+t193+t196-t198;
    const scalar t199 = t33*t48*t53;
    const scalar t200 = t33*t46*t55;
    const scalar t202 = t17*t46*t53*t195*2.0;
    const scalar t201 = t199+t200-t202;
    const scalar t203 = t33*t48*t62;
    const scalar t204 = t33*t46*t64;
    const scalar t206 = t17*t46*t62*t195*2.0;
    const scalar t205 = t203+t204-t206;
    const scalar t207 = t11*t33*t46;
    const scalar t208 = t33*t46*t48;
    const scalar t209 = t10+t192-t198+t207+t208;
    const scalar t210 = t13*t33*t46;
    const scalar t211 = t199-t202+t210;
    const scalar t212 = t15*t33*t46;
    const scalar t213 = t203-t206+t212;
    const scalar t214 = t2*t10*t23;
    const scalar t215 = t2*t23*t32*t57*2.0;
    const scalar t216 = t2*t3*t10*t23*t60*2.0;
    const scalar t217 = t2*t8*t10*t34*t57*2.0;
    const scalar t218 = t2*t4*t10*t21*t35*t57*2.0;
    const scalar t219 = t42*t60*2.0;
    const scalar t220 = t7*t10*t35*t57;
    const scalar t222 = t3*t33*t53;
    const scalar t221 = t10-t222;
    const scalar t223 = t3*t6*t10*t34*t57*2.0;
    const scalar t224 = t3*t8*t10*t34*t57*2.0;
    const scalar t225 = t2*t3*t5*t10*t35*t57*2.0;
    const scalar t226 = t3*t4*t10*t21*t35*t57*2.0;
    const scalar t227 = t4*t23*t44*t57*2.0;
    const scalar t228 = t3*t4*t10*t23*t60*2.0;
    const scalar t229 = t4*t6*t10*t34*t57*2.0;
    const scalar t230 = t2*t4*t5*t10*t35*t57*2.0;
    const scalar t231 = t23*t50;
    const scalar t232 = t4*t21*t23*t201*2.0;
    const scalar t233 = t2*t23*t57*t110*2.0;
    const scalar t234 = t3*t23*t50*t60*2.0;
    const scalar t235 = t8*t34*t50*t57*2.0;
    const scalar t236 = t4*t21*t35*t50*t57*2.0;
    const scalar t237 = t57*t57;
    const scalar t238 = t53*t53;
    const scalar t239 = t33*t53*t55*2.0;
    const scalar t241 = t17*t195*t238*2.0;
    const scalar t240 = t192+t193+t239-t241;
    const scalar t242 = t33*t55*t62;
    const scalar t243 = t33*t53*t64;
    const scalar t245 = t17*t53*t62*t195*2.0;
    const scalar t244 = t242+t243-t245;
    const scalar t246 = t11*t33*t53;
    const scalar t247 = t200-t202+t246;
    const scalar t248 = t13*t33*t53;
    const scalar t249 = t33*t53*t55;
    const scalar t250 = t10+t192-t241+t248+t249;
    const scalar t251 = t23*(t67-t86);
    const scalar t252 = t15*t33*t53;
    const scalar t253 = t242-t245+t252;
    const scalar t254 = t2*t23*t32*t66*2.0;
    const scalar t255 = t2*t4*t10*t23*t69*2.0;
    const scalar t256 = t2*t7*t10*t34*t66*2.0;
    const scalar t257 = t2*t3*t10*t20*t35*t66*2.0;
    const scalar t258 = t3*t10*t23;
    const scalar t259 = t3*t23*t42*t66*2.0;
    const scalar t260 = t3*t4*t10*t23*t69*2.0;
    const scalar t261 = t3*t6*t10*t34*t66*2.0;
    const scalar t262 = t2*t3*t5*t10*t35*t66*2.0;
    const scalar t263 = t44*t69*2.0;
    const scalar t264 = t8*t10*t35*t66;
    const scalar t266 = t4*t33*t62;
    const scalar t265 = t10-t266;
    const scalar t267 = t4*t6*t10*t34*t66*2.0;
    const scalar t268 = t4*t7*t10*t34*t66*2.0;
    const scalar t269 = t2*t4*t5*t10*t35*t66*2.0;
    const scalar t270 = t3*t4*t10*t20*t35*t66*2.0;
    const scalar t271 = t3*t20*t23*t205*2.0;
    const scalar t272 = t2*t23*t66*t110*2.0;
    const scalar t273 = t4*t23*t50*t69*2.0;
    const scalar t274 = t7*t34*t50*t66*2.0;
    const scalar t275 = t3*t20*t35*t50*t66*2.0;
    const scalar t276 = t23*t57;
    const scalar t277 = t2*t5*t23*t244*2.0;
    const scalar t278 = t3*t23*t60*t66*2.0;
    const scalar t279 = t4*t23*t57*t69*2.0;
    const scalar t280 = t6*t34*t57*t66*2.0;
    const scalar t281 = t2*t5*t35*t57*t66*2.0;
    const scalar t282 = t66*t66;
    const scalar t283 = t62*t62;
    const scalar t284 = t33*t62*t64*2.0;
    const scalar t286 = t17*t195*t283*2.0;
    const scalar t285 = t192+t193+t284-t286;
    const scalar t287 = t23*(t51-t72);
    const scalar t288 = t11*t33*t62;
    const scalar t289 = t204-t206+t288;
    const scalar t290 = t23*(t58-t79);
    const scalar t291 = t13*t33*t62;
    const scalar t292 = t243-t245+t291;
    const scalar t293 = t15*t33*t62;
    const scalar t294 = t33*t62*t64;
    const scalar t295 = t10+t192-t286+t293+t294;
    const scalar t296 = t2*t23*t177;
    const scalar t297 = t2*t23*(t51-t72);
    const scalar t298 = t19+t297;
    const scalar t299 = t2*t3*t23*t33*t46;
    const scalar t300 = t3*t23*t42*t71*2.0;
    const scalar t301 = t7*t20*t23*t33*t46*2.0;
    const scalar t302 = t3*t8*t10*t34*t71*2.0;
    const scalar t303 = t3*t4*t10*t21*t35*t71*2.0;
    const scalar t304 = t5*t298*2.0;
    const scalar t305 = t128+t129+t304;
    const scalar t306 = t4*t10*t23;
    const scalar t307 = t2*t4*t23*t33*t46;
    const scalar t308 = t4*t23*t44*t71*2.0;
    const scalar t309 = t8*t21*t23*t33*t46*2.0;
    const scalar t310 = t4*t7*t10*t34*t71*2.0;
    const scalar t311 = t3*t4*t10*t20*t35*t71*2.0;
    const scalar t312 = t2*t35*t50*t71;
    const scalar t313 = t231+t312-t23*t71-t2*t23*t209;
    const scalar t314 = t5*t313*2.0;
    const scalar t315 = t7*t34*t50*t71*2.0;
    const scalar t316 = t8*t34*t50*t71*2.0;
    const scalar t317 = t3*t20*t35*t50*t71*2.0;
    const scalar t318 = t4*t21*t35*t50*t71*2.0;
    const scalar t319 = t3*t23*t60*t71*2.0;
    const scalar t320 = t8*t34*t57*t71*2.0;
    const scalar t321 = t4*t21*t35*t57*t71*2.0;
    const scalar t322 = t23*t66;
    const scalar t323 = t4*t23*t69*t71*2.0;
    const scalar t324 = t7*t34*t66*t71*2.0;
    const scalar t325 = t3*t20*t35*t66*t71*2.0;
    const scalar t326 = t128+t129+t304;
    const scalar t327 = t51-t72;
    const scalar t328 = t51-t72;
    const scalar t329 = t51-t72;
    const scalar t330 = t11*t33*t46*2.0;
    const scalar t331 = t192-t198+t330;
    const scalar t332 = t51-t72;
    const scalar t333 = t51-t72;
    const scalar t334 = -t202+t210+t246;
    const scalar t335 = t3*t23*(t58-t79);
    const scalar t336 = t19+t335;
    const scalar t337 = -t206+t212+t288;
    const scalar t338 = t4*t23*(t67-t86);
    const scalar t339 = t19+t338;
    const scalar t340 = t2*t3*t23*t33*t53;
    const scalar t341 = t2*t23*t32*t78*2.0;
    const scalar t342 = t5*t6*t23*t33*t53*2.0;
    const scalar t343 = t2*t8*t10*t34*t78*2.0;
    const scalar t344 = t2*t4*t10*t21*t35*t78*2.0;
    const scalar t345 = t20*t336*2.0;
    const scalar t346 = t135+t136+t345;
    const scalar t347 = t3*t23*t221;
    const scalar t348 = t3*t4*t23*t33*t53;
    const scalar t349 = t4*t23*t44*t78*2.0;
    const scalar t350 = t8*t21*t23*t33*t53*2.0;
    const scalar t351 = t4*t6*t10*t34*t78*2.0;
    const scalar t352 = t2*t4*t5*t10*t35*t78*2.0;
    const scalar t353 = t2*t23*t78*t110*2.0;
    const scalar t354 = t8*t34*t50*t78*2.0;
    const scalar t355 = t4*t21*t35*t50*t78*2.0;
    const scalar t356 = t3*t35*t57*t78;
    const scalar t357 = t276+t356-t23*t78-t3*t23*t250;
    const scalar t358 = t20*t357*2.0;
    const scalar t359 = t6*t34*t57*t78*2.0;
    const scalar t360 = t8*t34*t57*t78*2.0;
    const scalar t361 = t2*t5*t35*t57*t78*2.0;
    const scalar t362 = t4*t21*t35*t57*t78*2.0;
    const scalar t363 = t4*t23*t69*t78*2.0;
    const scalar t364 = t6*t34*t66*t78*2.0;
    const scalar t365 = t2*t5*t35*t66*t78*2.0;
    const scalar t366 = t4*t21*t23*t334*2.0;
    const scalar t367 = t2*t23*t298*(t58-t79)*2.0;
    const scalar t368 = t3*t23*t336*(t51-t72)*2.0;
    const scalar t369 = t8*t34*(t51-t72)*(t58-t79)*2.0;
    const scalar t370 = t4*t21*t35*(t51-t72)*(t58-t79)*2.0;
    const scalar t371 = t135+t136+t345;
    const scalar t372 = t58-t79;
    const scalar t373 = t58-t79;
    const scalar t374 = t58-t79;
    const scalar t375 = t13*t33*t53*2.0;
    const scalar t376 = t192-t241+t375;
    const scalar t377 = t58-t79;
    const scalar t378 = t58-t79;
    const scalar t379 = -t245+t252+t291;
    const scalar t380 = t21*t339*2.0;
    const scalar t381 = t142+t143+t380;
    const scalar t382 = t2*t4*t23*t33*t62;
    const scalar t383 = t2*t23*t32*t85*2.0;
    const scalar t384 = t5*t6*t23*t33*t62*2.0;
    const scalar t385 = t2*t7*t10*t34*t85*2.0;
    const scalar t386 = t2*t3*t10*t20*t35*t85*2.0;
    const scalar t387 = t3*t4*t23*t33*t62;
    const scalar t388 = t3*t23*t42*t85*2.0;
    const scalar t389 = t7*t20*t23*t33*t62*2.0;
    const scalar t390 = t3*t6*t10*t34*t85*2.0;
    const scalar t391 = t2*t3*t5*t10*t35*t85*2.0;
    const scalar t392 = t4*t23*t265;
    const scalar t393 = t2*t23*t85*t110*2.0;
    const scalar t394 = t7*t34*t50*t85*2.0;
    const scalar t395 = t3*t20*t35*t50*t85*2.0;
    const scalar t396 = t3*t23*t60*t85*2.0;
    const scalar t397 = t6*t34*t57*t85*2.0;
    const scalar t398 = t2*t5*t35*t57*t85*2.0;
    const scalar t399 = t4*t35*t66*t85;
    const scalar t400 = t322+t399-t23*t85-t4*t23*t295;
    const scalar t401 = t21*t400*2.0;
    const scalar t402 = t6*t34*t66*t85*2.0;
    const scalar t403 = t7*t34*t66*t85*2.0;
    const scalar t404 = t2*t5*t35*t66*t85*2.0;
    const scalar t405 = t3*t20*t35*t66*t85*2.0;
    const scalar t406 = t3*t20*t23*t337*2.0;
    const scalar t407 = t2*t23*t298*(t67-t86)*2.0;
    const scalar t408 = t4*t23*t339*(t51-t72)*2.0;
    const scalar t409 = t7*t34*(t51-t72)*(t67-t86)*2.0;
    const scalar t410 = t3*t20*t35*(t51-t72)*(t67-t86)*2.0;
    const scalar t411 = t2*t5*t23*t379*2.0;
    const scalar t412 = t3*t23*t336*(t67-t86)*2.0;
    const scalar t413 = t4*t23*t339*(t58-t79)*2.0;
    const scalar t414 = t6*t34*(t58-t79)*(t67-t86)*2.0;
    const scalar t415 = t2*t5*t35*(t58-t79)*(t67-t86)*2.0;
    const scalar t416 = t142+t143+t380;
    const scalar t417 = t67-t86;
    const scalar t418 = t67-t86;
    const scalar t419 = t67-t86;
    const scalar t420 = t15*t33*t62*2.0;
    const scalar t421 = t192-t286+t420;
    const scalar t422 = t67-t86;
    const scalar t423 = t67-t86;
    hess(0, 0) = t36*(t104+t163+(t32*t32)*2.0+t2*t5*t6*t33*t35*2.0+t3*t6*t20*t33*t35*2.0+t4*t6*t21*t33*t35*2.0)*(1.0/2.0)-(t26*t26)*t37*(1.0/4.0);
    hess(0, 1) = t103;
    hess(0, 2) = t153;
    hess(0, 3) = t36*(t73+t74+t179+t180+t181+t182+t32*(-t19+t2*t23*(t10*(x0(0)-x1(0)*2.0+x2(0))-t17*t33*t46)+1.0)*2.0+t5*(t70+t176-t10*t23*t46)*2.0-t3*t10*t20*t23*2.0-t4*t10*t21*t23*2.0)*(1.0/2.0)-t26*t37*t114*(1.0/4.0);
    hess(0, 4) = t36*(t80+t81+t215+t216+t217+t218+t5*(t77+t6*t10*t35*(t10*(x0(1)-x1(1)*2.0+x2(1))-t17*t33*t53))*2.0-t2*t10*t20*t23*2.0+t2*t3*t10*t20*t35*t57*2.0)*(1.0/2.0)-t26*t37*t119*(1.0/4.0);
    hess(0, 5) = t36*(t87+t88+t254+t255+t256+t257+t5*(t84+t6*t10*t35*(t10*(x0(2)-x1(2)*2.0+x2(2))-t17*t33*t62))*2.0-t2*t10*t21*t23*2.0+t2*t4*t10*t21*t35*t66*2.0)*(1.0/2.0)-t26*t37*t123*(1.0/4.0);
    hess(0, 6) = t36*(-t73-t74+t115+t166+t32*t76*2.0+t5*(-t70+t10*t23*t46+t6*t10*t35*t71)*2.0+t2*t7*t10*t34*t71*2.0+t2*t8*t10*t34*t71*2.0+t2*t3*t10*t20*t35*t71*2.0+t2*t4*t10*t21*t35*t71*2.0)*(1.0/2.0)-t26*t37*t130*(1.0/4.0);
    hess(0, 7) = t36*(-t80-t81+t341+t343+t344-t5*(t77-t6*t10*t35*t78)*2.0+t2*t10*t20*t23*2.0+t2*t3*t10*t23*t83*2.0+t2*t3*t10*t20*t35*t78*2.0)*(1.0/2.0)-t26*t37*t137*(1.0/4.0);
    hess(0, 8) = t36*(-t87-t88+t383+t385+t386-t5*(t84-t6*t10*t35*t85)*2.0+t2*t10*t21*t23*2.0+t2*t4*t10*t23*t90*2.0+t2*t4*t10*t21*t35*t85*2.0)*(1.0/2.0)-t26*t37*t144*(1.0/4.0);
    hess(1, 0) = t103;
    hess(1, 1) = t36*(t104+t164+(t42*t42)*2.0+t2*t5*t7*t33*t35*2.0+t3*t7*t20*t33*t35*2.0+t4*t7*t21*t33*t35*2.0)*(1.0/2.0)-t37*(t102*t102)*(1.0/4.0);
    hess(1, 2) = t162;
    hess(1, 3) = t36*(t125+t126+t183+t184+t185+t186+t20*(t124+t7*t10*t35*t50)*2.0-t3*t5*t10*t23*2.0+t2*t3*t5*t10*t35*t50*2.0)*(1.0/2.0)-t37*t102*t114*(1.0/4.0);
    hess(1, 4) = t36*(-t115+t132+t133+t219+t223+t224+t225+t226+t20*(t131+t220-t10*t23*t53)*2.0-t2*t5*t10*t23*2.0)*(1.0/2.0)-t37*t102*t119*(1.0/4.0);
    hess(1, 5) = t36*(t139+t140+t259+t260+t261+t262+t20*(t138+t7*t10*t35*t66)*2.0-t3*t10*t21*t23*2.0+t3*t4*t10*t21*t35*t66*2.0)*(1.0/2.0)-t37*t102*t123*(1.0/4.0);
    hess(1, 6) = t36*(-t125-t126+t300+t302+t303-t20*(t124-t7*t10*t35*t71)*2.0+t3*t5*t10*t23*2.0+t2*t3*t10*t23*t76*2.0+t2*t3*t5*t10*t35*t71*2.0)*(1.0/2.0)-t37*t102*t130*(1.0/4.0);
    hess(1, 7) = t36*(t115-t132-t133+t165+t42*t83*2.0+t20*(-t131+t10*t23*t53+t7*t10*t35*t78)*2.0+t3*t6*t10*t34*t78*2.0+t3*t8*t10*t34*t78*2.0+t2*t3*t5*t10*t35*t78*2.0+t3*t4*t10*t21*t35*t78*2.0)*(1.0/2.0)-t37*t102*t137*(1.0/4.0);
    hess(1, 8) = t36*(-t139-t140+t388+t390+t391-t20*(t138-t7*t10*t35*t85)*2.0+t3*t10*t21*t23*2.0+t3*t4*t10*t23*t90*2.0+t3*t4*t10*t21*t35*t85*2.0)*(1.0/2.0)-t37*t102*t144*(1.0/4.0);
    hess(2, 0) = t153;
    hess(2, 1) = t162;
    hess(2, 2) = t36*(t163+t164+(t44*t44)*2.0+t2*t5*t8*t33*t35*2.0+t3*t8*t20*t33*t35*2.0+t4*t8*t21*t33*t35*2.0)*(1.0/2.0)-t37*(t108*t108)*(1.0/4.0);
    hess(2, 3) = t36*(t168+t169+t187+t188+t189+t190+t21*(t167+t8*t10*t35*t50)*2.0-t4*t5*t10*t23*2.0+t2*t4*t5*t10*t35*t50*2.0)*(1.0/2.0)-t37*t108*t114*(1.0/4.0);
    hess(2, 4) = t36*(t171+t172+t227+t228+t229+t230+t21*(t170+t8*t10*t35*t57)*2.0-t4*t10*t20*t23*2.0+t3*t4*t10*t20*t35*t57*2.0)*(1.0/2.0)-t37*t108*t119*(1.0/4.0);
    hess(2, 5) = t36*(-t165-t166+t174+t175+t263+t267+t268+t269+t270+t21*(t173+t264-t10*t23*t62)*2.0)*(1.0/2.0)-t37*t108*t123*(1.0/4.0);
    hess(2, 6) = t36*(-t168-t169+t308+t310+t311-t21*(t167-t8*t10*t35*t71)*2.0+t4*t5*t10*t23*2.0+t2*t4*t10*t23*t76*2.0+t2*t4*t5*t10*t35*t71*2.0)*(1.0/2.0)-t37*t108*t130*(1.0/4.0);
    hess(2, 7) = t36*(-t171-t172+t349+t351+t352-t21*(t170-t8*t10*t35*t78)*2.0+t4*t10*t20*t23*2.0+t3*t4*t10*t23*t83*2.0+t3*t4*t10*t20*t35*t78*2.0)*(1.0/2.0)-t37*t108*t137*(1.0/4.0);
    hess(2, 8) = t36*(t165+t166-t174-t175+t44*t90*2.0+t21*(-t173+t10*t23*t62+t8*t10*t35*t85)*2.0+t4*t6*t10*t34*t85*2.0+t4*t7*t10*t34*t85*2.0+t2*t4*t5*t10*t35*t85*2.0+t3*t4*t10*t20*t35*t85*2.0)*(1.0/2.0)-t37*t108*t144*(1.0/4.0);
    hess(3, 0) = t36*(t179+t180+t181+t182+t32*t110*2.0-t5*(-t176+t214+t296)*2.0-t3*t20*t23*t177*2.0-t4*t21*t23*t177*2.0)*(1.0/2.0)-t26*t37*t114*(1.0/4.0);
    hess(3, 1) = t36*(t126+t183+t184+t185+t186+t301+t5*(t299-t3*t10*t23+t2*t3*t10*t35*t50)*2.0+t7*t10*t20*t35*t50*2.0)*(1.0/2.0)-t37*t102*t114*(1.0/4.0);
    hess(3, 2) = t36*(t169+t187+t188+t189+t190+t309+t5*(t307-t4*t10*t23+t2*t4*t10*t35*t50)*2.0+t8*t10*t21*t35*t50*2.0)*(1.0/2.0)-t37*t108*t114*(1.0/4.0);
    hess(3, 3) = t36*((t110*t110)*2.0+t5*(t23*t50*-2.0+t2*t23*t197+t2*t35*t191)*2.0+t7*t34*t191*2.0+t8*t34*t191*2.0+t3*t20*t23*t197*2.0+t4*t21*t23*t197*2.0+t3*t20*t35*t191*2.0+t4*t21*t35*t191*2.0)*(1.0/2.0)-t37*(t114*t114)*(1.0/4.0);
    hess(3, 4) = t36*(t232+t233+t234+t235+t236+t5*(-t23*t57+t2*t23*t201+t2*t35*t50*t57)*2.0-t20*t23*t50*2.0+t3*t20*t23*t201*2.0+t3*t20*t35*t50*t57*2.0)*(1.0/2.0)-t37*t114*t119*(1.0/4.0);
    hess(3, 5) = t36*(t271+t272+t273+t274+t275+t5*(-t23*t66+t2*t23*t205+t2*t35*t50*t66)*2.0-t21*t23*t50*2.0+t4*t21*t23*t205*2.0+t4*t21*t35*t50*t66*2.0)*(1.0/2.0)-t37*t114*t123*(1.0/4.0);
    hess(3, 6) = t36*(t314+t315+t316+t317+t318+t76*t110*2.0-t3*t20*t23*t209*2.0-t4*t21*t23*t209*2.0)*(1.0/2.0)-t37*t114*t130*(1.0/4.0);
    hess(3, 7) = t36*(t353+t354+t355-t5*(t290+t2*t23*t211-t2*t35*t50*t78)*2.0+t20*t23*t50*2.0+t3*t23*t50*t83*2.0-t3*t20*t23*t211*2.0-t4*t21*t23*t211*2.0+t3*t20*t35*t50*t78*2.0)*(1.0/2.0)-t37*t114*t137*(1.0/4.0);
    hess(3, 8) = t36*(t393+t394+t395-t5*(t251+t2*t23*t213-t2*t35*t50*t85)*2.0+t21*t23*t50*2.0+t4*t23*t50*t90*2.0-t3*t20*t23*t213*2.0-t4*t21*t23*t213*2.0+t4*t21*t35*t50*t85*2.0)*(1.0/2.0)-t37*t114*t144*(1.0/4.0);
    hess(4, 0) = t36*(t81+t215+t216+t217+t218+t342+t20*(-t214+t340+t2*t3*t10*t35*t57)*2.0+t5*t6*t10*t35*t57*2.0)*(1.0/2.0)-t26*t37*t119*(1.0/4.0);
    hess(4, 1) = t36*(t219+t223+t224+t225+t226-t20*(-t220+t258+t347)*2.0-t2*t5*t23*t221*2.0-t4*t21*t23*t221*2.0)*(1.0/2.0)-t37*t102*t119*(1.0/4.0);
    hess(4, 2) = t36*(t171+t227+t228+t229+t230+t350+t20*(t348-t4*t10*t23+t3*t4*t10*t35*t57)*2.0+t8*t10*t21*t35*t57*2.0)*(1.0/2.0)-t37*t108*t119*(1.0/4.0);
    hess(4, 3) = t36*(t232+t233+t234+t235+t236+t20*(-t231+t3*t23*t201+t3*t35*t50*t57)*2.0-t5*t23*t57*2.0+t2*t5*t23*t201*2.0+t2*t5*t35*t50*t57*2.0)*(1.0/2.0)-t37*t114*t119*(1.0/4.0);
    hess(4, 4) = t36*((t60*t60)*2.0+t20*(t23*t57*-2.0+t3*t23*t240+t3*t35*t237)*2.0+t6*t34*t237*2.0+t8*t34*t237*2.0+t2*t5*t23*t240*2.0+t2*t5*t35*t237*2.0+t4*t21*t23*t240*2.0+t4*t21*t35*t237*2.0)*(1.0/2.0)-t37*(t119*t119)*(1.0/4.0);
    hess(4, 5) = t36*(t277+t278+t279+t280+t281+t20*(-t23*t66+t3*t23*t244+t3*t35*t57*t66)*2.0-t21*t23*t57*2.0+t4*t21*t23*t244*2.0+t4*t21*t35*t57*t66*2.0)*(1.0/2.0)-t37*t119*t123*(1.0/4.0);
    hess(4, 6) = t36*(t319+t320+t321-t20*(t287+t3*t23*t247-t3*t35*t57*t71)*2.0+t5*t23*t57*2.0+t2*t23*t57*t76*2.0-t2*t5*t23*t247*2.0-t4*t21*t23*t247*2.0+t2*t5*t35*t57*t71*2.0)*(1.0/2.0)-t37*t119*t130*(1.0/4.0);
    hess(4, 7) = t36*(t358+t359+t360+t361+t362+t60*t83*2.0-t2*t5*t23*t250*2.0-t4*t21*t23*t250*2.0)*(1.0/2.0)-t37*t119*t137*(1.0/4.0);
    hess(4, 8) = t36*(t396+t397+t398-t20*(t251+t3*t23*t253-t3*t35*t57*t85)*2.0+t21*t23*t57*2.0+t4*t23*t57*t90*2.0-t2*t5*t23*t253*2.0-t4*t21*t23*t253*2.0+t4*t21*t35*t57*t85*2.0)*(1.0/2.0)-t37*t119*t144*(1.0/4.0);
    hess(5, 0) = t36*(t87+t254+t255+t256+t257+t384+t21*(-t214+t382+t2*t4*t10*t35*t66)*2.0+t5*t6*t10*t35*t66*2.0)*(1.0/2.0)-t26*t37*t123*(1.0/4.0);
    hess(5, 1) = t36*(t139+t259+t260+t261+t262+t389+t21*(-t258+t387+t3*t4*t10*t35*t66)*2.0+t7*t10*t20*t35*t66*2.0)*(1.0/2.0)-t37*t102*t123*(1.0/4.0);
    hess(5, 2) = t36*(t263+t267+t268+t269+t270-t21*(-t264+t306+t392)*2.0-t2*t5*t23*t265*2.0-t3*t20*t23*t265*2.0)*(1.0/2.0)-t37*t108*t123*(1.0/4.0);
    hess(5, 3) = t36*(t271+t272+t273+t274+t275+t21*(-t231+t4*t23*t205+t4*t35*t50*t66)*2.0-t5*t23*t66*2.0+t2*t5*t23*t205*2.0+t2*t5*t35*t50*t66*2.0)*(1.0/2.0)-t37*t114*t123*(1.0/4.0);
    hess(5, 4) = t36*(t277+t278+t279+t280+t281+t21*(-t276+t4*t23*t244+t4*t35*t57*t66)*2.0-t20*t23*t66*2.0+t3*t20*t23*t244*2.0+t3*t20*t35*t57*t66*2.0)*(1.0/2.0)-t37*t119*t123*(1.0/4.0);
    hess(5, 5) = t36*((t69*t69)*2.0+t21*(t23*t66*-2.0+t4*t23*t285+t4*t35*t282)*2.0+t6*t34*t282*2.0+t7*t34*t282*2.0+t2*t5*t23*t285*2.0+t2*t5*t35*t282*2.0+t3*t20*t23*t285*2.0+t3*t20*t35*t282*2.0)*(1.0/2.0)-t37*(t123*t123)*(1.0/4.0);
    hess(5, 6) = t36*(t323+t324+t325-t21*(t287+t4*t23*t289-t4*t35*t66*t71)*2.0+t5*t23*t66*2.0+t2*t23*t66*t76*2.0-t2*t5*t23*t289*2.0-t3*t20*t23*t289*2.0+t2*t5*t35*t66*t71*2.0)*(1.0/2.0)-t37*t123*t130*(1.0/4.0);
    hess(5, 7) = t36*(t363+t364+t365-t21*(t290+t4*t23*t292-t4*t35*t66*t78)*2.0+t20*t23*t66*2.0+t3*t23*t66*t83*2.0-t2*t5*t23*t292*2.0-t3*t20*t23*t292*2.0+t3*t20*t35*t66*t78*2.0)*(1.0/2.0)-t37*t123*t137*(1.0/4.0);
    hess(5, 8) = t36*(t401+t402+t403+t404+t405+t69*t90*2.0-t2*t5*t23*t295*2.0-t3*t20*t23*t295*2.0)*(1.0/2.0)-t37*t123*t144*(1.0/4.0);
    hess(6, 0) = t36*(t32*t298*2.0+t5*(t214+t296+t6*t10*t35*(t51-t72))*2.0+t3*t20*t23*t177*2.0+t4*t21*t23*t177*2.0+t2*t7*t10*t34*(t51-t72)*2.0+t2*t8*t10*t34*(t51-t72)*2.0+t2*t3*t10*t20*t35*(t51-t72)*2.0+t2*t4*t10*t21*t35*(t51-t72)*2.0)*(1.0/2.0)-t26*t37*t305*(1.0/4.0);
    hess(6, 1) = t36*(-t126+t300-t301+t302+t303+t5*(t258-t299+t2*t3*t10*t35*t71)*2.0+t7*t10*t20*t35*t71*2.0+t2*t3*t10*t23*t298*2.0)*(1.0/2.0)-t37*t102*t305*(1.0/4.0);
    hess(6, 2) = t36*(-t169+t308-t309+t310+t311+t5*(t306-t307+t2*t4*t10*t35*t71)*2.0+t8*t10*t21*t35*t71*2.0+t2*t4*t10*t23*t298*2.0)*(1.0/2.0)-t37*t108*t305*(1.0/4.0);
    hess(6, 3) = t36*(t314+t315+t316+t317+t318+t110*t298*2.0-t3*t20*t23*t209*2.0-t4*t21*t23*t209*2.0)*(1.0/2.0)-t37*t114*t305*(1.0/4.0);
    hess(6, 4) = t36*(t319+t320+t321+t5*(t276-t2*t23*t247+t2*t35*t57*t71)*2.0-t20*t23*(t51-t72)*2.0-t3*t20*t23*t247*2.0-t4*t21*t23*t247*2.0+t2*t23*t57*t298*2.0+t3*t20*t35*t57*t71*2.0)*(1.0/2.0)-t37*t119*t305*(1.0/4.0);
    hess(6, 5) = t36*(t323+t324+t325+t5*(t322-t2*t23*t289+t2*t35*t66*t71)*2.0-t21*t23*(t51-t72)*2.0-t3*t20*t23*t289*2.0-t4*t21*t23*t289*2.0+t2*t23*t66*t298*2.0+t4*t21*t35*t66*t71*2.0)*(1.0/2.0)-t37*t123*t305*(1.0/4.0);
    hess(6, 6) = t36*(t5*(t23*(t51-t72)*2.0+t2*t35*(t327*t327)+t2*t23*t331)*2.0+(t298*t298)*2.0+t7*t34*(t328*t328)*2.0+t8*t34*(t329*t329)*2.0+t3*t20*t23*t331*2.0+t4*t21*t23*t331*2.0+t3*t20*t35*(t332*t332)*2.0+t4*t21*t35*(t333*t333)*2.0)*(1.0/2.0)-t37*(t326*t326)*(1.0/4.0);
    hess(6, 7) = t36*(t366+t367+t368+t369+t370+t5*(t290+t2*t23*t334+t2*t35*(t51-t72)*(t58-t79))*2.0+t20*t23*(t51-t72)*2.0+t3*t20*t23*t334*2.0+t3*t20*t35*(t51-t72)*(t58-t79)*2.0)*(1.0/2.0)-t37*t305*t346*(1.0/4.0);
    hess(6, 8) = t36*(t406+t407+t408+t409+t410+t5*(t251+t2*t23*t337+t2*t35*(t51-t72)*(t67-t86))*2.0+t21*t23*(t51-t72)*2.0+t4*t21*t23*t337*2.0+t4*t21*t35*(t51-t72)*(t67-t86)*2.0)*(1.0/2.0)-t37*t305*t381*(1.0/4.0);
    hess(7, 0) = t36*(-t81+t341-t342+t343+t344+t20*(t214-t340+t2*t3*t10*t35*t78)*2.0+t5*t6*t10*t35*t78*2.0+t2*t3*t10*t23*t336*2.0)*(1.0/2.0)-t26*t37*t346*(1.0/4.0);
    hess(7, 1) = t36*(t42*t336*2.0+t20*(t258+t347+t7*t10*t35*(t58-t79))*2.0+t2*t5*t23*t221*2.0+t4*t21*t23*t221*2.0+t3*t6*t10*t34*(t58-t79)*2.0+t3*t8*t10*t34*(t58-t79)*2.0+t2*t3*t5*t10*t35*(t58-t79)*2.0+t3*t4*t10*t21*t35*(t58-t79)*2.0)*(1.0/2.0)-t37*t102*t346*(1.0/4.0);
    hess(7, 2) = t36*(-t171+t349-t350+t351+t352+t20*(t306-t348+t3*t4*t10*t35*t78)*2.0+t8*t10*t21*t35*t78*2.0+t3*t4*t10*t23*t336*2.0)*(1.0/2.0)-t37*t108*t346*(1.0/4.0);
    hess(7, 3) = t36*(t353+t354+t355+t20*(t231-t3*t23*t211+t3*t35*t50*t78)*2.0-t5*t23*(t58-t79)*2.0-t2*t5*t23*t211*2.0-t4*t21*t23*t211*2.0+t3*t23*t50*t336*2.0+t2*t5*t35*t50*t78*2.0)*(1.0/2.0)-t37*t114*t346*(1.0/4.0);
    hess(7, 4) = t36*(t358+t359+t360+t361+t362+t60*t336*2.0-t2*t5*t23*t250*2.0-t4*t21*t23*t250*2.0)*(1.0/2.0)-t37*t119*t346*(1.0/4.0);
    hess(7, 5) = t36*(t363+t364+t365+t20*(t322-t3*t23*t292+t3*t35*t66*t78)*2.0-t21*t23*(t58-t79)*2.0-t2*t5*t23*t292*2.0-t4*t21*t23*t292*2.0+t3*t23*t66*t336*2.0+t4*t21*t35*t66*t78*2.0)*(1.0/2.0)-t37*t123*t346*(1.0/4.0);
    hess(7, 6) = t36*(t366+t367+t368+t369+t370+t20*(t287+t3*t23*t334+t3*t35*(t51-t72)*(t58-t79))*2.0+t5*t23*(t58-t79)*2.0+t2*t5*t23*t334*2.0+t2*t5*t35*(t51-t72)*(t58-t79)*2.0)*(1.0/2.0)-t37*t305*t346*(1.0/4.0);
    hess(7, 7) = t36*(t20*(t23*(t58-t79)*2.0+t3*t35*(t372*t372)+t3*t23*t376)*2.0+(t336*t336)*2.0+t6*t34*(t373*t373)*2.0+t8*t34*(t374*t374)*2.0+t2*t5*t23*t376*2.0+t4*t21*t23*t376*2.0+t2*t5*t35*(t377*t377)*2.0+t4*t21*t35*(t378*t378)*2.0)*(1.0/2.0)-t37*(t371*t371)*(1.0/4.0);
    hess(7, 8) = t36*(t411+t412+t413+t414+t415+t20*(t251+t3*t23*t379+t3*t35*(t58-t79)*(t67-t86))*2.0+t21*t23*(t58-t79)*2.0+t4*t21*t23*t379*2.0+t4*t21*t35*(t58-t79)*(t67-t86)*2.0)*(1.0/2.0)-t37*t346*t381*(1.0/4.0);
    hess(8, 0) = t36*(-t87+t383-t384+t385+t386+t21*(t214-t382+t2*t4*t10*t35*t85)*2.0+t5*t6*t10*t35*t85*2.0+t2*t4*t10*t23*t339*2.0)*(1.0/2.0)-t26*t37*t381*(1.0/4.0);
    hess(8, 1) = t36*(-t139+t388-t389+t390+t391+t21*(t258-t387+t3*t4*t10*t35*t85)*2.0+t7*t10*t20*t35*t85*2.0+t3*t4*t10*t23*t339*2.0)*(1.0/2.0)-t37*t102*t381*(1.0/4.0);
    hess(8, 2) = t36*(t44*t339*2.0+t21*(t306+t392+t8*t10*t35*(t67-t86))*2.0+t2*t5*t23*t265*2.0+t3*t20*t23*t265*2.0+t4*t6*t10*t34*(t67-t86)*2.0+t4*t7*t10*t34*(t67-t86)*2.0+t2*t4*t5*t10*t35*(t67-t86)*2.0+t3*t4*t10*t20*t35*(t67-t86)*2.0)*(1.0/2.0)-t37*t108*t381*(1.0/4.0);
    hess(8, 3) = t36*(t393+t394+t395+t21*(t231-t4*t23*t213+t4*t35*t50*t85)*2.0-t5*t23*(t67-t86)*2.0-t2*t5*t23*t213*2.0-t3*t20*t23*t213*2.0+t4*t23*t50*t339*2.0+t2*t5*t35*t50*t85*2.0)*(1.0/2.0)-t37*t114*t381*(1.0/4.0);
    hess(8, 4) = t36*(t396+t397+t398+t21*(t276-t4*t23*t253+t4*t35*t57*t85)*2.0-t20*t23*(t67-t86)*2.0-t2*t5*t23*t253*2.0-t3*t20*t23*t253*2.0+t4*t23*t57*t339*2.0+t3*t20*t35*t57*t85*2.0)*(1.0/2.0)-t37*t119*t381*(1.0/4.0);
    hess(8, 5) = t36*(t401+t402+t403+t404+t405+t69*t339*2.0-t2*t5*t23*t295*2.0-t3*t20*t23*t295*2.0)*(1.0/2.0)-t37*t123*t381*(1.0/4.0);
    hess(8, 6) = t36*(t406+t407+t408+t409+t410+t21*(t287+t4*t23*t337+t4*t35*(t51-t72)*(t67-t86))*2.0+t5*t23*(t67-t86)*2.0+t2*t5*t23*t337*2.0+t2*t5*t35*(t51-t72)*(t67-t86)*2.0)*(1.0/2.0)-t37*t305*t381*(1.0/4.0);
    hess(8, 7) = t36*(t411+t412+t413+t414+t415+t21*(t290+t4*t23*t379+t4*t35*(t58-t79)*(t67-t86))*2.0+t20*t23*(t67-t86)*2.0+t3*t20*t23*t379*2.0+t3*t20*t35*(t58-t79)*(t67-t86)*2.0)*(1.0/2.0)-t37*t346*t381*(1.0/4.0);
    hess(8, 8) = t36*(t21*(t23*(t67-t86)*2.0+t4*t35*(t417*t417)+t4*t23*t421)*2.0+(t339*t339)*2.0+t6*t34*(t418*t418)*2.0+t7*t34*(t419*t419)*2.0+t2*t5*t23*t421*2.0+t3*t20*t23*t421*2.0+t2*t5*t35*(t422*t422)*2.0+t3*t20*t35*(t423*t423)*2.0)*(1.0/2.0)-t37*(t416*t416)*(1.0/4.0);
  }
  
  template<typename T, int DIM>
  inline void pointedgevec(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, Eigen::Matrix<T, DIM, 1>& n, T& alpha)
  {
    alpha = std::min((T) 1.0, std::max((T) 0.0, (x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2)));
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha*(x3-x2);
    n = closest-x1;
  }
  
  template<typename T, int DIM>
  inline T pointlinealpha(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3)
  {
    return (x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2);
  }
  
  template<typename T, int DIM>
  inline void pointlinevec(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, Eigen::Matrix<T, DIM, 1>& n, T& alpha)
  {
    alpha = pointlinealpha<T, DIM>(x1, x2, x3);
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha*(x3-x2);
    n = closest-x1;
  }
  
  inline void grad_pointline_dist_fixed(const VectorXs& x0, const VectorXs& x1, 
	  const VectorXs& x2, const VectorXs& x3, const scalar& alpha_0, const scalar& alpha_1, 
	  VectorXs& grad)
  {
	const int DIM = x0.size();

	VectorXs cur = x0 + alpha_0 * (x1 - x0);
	VectorXs closest = x2 + alpha_1 * (x3 - x2);
	VectorXs n = closest - cur;
    
    scalar d = n.norm();
	VectorXs nhat = n / d;
    grad.segment(0, DIM) = -(1. - alpha_0) * nhat.segment(0, DIM);
    grad.segment(DIM, DIM) = -alpha_0 * nhat.segment(0, DIM);
    grad.segment(DIM * 2, DIM) = (1. - alpha_1) * nhat.segment(0, DIM);
    grad.segment(DIM * 3, DIM) = alpha_1 * nhat.segment(0, DIM);
  }
  
  template<typename T, int DIM>
  inline void grad_pointline_dist_fixed_constant(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, const scalar& alpha, const scalar& K, Eigen::Matrix<T, Eigen::Dynamic, 1>& grad)
  {
    Eigen::Matrix<T, DIM, 1> xc = x2 + alpha * (x3 - x2);
    
    grad.template segment<DIM>(0) = K * (x1 - xc);
    grad.template segment<DIM>(DIM) = -K * (x1 - xc) * (1.0 - alpha);
    grad.template segment<DIM>(DIM * 2) = -K * (x1 - xc) * alpha;
  }
  
  template<typename T, int DIM>
  inline void hess_pointline_dist_fixed(const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, const scalar& alpha, Eigen::Matrix<T, DIM * 3, DIM * 3>& hess)
  {
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha * (x3 - x2);
    Eigen::Matrix<T, DIM, 1> n = closest-x1;
    scalar d = n.norm();
    Eigen::Matrix<T, DIM, 1> nhat = n / d;
    Eigen::Matrix<T, DIM, DIM> innt = (Eigen::Matrix<T, DIM, DIM>::Identity() - nhat * nhat.transpose()) / d;
    hess.template block<DIM, DIM>(0, 0) = innt;
    hess.template block<DIM, DIM>(0, DIM) = innt * (alpha - 1.);
    hess.template block<DIM, DIM>(0, DIM * 2) = innt * (-alpha);
    
    hess.template block<DIM, DIM>(DIM, 0) = innt * (alpha - 1.);
    hess.template block<DIM, DIM>(DIM, DIM) = innt * (1. - alpha) * (1. - alpha);
    hess.template block<DIM, DIM>(DIM, DIM * 2) = innt * (1. - alpha) * alpha;
    
    hess.template block<DIM, DIM>(DIM * 2, 0) = innt * (-alpha);
    hess.template block<DIM, DIM>(DIM * 2, DIM) = innt * (1. - alpha) * alpha;
    hess.template block<DIM, DIM>(DIM * 2, DIM * 2) = innt * alpha * alpha;
  }
  
  template<typename T, int DIM>
  inline void hess_pointline_dist_fixed(const Eigen::Matrix<T, DIM, 1>& x0, const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, const scalar& alpha_0, const scalar& alpha_1, Eigen::Matrix<T, DIM * 4, DIM * 4>& hess)
  {
    Eigen::Matrix<T, DIM, 1> cur = x0 + alpha_0 * (x1 - x0);
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha_1 * (x3 - x2);
    Eigen::Matrix<T, DIM, 1> n = closest - cur;
    
    scalar d = n.norm();
    Eigen::Matrix<T, DIM, 1> nhat = n / d;
    Eigen::Matrix<T, DIM, DIM> innt = (Eigen::Matrix<T, DIM, DIM>::Identity() - nhat * nhat.transpose()) / d;

    hess.template block<DIM, DIM>(0, 0) = (1.0 - alpha_0) * (1.0 - alpha_0) * innt;
    hess.template block<DIM, DIM>(0, DIM) = (1.0 - alpha_0) * alpha_0 * innt;
    hess.template block<DIM, DIM>(0, DIM * 2) = -(1.0 - alpha_0) * (1.0 - alpha_1) * innt;
    hess.template block<DIM, DIM>(0, DIM * 3) = -(1.0 - alpha_0) * alpha_1 * innt;
    
    hess.template block<DIM, DIM>(DIM, 0) = (1.0 - alpha_0) * alpha_0 * innt;
    hess.template block<DIM, DIM>(DIM, DIM) = alpha_0 * alpha_0 * innt;
    hess.template block<DIM, DIM>(DIM, DIM * 2) = -alpha_0 * (1.0 - alpha_1) * innt;
    hess.template block<DIM, DIM>(DIM, DIM * 3) = -alpha_0 * alpha_1 * innt;
    
    hess.template block<DIM, DIM>(DIM * 2, 0) = -(1.0 - alpha_0) * (1.0 - alpha_1) * innt;
    hess.template block<DIM, DIM>(DIM * 2, DIM) = -alpha_0 * (1.0 - alpha_1) * innt;
    hess.template block<DIM, DIM>(DIM * 2, DIM * 2) = (1.0 - alpha_1) * (1.0 - alpha_1) * innt;
    hess.template block<DIM, DIM>(DIM * 2, DIM * 3) = alpha_1 * (1.0 - alpha_1) * innt;
    
    hess.template block<DIM, DIM>(DIM * 3, 0) = -(1.0 - alpha_0) * alpha_1 * innt;
    hess.template block<DIM, DIM>(DIM * 3, DIM) = -alpha_0 * alpha_1 * innt;
    hess.template block<DIM, DIM>(DIM * 3, DIM * 2) = alpha_1 * (1.0 - alpha_1) * innt;
    hess.template block<DIM, DIM>(DIM * 3, DIM * 3) = alpha_1 * alpha_1 * innt;
  }
  
  template<typename T, int DIM>
  inline void hess_pointline_dist_fixed(const Eigen::Matrix<T, DIM, 1>& x0, const Eigen::Matrix<T, DIM, 1>& x1, const Eigen::Matrix<T, DIM, 1>& x2, const Eigen::Matrix<T, DIM, 1>& x3, const scalar& alpha_0, const scalar& alpha_1, Eigen::Matrix<T, 1, DIM * 4>& hess, const int local_pidx, const int idir)
  {
    Eigen::Matrix<T, DIM, 1> cur = x0 + alpha_0 * (x1 - x0);
    Eigen::Matrix<T, DIM, 1> closest = x2 + alpha_1 * (x3 - x2);
    Eigen::Matrix<T, DIM, 1> n = closest - cur;
    
    scalar d = n.norm();
    Eigen::Matrix<T, DIM, 1> nhat = n / d;
    Eigen::Matrix<T, 1, DIM> innt = (-nhat(idir) * nhat.transpose()) / d;
    innt(0, idir) += 1.0 / d;
    
    if(local_pidx == 0) {
      hess.template block<1, DIM>(0, 0) = (1.0 - alpha_0) * (1.0 - alpha_0) * innt;
      hess.template block<1, DIM>(0, DIM) = (1.0 - alpha_0) * alpha_0 * innt;
      hess.template block<1, DIM>(0, DIM * 2) = -(1.0 - alpha_0) * (1.0 - alpha_1) * innt;
      hess.template block<1, DIM>(0, DIM * 3) = -(1.0 - alpha_0) * alpha_1 * innt;
    } else if (local_pidx == 1) {
      hess.template block<1, DIM>(0, 0) = (1.0 - alpha_0) * alpha_0 * innt;
      hess.template block<1, DIM>(0, DIM) = alpha_0 * alpha_0 * innt;
      hess.template block<1, DIM>(0, DIM * 2) = -alpha_0 * (1.0 - alpha_1) * innt;
      hess.template block<1, DIM>(0, DIM * 3) = -alpha_0 * alpha_1 * innt;
    } else if (local_pidx == 2) {
      hess.template block<1, DIM>(0, 0) = -(1.0 - alpha_0) * (1.0 - alpha_1) * innt;
      hess.template block<1, DIM>(0, DIM) = -alpha_0 * (1.0 - alpha_1) * innt;
      hess.template block<1, DIM>(0, DIM * 2) = (1.0 - alpha_1) * (1.0 - alpha_1) * innt;
      hess.template block<1, DIM>(0, DIM * 3) = alpha_1 * (1.0 - alpha_1) * innt;
    } else if (local_pidx == 3) {
      hess.template block<1, DIM>(0, 0) = -(1.0 - alpha_0) * alpha_1 * innt;
      hess.template block<1, DIM>(0, DIM) = -alpha_0 * alpha_1 * innt;
      hess.template block<1, DIM>(0, DIM * 2) = alpha_1 * (1.0 - alpha_1) * innt;
      hess.template block<1, DIM>(0, DIM * 3) = alpha_1 * alpha_1 * innt;
    }
  }
  
  template<typename T, int DIM>
  inline void hess_pointline_dist_fixed_constant(const scalar& alpha, const scalar& K, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& hess)
  {
    hess.template block<DIM, DIM>(0, 0) = K * Matrixs<DIM>::Identity();
    hess.template block<DIM, DIM>(0, DIM) = K * Matrixs<DIM>::Identity() * (alpha - 1.);
    hess.template block<DIM, DIM>(0, DIM * 2) = K * Matrixs<DIM>::Identity() * (-alpha);
    
    hess.template block<DIM, DIM>(DIM, 0) = K * Matrixs<DIM>::Identity() * (alpha - 1.);
    hess.template block<DIM, DIM>(DIM, DIM) = K * Matrixs<DIM>::Identity() * (1. - alpha) * (1. - alpha);
    hess.template block<DIM, DIM>(DIM, DIM * 2) = K * Matrixs<DIM>::Identity() * (1. - alpha) * alpha;
    
    hess.template block<DIM, DIM>(DIM * 2, 0) = K * Matrixs<DIM>::Identity() * (-alpha);
    hess.template block<DIM, DIM>(DIM * 2, DIM) = K * Matrixs<DIM>::Identity() * (1. - alpha) * alpha;
    hess.template block<DIM, DIM>(DIM * 2, DIM * 2) = K * Matrixs<DIM>::Identity() * alpha * alpha;
  }

  template<typename T>
  inline uint64 pphash64(T i, T j, T np)
  {
    return (uint64) std::max(i, j) * (uint64) np + (uint64) std::min(i, j);
  }
  
  template<typename T>
  inline void grad_pointline_dist(const Eigen::Matrix<T, 2, 1>& x0, const Eigen::Matrix<T, 2, 1>& x1, const Eigen::Matrix<T, 2, 1>& x2, Eigen::Matrix<T, 6, 1>& grad)
  {
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t5 = t2*t2;
    const scalar t6 = t3*t3;
    const scalar t7 = t5+t6;
    const scalar t8 = 1.0/t7;
    const scalar t9 = x0(0)-x1(0);
    const scalar t10 = t2*t9;
    const scalar t11 = x0(1)-x1(1);
    const scalar t12 = t3*t11;
    const scalar t13 = t10+t12;
    const scalar t15 = t2*t8*t13;
    const scalar t4 = t15-x0(0)+x1(0);
    const scalar t16 = t3*t8*t13;
    const scalar t14 = t16-x0(1)+x1(1);
    const scalar t17 = t4*t4;
    const scalar t18 = t14*t14;
    const scalar t19 = t17+t18;
    const scalar t20 = 1.0/sqrt(t19);
    const scalar t21 = x1(0)*2.0;
    const scalar t27 = x2(0)*2.0;
    const scalar t22 = t21-t27;
    const scalar t23 = 1.0/(t7*t7);
    const scalar t24 = t8*t13;
    const scalar t25 = x1(1)*2.0;
    const scalar t28 = x2(1)*2.0;
    const scalar t26 = t25-t28;
    grad(0) = t20*(t4*(t5*t8-1.0)*2.0+t2*t3*t8*t14*2.0)*(1.0/2.0);
    grad(1) = t20*(t14*(t6*t8-1.0)*2.0+t2*t3*t4*t8*2.0)*(1.0/2.0);
    grad(2) = t20*(t4*(t24+t2*t8*(-t21+x0(0)+x2(0))-t2*t13*t22*t23+1.0)*2.0+t14*(t3*t8*(x0(0)-x1(0)*2.0+x2(0))-t3*t13*t22*t23)*2.0)*(1.0/2.0);
    grad(3) = t20*(t14*(t24+t3*t8*(-t25+x0(1)+x2(1))-t3*t13*t23*t26+1.0)*2.0+t4*(t2*t8*(x0(1)-x1(1)*2.0+x2(1))-t2*t13*t23*t26)*2.0)*(1.0/2.0);
    grad(4) = t20*(t14*(t3*t8*t9-t3*t13*t22*t23)*2.0+t4*(t24+t2*t8*t9-t2*t13*t22*t23)*2.0)*(-1.0/2.0);
    grad(5) = t20*(t4*(t2*t8*t11-t2*t13*t23*t26)*2.0+t14*(t24+t3*t8*t11-t3*t13*t23*t26)*2.0)*(-1.0/2.0);
  }
  
  template<typename T>
  inline void grad_pointline_dist(const Eigen::Matrix<T, 3, 1>& x0, const Eigen::Matrix<T, 3, 1>& x1, const Eigen::Matrix<T, 3, 1>& x2, Eigen::Matrix<T, 9, 1>& gradE)
  {
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t4 = x1(2)-x2(2);
    const scalar t6 = t2*t2;
    const scalar t7 = t3*t3;
    const scalar t8 = t4*t4;
    const scalar t9 = t6+t7+t8;
    const scalar t10 = 1.0/t9;
    const scalar t11 = x0(0)-x1(0);
    const scalar t12 = t2*t11;
    const scalar t13 = x0(1)-x1(1);
    const scalar t14 = t3*t13;
    const scalar t15 = x0(2)-x1(2);
    const scalar t16 = t4*t15;
    const scalar t17 = t12+t14+t16;
    const scalar t20 = t2*t10*t17;
    const scalar t5 = t20-x0(0)+x1(0);
    const scalar t21 = t3*t10*t17;
    const scalar t18 = t21-x0(1)+x1(1);
    const scalar t22 = t4*t10*t17;
    const scalar t19 = t22-x0(2)+x1(2);
    const scalar t23 = t5*t5;
    const scalar t24 = t18*t18;
    const scalar t25 = t19*t19;
    const scalar t26 = t23+t24+t25;
    const scalar t27 = 1.0/sqrt(t26);
    const scalar t28 = x1(0)*2.0;
    const scalar t32 = x2(0)*2.0;
    const scalar t29 = t28-t32;
    const scalar t30 = 1.0/(t9*t9);
    const scalar t31 = -t28+x0(0)+x2(0);
    const scalar t33 = t10*t17;
    const scalar t34 = x1(1)*2.0;
    const scalar t37 = x2(1)*2.0;
    const scalar t35 = t34-t37;
    const scalar t36 = -t34+x0(1)+x2(1);
    const scalar t38 = x1(2)*2.0;
    const scalar t41 = x2(2)*2.0;
    const scalar t39 = t38-t41;
    const scalar t40 = -t38+x0(2)+x2(2);
    gradE(0) = t27*(t5*(t6*t10-1.0)*2.0+t2*t3*t10*t18*2.0+t2*t4*t10*t19*2.0)*(1.0/2.0);
    gradE(1) = t27*(t18*(t7*t10-1.0)*2.0+t2*t3*t5*t10*2.0+t3*t4*t10*t19*2.0)*(1.0/2.0);
    gradE(2) = t27*(t19*(t8*t10-1.0)*2.0+t2*t4*t5*t10*2.0+t3*t4*t10*t18*2.0)*(1.0/2.0);
    gradE(3) = t27*(t18*(t3*t10*t31-t3*t17*t29*t30)*2.0+t19*(t4*t10*t31-t4*t17*t29*t30)*2.0+t5*(t33+t2*t10*(x0(0)-x1(0)*2.0+x2(0))-t2*t17*t29*t30+1.0)*2.0)*(1.0/2.0);
    gradE(4) = t27*(t5*(t2*t10*t36-t2*t17*t30*t35)*2.0+t19*(t4*t10*t36-t4*t17*t30*t35)*2.0+t18*(t33+t3*t10*(x0(1)-x1(1)*2.0+x2(1))-t3*t17*t30*t35+1.0)*2.0)*(1.0/2.0);
    gradE(5) = t27*(t5*(t2*t10*t40-t2*t17*t30*t39)*2.0+t18*(t3*t10*t40-t3*t17*t30*t39)*2.0+t19*(t33+t4*t10*(x0(2)-x1(2)*2.0+x2(2))-t4*t17*t30*t39+1.0)*2.0)*(1.0/2.0);
    gradE(6) = t27*(t18*(t3*t10*t11-t3*t17*t29*t30)*2.0+t19*(t4*t10*t11-t4*t17*t29*t30)*2.0+t5*(t33+t2*t10*t11-t2*t17*t29*t30)*2.0)*(-1.0/2.0);
    gradE(7) = t27*(t5*(t2*t10*t13-t2*t17*t30*t35)*2.0+t19*(t4*t10*t13-t4*t17*t30*t35)*2.0+t18*(t33+t3*t10*t13-t3*t17*t30*t35)*2.0)*(-1.0/2.0);
    gradE(8) = t27*(t5*(t2*t10*t15-t2*t17*t30*t39)*2.0+t18*(t3*t10*t15-t3*t17*t30*t39)*2.0+t19*(t33+t4*t10*t15-t4*t17*t30*t39)*2.0)*(-1.0/2.0);
  }
  
  template<typename T>
  inline void hess_pointline_dist(const Eigen::Matrix<T, 2, 1>& x0, const Eigen::Matrix<T, 2, 1>& x1, const Eigen::Matrix<T, 2, 1>& x2, Eigen::Matrix<T, 6, 6>& hessE)
  {
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t5 = t2*t2;
    const scalar t6 = t3*t3;
    const scalar t7 = t5+t6;
    const scalar t8 = 1.0/t7;
    const scalar t9 = x0(0)-x1(0);
    const scalar t10 = t2*t9;
    const scalar t11 = x0(1)-x1(1);
    const scalar t12 = t3*t11;
    const scalar t13 = t10+t12;
    const scalar t15 = t2*t8*t13;
    const scalar t4 = t15-x0(0)+x1(0);
    const scalar t16 = t3*t8*t13;
    const scalar t14 = t16-x0(1)+x1(1);
    const scalar t21 = t5*t8;
    const scalar t22 = t21-1.0;
    const scalar t25 = t4*t22*2.0;
    const scalar t26 = t2*t3*t8*t14*2.0;
    const scalar t17 = t25+t26;
    const scalar t18 = t4*t4;
    const scalar t19 = t14*t14;
    const scalar t20 = t18+t19;
    const scalar t23 = 1.0/sqrt(t20);
    const scalar t24 = 1.0/pow(t20,3.0/2.0);
    const scalar t27 = t6*t8;
    const scalar t28 = t27-1.0;
    const scalar t29 = 1.0/(t7*t7);
    const scalar t30 = x1(0)*2.0;
    const scalar t32 = x2(0)*2.0;
    const scalar t31 = t30-t32;
    const scalar t33 = -t30+x0(0)+x2(0);
    const scalar t34 = t3*t8*t33;
    const scalar t47 = t3*t13*t29*t31;
    const scalar t35 = t34-t47;
    const scalar t36 = t8*t13;
    const scalar t37 = x1(1)*2.0;
    const scalar t39 = x2(1)*2.0;
    const scalar t38 = t37-t39;
    const scalar t40 = -t37+x0(1)+x2(1);
    const scalar t41 = t3*t8*t40;
    const scalar t52 = t3*t13*t29*t38;
    const scalar t42 = t36+t41-t52+1.0;
    const scalar t43 = t8*t31;
    const scalar t85 = t5*t29*t31;
    const scalar t44 = t43-t85;
    const scalar t45 = t4*t44*2.0;
    const scalar t46 = t3*t8*t14*2.0;
    const scalar t74 = t3*t8*t9;
    const scalar t48 = t47-t74;
    const scalar t49 = t2*t8*t9;
    const scalar t67 = t2*t13*t29*t31;
    const scalar t50 = t36+t49-t67;
    const scalar t51 = t2*t8*t14*2.0;
    const scalar t53 = t2*t8*t11;
    const scalar t71 = t2*t13*t29*t38;
    const scalar t54 = t53-t71;
    const scalar t55 = t3*t8*t11;
    const scalar t56 = t36-t52+t55;
    const scalar t57 = t2*t3*t8*t22*2.0;
    const scalar t58 = t2*t3*t8*t28*2.0;
    const scalar t59 = t57+t58;
    const scalar t60 = t23*t59*(1.0/2.0);
    const scalar t61 = t14*t28*2.0;
    const scalar t62 = t2*t3*t4*t8*2.0;
    const scalar t63 = t61+t62;
    const scalar t64 = t60-t17*t24*t63*(1.0/4.0);
    const scalar t65 = t5*t6*t29*2.0;
    const scalar t66 = t2*t8*t33;
    const scalar t68 = t14*t35*2.0;
    const scalar t69 = t36+t66-t67+1.0;
    const scalar t70 = t2*t8*t40;
    const scalar t72 = t70-t71;
    const scalar t73 = t14*t42*2.0;
    const scalar t75 = t3*t4*t8*2.0;
    const scalar t76 = t4*t50*2.0;
    const scalar t101 = t14*t48*2.0;
    const scalar t77 = t76-t101;
    const scalar t78 = t8*t38;
    const scalar t110 = t6*t29*t38;
    const scalar t79 = t78-t110;
    const scalar t80 = t14*t79*2.0;
    const scalar t81 = t2*t4*t8*2.0;
    const scalar t82 = t4*t54*2.0;
    const scalar t83 = t14*t56*2.0;
    const scalar t84 = t82+t83;
    const scalar t86 = t2*t3*t8*t35*2.0;
    const scalar t87 = t4*t69*2.0;
    const scalar t88 = t68+t87;
    const scalar t89 = t28*t35*2.0;
    const scalar t90 = t3*t8;
    const scalar t127 = t2*t3*t29*t31;
    const scalar t91 = t90-t127;
    const scalar t92 = t6*t14*t29*t31*2.0;
    const scalar t93 = t2*t3*t8*t69*2.0;
    const scalar t94 = t31*t31;
    const scalar t95 = 1.0/(t7*t7*t7);
    const scalar t96 = t4*t72*2.0;
    const scalar t97 = t73+t96;
    const scalar t98 = t13*t29*t31*2.0;
    const scalar t99 = t2*t13*t29*2.0;
    const scalar t100 = t3*t13*t29*2.0;
    const scalar t102 = t13*t29*t31;
    const scalar t103 = t3*t29*t33*t38;
    const scalar t104 = t13*t29*t38;
    const scalar t105 = t2*t29*t33*t38;
    const scalar t106 = t2*t8;
    const scalar t107 = t2*t3*t8*t42*2.0;
    const scalar t108 = t28*t42*2.0;
    const scalar t109 = t3*t8*2.0;
    const scalar t159 = t2*t3*t29*t38;
    const scalar t111 = t106-t159;
    const scalar t112 = t2*t3*t8*t72*2.0;
    const scalar t113 = t2*t29*t31*t40;
    const scalar t123 = t8*t40;
    const scalar t124 = t2*t13*t31*t38*t95*2.0;
    const scalar t114 = t104+t105+t113-t123-t124;
    const scalar t115 = t4*t114*2.0;
    const scalar t116 = t3*t29*t31*t40;
    const scalar t125 = t3*t13*t31*t38*t95*2.0;
    const scalar t134 = t8*t33;
    const scalar t117 = t102+t103+t116-t125-t134;
    const scalar t118 = t14*t117*2.0;
    const scalar t119 = t115+t118-t35*t42*2.0-t69*t72*2.0;
    const scalar t120 = t23*t119*(-1.0/2.0)-t24*t88*t97*(1.0/4.0);
    const scalar t121 = t2*t8*2.0;
    const scalar t122 = t38*t38;
    const scalar t126 = t13*t29*t38*2.0;
    const scalar t128 = t14*t91*2.0;
    const scalar t129 = t22*t50*2.0;
    const scalar t130 = t17*t24*t77*(1.0/4.0);
    const scalar t131 = t4*t91*2.0;
    const scalar t132 = t28*t48*2.0;
    const scalar t133 = t24*t63*t77*(1.0/4.0);
    const scalar t135 = t2*t29*t31*t33;
    const scalar t136 = t2*t9*t29*t31;
    const scalar t137 = t35*(t47-t74)*2.0;
    const scalar t138 = t3*t29*t31*t33;
    const scalar t139 = t3*t9*t29*t31;
    const scalar t155 = t3*t13*t94*t95*2.0;
    const scalar t140 = t90+t100+t138+t139-t155;
    const scalar t141 = t14*t140*2.0;
    const scalar t142 = t24*t77*t88*(1.0/4.0);
    const scalar t143 = t2*t9*t29*t38;
    const scalar t144 = t104+t113-t123-t124+t143;
    const scalar t145 = t4*t144*2.0;
    const scalar t146 = t42*(t47-t74)*2.0;
    const scalar t147 = t3*t9*t29*t38;
    const scalar t157 = t8*t9;
    const scalar t148 = t102+t116-t125+t147-t157;
    const scalar t149 = t14*t148*2.0;
    const scalar t150 = t145+t146+t149-t50*t72*2.0;
    const scalar t151 = t23*t150*(1.0/2.0);
    const scalar t152 = t24*t77*t97*(1.0/4.0);
    const scalar t153 = t151+t152;
    const scalar t154 = t47-t74;
    const scalar t156 = t2*t11*t29*t31;
    const scalar t158 = t3*t11*t29*t31;
    const scalar t160 = t14*t111*2.0;
    const scalar t161 = t22*t54*2.0;
    const scalar t162 = t2*t3*t8*t56*2.0;
    const scalar t163 = t17*t24*t84*(1.0/4.0);
    const scalar t164 = t109-t110;
    const scalar t165 = t14*t164*2.0;
    const scalar t166 = t4*t111*2.0;
    const scalar t167 = t28*t56*2.0;
    const scalar t168 = t2*t3*t8*t54*2.0;
    const scalar t169 = t24*t63*t84*(1.0/4.0);
    const scalar t170 = t54*t69*2.0;
    const scalar t171 = t35*t56*2.0;
    const scalar t172 = t24*t84*t88*(1.0/4.0);
    const scalar t173 = t3*t29*t38*t40;
    const scalar t174 = t3*t11*t29*t38;
    const scalar t185 = t8*t11;
    const scalar t189 = t3*t13*t95*t122*2.0;
    const scalar t175 = t90+t100-t123+t126+t173+t174-t185-t189;
    const scalar t176 = t54*t72*2.0;
    const scalar t177 = t42*t56*2.0;
    const scalar t178 = t2*t29*t38*t40;
    const scalar t179 = t2*t11*t29*t38;
    const scalar t188 = t2*t13*t95*t122*2.0;
    const scalar t180 = t99+t106+t178+t179-t188;
    const scalar t181 = t176+t177-t4*t180*2.0-t14*t175*2.0;
    const scalar t182 = t24*t84*t97*(1.0/4.0);
    const scalar t183 = t182-t23*t181*(1.0/2.0);
    const scalar t184 = t48*t56*2.0;
    const scalar t186 = t102-t125+t147-t157+t158;
    const scalar t187 = t14*t186*2.0;
    hessE(0, 0) = (t17*t17)*t24*(-1.0/4.0)+t23*(t65+(t22*t22)*2.0)*(1.0/2.0);
    hessE(0, 1) = t64;
    hessE(0, 2) = t23*(t45+t46+t86+t22*(t36+t2*t8*(x0(0)-x1(0)*2.0+x2(0))-t2*t13*t29*t31+1.0)*2.0-t2*t3*t14*t29*t31*2.0)*(1.0/2.0)-t17*t24*(t68+t4*(t36+t66-t2*t13*t29*t31+1.0)*2.0)*(1.0/4.0);
    hessE(0, 3) = t23*(t51+t107+t22*(t2*t8*(x0(1)-x1(1)*2.0+x2(1))-t2*t13*t29*t38)*2.0-t4*t5*t29*t38*2.0-t2*t3*t14*t29*t38*2.0)*(1.0/2.0)-t17*t24*(t73+t4*(t70-t2*t13*t29*t38)*2.0)*(1.0/4.0);
    hessE(0, 4) = t130-t23*(t45+t46+t129-t2*t3*t8*t48*2.0-t2*t3*t14*t29*t31*2.0)*(1.0/2.0);
    hessE(0, 5) = t163-t23*(t51+t161+t162-t4*t5*t29*t38*2.0-t2*t3*t14*t29*t38*2.0)*(1.0/2.0);
    hessE(1, 0) = t64;
    hessE(1, 1) = t24*(t63*t63)*(-1.0/4.0)+t23*(t65+(t28*t28)*2.0)*(1.0/2.0);
    hessE(1, 2) = t23*(t75+t89+t93-t6*t14*t29*t31*2.0-t2*t3*t4*t29*t31*2.0)*(1.0/2.0)-t24*t63*t88*(1.0/4.0);
    hessE(1, 3) = t23*(t80+t81+t108+t112-t2*t3*t4*t29*t38*2.0)*(1.0/2.0)-t24*t63*t97*(1.0/4.0);
    hessE(1, 4) = t133+t23*(-t75+t92+t132-t2*t3*t8*t50*2.0+t2*t3*t4*t29*t31*2.0)*(1.0/2.0);
    hessE(1, 5) = t169-t23*(t80+t81+t167+t168-t2*t3*t4*t29*t38*2.0)*(1.0/2.0);
    hessE(2, 0) = t23*(t86+t128+t22*t69*2.0-t4*(t85-t2*t8*2.0)*2.0)*(1.0/2.0)-t17*t24*t88*(1.0/4.0);
    hessE(2, 1) = t23*(t89-t92+t93+t131)*(1.0/2.0)-t24*t63*t88*(1.0/4.0);
    hessE(2, 2) = t23*(t4*(t98+t99+t121-t8*t33*2.0+t2*t29*t31*t33*2.0-t2*t13*t94*t95*2.0)*2.0+t14*(t100+t109+t3*t29*t31*t33*2.0-t3*t13*t94*t95*2.0)*2.0-(t35*t35)*2.0-(t69*t69)*2.0)*(-1.0/2.0)-t24*(t88*t88)*(1.0/4.0);
    hessE(2, 3) = t120;
    hessE(2, 4) = t142+t23*(t137+t141+t4*(t98+t99+t106+t135+t136-t8*t9-t8*t33-t2*t13*t94*t95*2.0)*2.0-t50*t69*2.0)*(1.0/2.0);
    hessE(2, 5) = t172-t23*(t170+t171-t4*(t104+t105+t156-t8*t11-t2*t13*t31*t38*t95*2.0)*2.0-t14*(t102+t103+t158-t8*t33-t3*t13*t31*t38*t95*2.0)*2.0)*(1.0/2.0);
    hessE(3, 0) = t23*(t107+t160+t22*t72*2.0-t4*t5*t29*t38*2.0)*(1.0/2.0)-t17*t24*t97*(1.0/4.0);
    hessE(3, 1) = t23*(t108+t112+t165+t166)*(1.0/2.0)-t24*t63*t97*(1.0/4.0);
    hessE(3, 2) = t120;
    hessE(3, 3) = t23*(t14*(t100+t109+t126-t8*t40*2.0+t3*t29*t38*t40*2.0-t3*t13*t95*t122*2.0)*2.0+t4*(t99+t121+t2*t29*t38*t40*2.0-t2*t13*t95*t122*2.0)*2.0-(t42*t42)*2.0-(t72*t72)*2.0)*(-1.0/2.0)-t24*(t97*t97)*(1.0/4.0);
    hessE(3, 4) = t153;
    hessE(3, 5) = t183;
    hessE(4, 0) = t130-t23*(t128+t129-t4*(t85-t121)*2.0-t2*t3*t8*t48*2.0)*(1.0/2.0);
    hessE(4, 1) = t133+t23*(t92-t131+t132-t2*t3*t8*t50*2.0)*(1.0/2.0);
    hessE(4, 2) = t142+t23*(t137+t141-t50*t69*2.0+t4*(t98+t99+t106-t134+t135+t136-t8*t9-t2*t13*t94*t95*2.0)*2.0)*(1.0/2.0);
    hessE(4, 3) = t153;
    hessE(4, 4) = t23*(t14*(t100-t155+t3*t9*t29*t31*2.0)*2.0+t4*(t98+t99-t8*t9*2.0+t2*t9*t29*t31*2.0-t2*t13*t94*t95*2.0)*2.0-(t50*t50)*2.0-(t154*t154)*2.0)*(-1.0/2.0)-t24*(t77*t77)*(1.0/4.0);
    hessE(4, 5) = t23*(t184+t187-t50*t54*2.0+t4*(t104-t124+t143+t156-t8*t11)*2.0)*(-1.0/2.0)-t24*t77*t84*(1.0/4.0);
    hessE(5, 0) = t163-t23*(t160+t161+t162-t4*t5*t29*t38*2.0)*(1.0/2.0);
    hessE(5, 1) = t169-t23*(t165+t166+t167+t168)*(1.0/2.0);
    hessE(5, 2) = t172-t23*(t170+t171-t14*(t102+t103-t125-t134+t158)*2.0-t4*(t104+t105-t124+t156-t8*t11)*2.0)*(1.0/2.0);
    hessE(5, 3) = t183;
    hessE(5, 4) = t23*(t184+t187+t4*(t104-t124+t143+t156-t185)*2.0-t50*t54*2.0)*(-1.0/2.0)-t24*t77*t84*(1.0/4.0);
    hessE(5, 5) = t24*(t84*t84)*(-1.0/4.0)-t23*(t14*(t100+t126-t189-t8*t11*2.0+t3*t11*t29*t38*2.0)*2.0+t4*(t99-t188+t2*t11*t29*t38*2.0)*2.0-(t54*t54)*2.0-(t56*t56)*2.0)*(1.0/2.0);
  }
  
  template<typename T>
  inline void hess_pointline_dist(const Eigen::Matrix<T, 3, 1>& x0, const Eigen::Matrix<T, 3, 1>& x1, const Eigen::Matrix<T, 3, 1>& x2, Eigen::Matrix<T, 9, 9>& hessE)
  {
    const scalar t2 = x1(0)-x2(0);
    const scalar t3 = x1(1)-x2(1);
    const scalar t4 = x1(2)-x2(2);
    const scalar t6 = t2*t2;
    const scalar t7 = t3*t3;
    const scalar t8 = t4*t4;
    const scalar t9 = t6+t7+t8;
    const scalar t10 = 1.0/t9;
    const scalar t11 = x0(0)-x1(0);
    const scalar t12 = t2*t11;
    const scalar t13 = x0(1)-x1(1);
    const scalar t14 = t3*t13;
    const scalar t15 = x0(2)-x1(2);
    const scalar t16 = t4*t15;
    const scalar t17 = t12+t14+t16;
    const scalar t22 = t2*t10*t17;
    const scalar t5 = t22-x0(0)+x1(0);
    const scalar t24 = t3*t10*t17;
    const scalar t18 = t24-x0(1)+x1(1);
    const scalar t26 = t4*t10*t17;
    const scalar t19 = t26-x0(2)+x1(2);
    const scalar t29 = t6*t10;
    const scalar t20 = t29-1.0;
    const scalar t21 = 1.0/(t9*t9);
    const scalar t23 = t5*t5;
    const scalar t25 = t18*t18;
    const scalar t27 = t19*t19;
    const scalar t28 = t23+t25+t27;
    const scalar t33 = t5*t20*2.0;
    const scalar t34 = t2*t3*t10*t18*2.0;
    const scalar t35 = t2*t4*t10*t19*2.0;
    const scalar t30 = t33+t34+t35;
    const scalar t31 = 1.0/sqrt(t28);
    const scalar t32 = 1.0/pow(t28,3.0/2.0);
    const scalar t36 = t7*t10;
    const scalar t37 = t36-1.0;
    const scalar t38 = t8*t10;
    const scalar t39 = t38-1.0;
    const scalar t40 = x1(0)*2.0;
    const scalar t42 = x2(0)*2.0;
    const scalar t41 = t40-t42;
    const scalar t43 = -t40+x0(0)+x2(0);
    const scalar t44 = t10*t17;
    const scalar t45 = t3*t10*t43;
    const scalar t70 = t3*t17*t21*t41;
    const scalar t46 = t45-t70;
    const scalar t47 = t4*t10*t43;
    const scalar t71 = t4*t17*t21*t41;
    const scalar t48 = t47-t71;
    const scalar t49 = x1(1)*2.0;
    const scalar t51 = x2(1)*2.0;
    const scalar t50 = t49-t51;
    const scalar t52 = -t49+x0(1)+x2(1);
    const scalar t53 = t3*t10*t52;
    const scalar t77 = t3*t17*t21*t50;
    const scalar t54 = t44+t53-t77+1.0;
    const scalar t55 = t4*t10*t52;
    const scalar t78 = t4*t17*t21*t50;
    const scalar t56 = t55-t78;
    const scalar t57 = x1(2)*2.0;
    const scalar t59 = x2(2)*2.0;
    const scalar t58 = t57-t59;
    const scalar t60 = -t57+x0(2)+x2(2);
    const scalar t61 = t4*t10*t60;
    const scalar t85 = t4*t17*t21*t58;
    const scalar t62 = t44+t61-t85+1.0;
    const scalar t63 = t3*t10*t60;
    const scalar t86 = t3*t17*t21*t58;
    const scalar t64 = t63-t86;
    const scalar t65 = t10*t41;
    const scalar t169 = t6*t21*t41;
    const scalar t66 = t65-t169;
    const scalar t67 = t5*t66*2.0;
    const scalar t68 = t3*t10*t18*2.0;
    const scalar t69 = t4*t10*t19*2.0;
    const scalar t72 = t2*t10*t11;
    const scalar t108 = t2*t17*t21*t41;
    const scalar t73 = t44+t72-t108;
    const scalar t122 = t3*t10*t11;
    const scalar t74 = t70-t122;
    const scalar t124 = t4*t10*t11;
    const scalar t75 = t71-t124;
    const scalar t76 = t2*t10*t18*2.0;
    const scalar t79 = t3*t10*t13;
    const scalar t80 = t44-t77+t79;
    const scalar t81 = t2*t10*t13;
    const scalar t113 = t2*t17*t21*t50;
    const scalar t82 = t81-t113;
    const scalar t132 = t4*t10*t13;
    const scalar t83 = t78-t132;
    const scalar t84 = t2*t10*t19*2.0;
    const scalar t87 = t4*t10*t15;
    const scalar t88 = t44-t85+t87;
    const scalar t89 = t2*t10*t15;
    const scalar t118 = t2*t17*t21*t58;
    const scalar t90 = t89-t118;
    const scalar t136 = t3*t10*t15;
    const scalar t91 = t86-t136;
    const scalar t92 = t2*t3*t10*t20*2.0;
    const scalar t93 = t2*t3*t10*t37*2.0;
    const scalar t94 = t2*t3*t8*t21*2.0;
    const scalar t95 = t92+t93+t94;
    const scalar t96 = t31*t95*(1.0/2.0);
    const scalar t97 = t18*t37*2.0;
    const scalar t98 = t2*t3*t5*t10*2.0;
    const scalar t99 = t3*t4*t10*t19*2.0;
    const scalar t100 = t97+t98+t99;
    const scalar t101 = t96-t30*t32*t100*(1.0/4.0);
    const scalar t102 = t6*t7*t21*2.0;
    const scalar t103 = t19*t39*2.0;
    const scalar t104 = t2*t4*t5*t10*2.0;
    const scalar t105 = t3*t4*t10*t18*2.0;
    const scalar t106 = t103+t104+t105;
    const scalar t107 = t2*t10*t43;
    const scalar t109 = t44+t107-t108+1.0;
    const scalar t110 = t18*t46*2.0;
    const scalar t111 = t19*t48*2.0;
    const scalar t112 = t2*t10*t52;
    const scalar t114 = t18*t54*2.0;
    const scalar t115 = t112-t113;
    const scalar t116 = t19*t56*2.0;
    const scalar t117 = t2*t10*t60;
    const scalar t119 = t19*t62*2.0;
    const scalar t120 = t117-t118;
    const scalar t121 = t18*t64*2.0;
    const scalar t123 = t3*t5*t10*2.0;
    const scalar t125 = t18*t74*2.0;
    const scalar t126 = t19*t75*2.0;
    const scalar t162 = t5*t73*2.0;
    const scalar t127 = t125+t126-t162;
    const scalar t128 = t10*t50;
    const scalar t208 = t7*t21*t50;
    const scalar t129 = t128-t208;
    const scalar t130 = t18*t129*2.0;
    const scalar t131 = t2*t5*t10*2.0;
    const scalar t133 = t18*t80*2.0;
    const scalar t134 = t5*t82*2.0;
    const scalar t164 = t19*t83*2.0;
    const scalar t135 = t133+t134-t164;
    const scalar t137 = t3*t10*t19*2.0;
    const scalar t138 = t19*t88*2.0;
    const scalar t139 = t5*t90*2.0;
    const scalar t168 = t18*t91*2.0;
    const scalar t140 = t138+t139-t168;
    const scalar t141 = t2*t4*t10*t20*2.0;
    const scalar t142 = t2*t4*t10*t39*2.0;
    const scalar t143 = t2*t4*t7*t21*2.0;
    const scalar t144 = t141+t142+t143;
    const scalar t145 = t31*t144*(1.0/2.0);
    const scalar t146 = t145-t30*t32*t106*(1.0/4.0);
    const scalar t147 = t3*t4*t10*t37*2.0;
    const scalar t148 = t3*t4*t10*t39*2.0;
    const scalar t149 = t3*t4*t6*t21*2.0;
    const scalar t150 = t147+t148+t149;
    const scalar t151 = t31*t150*(1.0/2.0);
    const scalar t152 = t151-t32*t100*t106*(1.0/4.0);
    const scalar t153 = t6*t8*t21*2.0;
    const scalar t154 = t7*t8*t21*2.0;
    const scalar t155 = t5*t109*2.0;
    const scalar t156 = t110+t111+t155;
    const scalar t157 = t5*t115*2.0;
    const scalar t158 = t114+t116+t157;
    const scalar t159 = t5*t120*2.0;
    const scalar t160 = t119+t121+t159;
    const scalar t161 = t4*t5*t10*2.0;
    const scalar t163 = t4*t10*t18*2.0;
    const scalar t165 = t10*t58;
    const scalar t250 = t8*t21*t58;
    const scalar t166 = t165-t250;
    const scalar t167 = t19*t166*2.0;
    const scalar t170 = t2*t3*t10*t46*2.0;
    const scalar t171 = t2*t4*t10*t48*2.0;
    const scalar t172 = t37*t46*2.0;
    const scalar t173 = t3*t10;
    const scalar t278 = t2*t3*t21*t41;
    const scalar t174 = t173-t278;
    const scalar t175 = t2*t3*t10*t109*2.0;
    const scalar t176 = t7*t18*t21*t41*2.0;
    const scalar t177 = t3*t4*t10*t48*2.0;
    const scalar t178 = t3*t4*t19*t21*t41*2.0;
    const scalar t179 = t39*t48*2.0;
    const scalar t180 = t4*t10;
    const scalar t280 = t2*t4*t21*t41;
    const scalar t181 = t180-t280;
    const scalar t182 = t2*t4*t10*t109*2.0;
    const scalar t183 = t8*t19*t21*t41*2.0;
    const scalar t184 = t3*t4*t10*t46*2.0;
    const scalar t185 = t3*t4*t18*t21*t41*2.0;
    const scalar t186 = t41*t41;
    const scalar t187 = 1.0/(t9*t9*t9);
    const scalar t188 = t17*t21*t41;
    const scalar t189 = t3*t17*t21*2.0;
    const scalar t190 = t4*t17*t21*2.0;
    const scalar t191 = t2*t17*t21*2.0;
    const scalar t192 = t17*t21*t41*2.0;
    const scalar t193 = t17*t21*t50;
    const scalar t194 = t2*t21*t43*t50;
    const scalar t195 = t4*t21*t43*t50;
    const scalar t196 = t3*t21*t43*t50;
    const scalar t197 = t17*t21*t58;
    const scalar t198 = t2*t21*t43*t58;
    const scalar t199 = t3*t21*t43*t58;
    const scalar t200 = t4*t21*t43*t58;
    const scalar t201 = t2*t10;
    const scalar t202 = t2*t3*t10*t54*2.0;
    const scalar t203 = t5*t6*t21*t50*2.0;
    const scalar t204 = t2*t4*t10*t56*2.0;
    const scalar t205 = t2*t4*t19*t21*t50*2.0;
    const scalar t206 = t37*t54*2.0;
    const scalar t207 = t3*t10*2.0;
    const scalar t334 = t2*t3*t21*t50;
    const scalar t209 = t201-t334;
    const scalar t210 = t2*t3*t10*t115*2.0;
    const scalar t211 = t3*t4*t10*t56*2.0;
    const scalar t212 = t39*t56*2.0;
    const scalar t342 = t3*t4*t21*t50;
    const scalar t213 = t180-t342;
    const scalar t214 = t3*t4*t10*t54*2.0;
    const scalar t215 = t8*t19*t21*t50*2.0;
    const scalar t216 = t2*t4*t10*t115*2.0;
    const scalar t217 = t2*t4*t5*t21*t50*2.0;
    const scalar t218 = t4*t21*t41*t52;
    const scalar t234 = t4*t17*t41*t50*t187*2.0;
    const scalar t219 = t195+t218-t234;
    const scalar t220 = t19*t219*2.0;
    const scalar t221 = t2*t21*t41*t52;
    const scalar t232 = t10*t52;
    const scalar t235 = t2*t17*t41*t50*t187*2.0;
    const scalar t222 = t193+t194+t221-t232-t235;
    const scalar t223 = t5*t222*2.0;
    const scalar t224 = t3*t21*t41*t52;
    const scalar t233 = t3*t17*t41*t50*t187*2.0;
    const scalar t261 = t10*t43;
    const scalar t225 = t188+t196+t224-t233-t261;
    const scalar t226 = t18*t225*2.0;
    const scalar t227 = t220+t223+t226-t46*t54*2.0-t48*t56*2.0-t109*t115*2.0;
    const scalar t228 = t31*t227*(-1.0/2.0)-t32*t156*t158*(1.0/4.0);
    const scalar t229 = t2*t10*2.0;
    const scalar t230 = t50*t50;
    const scalar t231 = t4*t10*2.0;
    const scalar t236 = t17*t21*t50*2.0;
    const scalar t237 = t3*t21*t52*t58;
    const scalar t238 = t2*t21*t52*t58;
    const scalar t239 = t4*t21*t52*t58;
    const scalar t240 = t2*t4*t10*t62*2.0;
    const scalar t241 = t5*t6*t21*t58*2.0;
    const scalar t242 = t2*t3*t10*t64*2.0;
    const scalar t243 = t2*t3*t18*t21*t58*2.0;
    const scalar t244 = t37*t64*2.0;
    const scalar t245 = t3*t4*t10*t62*2.0;
    const scalar t246 = t7*t18*t21*t58*2.0;
    const scalar t247 = t2*t3*t10*t120*2.0;
    const scalar t248 = t2*t3*t5*t21*t58*2.0;
    const scalar t249 = t39*t62*2.0;
    const scalar t391 = t2*t4*t21*t58;
    const scalar t251 = t201-t391;
    const scalar t395 = t3*t4*t21*t58;
    const scalar t252 = t173-t395;
    const scalar t253 = t2*t4*t10*t120*2.0;
    const scalar t254 = t3*t4*t10*t64*2.0;
    const scalar t255 = t3*t21*t41*t60;
    const scalar t272 = t3*t17*t41*t58*t187*2.0;
    const scalar t256 = t199+t255-t272;
    const scalar t257 = t18*t256*2.0;
    const scalar t258 = t2*t21*t41*t60;
    const scalar t266 = t10*t60;
    const scalar t273 = t2*t17*t41*t58*t187*2.0;
    const scalar t259 = t197+t198+t258-t266-t273;
    const scalar t260 = t5*t259*2.0;
    const scalar t262 = t4*t21*t41*t60;
    const scalar t263 = t2*t21*t50*t60;
    const scalar t275 = t2*t17*t50*t58*t187*2.0;
    const scalar t264 = t238+t263-t275;
    const scalar t265 = t5*t264*2.0;
    const scalar t267 = t3*t21*t50*t60;
    const scalar t268 = t4*t21*t50*t60;
    const scalar t274 = t4*t17*t50*t58*t187*2.0;
    const scalar t269 = t193-t232+t239+t268-t274;
    const scalar t270 = t19*t269*2.0;
    const scalar t271 = t58*t58;
    const scalar t276 = t17*t21*t58*2.0;
    const scalar t277 = t20*t73*2.0;
    const scalar t279 = t18*t174*2.0;
    const scalar t281 = t19*t181*2.0;
    const scalar t282 = t5*t174*2.0;
    const scalar t283 = t37*t74*2.0;
    const scalar t284 = t3*t4*t10*t75*2.0;
    const scalar t285 = t5*t181*2.0;
    const scalar t286 = t39*t75*2.0;
    const scalar t287 = t3*t4*t10*t74*2.0;
    const scalar t288 = t46*t74*2.0;
    const scalar t289 = t48*t75*2.0;
    const scalar t290 = t3*t21*t41*t43;
    const scalar t291 = t3*t11*t21*t41;
    const scalar t326 = t3*t17*t186*t187*2.0;
    const scalar t292 = t173+t189+t290+t291-t326;
    const scalar t293 = t18*t292*2.0;
    const scalar t294 = t4*t21*t41*t43;
    const scalar t295 = t4*t11*t21*t41;
    const scalar t327 = t4*t17*t186*t187*2.0;
    const scalar t296 = t180+t190+t294+t295-t327;
    const scalar t297 = t19*t296*2.0;
    const scalar t298 = t2*t21*t41*t43;
    const scalar t299 = t2*t11*t21*t41;
    const scalar t300 = t56*t75*2.0;
    const scalar t301 = t3*t11*t21*t50;
    const scalar t314 = t10*t11;
    const scalar t302 = t188+t224-t233+t301-t314;
    const scalar t303 = t18*t302*2.0;
    const scalar t304 = t4*t11*t21*t50;
    const scalar t305 = t218-t234+t304;
    const scalar t306 = t19*t305*2.0;
    const scalar t307 = t2*t11*t21*t50;
    const scalar t308 = t193+t221-t232-t235+t307;
    const scalar t309 = t5*t308*2.0;
    const scalar t310 = t54*t74*2.0;
    const scalar t311 = t31*(t300+t303+t306+t309+t310-t73*t115*2.0)*(1.0/2.0);
    const scalar t312 = t311-t32*t127*t158*(1.0/4.0);
    const scalar t313 = t64*t74*2.0;
    const scalar t315 = t4*t11*t21*t58;
    const scalar t316 = t3*t11*t21*t58;
    const scalar t317 = t255-t272+t316;
    const scalar t318 = t18*t317*2.0;
    const scalar t319 = t2*t11*t21*t58;
    const scalar t320 = t197+t258-t266-t273+t319;
    const scalar t321 = t5*t320*2.0;
    const scalar t322 = t62*t75*2.0;
    const scalar t323 = t125+t126-t162;
    const scalar t324 = t70-t122;
    const scalar t325 = t71-t124;
    const scalar t328 = t4*t13*t21*t41;
    const scalar t329 = t2*t13*t21*t41;
    const scalar t330 = t3*t13*t21*t41;
    const scalar t331 = t3*t15*t21*t41;
    const scalar t332 = t2*t15*t21*t41;
    const scalar t333 = t4*t15*t21*t41;
    const scalar t335 = t18*t209*2.0;
    const scalar t336 = t2*t4*t10*t83*2.0;
    const scalar t337 = t30*t32*t135*(1.0/4.0);
    const scalar t338 = t207-t208;
    const scalar t339 = t18*t338*2.0;
    const scalar t340 = t37*t80*2.0;
    const scalar t341 = t5*t209*2.0;
    const scalar t343 = t19*t213*2.0;
    const scalar t344 = t2*t3*t10*t82*2.0;
    const scalar t345 = t32*t100*t135*(1.0/4.0);
    const scalar t346 = t18*t213*2.0;
    const scalar t347 = t39*t83*2.0;
    const scalar t348 = t32*t106*t135*(1.0/4.0);
    const scalar t349 = t48*t83*2.0;
    const scalar t350 = t32*t135*t156*(1.0/4.0);
    const scalar t351 = t56*t83*2.0;
    const scalar t352 = t2*t21*t50*t52;
    const scalar t353 = t2*t13*t21*t50;
    const scalar t385 = t2*t17*t187*t230*2.0;
    const scalar t354 = t191+t201+t352+t353-t385;
    const scalar t355 = t5*t354*2.0;
    const scalar t356 = t4*t21*t50*t52;
    const scalar t357 = t4*t13*t21*t50;
    const scalar t386 = t4*t17*t187*t230*2.0;
    const scalar t358 = t180+t190+t356+t357-t386;
    const scalar t359 = t19*t358*2.0;
    const scalar t360 = t3*t21*t50*t52;
    const scalar t361 = t3*t13*t21*t50;
    const scalar t367 = t10*t13;
    const scalar t384 = t3*t17*t187*t230*2.0;
    const scalar t362 = t173+t189-t232+t236+t360+t361-t367-t384;
    const scalar t363 = t18*t362*2.0;
    const scalar t364 = t31*(t351+t355+t359+t363-t54*t80*2.0-t82*t115*2.0)*(1.0/2.0);
    const scalar t365 = t32*t135*t158*(1.0/4.0);
    const scalar t366 = t364+t365;
    const scalar t368 = t4*t13*t21*t58;
    const scalar t369 = t2*t13*t21*t58;
    const scalar t370 = t263-t275+t369;
    const scalar t371 = t5*t370*2.0;
    const scalar t372 = t3*t13*t21*t58;
    const scalar t389 = t3*t17*t50*t58*t187*2.0;
    const scalar t373 = t197-t266+t267+t372-t389;
    const scalar t374 = t18*t373*2.0;
    const scalar t375 = t62*t83*2.0;
    const scalar t376 = t32*t135*t160*(1.0/4.0);
    const scalar t377 = -t234+t304+t328;
    const scalar t378 = t19*t377*2.0;
    const scalar t379 = t188-t233+t301-t314+t330;
    const scalar t380 = t18*t379*2.0;
    const scalar t381 = t74*t80*2.0;
    const scalar t382 = t32*t135*(t125+t126-t162)*(1.0/4.0);
    const scalar t383 = t78-t132;
    const scalar t387 = t2*t15*t21*t50;
    const scalar t388 = t3*t15*t21*t50;
    const scalar t390 = t4*t15*t21*t50;
    const scalar t392 = t19*t251*2.0;
    const scalar t393 = t2*t3*t10*t91*2.0;
    const scalar t394 = t30*t32*t140*(1.0/4.0);
    const scalar t396 = t19*t252*2.0;
    const scalar t397 = t37*t91*2.0;
    const scalar t398 = t32*t100*t140*(1.0/4.0);
    const scalar t399 = t231-t250;
    const scalar t400 = t19*t399*2.0;
    const scalar t401 = t39*t88*2.0;
    const scalar t402 = t5*t251*2.0;
    const scalar t403 = t18*t252*2.0;
    const scalar t404 = t2*t4*t10*t90*2.0;
    const scalar t405 = t32*t106*t140*(1.0/4.0);
    const scalar t406 = t46*t91*2.0;
    const scalar t407 = t32*t140*t156*(1.0/4.0);
    const scalar t408 = t54*t91*2.0;
    const scalar t409 = t32*t140*t158*(1.0/4.0);
    const scalar t410 = t64*t91*2.0;
    const scalar t411 = t2*t21*t58*t60;
    const scalar t412 = t2*t15*t21*t58;
    const scalar t441 = t2*t17*t187*t271*2.0;
    const scalar t413 = t191+t201+t411+t412-t441;
    const scalar t414 = t5*t413*2.0;
    const scalar t415 = t3*t21*t58*t60;
    const scalar t416 = t3*t15*t21*t58;
    const scalar t442 = t3*t17*t187*t271*2.0;
    const scalar t417 = t173+t189+t415+t416-t442;
    const scalar t418 = t18*t417*2.0;
    const scalar t419 = t4*t21*t58*t60;
    const scalar t420 = t4*t15*t21*t58;
    const scalar t428 = t10*t15;
    const scalar t440 = t4*t17*t187*t271*2.0;
    const scalar t421 = t180+t190-t266+t276+t419+t420-t428-t440;
    const scalar t422 = t19*t421*2.0;
    const scalar t423 = t31*(t410+t414+t418+t422-t62*t88*2.0-t90*t120*2.0)*(1.0/2.0);
    const scalar t424 = t32*t140*t160*(1.0/4.0);
    const scalar t425 = t423+t424;
    const scalar t426 = -t272+t316+t331;
    const scalar t427 = t18*t426*2.0;
    const scalar t429 = t188-t314+t315+t333-t4*t17*t41*t58*t187*2.0;
    const scalar t430 = t19*t429*2.0;
    const scalar t431 = t75*t88*2.0;
    const scalar t432 = t32*t140*(t125+t126-t162)*(1.0/4.0);
    const scalar t433 = -t275+t369+t387;
    const scalar t434 = t5*t433*2.0;
    const scalar t435 = t193-t274-t367+t368+t390;
    const scalar t436 = t19*t435*2.0;
    const scalar t437 = t83*t88*2.0;
    const scalar t438 = t80*t91*2.0;
    const scalar t439 = t86-t136;
    hessE(0, 0) = t31*(t102+t153+(t20*t20)*2.0)*(1.0/2.0)-(t30*t30)*t32*(1.0/4.0);
    hessE(0, 1) = t101;
    hessE(0, 2) = t146;
    hessE(0, 3) = t31*(t67+t68+t69+t170+t171+t20*(t44+t2*t10*(x0(0)-x1(0)*2.0+x2(0))-t2*t17*t21*t41+1.0)*2.0-t2*t3*t18*t21*t41*2.0-t2*t4*t19*t21*t41*2.0)*(1.0/2.0)-t30*t32*(t110+t111+t5*(t44+t107-t2*t17*t21*t41+1.0)*2.0)*(1.0/4.0);
    hessE(0, 4) = t31*(t76+t202+t204+t20*(t2*t10*(x0(1)-x1(1)*2.0+x2(1))-t2*t17*t21*t50)*2.0-t5*t6*t21*t50*2.0-t2*t3*t18*t21*t50*2.0-t2*t4*t19*t21*t50*2.0)*(1.0/2.0)-t30*t32*(t114+t116+t5*(t112-t2*t17*t21*t50)*2.0)*(1.0/4.0);
    hessE(0, 5) = t31*(t84+t240+t242+t20*(t2*t10*(x0(2)-x1(2)*2.0+x2(2))-t2*t17*t21*t58)*2.0-t5*t6*t21*t58*2.0-t2*t3*t18*t21*t58*2.0-t2*t4*t19*t21*t58*2.0)*(1.0/2.0)-t30*t32*(t119+t121+t5*(t117-t2*t17*t21*t58)*2.0)*(1.0/4.0);
    hessE(0, 6) = t31*(t67+t68+t69+t277-t2*t3*t10*t74*2.0-t2*t4*t10*t75*2.0-t2*t3*t18*t21*t41*2.0-t2*t4*t19*t21*t41*2.0)*(-1.0/2.0)-t30*t32*t127*(1.0/4.0);
    hessE(0, 7) = t337+t31*(-t76+t203+t205+t336-t20*t82*2.0-t2*t3*t10*t80*2.0+t2*t3*t18*t21*t50*2.0)*(1.0/2.0);
    hessE(0, 8) = t394+t31*(-t84+t241+t243+t393-t20*t90*2.0-t2*t4*t10*t88*2.0+t2*t4*t19*t21*t58*2.0)*(1.0/2.0);
    hessE(1, 0) = t101;
    hessE(1, 1) = t31*(t102+t154+(t37*t37)*2.0)*(1.0/2.0)-t32*(t100*t100)*(1.0/4.0);
    hessE(1, 2) = t152;
    hessE(1, 3) = t31*(t123+t172+t175+t177-t7*t18*t21*t41*2.0-t2*t3*t5*t21*t41*2.0-t3*t4*t19*t21*t41*2.0)*(1.0/2.0)-t32*t100*t156*(1.0/4.0);
    hessE(1, 4) = t31*(t69+t130+t131+t206+t210+t211-t2*t3*t5*t21*t50*2.0-t3*t4*t19*t21*t50*2.0)*(1.0/2.0)-t32*t100*t158*(1.0/4.0);
    hessE(1, 5) = t31*(t137+t244+t245+t247-t7*t18*t21*t58*2.0-t2*t3*t5*t21*t58*2.0-t3*t4*t19*t21*t58*2.0)*(1.0/2.0)-t32*t100*t160*(1.0/4.0);
    hessE(1, 6) = t31*(-t123+t176+t178+t283+t284-t2*t3*t10*t73*2.0+t2*t3*t5*t21*t41*2.0)*(1.0/2.0)-t32*t100*t127*(1.0/4.0);
    hessE(1, 7) = t345-t31*(t69+t130+t131+t340+t344-t3*t4*t10*t83*2.0-t2*t3*t5*t21*t50*2.0-t3*t4*t19*t21*t50*2.0)*(1.0/2.0);
    hessE(1, 8) = t398+t31*(-t137+t246+t248+t397-t2*t3*t10*t90*2.0-t3*t4*t10*t88*2.0+t3*t4*t19*t21*t58*2.0)*(1.0/2.0);
    hessE(2, 0) = t146;
    hessE(2, 1) = t152;
    hessE(2, 2) = t31*(t153+t154+(t39*t39)*2.0)*(1.0/2.0)-t32*(t106*t106)*(1.0/4.0);
    hessE(2, 3) = t31*(t161+t179+t182+t184-t8*t19*t21*t41*2.0-t2*t4*t5*t21*t41*2.0-t3*t4*t18*t21*t41*2.0)*(1.0/2.0)-t32*t106*t156*(1.0/4.0);
    hessE(2, 4) = t31*(t163+t212+t214+t216-t8*t19*t21*t50*2.0-t2*t4*t5*t21*t50*2.0-t3*t4*t18*t21*t50*2.0)*(1.0/2.0)-t32*t106*t158*(1.0/4.0);
    hessE(2, 5) = t31*(t68+t131+t167+t249+t253+t254-t2*t4*t5*t21*t58*2.0-t3*t4*t18*t21*t58*2.0)*(1.0/2.0)-t32*t106*t160*(1.0/4.0);
    hessE(2, 6) = t31*(-t161+t183+t185+t286+t287-t2*t4*t10*t73*2.0+t2*t4*t5*t21*t41*2.0)*(1.0/2.0)-t32*t106*t127*(1.0/4.0);
    hessE(2, 7) = t348+t31*(-t163+t215+t217+t347-t3*t4*t10*t80*2.0-t2*t4*t10*t82*2.0+t3*t4*t18*t21*t50*2.0)*(1.0/2.0);
    hessE(2, 8) = t405-t31*(t68+t131+t167+t401+t404-t3*t4*t10*t91*2.0-t2*t4*t5*t21*t58*2.0-t3*t4*t18*t21*t58*2.0)*(1.0/2.0);
    hessE(3, 0) = t31*(t170+t171+t279+t281+t20*t109*2.0-t5*(t169-t2*t10*2.0)*2.0)*(1.0/2.0)-t30*t32*t156*(1.0/4.0);
    hessE(3, 1) = t31*(t172+t175-t176+t177-t178+t282)*(1.0/2.0)-t32*t100*t156*(1.0/4.0);
    hessE(3, 2) = t31*(t179+t182-t183+t184-t185+t285)*(1.0/2.0)-t32*t106*t156*(1.0/4.0);
    hessE(3, 3) = t32*(t156*t156)*(-1.0/4.0)-t31*(t5*(t191+t192+t229-t10*t43*2.0+t2*t21*t41*t43*2.0-t2*t17*t186*t187*2.0)*2.0+t18*(t189+t207+t3*t21*t41*t43*2.0-t3*t17*t186*t187*2.0)*2.0+t19*(t190+t231+t4*t21*t41*t43*2.0-t4*t17*t186*t187*2.0)*2.0-(t46*t46)*2.0-(t48*t48)*2.0-(t109*t109)*2.0)*(1.0/2.0);
    hessE(3, 4) = t228;
    hessE(3, 5) = t31*(t257+t260-t46*t64*2.0-t48*t62*2.0-t109*t120*2.0+t19*(t188+t200+t262-t10*t43-t4*t17*t41*t58*t187*2.0)*2.0)*(-1.0/2.0)-t32*t156*t160*(1.0/4.0);
    hessE(3, 6) = t31*(t288+t289+t293+t297+t5*(t191+t192+t201+t298+t299-t10*t11-t10*t43-t2*t17*t186*t187*2.0)*2.0-t73*t109*2.0)*(1.0/2.0)-t32*t127*t156*(1.0/4.0);
    hessE(3, 7) = t350+t31*(t349-t46*t80*2.0-t82*t109*2.0+t5*(t193+t194+t329-t10*t13-t2*t17*t41*t50*t187*2.0)*2.0+t18*(t188+t196+t330-t10*t43-t3*t17*t41*t50*t187*2.0)*2.0+t19*(t195+t328-t4*t17*t41*t50*t187*2.0)*2.0)*(1.0/2.0);
    hessE(3, 8) = t407+t31*(t406-t48*t88*2.0-t90*t109*2.0+t5*(t197+t198+t332-t10*t15-t2*t17*t41*t58*t187*2.0)*2.0+t19*(t188+t200+t333-t10*t43-t4*t17*t41*t58*t187*2.0)*2.0+t18*(t199+t331-t3*t17*t41*t58*t187*2.0)*2.0)*(1.0/2.0);
    hessE(4, 0) = t31*(t202-t203+t204-t205+t335+t20*t115*2.0)*(1.0/2.0)-t30*t32*t158*(1.0/4.0);
    hessE(4, 1) = t31*(t206+t210+t211+t339+t341+t343)*(1.0/2.0)-t32*t100*t158*(1.0/4.0);
    hessE(4, 2) = t31*(t212+t214-t215+t216-t217+t346)*(1.0/2.0)-t32*t106*t158*(1.0/4.0);
    hessE(4, 3) = t228;
    hessE(4, 4) = t32*(t158*t158)*(-1.0/4.0)-t31*(t18*(t189+t207+t236-t10*t52*2.0+t3*t21*t50*t52*2.0-t3*t17*t187*t230*2.0)*2.0+t5*(t191+t229+t2*t21*t50*t52*2.0-t2*t17*t187*t230*2.0)*2.0+t19*(t190+t231+t4*t21*t50*t52*2.0-t4*t17*t187*t230*2.0)*2.0-(t54*t54)*2.0-(t56*t56)*2.0-(t115*t115)*2.0)*(1.0/2.0);
    hessE(4, 5) = t31*(t265+t270-t54*t64*2.0-t56*t62*2.0-t115*t120*2.0+t18*(t197+t237+t267-t10*t60-t3*t17*t50*t58*t187*2.0)*2.0)*(-1.0/2.0)-t32*t158*t160*(1.0/4.0);
    hessE(4, 6) = t312;
    hessE(4, 7) = t366;
    hessE(4, 8) = t409+t31*(t408+t19*(t193-t232+t239+t390-t4*t17*t50*t58*t187*2.0)*2.0-t56*t88*2.0-t90*t115*2.0+t18*(t197+t237+t388-t10*t15-t3*t17*t50*t58*t187*2.0)*2.0+t5*(t238+t387-t2*t17*t50*t58*t187*2.0)*2.0)*(1.0/2.0);
    hessE(5, 0) = t31*(t240-t241+t242-t243+t392+t20*t120*2.0)*(1.0/2.0)-t30*t32*t160*(1.0/4.0);
    hessE(5, 1) = t31*(t244+t245-t246+t247-t248+t396)*(1.0/2.0)-t32*t100*t160*(1.0/4.0);
    hessE(5, 2) = t31*(t249+t253+t254+t400+t402+t403)*(1.0/2.0)-t32*t106*t160*(1.0/4.0);
    hessE(5, 3) = t31*(t257+t260+t19*(t188+t200-t261+t262-t4*t17*t41*t58*t187*2.0)*2.0-t46*t64*2.0-t48*t62*2.0-t109*t120*2.0)*(-1.0/2.0)-t32*t156*t160*(1.0/4.0);
    hessE(5, 4) = t31*(t265+t270+t18*(t197+t237-t266+t267-t3*t17*t50*t58*t187*2.0)*2.0-t54*t64*2.0-t56*t62*2.0-t115*t120*2.0)*(-1.0/2.0)-t32*t158*t160*(1.0/4.0);
    hessE(5, 5) = t32*(t160*t160)*(-1.0/4.0)-t31*(t19*(t190+t231+t276-t10*t60*2.0+t4*t21*t58*t60*2.0-t4*t17*t187*t271*2.0)*2.0+t18*(t189+t207+t3*t21*t58*t60*2.0-t3*t17*t187*t271*2.0)*2.0+t5*(t191+t229+t2*t21*t58*t60*2.0-t2*t17*t187*t271*2.0)*2.0-(t62*t62)*2.0-(t64*t64)*2.0-(t120*t120)*2.0)*(1.0/2.0);
    hessE(5, 6) = t31*(t313+t318+t321+t322-t73*t120*2.0+t19*(t188+t262+t315-t10*t11-t4*t17*t41*t58*t187*2.0)*2.0)*(1.0/2.0)-t32*t127*t160*(1.0/4.0);
    hessE(5, 7) = t376+t31*(t371+t374+t375-t64*t80*2.0-t82*t120*2.0+t19*(t193+t268-t274+t368-t10*t13)*2.0)*(1.0/2.0);
    hessE(5, 8) = t425;
    hessE(6, 0) = t31*(t277+t279+t281-t5*(t169-t229)*2.0-t2*t3*t10*t74*2.0-t2*t4*t10*t75*2.0)*(-1.0/2.0)-t30*t32*t127*(1.0/4.0);
    hessE(6, 1) = t31*(t176+t178-t282+t283+t284-t2*t3*t10*t73*2.0)*(1.0/2.0)-t32*t100*t127*(1.0/4.0);
    hessE(6, 2) = t31*(t183+t185-t285+t286+t287-t2*t4*t10*t73*2.0)*(1.0/2.0)-t32*t106*t127*(1.0/4.0);
    hessE(6, 3) = t31*(t288+t289+t293+t297-t73*t109*2.0+t5*(t191+t192+t201-t261+t298+t299-t10*t11-t2*t17*t186*t187*2.0)*2.0)*(1.0/2.0)-t32*t127*t156*(1.0/4.0);
    hessE(6, 4) = t312;
    hessE(6, 5) = t31*(t313+t318+t321+t322+t19*(t188+t262-t314+t315-t4*t17*t41*t58*t187*2.0)*2.0-t73*t120*2.0)*(1.0/2.0)-t32*t127*t160*(1.0/4.0);
    hessE(6, 6) = t32*(t323*t323)*(-1.0/4.0)-t31*(t18*(t189-t326+t3*t11*t21*t41*2.0)*2.0+t19*(t190-t327+t4*t11*t21*t41*2.0)*2.0+t5*(t191+t192-t10*t11*2.0+t2*t11*t21*t41*2.0-t2*t17*t186*t187*2.0)*2.0-(t73*t73)*2.0-(t324*t324)*2.0-(t325*t325)*2.0)*(1.0/2.0);
    hessE(6, 7) = t382-t31*(t378+t380+t381-t73*t82*2.0-t75*t83*2.0+t5*(t193-t235+t307+t329-t10*t13)*2.0)*(1.0/2.0);
    hessE(6, 8) = t432-t31*(t427+t430+t431-t73*t90*2.0-t74*t91*2.0+t5*(t197-t273+t319+t332-t10*t15)*2.0)*(1.0/2.0);
    hessE(7, 0) = t337+t31*(t203+t205-t335+t336-t20*t82*2.0-t2*t3*t10*t80*2.0)*(1.0/2.0);
    hessE(7, 1) = t345-t31*(t339+t340+t341+t343+t344-t3*t4*t10*t83*2.0)*(1.0/2.0);
    hessE(7, 2) = t348+t31*(t215+t217-t346+t347-t3*t4*t10*t80*2.0-t2*t4*t10*t82*2.0)*(1.0/2.0);
    hessE(7, 3) = t350+t31*(t349+t18*(t188+t196-t233-t261+t330)*2.0-t46*t80*2.0-t82*t109*2.0+t5*(t193+t194-t235+t329-t10*t13)*2.0+t19*(t195-t234+t328)*2.0)*(1.0/2.0);
    hessE(7, 4) = t366;
    hessE(7, 5) = t376+t31*(t371+t374+t375+t19*(t193+t268-t274-t367+t368)*2.0-t64*t80*2.0-t82*t120*2.0)*(1.0/2.0);
    hessE(7, 6) = t382-t31*(t378+t380+t381+t5*(t193-t235+t307+t329-t367)*2.0-t73*t82*2.0-t75*t83*2.0)*(1.0/2.0);
    hessE(7, 7) = t31*(t18*(t189+t236-t384-t10*t13*2.0+t3*t13*t21*t50*2.0)*2.0+t5*(t191-t385+t2*t13*t21*t50*2.0)*2.0+t19*(t190-t386+t4*t13*t21*t50*2.0)*2.0-(t80*t80)*2.0-(t82*t82)*2.0-(t383*t383)*2.0)*(-1.0/2.0)-t32*(t135*t135)*(1.0/4.0);
    hessE(7, 8) = t31*(t434+t436+t437+t438-t82*t90*2.0+t18*(t197+t372+t388-t389-t10*t15)*2.0)*(-1.0/2.0)-t32*t135*t140*(1.0/4.0);
    hessE(8, 0) = t394+t31*(t241+t243-t392+t393-t20*t90*2.0-t2*t4*t10*t88*2.0)*(1.0/2.0);
    hessE(8, 1) = t398+t31*(t246+t248-t396+t397-t2*t3*t10*t90*2.0-t3*t4*t10*t88*2.0)*(1.0/2.0);
    hessE(8, 2) = t405-t31*(t400+t401+t402+t403+t404-t3*t4*t10*t91*2.0)*(1.0/2.0);
    hessE(8, 3) = t407+t31*(t406+t19*(t188+t200-t261+t333-t4*t17*t41*t58*t187*2.0)*2.0-t48*t88*2.0-t90*t109*2.0+t5*(t197+t198-t273+t332-t10*t15)*2.0+t18*(t199-t272+t331)*2.0)*(1.0/2.0);
    hessE(8, 4) = t409+t31*(t408+t19*(t193-t232+t239-t274+t390)*2.0-t56*t88*2.0-t90*t115*2.0+t18*(t197+t237+t388-t389-t10*t15)*2.0+t5*(t238-t275+t387)*2.0)*(1.0/2.0);
    hessE(8, 5) = t425;
    hessE(8, 6) = t432-t31*(t427+t430+t431+t5*(t197-t273+t319+t332-t428)*2.0-t73*t90*2.0-t74*t91*2.0)*(1.0/2.0);
    hessE(8, 7) = t31*(t434+t436+t437+t438+t18*(t197+t372+t388-t389-t428)*2.0-t82*t90*2.0)*(-1.0/2.0)-t32*t135*t140*(1.0/4.0);
    hessE(8, 8) = t31*(t19*(t190+t276-t440-t10*t15*2.0+t4*t15*t21*t58*2.0)*2.0+t5*(t191-t441+t2*t15*t21*t58*2.0)*2.0+t18*(t189-t442+t3*t15*t21*t58*2.0)*2.0-(t88*t88)*2.0-(t90*t90)*2.0-(t439*t439)*2.0)*(-1.0/2.0)-t32*(t140*t140)*(1.0/4.0);
  }
  
  template<typename T>
  inline Eigen::Matrix<T, 1, 2> cross2Mat(const Eigen::Matrix<T, 2, 1>& r)
  {
    Eigen::Matrix<T, 1, 2> ret;
    ret(0) = -r(1);
    ret(1) = r(0);
    return ret;
  }
  
  template<typename T>
  inline T det3(const Eigen::Matrix<T, 3, 1>& a,
                    const Eigen::Matrix<T, 3, 1>& v1,
                    const Eigen::Matrix<T, 3, 1>& v2)
  {
    return a[0] * (v1[1] * v2[2] - v1[2] * v2[1]) + a[1] * (v1[2] * v2[0] - v1[0] * v2[2]) + a[2] * (v1[0] * v2[1] - v1[1] * v2[0]);
  }
  
  template<typename T>
  T fraction_inside(const T& phi_left, const T& phi_right)
  {
    if(phi_left < (T) 0.0 && phi_right < (T) 0.0)
      return (T) 1.0;
    if (phi_left < (T) 0.0 && phi_right >= (T) 0.0)
      return phi_left / (phi_left - phi_right);
    if(phi_left >= (T) 0.0 && phi_right < (T) 0.0)
      return phi_right / (phi_right - phi_left);
    else
      return (T) 0.0;
  }
    
    inline void cycle_array(scalar* arr, int size) {
        scalar t = arr[0];
        for(int i = 0; i < size-1; ++i)
            arr[i] = arr[i+1];
        arr[size-1] = t;
    }
    
    inline scalar fraction_inside(const scalar& phi_bl, const scalar& phi_br,
                                  const scalar& phi_tl, const scalar& phi_tr) {
        
        int inside_count = (phi_bl<0?1:0) + (phi_tl<0?1:0) + (phi_br<0?1:0) + (phi_tr<0?1:0);
        scalar list[] = { phi_bl, phi_br, phi_tr, phi_tl };
        
        if(inside_count == 4)
            return 1.;
        else if (inside_count == 3) {
            //rotate until the positive value is in the first position
            while(list[0] < 0) {
                cycle_array(list,4);
            }
            
            //Work out the area of the exterior triangle
            scalar side0 = 1.-fraction_inside(list[0], list[3]);
            scalar side1 = 1.-fraction_inside(list[0], list[1]);
            return 1. - 0.5*side0*side1;
        }
        else if(inside_count == 2) {
            
            //rotate until a negative value is in the first position, and the next negative is in either slot 1 or 2.
            while(list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) {
                cycle_array(list,4);
            }
            
            if(list[1] < 0) { //the matching signs are adjacent
                scalar side_left = fraction_inside(list[0], list[3]);
                scalar side_right = fraction_inside(list[1], list[2]);
                return  0.5*(side_left + side_right);
            }
            else { //matching signs are diagonally opposite
                //determine the centre point's sign to disambiguate this case
                scalar middle_point = 0.25*(list[0] + list[1] + list[2] + list[3]);
                if(middle_point < 0) {
                    scalar area = 0.;
                    
                    //first triangle (top left)
                    scalar side1 = 1.-fraction_inside(list[0], list[3]);
                    scalar side3 = 1.-fraction_inside(list[2], list[3]);
                    
                    area += 0.5*side1*side3;
                    
                    //second triangle (top right)
                    scalar side2 = 1-fraction_inside(list[2], list[1]);
                    scalar side0 = 1-fraction_inside(list[0], list[1]);
                    area += 0.5*side0*side2;
                    
                    return 1.-area;
                }
                else {
                    scalar area = 0.;
                    
                    //first triangle (bottom left)
                    scalar side0 = fraction_inside(list[0], list[1]);
                    scalar side1 = fraction_inside(list[0], list[3]);
                    area += 0.5*side0*side1;
                    
                    //second triangle (top right)
                    scalar side2 = fraction_inside(list[2], list[1]);
                    scalar side3 = fraction_inside(list[2], list[3]);
                    area += 0.5*side2*side3;
                    return area;
                }
                
            }
        }
        else if(inside_count == 1) {
            //rotate until the negative value is in the first position
            while(list[0] >= 0) {
                cycle_array(list,4);
            }
            
            //Work out the area of the interior triangle, and subtract from 1.
            scalar side0 = fraction_inside(list[0], list[3]);
            scalar side1 = fraction_inside(list[0], list[1]);
            return 0.5*side0*side1;
        }
        else
            return 0;
        
    }
  
  template<typename T>
  inline void closestDistanceBetweenLines(const Eigen::Matrix<T, 3, 1>& a0,
                                          const Eigen::Matrix<T, 3, 1>& a1,
                                          const Eigen::Matrix<T, 3, 1>& b0,
                                          const Eigen::Matrix<T, 3, 1>& b1,
                                          Eigen::Matrix<T, 3, 1>& pA,
                                          Eigen::Matrix<T, 3, 1>& pB,
                                          T& d,
                                          bool clampAll = true, bool clampA0 = true, bool clampA1 = true, bool clampB0 = true, bool clampB1 = true)
  {
    //Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
    //Return distance, the two closest points, and their average
    if(clampAll){
      clampA0 = true;
      clampA1 = true;
      clampB0 = true;
      clampB1 = true;
    }
    
    //Calculate denomitator
    auto A = a1 - a0;
    auto B = b1 - b0;
    auto _A = A.normalized();
    auto _B = B.normalized();
    auto cross = _A.cross(_B);
    auto denom = cross.squaredNorm();
    auto nA = A.norm();
    auto nB = B.norm();
    
    //If denominator is 0, lines are parallel: Calculate distance with a projection and evaluate clamp edge cases
    if (denom == 0){
      auto d0 = _A.dot(b0 - a0);
      
      d = (_A * d0 + a0 - b0).norm();
      

      
      //If clamping: the only time we'll get closest points will be when lines don't overlap at all. Find if segments overlap using dot products.
      if(clampA0 || clampA1 || clampB0 || clampB1){
        auto d1 = _A.dot(b1 - a0);
        
        //Is segment B before A?
        if(d0 <= 0 && 0 >= d1){
          if(clampA0 && clampB1){
            if(fabs(d0) < fabs(d1)){
              pA = b0, pB = a0, d = (b0 - a0).norm();
              return;
            }
            pA = b1, pB = a0, d = (b1 - a0).norm();
            return;
          }
        }
        //Is segment B after A?
        else if(d0 >= nA && nA <= d1){
          if(clampA1 && clampB0){
            if(fabs(d0) < fabs(d1)){
              pA = b0, pB = a1, d = (b0 - a1).norm();
              return;
            }
            pA = b1, pB = a1, d = (b1 - a1).norm();
          }
        }
      }
      
      //If clamping is off, or segments overlapped, we have infinite results, just return position.
      return;
    }
    
    //Lines criss-cross: Calculate the dereminent and return points
    auto t = b0 - a0;
    auto det0 = det3<T>(t, _B, cross);
    auto det1 = det3<T>(t, _A, cross);
    
    auto t0 = det0 / denom;
    auto t1 = det1 / denom;
    
    pA = a0 + _A * t0;
    pB = b0 + _B * t1;
    
    //Clamp results to line segments if needed
    if(clampA0 || clampA1 || clampB0 || clampB1){
      if(t0 < 0 && clampA0)
        pA = a0;
      else if(t0 > nA && clampA1)
        pA = a1;
      
      if(t1 < 0 && clampB0)
        pB = b0;
      else if(t1 > nB && clampB1)
        pB = b1;
    }
    
    d = (pA - pB).norm();
  }

  template<typename T>
  inline void closestDistanceBetweenLines(const Eigen::Matrix<T, 3, 1>& a0,
                                          const Eigen::Matrix<T, 3, 1>& a1,
                                          const Eigen::Matrix<T, 3, 1>& b0,
                                          const Eigen::Matrix<T, 3, 1>& b1,
                                          T& d)
  {
    Eigen::Matrix<T, 3, 1> pA, pB;
    closestDistanceBetweenLines<T>(a0, a1, b0, b1, pA, pB, d);
  }
  
  template<typename T>
  inline void closestDistanceBetweenLines(const Eigen::Matrix<T, 2, 1>& a0,
                                          const Eigen::Matrix<T, 2, 1>& a1,
                                          const Eigen::Matrix<T, 2, 1>& b0,
                                          const Eigen::Matrix<T, 2, 1>& b1,
                                          T& d)
  {
    Eigen::Matrix<T, 3, 1> pA, pB;
    Eigen::Matrix<T, 3, 1> a0_(a0(0), a0(1), 0.0);
    Eigen::Matrix<T, 3, 1> a1_(a1(0), a1(1), 0.0);
    Eigen::Matrix<T, 3, 1> b0_(b0(0), b0(1), 0.0);
    Eigen::Matrix<T, 3, 1> b1_(b1(0), b1(1), 0.0);
    
    closestDistanceBetweenLines(a0_, a1_, b0_, b1_, pA, pB, d);
  }
  
  template<typename S>
  inline S gamma_curved(const S& x, const S& power)
  {
    return 1.0 - exp(-x * power);
  }
  
  inline scalar compute_avg( const VectorXs& x )
  {
    return x.lpNorm<1>() / (scalar) x.size();
  }
  
  // DDG operators
  template<int DIM>
  inline scalar compute_edge_area( const Vectors<DIM>& e0, const Vectors<DIM>& e1 )
  {
    return (e1 - e0).norm();
  }
  
  template<int DIM>
  inline void compute_edge_area( const VectorXs& x, const std::vector<std::pair<int, int> >& edges, VectorXs& area, TwoDScene<DIM>* scene )
  {
    int ne = edges.size();
    for(int i = 0; i < ne; ++i)
    {
      area(i) = compute_edge_area<DIM>(x.segment<DIM>( scene->getDof( edges[i].first ) ), x.segment<DIM>( scene->getDof(edges[i].second) ));
    }
  }
  
  template<int DIM>
  inline void compute_edge_dir( const VectorXs& x, const std::vector<std::pair<int, int> >& edges, MatrixXs& dir, TwoDScene<DIM>* scene )
  {
    int ne = edges.size();
    for(int i = 0; i < ne; ++i)
    {
      const Vectors<DIM>& x0 = x.segment<DIM>( scene->getDof(edges[i].first) );
      const Vectors<DIM>& x1 = x.segment<DIM>( scene->getDof(edges[i].second) );
      dir.row(i) = (x1 - x0).normalized().transpose();
    }
  }
  
  inline scalar compute_vertex_area( const scalar& a0, const scalar& a1 )
  {
    return 0.5 * (a0 + a1);
  }
  
  inline void compute_radius_vector( const VectorXs& radii, const std::vector< int >& particles, VectorXs& rad_vec)
  {
    int np = rad_vec.size();
    for(int i = 0; i < np; ++i)
    {
      rad_vec(i) = radii(particles[i]);
    }
  }
  
  template<typename T, typename S>
  inline void compute_vertex_val(const T& u_k, const S& W_fv, T& u_vert)
  {
    u_vert = W_fv * u_k;
  }
  
  template<int DIM>
  inline void compute_vertex_area( const std::vector<std::pair<int, int> >& edges, const VectorXs& edge_area, VectorXs& vertex_area )
  {
    int ne = edges.size();
    vertex_area.setZero();
    
    for(int i = 0; i < ne; ++i)
    {
      int i0 = edges[i].first;
      int i1 = edges[i].second;
      vertex_area(i0) += edge_area(i) * 0.5;
      vertex_area(i1) += edge_area(i) * 0.5;
    }
  }
  
  template<int DIM>
  inline void compute_normalized_vector( MatrixXs& v )
  {
    int np = v.rows();
    for(int i = 0; i < np; ++i)
    {
      v.row(i).normalize();
    }
  }
  
  inline void compute_bracket_matrix( const VectorXs& vertex_value, SparseXs& Gv )
  {
    int np = vertex_value.size();
    
    TripletXs tri;
    
    tri.reserve(np);
    
    for(int i = 0; i < np; ++i)
    {
      tri.push_back(Triplets(i, i, vertex_value[i]));
    }
    
    Gv.setFromTriplets(tri.begin(), tri.end());
  }
  
  inline void compute_flat_bracket_matrix( const VectorXs& vertex_value, SparseXs& Gv )
  {
    int np = vertex_value.size();
    
    TripletXs tri;
    
    tri.reserve(np);
    
    for(int i = 0; i < np; ++i)
    {
      tri.push_back(Triplets(0, i, vertex_value[i]));
    }
    
    Gv.setFromTriplets(tri.begin(), tri.end());
  }
  
  inline void compute_bracket_matrix( const VectorXs& vertex_value, SparseXs& Gf, int dim )
  {
    int ne = vertex_value.size();
    
    TripletXs tri;
    
    tri.reserve(ne * dim);
    
    for(int i = 0; i < ne; ++i)
    {
      for(int r = 0; r < dim; ++r)
        tri.push_back(Triplets(i * dim + r, i * dim + r, vertex_value[i]));
    }
    
    Gf.setFromTriplets(tri.begin(), tri.end());
  }
  
  
  inline void compute_inv_bracket_matrix( const VectorXs& vertex_value, SparseXs& iGv, scalar min_val = -1e+20 )
  {
    int np = vertex_value.size();
    
    TripletXs tri;
    
    tri.resize(np);
    
    for(int i = 0; i < np; ++i)
    {
      tri[i] = Triplets(i, i, 1.0 / std::max(min_val, vertex_value[i]));
    }
    
    iGv.setFromTriplets(tri.begin(), tri.end());
  }
  
  inline void parallel_compute_inv_bracket_matrix( const VectorXs& vertex_value, SparseXs& iGv, scalar min_val = -1e+20 )
  {
    int np = vertex_value.size();
    
    TripletXs tri;
    
    tri.resize(np);
    
    threadutils::thread_pool::ParallelFor(0, np, [&] (int i) {
      tri[i] = Triplets(i, i, 1.0 / std::max(min_val, vertex_value[i]));
    });

    iGv.setFromTriplets(tri.begin(), tri.end());
  }
  
  
  inline void compute_vertex_divergence_matrix( const SparseXs& iGv, const SparseXs& Gf, const SparseXs& grad, SparseXs& div )
  {
    div = -iGv * grad.transpose() * Gf;
  }
  
  inline void compute_vertex_laplacian_matrix( const SparseXs& grad, const SparseXs& div, SparseXs& L)
  {
    L = div * grad;
  }
  
  inline void compute_edge_to_vertex_matrix( const std::vector<std::pair<int, int> >& edges, const VectorXs& vertex_area, const VectorXs& edge_area, SparseXs& W )
  {
    // W.resize(vertex_area.size(), edge_area.size());
    TripletXs tri;
    
    int ne = edge_area.size();
    for(int j = 0; j < ne; ++j)
    {
      const scalar& A_Fj = edge_area[j];
      const int& i0 = edges[j].first;
      const int& i1 = edges[j].second;
      const scalar& A_Vi0 = vertex_area[i0];
      const scalar& A_Vi1 = vertex_area[i1];
      
      tri.push_back(Triplets(i0, j, A_Fj / A_Vi0 * 0.5));
      tri.push_back(Triplets(i1, j, A_Fj / A_Vi1 * 0.5));
    }
    
    W.setFromTriplets(tri.begin(), tri.end());
  }
  
  inline void compute_edge_to_vertex_matrix( const std::vector<std::pair<int, int> >& edges, const VectorXs& vertex_area, const VectorXs& edge_area, const VectorXi& ppp_count, SparseXs& W )
  {
    // W.resize(vertex_area.size(), edge_area.size());
    TripletXs tri;
    
    int ne = edge_area.size();
    for(int j = 0; j < ne; ++j)
    {
      const scalar& A_Fj = edge_area[j];
      const int& i0 = edges[j].first;
      const int& i1 = edges[j].second;
      const scalar& A_Vi0 = vertex_area[i0];
      const scalar& A_Vi1 = vertex_area[i1];
      const int& N_Vi0 = ppp_count[i0];
      const int& N_Vi1 = ppp_count[i1];
      
      if(N_Vi0 > 0) tri.push_back(Triplets(i0, j, A_Fj / (A_Vi0 * (scalar) N_Vi0)));
      if(N_Vi1 > 0) tri.push_back(Triplets(i1, j, A_Fj / (A_Vi1 * (scalar) N_Vi1)));
    }
    
    W.setFromTriplets(tri.begin(), tri.end());
  }
  
  inline void inverse_diagonal( MatrixXs& A )
  {
    int np = A.rows();
    for(int i = 0; i < np; ++i)
    {
      A(i, i) = 1.0 / A(i, i);
    }
  }
  
  template<int DIM>
  inline void compute_mass_vector_solid( const VectorXs& rho_hair,
                                        const VectorXs& area_v,
                                        const std::vector<int>& indices,
                                        const VectorXs& radii,
                                        VectorXs& m,
                                        TwoDScene<DIM>* scene)
  {
    int np = indices.size();
    for(int i = 0; i < np; ++i)
    {
      int idx = indices[i];
      scalar masshair = rho_hair(idx) * M_PI * radii[idx] * radii[idx] * area_v(i);
      scalar total_mass = masshair;
      m.segment<DIM>( scene->getDof(idx) ).setConstant(total_mass);
    }
  }
  
  
  inline scalar compute_integrated_area( const scalar& r0, const scalar& r1 )
  {
    return M_PI / 8.0 * (3.0 * r0 * r0 + 2.0 * r0 * r1 + 3.0 * r1 * r1);
  }
  
  inline void compute_next_pressure(const SparseXs& L, const VectorXs& eta, const scalar& sigma, VectorXs& pressure)
  {
    pressure = -sigma * (L * eta);
  }
  
  
  inline void compute_D_matrix(const SparseXs& W_fv, const VectorXs& u_next, const SparseXs& gradF, const SparseXs& dir_f_exp, SparseXs& D)
  {
    SparseXs u_mat(u_next.size(), u_next.size());
    
    compute_bracket_matrix(u_next, u_mat);
    D = W_fv * u_mat * dir_f_exp * gradF;
  }
  
  
  inline void compute_A_matrix(const SparseXs& D, const SparseXs& divV, const VectorXs& u_next, const SparseXs& dir_f_exp, const SparseXs& L, const scalar& mu, SparseXs& A)
  {
    VectorXs divu = divV * dir_f_exp.transpose() * u_next;
    
    SparseXs divu_mat(D.rows(), D.cols());
    
    compute_bracket_matrix(divu, divu_mat);
    
    A = D + divu_mat - mu * L;
  }

  template<int DIM>
  inline void compute_mass_vector_totalliquid(
                                              const VectorXs& eta,
                                              const std::vector<int>& indices,
                                              const VectorXs& area_v,
                                              const VectorXs& radii,
                                              const scalar& rho,
                                              const scalar& H,
                                              VectorXs& m, 
                                              TwoDScene<DIM>* scene)
  {
    int np = indices.size();
    for(int i = 0; i < np; ++i)
    {
      int idx = indices[i];
      scalar realeta = eta(i) * H;
      scalar massliquid = rho * M_PI * std::max(0.0, realeta * realeta - radii[idx] * radii[idx]) * area_v(i);
      scalar total_mass = massliquid;
      m.segment<DIM>( scene->getDof(idx) ) += Vectors<DIM>::Ones * total_mass;
    }
  }
  
  template<int DIM>
  inline scalar compute_hair_length_geodesic( const std::vector<int>& hair, const VectorXs& x, int start_i, const scalar& start_alpha, int end_i, const scalar& end_alpha, TwoDScene<DIM>* scene )
  {
    start_i = std::min((int) hair.size() - 2, std::max(0, start_i));
    end_i = std::min((int) hair.size() - 1, std::max(0, end_i));
    
    int ipstart = hair[start_i];
    int ipend = hair[end_i];
    int ipnext = hair[start_i + 1];
    scalar cut_start_len = (x.segment<DIM>( scene->getDof(ipnext) ) - x.segment<DIM>( scene->getDof(ipstart) )).norm() * start_alpha;
    
    scalar cut_end_len = 0.0;
    if(end_i < (int) hair.size() - 1) {
      int inext = end_i + 1;
      cut_end_len = (x.segment<DIM>( scene->getDof(hair[inext]) ) - x.segment<DIM>( scene->getDof(ipend) )).norm() * (1.0 - end_alpha);
    }
    
    scalar accu_len = 0.0;
    
    for(int i = start_i; i <= end_i; ++i)
    {
      int inext = i + 1;
      if(inext >= (int) hair.size()) break;
      
      accu_len += (x.segment<DIM>( scene->getDof(hair[inext]) ) - x.segment<DIM>( scene->getDof(hair[i]) )).norm();
    }
    
    return accu_len - cut_start_len - cut_end_len;
  }
  
  
  inline scalar compute_point_edge_closest_alpha( const Vector2s& x1, const Vector2s& x2, const Vector2s& x3 )
  {
    scalar alpha = (x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2);
    return alpha;
  }

  
  inline Vector2s compute_point_edge_closest_vector( const Vector2s& x1, const Vector2s& x2, const Vector2s& x3 )
  {
    scalar alpha = (x1-x2).dot(x3-x2)/(x3-x2).dot(x3-x2);
    alpha = std::min(1.0, std::max(0.0, alpha));
    Vector2s closest = x2 + alpha*(x3-x2);
    Vector2s n = closest-x1;
    
    return n;
  }
  
  inline uint64 make_uint64(int i0, int i1)
  {
    return (((uint64) i0) << 32UL) | (uint64) i1;
  }
  
  inline void extract_uint64(const uint64& i, int& i0, int& i1)
  {
    i1 = (int)(i & 0xFFFFFFFFUL);
    i0 = (int)((i >> 32UL) & 0xFFFFFFFFUL);
  }
  
  
  inline void compute_cwiseProduct( VectorXs& b, const VectorXs& x, const VectorXs& y )
  {
    scalar* b_ptr = b.data();
    const scalar* x_ptr = x.data();
    const scalar* y_ptr = y.data();
    
    
    threadutils::thread_pool::ParallelFor(0, (int) x.size(), [&] (int k) {
      b_ptr[k] = x_ptr[k] * y_ptr[k];
    });
  }
  
  inline void accumulate_cwiseProduct( VectorXs& b, const VectorXs& x )
  {
    scalar* b_ptr = b.data();
    const scalar* x_ptr = x.data();
    
    threadutils::thread_pool::ParallelFor(0, (int) x.size(), [&] (int k) {
      b_ptr[k] *= x_ptr[k];
    });
  }
  
  inline void computeJTPhi_coeff_add( VectorXs& b, const VectorXs& a, const scalar& c, const VectorXs& x, const SparseXs& JT )
  {
    const scalar* a_ptr = a.data();
    scalar* b_ptr = b.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum * c + a_ptr[k];
    });
  }
  
  inline void computeJTPhi_coeff_add( VectorXs& b, const VectorXs& a, const VectorXs& c, const VectorXs& x, const SparseXs& JT )
  {
    const scalar* a_ptr = a.data();
    const scalar* c_ptr = c.data();
    scalar* b_ptr = b.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum * c_ptr[k] + a_ptr[k];
    });
  }
  
  inline void computeJTPhi_coeff( VectorXs& b, const VectorXs& c, const VectorXs& x, const SparseXs& JT )
  {
    scalar* b_ptr = b.data();
    const scalar* c_ptr = c.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum * c_ptr[k];
    });
  }
  
  inline void computeJTPhi_coeff( MatrixXs& b, const VectorXs& c, const VectorXs& x, const SparseXs& JT, int col_idx )
  {
    scalar* b_ptr = b.data() + col_idx * b.rows();
    const scalar* c_ptr = c.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum * c_ptr[k];
    });
  }
  
  inline void computeJTPhi( VectorXs& b, const VectorXs& x, const SparseXs& JT )
  {
    scalar* b_ptr = b.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum;
    });
  }
  
  inline void accumulateJTPhi_coeff( VectorXs& b, const VectorXs& c, const VectorXs& x, const SparseXs& JT )
  {
    scalar* b_ptr = b.data();
    const scalar* c_ptr = c.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] += sum * c_ptr[k];
    });
  }
  
  inline void computeJTPhi_coeff( VectorXs& b, const scalar& c, const VectorXs& x, const SparseXs& JT )
  {
    scalar* b_ptr = b.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum * c;
    });
  }
  
  inline void computeJTPhi_coeff( MatrixXs& b, const scalar& c, const VectorXs& x, const SparseXs& JT, int col_idx )
  {
    scalar* b_ptr = b.data() + col_idx * b.rows();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] = sum * c;
    });
  }
  
  inline void accumulateJTPhi_coeff( VectorXs& b, const scalar& c, const VectorXs& x, const SparseXs& JT )
  {
    scalar* b_ptr = b.data();
    threadutils::thread_pool::ParallelFor(0, (int) JT.outerSize(), [&] (int k) {
      scalar sum = 0.0;
      for (SparseXs::InnerIterator it(JT, k); it; ++it)
      {
        sum += x(it.row()) * it.value();
      }
      b_ptr[k] += sum * c;
    });
  }
  

}

#endif

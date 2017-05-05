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

#include "MathDefs.h"
#include <vector>
#include <tbb/tbb.h>

#ifndef _SORTER_H
#define _SORTER_H

class Sorter {
public:
  Sorter( int ni_, int nj_, int nk_ );
  void resize( int ni_, int nj_, int nk_ );
  ~Sorter();
  
  inline uint64_t hash(int i, int j, int k, int pidx)
  {
    int high_part = k * ni * nj + j * ni + i;
    return (uint64_t) high_part << 32UL | (uint64_t) pidx;
  }
  
  template<typename Callable>
  void sort( size_t total_size, Callable func )
  {
    if(array_idx.size() != total_size) {
      array_idx.resize(total_size);
    }
    
    memset(&array_sup[0], 0, array_sup.size() * sizeof(std::pair<int, int>));
    
    const int np = (int) total_size;
    
    tbb::parallel_for(0, np, 1, [&] (int pidx) {
      
      int i, j, k;
      func(pidx, i, j, k);
      array_idx[pidx] = hash(i, j, k, pidx);
    });

    tbb::parallel_sort(array_idx.begin(), array_idx.end());

    tbb::parallel_for(0, np, 1, [&] (int pidx) {
      unsigned int G_ID = (unsigned int) pidx;
      unsigned int G_ID_PREV = (G_ID == 0)? (unsigned int) np : G_ID; G_ID_PREV--;
      unsigned int G_ID_NEXT = G_ID + 1; if (G_ID_NEXT == (unsigned int) np) { G_ID_NEXT = 0; }
      
      unsigned int cell = (unsigned int) ( array_idx[G_ID] >> 32UL );
      unsigned int cell_prev = (unsigned int) ( array_idx[G_ID_PREV] >> 32UL );
      unsigned int cell_next = (unsigned int) ( array_idx[G_ID_NEXT] >> 32UL );
      if (cell != cell_prev)
      {
        // I'm the start of a cell
        array_sup[cell].first = G_ID;
      }
      if (cell != cell_next)
      {
        // I'm the end of a cell
        array_sup[cell].second = G_ID + 1;
      }
    });
    
    
  }
  
  template<typename Callable>
  void getNeigboringParticles_cell( int i, int j, int k, int wl, int wh, int hl, int hh, int dl, int dh, Callable func )
  {
    for(int sk=k+dl; sk <= k+dh; ++sk) for( int si=i+wl; si<=i+wh; si++ ) for( int sj=j+hl; sj<=j+hh; sj++ ) {
      getCellAt(si, sj, sk, func);
    }
  }
  
  template<typename Callable>
  void getCellAt( int i, int j, int k, Callable func )
  {
    if( i < 0 || i > ni-1 || j < 0 || j > nj-1 || k < 0 || k > nk-1 ) return;
    unsigned int G_CELL = (unsigned int) (k * (ni * nj) + j * ni + i);
    const std::pair<int, int>& G_START_END = array_sup[G_CELL];
    for(int N_ID = G_START_END.first; N_ID < G_START_END.second; ++N_ID)
    {
      func((int) (array_idx[N_ID] & 0xFFFFFFFFUL));
    }
  }
  
  std::vector<uint64_t> array_idx;
  std::vector< std::pair<int, int> > array_sup;
  
  int ni;
  int nj;
  int nk;
};

#endif

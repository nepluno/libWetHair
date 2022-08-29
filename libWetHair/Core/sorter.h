//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2017 Yun (Raymond) Fei, Henrique Teles Maia, Christopher Batty,
// Changxi Zheng, and Eitan Grinspun
//

#ifndef LIBWETHAIR_CORE_SORTER_H_
#define LIBWETHAIR_CORE_SORTER_H_

#include <tbb/tbb.h>

#include <vector>

#include "MathDefs.h"

namespace libwethair {

class Sorter {
 public:
  Sorter(int ni_, int nj_, int nk_);
  void resize(int ni_, int nj_, int nk_);
  ~Sorter();

  inline uint64_t hash(int i, int j, int k, int pidx) {
    int high_part = k * ni * nj + j * ni + i;
    return (uint64_t)high_part << 32UL | (uint64_t)pidx;
  }

  template <typename Callable>
  void sort(size_t total_size, Callable func) {
    if (array_idx.size() != total_size) {
      array_idx.resize(total_size);
    }

    memset(&array_sup[0], 0, array_sup.size() * sizeof(std::pair<int, int>));

    const int np = (int)total_size;

    tbb::parallel_for(0, np, 1, [&](int pidx) {
      int i, j, k;
      func(pidx, i, j, k);
      array_idx[pidx] = hash(i, j, k, pidx);
    });

    tbb::parallel_sort(array_idx.begin(), array_idx.end());

    tbb::parallel_for(0, np, 1, [&](int pidx) {
      int G_ID = pidx;
      int G_ID_PREV = G_ID - 1;
      int G_ID_NEXT = G_ID + 1;

      unsigned int cell = (unsigned int)(array_idx[G_ID] >> 32UL);
      unsigned int cell_prev =
          G_ID_PREV < 0 ? -1U : (unsigned int)(array_idx[G_ID_PREV] >> 32UL);
      unsigned int cell_next =
          G_ID_NEXT >= np ? -1U : (unsigned int)(array_idx[G_ID_NEXT] >> 32UL);
      if (cell != cell_prev) {
        // I'm the start of a cell
        array_sup[cell].first = G_ID;
      }
      if (cell != cell_next) {
        // I'm the end of a cell
        array_sup[cell].second = G_ID + 1;
      }
    });
  }

  template <typename Callable>
  void getNeigboringParticles_cell(int i, int j, int k, int wl, int wh, int hl,
                                   int hh, int dl, int dh, Callable func) {
    for (int sk = k + dl; sk <= k + dh; ++sk)
      for (int si = i + wl; si <= i + wh; si++)
        for (int sj = j + hl; sj <= j + hh; sj++) {
          getCellAt(si, sj, sk, func);
        }
  }

  template <typename Callable>
  void getCellAt(int i, int j, int k, Callable func) {
    if (i < 0 || i > ni - 1 || j < 0 || j > nj - 1 || k < 0 || k > nk - 1)
      return;
    unsigned int G_CELL = (unsigned int)(k * (ni * nj) + j * ni + i);
    const std::pair<int, int>& G_START_END = array_sup[G_CELL];
    for (int N_ID = G_START_END.first; N_ID < G_START_END.second; ++N_ID) {
      func((int)(array_idx[N_ID] & 0xFFFFFFFFUL));
    }
  }

  std::vector<uint64_t> array_idx;
  std::vector<std::pair<int, int> > array_sup;

  int ni;
  int nj;
  int nk;
};

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_SORTER_H_

//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2015 Xinxin Zhang
//

#include "GeometricLevelGen.h"

//#define AMG_VERBOSE

template <class T>
void levelGen<T>::generateRP(const robertbridson::FixedSparseMatrix<T> &A,
                             robertbridson::FixedSparseMatrix<T> &R,
                             robertbridson::FixedSparseMatrix<T> &P,
                             int ni, int nj, int nk) {
  int nni = ceil((float)ni / 2.0);
  int nnj = ceil((float)nj / 2.0);
  int nnk = ceil((float)nk / 2.0);
  robertbridson::SparseMatrix<T> r;
  robertbridson::SparseMatrix<T> p;
  p.resize(ni * nj * nk);
  p.zero();
  r.resize(nni * nnj * nnk);
  r.zero();

  for (int k = 0; k < nnk; k++)
    for (int j = 0; j < nnj; j++)
      for (int i = 0; i < nni; i++) {
        unsigned int index = (k * nnj + j) * nni + i;
        for (int kk = 0; kk <= 1; kk++)
          for (int jj = 0; jj <= 1; jj++)
            for (int ii = 0; ii <= 1; ii++) {
              int iii = i * 2 + ii;
              int jjj = j * 2 + jj;
              int kkk = k * 2 + kk;
              if (iii < ni && jjj < nj && kkk < nk) {
                unsigned int index2 = (kkk * nj + jjj) * ni + iii;
                r.set_element(index, index2, (T)0.125);
                p.set_element(index2, index, 1.0);
              }
            }
      }

  R.construct_from_matrix(r);
  P.construct_from_matrix(p);
  r.clear();
  p.clear();

  // transposeMat(R,P,(T)8.0);
}

template struct levelGen<double>;

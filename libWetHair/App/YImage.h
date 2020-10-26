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

#ifndef LIBWETHAIR_APP_YIMAGE_H_
#define LIBWETHAIR_APP_YIMAGE_H_

// file loading/saving automatically picks up changes to this struct.
// The possibilities are: ARGB, ABGR, RGBA, BGRA.

struct YPixel {
  unsigned char r;
  unsigned char g;
  unsigned char b;
  unsigned char a;
};

class YImage {
 public:
  YImage();

  YImage(const YImage&);

  virtual ~YImage();

  YImage& operator=(const YImage&);

  bool save(const char* fname) const;

  bool load(const char* fname);

  YPixel* data();

  const YPixel* data() const;

  YPixel& at(int i, int j);

  const YPixel& at(int i, int j) const;

  int width() const;

  int height() const;

  void resize(int width, int height);

  // flip vertically
  void flip();

  // flip horizontally
  void mirror();

  // average rgb
  void greyscale();

 protected:
  int m_width;
  int m_height;
  YPixel* m_data;  // raw image data
};

#endif  // LIBWETHAIR_APP_YIMAGE_H_

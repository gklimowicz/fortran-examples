// Copyright (c) 2012-2013, 2015 Harald Klimach <harald@klimachs.de>
//
// Parts of this file were written by Harald Klimach for
// German Research School of Simulation Sciences and University of
// Siegen.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// OR OTHER DEALINGS IN THE SOFTWARE.
// *************************************************************************** !

#include <stdlib.h>
#include "lua.h"

typedef struct
{
  int length;
  int space;
  char *container;
} charbuf;

// Writer to use during lua_dump.
static int buf_writer(lua_State *L, const void* p, size_t sz, void* ud)
{
  charbuf *dat;
  const char *buf;
  int i;

  dat = ud;
  buf = p;

  if ( sz + dat->length > dat->space ) {
    // Increase the size of the buffer, if needed.
    dat->space = ((dat->space*2) > (sz + dat->length))
               ? (dat->space*2) : (sz + dat->length);
    dat->container = realloc(dat->container, dat->space);
    if (!dat->container) return -10;
  }

  // Append the data to write into the buffer.
  for (i=0; i<sz; i++) {
    dat->container[dat->length + i] = buf[i];
  }
  dat->length = dat->length + sz;
  return 0;
}


// Wrapper around lua_dump to write into a memory buffer.
// Return Fortran friendly arguments.
const char* dump_lua_toBuf(lua_State *L, int *length, int *ierr)
{
  charbuf dat;
  char *buf;
  int i;
  int errcode;
  size_t sz;

  dat.length = 0;
  dat.space = 1024;
  dat.container = malloc(dat.space);

  errcode = lua_dump(L, buf_writer, &dat, 0);

  (*ierr) = errcode;
  (*length) = dat.length;
  sz = dat.length;
  buf = malloc(dat.length);
  for (i=0; i<dat.length; i++) {
    buf[i] = dat.container[i];
  }
  free(dat.container);
  return buf;
}

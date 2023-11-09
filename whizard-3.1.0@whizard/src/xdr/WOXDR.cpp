//// WOXDR.cpp
//
// Simple XDR class, see header
//
// WGL, 24 May 2002
//
////
#include "WOXDR.hh"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__APPLE_CC__)
#include "sys/types.h"
#endif

#if defined(__linux) || defined(__CYGWIN__) || defined(__APPLE_CC__)
#include <netinet/in.h>
#endif
#ifdef _MSC_VER
#include <winsock.h>
#else
#include <sys/socket.h>
#endif

namespace WOUTIL{
////
//
// Constructor, destructor
//
////
WOXDR::~WOXDR()
{
   if (_fp) {
      fclose(_fp);
      _fp = 0;
   }
   if (_fileName) {
      delete [] _fileName;
      _fileName = 0;
   }
   return;
}

WOXDR::WOXDR(const char *filename, bool open_for_write) : _fileName(0), _fp(0) 
{
   setFileName(filename, open_for_write);
   if (htonl(1L) == 1L) _hasNetworkOrder = true;
   else                 _hasNetworkOrder = false;
   return;
}

void WOXDR::setFileName(const char *filename, bool open_for_write)
{
//
// First check if we can open this file
//
   if (filename == 0) {
      _error = WOXDR_OPENFAILURE;
      return;
   }
#ifdef _MSC_VER
   FILE *fp = fopen(filename, open_for_write ? "wb" : "rb");
#else
   FILE *fp = fopen(filename, open_for_write ? "w" : "r");
#endif
   if (fp == 0) {
      _error = WOXDR_OPENFAILURE;
      return;
   }

   if (_fp) fclose(_fp);
   _fp = fp;

   if (_fileName) {
      delete [] _fileName;
      _fileName = 0;
   }

   int n = strlen(filename);
   _fileName = new char [n + 1];
   strncpy(_fileName, filename, n);
   _fileName[n] = '\0';

   _openForWrite = open_for_write;

   _error = WOXDR_SUCCESS;
   return;
}

double WOXDR::ntohd(double d) const
{
//
// If we already have network order, we don't swap
//
   if (_hasNetworkOrder == false) {
      union {
         double        d;
         unsigned char b[8];
      } dd;
      int i;

      dd.d = d;
      for (i = 0; i < 4; i++) {
         unsigned char c = dd.b[i];
         dd.b[i]         = dd.b[7 - i];
         dd.b[7 - i]     = c;
      }
      d = dd.d;
   }
   return(d);
}

long WOXDR::checkRead(long *l)
{
   if (_openForWrite) return(_error = WOXDR_READONLY);
   if (_fp == 0)      return(_error = WOXDR_NOFILE);
   if (l) {
      // je: in architectures where long isn't 4 byte long this code crashes
      //long nr;
      //if ((nr = fread(l, 4, 1, _fp)) != 1) return(_error = WOXDR_READERROR);
      //*l = ntohl(*l);

      int32_t buf;
      if (fread(&buf, 4, 1, _fp) != 1) return(_error = WOXDR_READERROR);
      *l = ((int32_t)ntohl(buf));
   }
   return(WOXDR_SUCCESS);
}

long WOXDR::checkRead(double *d)
{
   if (_openForWrite) return(_error = WOXDR_READONLY);
   if (_fp == 0)      return(_error = WOXDR_NOFILE);
   if (d) {
      if (fread(d, 8, 1, _fp) != 1) return(_error = WOXDR_READERROR);
      *d = ntohd(*d);
   }
   return(WOXDR_SUCCESS);
}

long WOXDR::checkRead(float *f)
{
   if (_openForWrite) return(_error = WOXDR_READONLY);
   if (_fp == 0)      return(_error = WOXDR_NOFILE);
   if (f) {
      if (fread(f, 4, 1, _fp) != 1) return(_error = WOXDR_READERROR);
      // je: in architectures where long isn't 4 byte long this code crashes
      //*((long *) f) = ntohl(*((long *) f));

      *((int32_t *) f) = ntohl(*((int32_t *) f));
   }
   return(WOXDR_SUCCESS);
}

long WOXDR::readLong(void)
{
   long l = 0;
   checkRead(&l);
   return(l);
}

double WOXDR::readDouble(void)
{
   double d = 0.0;
   checkRead(&d);
   return(d);
}

double WOXDR::readFloat(void)
{
   float f = 0.0;
   checkRead(&f);
   return((double) f);
}

const char *WOXDR::readString(long &length)
{
   if (checkRead(&length)) return(0);
   long rl = (length + 3) & 0xFFFFFFFC;
   char *s = new char[rl + 1];
   if (fread(s, 1, rl, _fp) != (unsigned long) rl) {
      _error = WOXDR_READERROR;
      delete [] s;
      return(0);
   }
   s[rl] = '\0';
   _error = WOXDR_SUCCESS;
   return(s);
}

long *WOXDR::readLongArray(long &length)
{
   if (checkRead(&length)) return(0);
   long *s = new long[length];
   // je: in architectures where long isn't 4 byte long this code crashes
   //if (fread(s, 4, length, _fp) != (unsigned long) length) {
   //   _error = WOXDR_READERROR;
   //   delete [] s;
   //   return(0);
   //}
   //if (_hasNetworkOrder == false) for (long i = 0; i < length; i++) s[i] = ntohl(s[i]);

   int32_t *buf = new int32_t[length];
   if (fread(buf, 4, length, _fp) != (unsigned long) length) {
       _error = WOXDR_READERROR;
       delete [] buf;
       delete [] s;
       return(0);
   }
   for (long i = 0; i < length; i++){
       if (_hasNetworkOrder == false){
           s[i] = ((int32_t)ntohl(buf[i]));
       }
       else{
            s[i] = (long)buf[i];
       }
   }
   delete [] buf;
   _error = WOXDR_SUCCESS;
   return(s);
}

double *WOXDR::readDoubleArray(long &length)
{
   if (checkRead(&length)) return(0);
   double *s = new double[length];
   if (fread(s, 8, length, _fp) != (unsigned long) length) {
      _error = WOXDR_READERROR;
      delete [] s;
      return(0);
   }
   if (_hasNetworkOrder == false) for (long i = 0; i < length; i++) s[i] = ntohd(s[i]);
   _error = WOXDR_SUCCESS;
   return(s);
}

double *WOXDR::readFloatArray(long &length)
{
   if (checkRead(&length)) return(0);
   long *st = new long[length];
   // je: FIXME this will cause problems in architectures where long isn't 4 byte long
   if (fread(st, 4, length, _fp) != (unsigned long) length) {
      _error = WOXDR_READERROR;
      delete [] st;
      return(0);
   }
   double *s = new double[length];
   // je: FIXME what happens if _hasNetworkOrder == true?!
   if (_hasNetworkOrder == false) {
      for (long i = 0; i < length; i++) {
         long l = ntohl(st[i]);
         s[i] = (double) (*((float *) &l));
      }
   }
   _error = WOXDR_SUCCESS;
   delete [] st;
   return(s);
}

long WOXDR::checkWrite(long *l)
{
  if (_openForWrite == false) return(_error = WOXDR_WRITEONLY);
  if (_fp == 0)               return(_error = WOXDR_NOFILE);
  if (l) {
    long ll = htonl(*l);
    // je: FIXME this will cause problems in architectures where long isn't 4 byte long
    if (fwrite(&ll, 4, 1, _fp) != 4) return(_error = WOXDR_WRITEERROR);
   }
   return(WOXDR_SUCCESS);
}

long WOXDR::checkWrite(double *d)
{
  if (_openForWrite == false) return(_error = WOXDR_WRITEONLY);
  if (_fp == 0)               return(_error = WOXDR_NOFILE);
  if (d) {
    double dd = htond(*d);
    if (fwrite(&dd, 8, 1, _fp) != 8) return(_error = WOXDR_WRITEERROR);
   }
   return(WOXDR_SUCCESS);
}

long WOXDR::writeLong(long data)
{
   return(checkWrite(&data));
}

long WOXDR::writeDouble(double data)
{
   return(checkWrite(&data));
}

long WOXDR::writeString(const char *data)
{
   return(writeString(data, strlen(data)));
}

long WOXDR::writeString(const char *data, long length)
{
   if (checkWrite(&length)) return(_error);
   if (fwrite(data, 1, length, _fp) != (unsigned long) length) return(_error = WOXDR_WRITEERROR);
   long l = ((length + 3) & 0xFFFFFFFC) - length;
   if (fwrite(&l, 1, l, _fp) != (unsigned long) l) return(_error = WOXDR_WRITEERROR);
   return(_error = WOXDR_SUCCESS);
}

long WOXDR::writeLongArray(const long *data, long length)
{
   if (checkWrite(&length)) return(_error);
   long *s = (long *) data;
   if (_hasNetworkOrder == false) {
      s = new long[length];
      for (long i = 0; i < length; i++) s[i] = htonl(data[i]);
   }
   // je: FIXME this will cause problems in architectures where long isn't 4 byte long
   long l = fwrite(s, 4, length, _fp);
   if (_hasNetworkOrder == false) delete [] s;
   if (l != length) return(_error = WOXDR_WRITEERROR);
   return(_error = WOXDR_SUCCESS);
}

long WOXDR::writeDoubleArray(const double *data, long length)
{
   if (checkWrite(&length)) return(_error);
   double *s = (double *) data;
   if (_hasNetworkOrder == false) {
      s = new double[length];
      for (long i = 0; i < length; i++) s[i] = htond(data[i]);
   }
   long l = fwrite(s, 8, length, _fp);
   if (_hasNetworkOrder == false) delete [] s;
   if (l != length) return(_error = WOXDR_WRITEERROR);
   return(_error = WOXDR_SUCCESS);
}


long WOXDR::filePosition(long pos)
{
   if (_fp == 0) {
      _error = WOXDR_NOFILE;
      return(-1);
   }
   if (pos == -1) return(ftell(_fp));
   if (fseek(_fp, pos, SEEK_SET)) {
      _error = WOXDR_SEEKERROR;
      return(-1);
   }
   return(pos);
}

}// end namespace

/*
! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'
*/
#include <signal.h>
#include <stdlib.h>

extern int wo_sigint;
extern int wo_sigterm;
extern int wo_sigxcpu;
extern int wo_sigxfsz;

static void wo_handler_sigint (int sig) {
  wo_sigint = sig;
}

static void wo_handler_sigterm (int sig) {
  wo_sigterm = sig;
}

static void wo_handler_sigxcpu (int sig) {
  wo_sigxcpu = sig;
}

static void wo_handler_sigxfsz (int sig) {
  wo_sigxfsz = sig;
}

int wo_mask_sigint () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = wo_handler_sigint;
  return sigaction(SIGINT, &sa, NULL);
}

int wo_mask_sigterm () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = wo_handler_sigterm;
  return sigaction(SIGTERM, &sa, NULL);
}

int wo_mask_sigxcpu () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = wo_handler_sigxcpu;
  return sigaction(SIGXCPU, &sa, NULL);
}

int wo_mask_sigxfsz () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = wo_handler_sigxfsz;
  return sigaction(SIGXFSZ, &sa, NULL);
}

int wo_release_sigint () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = SIG_DFL;
  return sigaction(SIGINT, &sa, NULL);
}

int wo_release_sigterm () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = SIG_DFL;
  return sigaction(SIGTERM, &sa, NULL);
}

int wo_release_sigxcpu () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = SIG_DFL;
  return sigaction(SIGXCPU, &sa, NULL);
}

int wo_release_sigxfsz () {
  struct sigaction sa;
  sigset_t blocks;
  sigfillset (&blocks);
  sa.sa_flags = 0;
  sa.sa_mask = blocks;
  sa.sa_handler = SIG_DFL;
  return sigaction(SIGXFSZ, &sa, NULL);
}

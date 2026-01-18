
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin
! and Lars Nordstrom. This file is distributed under the terms of the GNU
! General Public License. See the file COPYING for license details.

subroutine writew90
use modmain
use modw90
implicit none
! initialise universal and Wannier90 variables
call initw90
! write the .win file
call writew90win
! write the .eig file
call writew90eig
! call the Wannier90 setup routine
call setupw90
! write the .mmn file
call writew90mmn
! write the .spn file
call writew90spn
end subroutine


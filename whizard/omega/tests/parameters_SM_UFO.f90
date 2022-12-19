module parameters_sm_ufo
  use kinds
  use constants
  implicit none
  private
  public :: setup_parameters, import_from_whizard, model_update_alpha_s, &
    print_parameters
  real(kind=default), public, save :: aEWM1 = 127.9_default, &
    Gf = 1.16637e-05_default, aS = 0.1184_default, ymdo = 0.00504_default, &
    ymup = 0.00255_default, yms = 0.101_default, ymc = 1.27_default, &
    ymb = 4.7_default, ymt = 172._default, yme = 0.000511_default, &
    ymm = 0.10566_default, ymtau = 1.777_default, MZ = 91.1876_default, &
    Me = 0.000511_default, MMU = 0.10566_default, MTA = 1.777_default, &
    MU = 0.00255_default, MC = 1.27_default, MT = 172._default, &
    MD = 0.00504_default, MS = 0.101_default, MB = 4.7_default, &
    MH = 125._default, WZ = 2.4952_default, WW = 2.085_default, &
    WT = 1.50833649_default, WH = 0.00407_default
  real(kind=default), public, save :: muH, yup, ytau, yt, ys, ym, ye, ydo, &
    yc, yb, lam, vev, gw, g1, sw, cw, sw2, ee, MW, G, aEW
  complex(kind=default), public, save :: GC_1, GC_10, GC_11, GC_12, GC_13, &
    GC_14, GC_15, GC_16, GC_17, GC_18, GC_19, GC_2, GC_20, GC_21, GC_22, &
    GC_23, GC_24, GC_25, GC_26, GC_27, GC_28, GC_29, GC_3, GC_30, GC_31, &
    GC_32, GC_33, GC_34, GC_35, GC_36, GC_37, GC_38, GC_39, GC_4, GC_40, &
    GC_41, GC_42, GC_43, GC_44, GC_45, GC_46, GC_47, GC_48, GC_49, GC_5, &
    GC_50, GC_51, GC_52, GC_53, GC_54, GC_55, GC_56, GC_57, GC_58, GC_59, &
    GC_6, GC_60, GC_61, GC_62, GC_63, GC_64, GC_65, GC_66, GC_67, GC_68, &
    GC_69, GC_7, GC_70, GC_71
  complex(kind=default), public, save :: GC_72, GC_73, GC_74, GC_75, GC_76, &
    GC_77, GC_78, GC_79, GC_8, GC_80, GC_81, GC_82, GC_83, GC_84, GC_85, &
    GC_86, GC_87, GC_88, GC_89, GC_9, GC_90, GC_91, I4a33, I4a22, I4a11, &
    I3a33, I3a22, I3a11, I2a33, I2a22, I2a11, I1a33, I1a22, I1a11
  interface cconjg
    module procedure cconjg_real, cconjg_complex
  end interface
  private :: cconjg_real, cconjg_complex
contains
  function cconjg_real (x) result (xc)
    real(kind=default), intent(in) :: x
    real(kind=default) :: xc
    xc = x
  end function cconjg_real
  function cconjg_complex (z) result (zc)
    complex(kind=default), intent(in) :: z
    complex(kind=default) :: zc
    zc = conjg (z)
  end function cconjg_complex
  ! derived parameters:
  subroutine setup_parameters_001 ()
    aEW = (1.0_default / aEWM1)
    G = (sqrt (pi) * 2.0_default * sqrt ((aS / 2.0_default)))
    MW = sqrt ((((MZ**2) / 2.0_default) + sqrt ((((MZ**4) / 4.0_default) - (( &
      (MZ**2) * aEW * pi) / (Gf * sqrt (2.0_default)))))))
    ee = (sqrt (pi) * 2.0_default * sqrt (aEW))
    sw2 = (1.0_default - ((MW**2) / (MZ**2)))
    cw = sqrt ((1.0_default - sw2))
    sw = sqrt (sw2)
    g1 = (ee / cw)
    gw = (ee / sw)
    vev = ((sw * 2.0_default * MW) / ee)
    lam = ((MH**2) / (2.0_default * (vev**2)))
    yb = ((ymb * sqrt (2.0_default)) / vev)
    yc = ((ymc * sqrt (2.0_default)) / vev)
    ydo = ((ymdo * sqrt (2.0_default)) / vev)
    ye = ((yme * sqrt (2.0_default)) / vev)
    ym = ((ymm * sqrt (2.0_default)) / vev)
    ys = ((yms * sqrt (2.0_default)) / vev)
    yt = ((ymt * sqrt (2.0_default)) / vev)
    ytau = ((ymtau * sqrt (2.0_default)) / vev)
    yup = ((ymup * sqrt (2.0_default)) / vev)
    muH = sqrt ((lam * (vev**2)))
    I1a11 = ydo
    I1a22 = ys
    I1a33 = yb
    I2a11 = yup
    I2a22 = yc
    I2a33 = yt
    I3a11 = yup
    I3a22 = yc
    I3a33 = yt
    I4a11 = ydo
    I4a22 = ys
    I4a33 = yb
    GC_91 = ((yup / sqrt (2.0_default)) / (0,1))
    GC_90 = (((-1.0_default) * (((0,1) * yup) / sqrt (2.0_default))) / (0,1))
    GC_9 = (((ee**2) / (2.0_default * cw)) / (0,1))
    GC_89 = (((-1.0_default) * (((0,1) * ytau) /  &
      sqrt (2.0_default))) / (0,1))
    GC_88 = (((-1.0_default) * (ytau / sqrt (2.0_default))) / (0,1))
    GC_87 = (ytau / (0,1))
    GC_86 = (((-1.0_default) * ytau) / (0,1))
    GC_85 = ((yt / sqrt (2.0_default)) / (0,1))
    GC_84 = (((-1.0_default) * (((0,1) * yt) / sqrt (2.0_default))) / (0,1))
    GC_83 = (((-1.0_default) * (((0,1) * ys) / sqrt (2.0_default))) / (0,1))
    GC_82 = (((-1.0_default) * (ys / sqrt (2.0_default))) / (0,1))
    GC_81 = (((-1.0_default) * (((0,1) * ym) / sqrt (2.0_default))) / (0,1))
    GC_80 = (((-1.0_default) * (ym / sqrt (2.0_default))) / (0,1))
    GC_8 = ((((ee**2) * (0,1)) / (2.0_default * cw)) / (0,1))
    GC_79 = (ym / (0,1))
    GC_78 = (((-1.0_default) * ym) / (0,1))
    GC_77 = (((-1.0_default) * (((0,1) * ye) / sqrt (2.0_default))) / (0,1))
    GC_76 = (((-1.0_default) * (ye / sqrt (2.0_default))) / (0,1))
    GC_75 = (ye / (0,1))
    GC_74 = (((-1.0_default) * ye) / (0,1))
    GC_73 = (((-1.0_default) * (((0,1) * ydo) / sqrt (2.0_default))) / (0,1))
    GC_72 = (((-1.0_default) * (ydo / sqrt (2.0_default))) / (0,1))
    GC_71 = ((yc / sqrt (2.0_default)) / (0,1))
    GC_70 = (((-1.0_default) * (((0,1) * yc) / sqrt (2.0_default))) / (0,1))
    GC_7 = (((((-1.0_default) * ee)**2) / (2.0_default * cw)) / (0,1))
    GC_69 = (((-1.0_default) * (((0,1) * yb) / sqrt (2.0_default))) / (0,1))
    GC_68 = (((-1.0_default) * (yb / sqrt (2.0_default))) / (0,1))
    GC_67 = ((((vev * (sw**2) * (ee**2) * (0,1)) / (2.0_default * (cw**2))) &
       + (vev * (ee**2) * (0,1)) + ((vev * (0,1) * (cw**2) * (ee**2)) /  &
      (2.0_default * (sw**2)))) / (0,1))
    GC_66 = ((((((-1.0_default) * vev * (ee**2) * (0,1)) / 2.0_default) - ( &
      (vev * (0,1) * (cw**2) * (ee**2)) / (4.0_default * (sw**2)))) - ( &
      (vev * (sw**2) * (ee**2) * (0,1)) / (4.0_default * (cw**2)))) / (0,1))
    GC_65 = (((((ee**2) * vev) / (4.0_default * cw)) + ((vev * cw *  &
      (ee**2)) / (4.0_default * (sw**2)))) / (0,1))
    GC_64 = (((((-1.0_default) * (ee**2) * vev) / (4.0_default * cw)) + ( &
      (vev * cw * (ee**2)) / (4.0_default * (sw**2)))) / (0,1))
    GC_63 = (((((ee**2) * vev) / (4.0_default * cw)) - ((vev * cw *  &
      (ee**2)) / (4.0_default * (sw**2)))) / (0,1))
    GC_62 = (((((-1.0_default) * (ee**2) * vev) / (4.0_default * cw)) - ( &
      (vev * cw * (ee**2)) / (4.0_default * (sw**2)))) / (0,1))
    GC_61 = ((((ee**2) * vev) / (2.0_default * sw)) / (0,1))
    GC_60 = ((((-1.0_default) * (ee**2) * vev) / (2.0_default * sw)) / (0,1))
    GC_6 = (((0,1) * 2.0_default * (ee**2)) / (0,1))
    GC_59 = ((((ee**2) * vev) / (4.0_default * (sw**2))) / (0,1))
    GC_58 = (((vev * (ee**2) * (0,1)) / (2.0_default * (sw**2))) / (0,1))
    GC_57 = ((((-1.0_default) * vev * (ee**2) * (0,1)) / (4.0_default *  &
      (sw**2))) / (0,1))
    GC_56 = ((((-1.0_default) * (ee**2) * vev) / (4.0_default *  &
      (sw**2))) / (0,1))
    GC_55 = ((vev * lam * (-6.0_default) * (0,1)) / (0,1))
    GC_54 = ((vev * lam * (-2.0_default) * (0,1)) / (0,1))
    GC_53 = ((((ee**2) * vev) / (2.0_default * cw)) / (0,1))
    GC_52 = ((((-1.0_default) * (ee**2) * vev) / (2.0_default * cw)) / (0,1))
    GC_51 = (((((sw**2) * (ee**2) * (0,1)) / (2.0_default * (cw**2))) + ( &
      (ee**2) * (0,1)) + (((0,1) * (cw**2) * (ee**2)) / (2.0_default *  &
      (sw**2)))) / (0,1))
    GC_50 = (((((sw**2) * (ee**2) * (0,1)) / (2.0_default * (cw**2))) +  &
      ((-1.0_default) * (ee**2) * (0,1)) + (((0,1) * (cw**2) * (ee**2)) /  &
      (2.0_default * (sw**2)))) / (0,1))
    GC_5 = (((ee**2) * (0,1)) / (0,1))
    GC_49 = (((((0,1) * cw * (ee**2)) / sw) - ((sw *  &
      (ee**2) * (0,1)) / cw)) / (0,1))
    GC_48 = (((((0,1) * cw * ee) / (2.0_default * sw)) + ( &
      (sw * ee * (0,1)) / (2.0_default * cw))) / (0,1))
    GC_47 = (((((-1.0_default) * (0,1) * cw * ee) / (2.0_default * sw)) + ( &
      (sw * ee * (0,1)) / (2.0_default * cw))) / (0,1))
    GC_46 = (((((0,1) * cw * ee) / (2.0_default * sw)) - ( &
      (sw * ee * (0,1)) / (6.0_default * cw))) / (0,1))
    GC_45 = (((((-1.0_default) * (0,1) * cw * ee) / (2.0_default * sw)) - ( &
      (sw * ee * (0,1)) / (6.0_default * cw))) / (0,1))
    GC_44 = (((((-1.0_default) * cw * ee) / (2.0_default * sw)) - ( &
      (ee * sw) / (2.0_default * cw))) / (0,1))
    GC_43 = (((sw * ee * (0,1)) / cw) / (0,1))
    GC_42 = (((sw * (0,1) * (-2.0_default) * ee) /  &
      (3.0_default * cw)) / (0,1))
    GC_41 = (((sw * ee * (0,1)) / (3.0_default * cw)) / (0,1))
    GC_40 = ((((0,1) * (ee**2) * (-2.0_default) * cw) / sw) / (0,1))
    GC_4 = ((ee * (0,1)) / (0,1))
    GC_39 = (((ee**2) / (2.0_default * sw)) / (0,1))
    GC_38 = ((((-1.0_default) * (ee**2) * (0,1)) /  &
      (2.0_default * sw)) / (0,1))
    GC_37 = (((((-1.0_default) * ee)**2) / (2.0_default * sw)) / (0,1))
    GC_36 = ((((0,1) * cw * ee) / sw) / (0,1))
    GC_35 = (((-1.0_default) * (((0,1) * cw * ee) / sw)) / (0,1))
    GC_34 = (((ee * (0,1)) / (sw * sqrt (2.0_default))) / (0,1))
    GC_33 = (((ee * (0,1)) / (2.0_default * sw)) / (0,1))
    GC_32 = ((((-1.0_default) * ee * (0,1)) / (2.0_default * sw)) / (0,1))
    GC_31 = ((((-1.0_default) * ee) / (2.0_default * sw)) / (0,1))
    GC_30 = ((((0,1) * (cw**2) * (ee**2)) / (sw**2)) / (0,1))
    GC_3 = (((-1.0_default) * ee * (0,1)) / (0,1))
    GC_29 = (((-1.0_default) * (((ee**2) * (0,1)) / (sw**2))) / (0,1))
    GC_28 = ((((ee**2) * (0,1)) / (2.0_default * (sw**2))) / (0,1))
    GC_27 = ((lam * (-6.0_default) * (0,1)) / (0,1))
    GC_26 = ((lam * (-4.0_default) * (0,1)) / (0,1))
    GC_25 = ((lam * (-2.0_default) * (0,1)) / (0,1))
    GC_24 = (((-1.0_default) * I4a33) / (0,1))
    GC_23 = (((-1.0_default) * I4a22) / (0,1))
    GC_22 = (((-1.0_default) * I4a11) / (0,1))
    GC_21 = (I3a33 / (0,1))
    GC_20 = (I3a22 / (0,1))
    GC_2 = ((((0,1) * 2.0_default * ee) / 3.0_default) / (0,1))
    GC_19 = (I3a11 / (0,1))
    GC_18 = (((-1.0_default) * I2a33) / (0,1))
    GC_17 = (((-1.0_default) * I2a22) / (0,1))
    GC_16 = (((-1.0_default) * I2a11) / (0,1))
    GC_15 = (I1a33 / (0,1))
    GC_14 = (I1a22 / (0,1))
    GC_13 = (I1a11 / (0,1))
  end subroutine setup_parameters_001
  subroutine setup_parameters_002 ()
    GC_12 = (((0,1) * (G**2)) / (0,1))
    GC_11 = (((0,1) * G) / (0,1))
    GC_10 = (((-1.0_default) * G) / (0,1))
    GC_1 = ((((-1.0_default) * ee * (0,1)) / 3.0_default) / (0,1))
  end subroutine setup_parameters_002
  subroutine setup_parameters ()
    call setup_parameters_001 ()
    call setup_parameters_002 ()
  end subroutine setup_parameters
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(27), intent(in) :: par_array
    integer, intent(in) :: scheme
    aEWM1 = par_array(1)
    Gf = par_array(2)
    aS = par_array(3)
    ymdo = par_array(4)
    ymup = par_array(5)
    yms = par_array(6)
    ymc = par_array(7)
    ymb = par_array(8)
    ymt = par_array(9)
    yme = par_array(10)
    ymm = par_array(11)
    ymtau = par_array(12)
    MZ = par_array(13)
    Me = par_array(14)
    MMU = par_array(15)
    MTA = par_array(16)
    MU = par_array(17)
    MC = par_array(18)
    MT = par_array(19)
    MD = par_array(20)
    MS = par_array(21)
    MB = par_array(22)
    MH = par_array(23)
    WZ = par_array(24)
    WW = par_array(25)
    WT = par_array(26)
    WH = par_array(27)
    call setup_parameters ()
  end subroutine import_from_whizard
  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    aS = alpha_s
    G = (sqrt (pi) * 2.0_default * sqrt ((aS / 2.0_default)))
    GC_12 = (((0,1) * (G**2)) / (0,1))
    GC_11 = (((0,1) * G) / (0,1))
    GC_10 = (((-1.0_default) * G) / (0,1))
  end subroutine model_update_alpha_s
  subroutine print_parameters ()
    character(len=*), parameter :: fmt_real = "(A12,4X,' = ',E25.18)", &
      fmt_complex = "(A12,4X,' = ',E25.18,' + i*',E25.18)", &
      fmt_real_array = "(A12,'(',I2.2,')',' = ',E25.18)", &
      fmt_complex_array = "(A12,'(',I2.2,')',' = ',E25.18,' + i*',E25.18)"
    write (unit = *, fmt = "(A)") "default values for the input parameters:"
    write (unit = *, fmt = fmt_real) "aEWM1", aEWM1
    write (unit = *, fmt = fmt_real) "Gf", Gf
    write (unit = *, fmt = fmt_real) "aS", aS
    write (unit = *, fmt = fmt_real) "ymdo", ymdo
    write (unit = *, fmt = fmt_real) "ymup", ymup
    write (unit = *, fmt = fmt_real) "yms", yms
    write (unit = *, fmt = fmt_real) "ymc", ymc
    write (unit = *, fmt = fmt_real) "ymb", ymb
    write (unit = *, fmt = fmt_real) "ymt", ymt
    write (unit = *, fmt = fmt_real) "yme", yme
    write (unit = *, fmt = fmt_real) "ymm", ymm
    write (unit = *, fmt = fmt_real) "ymtau", ymtau
    write (unit = *, fmt = fmt_real) "MZ", MZ
    write (unit = *, fmt = fmt_real) "Me", Me
    write (unit = *, fmt = fmt_real) "MMU", MMU
    write (unit = *, fmt = fmt_real) "MTA", MTA
    write (unit = *, fmt = fmt_real) "MU", MU
    write (unit = *, fmt = fmt_real) "MC", MC
    write (unit = *, fmt = fmt_real) "MT", MT
    write (unit = *, fmt = fmt_real) "MD", MD
    write (unit = *, fmt = fmt_real) "MS", MS
    write (unit = *, fmt = fmt_real) "MB", MB
    write (unit = *, fmt = fmt_real) "MH", MH
    write (unit = *, fmt = fmt_real) "WZ", WZ
    write (unit = *, fmt = fmt_real) "WW", WW
    write (unit = *, fmt = fmt_real) "WT", WT
    write (unit = *, fmt = fmt_real) "WH", WH
    write (unit = *, fmt = "(A)") "derived parameters:"
    write (unit = *, fmt = fmt_real) "muH", muH
    write (unit = *, fmt = fmt_real) "yup", yup
    write (unit = *, fmt = fmt_real) "ytau", ytau
    write (unit = *, fmt = fmt_real) "yt", yt
    write (unit = *, fmt = fmt_real) "ys", ys
    write (unit = *, fmt = fmt_real) "ym", ym
    write (unit = *, fmt = fmt_real) "ye", ye
    write (unit = *, fmt = fmt_real) "ydo", ydo
    write (unit = *, fmt = fmt_real) "yc", yc
    write (unit = *, fmt = fmt_real) "yb", yb
    write (unit = *, fmt = fmt_real) "lam", lam
    write (unit = *, fmt = fmt_real) "vev", vev
    write (unit = *, fmt = fmt_real) "gw", gw
    write (unit = *, fmt = fmt_real) "g1", g1
    write (unit = *, fmt = fmt_real) "sw", sw
    write (unit = *, fmt = fmt_real) "cw", cw
    write (unit = *, fmt = fmt_real) "sw2", sw2
    write (unit = *, fmt = fmt_real) "ee", ee
    write (unit = *, fmt = fmt_real) "MW", MW
    write (unit = *, fmt = fmt_real) "G", G
    write (unit = *, fmt = fmt_real) "aEW", aEW
    write (unit = *, fmt = fmt_complex) "GC_1", GC_1
    write (unit = *, fmt = fmt_complex) "GC_10", GC_10
    write (unit = *, fmt = fmt_complex) "GC_11", GC_11
    write (unit = *, fmt = fmt_complex) "GC_12", GC_12
    write (unit = *, fmt = fmt_complex) "GC_13", GC_13
    write (unit = *, fmt = fmt_complex) "GC_14", GC_14
    write (unit = *, fmt = fmt_complex) "GC_15", GC_15
    write (unit = *, fmt = fmt_complex) "GC_16", GC_16
    write (unit = *, fmt = fmt_complex) "GC_17", GC_17
    write (unit = *, fmt = fmt_complex) "GC_18", GC_18
    write (unit = *, fmt = fmt_complex) "GC_19", GC_19
    write (unit = *, fmt = fmt_complex) "GC_2", GC_2
    write (unit = *, fmt = fmt_complex) "GC_20", GC_20
    write (unit = *, fmt = fmt_complex) "GC_21", GC_21
    write (unit = *, fmt = fmt_complex) "GC_22", GC_22
    write (unit = *, fmt = fmt_complex) "GC_23", GC_23
    write (unit = *, fmt = fmt_complex) "GC_24", GC_24
    write (unit = *, fmt = fmt_complex) "GC_25", GC_25
    write (unit = *, fmt = fmt_complex) "GC_26", GC_26
    write (unit = *, fmt = fmt_complex) "GC_27", GC_27
    write (unit = *, fmt = fmt_complex) "GC_28", GC_28
    write (unit = *, fmt = fmt_complex) "GC_29", GC_29
    write (unit = *, fmt = fmt_complex) "GC_3", GC_3
    write (unit = *, fmt = fmt_complex) "GC_30", GC_30
    write (unit = *, fmt = fmt_complex) "GC_31", GC_31
    write (unit = *, fmt = fmt_complex) "GC_32", GC_32
    write (unit = *, fmt = fmt_complex) "GC_33", GC_33
    write (unit = *, fmt = fmt_complex) "GC_34", GC_34
    write (unit = *, fmt = fmt_complex) "GC_35", GC_35
    write (unit = *, fmt = fmt_complex) "GC_36", GC_36
    write (unit = *, fmt = fmt_complex) "GC_37", GC_37
    write (unit = *, fmt = fmt_complex) "GC_38", GC_38
    write (unit = *, fmt = fmt_complex) "GC_39", GC_39
    write (unit = *, fmt = fmt_complex) "GC_4", GC_4
    write (unit = *, fmt = fmt_complex) "GC_40", GC_40
    write (unit = *, fmt = fmt_complex) "GC_41", GC_41
    write (unit = *, fmt = fmt_complex) "GC_42", GC_42
    write (unit = *, fmt = fmt_complex) "GC_43", GC_43
    write (unit = *, fmt = fmt_complex) "GC_44", GC_44
    write (unit = *, fmt = fmt_complex) "GC_45", GC_45
    write (unit = *, fmt = fmt_complex) "GC_46", GC_46
    write (unit = *, fmt = fmt_complex) "GC_47", GC_47
    write (unit = *, fmt = fmt_complex) "GC_48", GC_48
    write (unit = *, fmt = fmt_complex) "GC_49", GC_49
    write (unit = *, fmt = fmt_complex) "GC_5", GC_5
    write (unit = *, fmt = fmt_complex) "GC_50", GC_50
    write (unit = *, fmt = fmt_complex) "GC_51", GC_51
    write (unit = *, fmt = fmt_complex) "GC_52", GC_52
    write (unit = *, fmt = fmt_complex) "GC_53", GC_53
    write (unit = *, fmt = fmt_complex) "GC_54", GC_54
    write (unit = *, fmt = fmt_complex) "GC_55", GC_55
    write (unit = *, fmt = fmt_complex) "GC_56", GC_56
    write (unit = *, fmt = fmt_complex) "GC_57", GC_57
    write (unit = *, fmt = fmt_complex) "GC_58", GC_58
    write (unit = *, fmt = fmt_complex) "GC_59", GC_59
    write (unit = *, fmt = fmt_complex) "GC_6", GC_6
    write (unit = *, fmt = fmt_complex) "GC_60", GC_60
    write (unit = *, fmt = fmt_complex) "GC_61", GC_61
    write (unit = *, fmt = fmt_complex) "GC_62", GC_62
    write (unit = *, fmt = fmt_complex) "GC_63", GC_63
    write (unit = *, fmt = fmt_complex) "GC_64", GC_64
    write (unit = *, fmt = fmt_complex) "GC_65", GC_65
    write (unit = *, fmt = fmt_complex) "GC_66", GC_66
    write (unit = *, fmt = fmt_complex) "GC_67", GC_67
    write (unit = *, fmt = fmt_complex) "GC_68", GC_68
    write (unit = *, fmt = fmt_complex) "GC_69", GC_69
    write (unit = *, fmt = fmt_complex) "GC_7", GC_7
    write (unit = *, fmt = fmt_complex) "GC_70", GC_70
    write (unit = *, fmt = fmt_complex) "GC_71", GC_71
    write (unit = *, fmt = fmt_complex) "GC_72", GC_72
    write (unit = *, fmt = fmt_complex) "GC_73", GC_73
    write (unit = *, fmt = fmt_complex) "GC_74", GC_74
    write (unit = *, fmt = fmt_complex) "GC_75", GC_75
    write (unit = *, fmt = fmt_complex) "GC_76", GC_76
    write (unit = *, fmt = fmt_complex) "GC_77", GC_77
    write (unit = *, fmt = fmt_complex) "GC_78", GC_78
    write (unit = *, fmt = fmt_complex) "GC_79", GC_79
    write (unit = *, fmt = fmt_complex) "GC_8", GC_8
    write (unit = *, fmt = fmt_complex) "GC_80", GC_80
    write (unit = *, fmt = fmt_complex) "GC_81", GC_81
    write (unit = *, fmt = fmt_complex) "GC_82", GC_82
    write (unit = *, fmt = fmt_complex) "GC_83", GC_83
    write (unit = *, fmt = fmt_complex) "GC_84", GC_84
    write (unit = *, fmt = fmt_complex) "GC_85", GC_85
    write (unit = *, fmt = fmt_complex) "GC_86", GC_86
    write (unit = *, fmt = fmt_complex) "GC_87", GC_87
    write (unit = *, fmt = fmt_complex) "GC_88", GC_88
    write (unit = *, fmt = fmt_complex) "GC_89", GC_89
    write (unit = *, fmt = fmt_complex) "GC_9", GC_9
    write (unit = *, fmt = fmt_complex) "GC_90", GC_90
    write (unit = *, fmt = fmt_complex) "GC_91", GC_91
    write (unit = *, fmt = fmt_complex) "I4a33", I4a33
    write (unit = *, fmt = fmt_complex) "I4a22", I4a22
    write (unit = *, fmt = fmt_complex) "I4a11", I4a11
    write (unit = *, fmt = fmt_complex) "I3a33", I3a33
    write (unit = *, fmt = fmt_complex) "I3a22", I3a22
    write (unit = *, fmt = fmt_complex) "I3a11", I3a11
    write (unit = *, fmt = fmt_complex) "I2a33", I2a33
    write (unit = *, fmt = fmt_complex) "I2a22", I2a22
    write (unit = *, fmt = fmt_complex) "I2a11", I2a11
    write (unit = *, fmt = fmt_complex) "I1a33", I1a33
    write (unit = *, fmt = fmt_complex) "I1a22", I1a22
    write (unit = *, fmt = fmt_complex) "I1a11", I1a11
  end subroutine print_parameters
end module parameters_sm_ufo

! ----------------------------------------------------------------------
! This is new PDE_3D_src that works.
! ----------------------------------------------------------------------

! Must set "NVAR, IPDE_SYSTEM" for this to work.







      SUBROUTINE PRIM_TO_WAVSPD ( idirection,
     1         gamma, B0, Vpt, S_min, S_max)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara, Deepak Bhoriya

! In 3D, this subroutine uses the primitive variables "Vpt" to build the
! extremal wave speeds "S_min, S_max" in the "idirection", for the PDE
! system given by "IPDE_SYSTEM".

! "Vpt ( NVAR + NFLUID)" is input, "S_min, S_max" are output.

! "IPDE_SYSTEM == 1" for MHD
! "IPDE_SYSTEM == 2" for CGL System
! "IPDE_SYSTEM == 3" for CGL System with p_elec

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER idirection

      REAL gamma, Vpt ( NVAR + NFLUID), S_min, S_max

      REAL, DIMENSION ( 3) :: B0

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      REAL rho, prs, vlx, vly, vlz,
     1     tempa,
     1     bfx, bfy, bfz,  b2_4pirho, cs, aa, bb, cc, m_f,

     1     ob4pi, ob8pi, ubx, uby, ubz,
     1     Bsq, p_perp, p_par, p_elec

      REAL, PARAMETER :: pi = 4.0 * ATAN ( 1.0), smallnum = 1.0e-12

! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 1)

      rho = Vpt ( 1)

      vlx = Vpt ( 2)
      vly = Vpt ( 3)
      vlz = Vpt ( 4)

      prs = Vpt ( 5)

      bfx = Vpt ( 6) + B0 ( 1)
      bfy = Vpt ( 7) + B0 ( 2)
      bfz = Vpt ( 8) + B0 ( 3)

      cs = SQRT ( gamma*prs/rho)

      b2_4pirho = ( bfx**2 + bfy**2 + bfz**2) / (4.0*pi*rho)

      aa = 1.0
      bb = - (cs**2 + b2_4pirho)

! ----------------------------------------------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the x-direction

      IF ( idirection .EQ. 1) THEN

         cc = cs**2 * bfx**2 / (4.0*pi*rho)
         m_f = SQRT( ( -bb + SQRT(bb**2 - 4.0 * aa * cc)) / (2.0*aa))

         S_min = vlx - m_f

         S_max = vlx + m_f

      END IF

! ------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the y-direction

      IF ( idirection .EQ. 2) THEN

         cc = cs**2 * bfy**2 / (4.0*pi*rho)
         m_f = SQRT( ( -bb + SQRT(bb**2 - 4.0 * aa * cc)) / (2.0*aa))

         S_min = vly - m_f

         S_max = vly + m_f

      END IF

! ------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the z-direction

      IF ( idirection .EQ. 3) THEN

         cc = cs**2 * bfz**2 / (4.0*pi*rho)
         m_f = SQRT( ( -bb + SQRT(bb**2 - 4.0 * aa * cc)) / (2.0*aa))

         S_min = vlz - m_f

         S_max = vlz + m_f

      END IF

#endif
! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 2)

      ob4pi = 0.25 / pi  ;  ob8pi = 0.125 / pi

      rho = Vpt ( 1)

      vlx = Vpt ( 2)
      vly = Vpt ( 3)
      vlz = Vpt ( 4)

      p_par  = Vpt ( 5)

      bfx = Vpt ( 6) + B0 ( 1)
      bfy = Vpt ( 7) + B0 ( 2)
      bfz = Vpt ( 8) + B0 ( 3)

      p_perp = Vpt ( 9)

! Some temporary variables
      Bsq = ( bfx**2 + bfy**2 + bfz**2)
      ubx = bfx / (SQRT(Bsq) + smallnum )
      uby = bfy / (SQRT(Bsq) + smallnum )
      ubz = bfz / (SQRT(Bsq) + smallnum )

! ----------------------------------------------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the x-direction

      IF ( idirection .EQ. 1) THEN

        aa = 2.0*rho
        bb = -(Bsq*ob4pi + 2.0*p_perp + ubx**2*(2.0*p_par-p_perp))
        cc = - ( p_perp**2 * ubx**2 * (1.0-ubx**2)
     1         - 3.0*p_par*p_perp*ubx**2 *( 2.0 - ubx**2)
     1         + 3.0*p_par**2*ubx**4 - 3.0*bfx**2*p_par*ob4pi
     1         ) / (2.0*rho)

! fast magneto
        m_f = SQRT( ( - bb + SQRT (bb**2 - 4.0*aa*cc)) / (2.0*rho) )

! alfven
        tempa = SQRT ( (bfx**2 * ob4pi / rho)
     1                  - (p_par - p_perp) * ubx**2 / rho )

        S_min = AMIN1 (vlx - tempa, vlx - m_f)

        S_max = AMAX1 (vlx + tempa, vlx + m_f)

      END IF

! ------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the y-direction

      IF ( idirection .EQ. 2) THEN

        aa = 2.0*rho
        bb = -(Bsq*ob4pi + 2.0*p_perp + uby**2*(2.0*p_par-p_perp))
        cc = - ( p_perp**2 * uby**2 * (1.0-uby**2)
     1         - 3.0*p_par*p_perp*uby**2 *( 2.0 - uby**2)
     1         + 3.0*p_par**2*uby**4 - 3.0*bfy**2*p_par*ob4pi
     1         ) / (2.0*rho)

! fast magneto
        m_f = SQRT( ( - bb + SQRT (bb**2 - 4.0*aa*cc)) / (2.0*rho) )

! alfven
        tempa = SQRT ( (bfy**2 * ob4pi / rho)
     1        - (p_par - p_perp) * uby**2 / rho)

        S_min = AMIN1 (vly - tempa, vly - m_f)

        S_max = AMAX1 (vly + tempa, vly + m_f)

      END IF

! ------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the z-direction

      IF ( idirection .EQ. 3) THEN

        aa = 2.0*rho
        bb = -(Bsq*ob4pi + 2.0*p_perp + ubz**2*(2.0*p_par-p_perp))
        cc = - ( p_perp**2 * ubz**2 * (1.0-ubz**2)
     1         - 3.0*p_par*p_perp*ubz**2 *( 2.0 - ubz**2)
     1         + 3.0*p_par**2*ubz**4 - 3.0*bfz**2*p_par*ob4pi
     1         ) / (2.0*rho)

! fast magneto
        m_f = SQRT( ( - bb + SQRT (bb**2 - 4.0*aa*cc)) / (2.0*rho) )

! alfven
        tempa = SQRT( (bfz**2 * ob4pi / rho)
     1        - (p_par - p_perp) * ubz**2 / rho)

        S_min = AMIN1 (vlz - tempa, vlz - m_f)

        S_max = AMAX1 (vlz + tempa, vlz + m_f)

      END IF

#endif
! End of "#if ( IPDE_SYSTEM == 2)"
! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 3)

      ob4pi = 0.25 / pi  ;  ob8pi = 0.125 / pi

      rho = Vpt ( 1)

      vlx = Vpt ( 2)
      vly = Vpt ( 3)
      vlz = Vpt ( 4)

      p_par  = Vpt ( 5)

      bfx = Vpt ( 6) + B0 ( 1)
      bfy = Vpt ( 7) + B0 ( 2)
      bfz = Vpt ( 8) + B0 ( 3)

      p_perp = Vpt ( 9)
      p_elec = Vpt (10)

! Some temporary variables
      Bsq = ( bfx**2 + bfy**2 + bfz**2)
      ubx = bfx / (SQRT(Bsq) + smallnum )
      uby = bfy / (SQRT(Bsq) + smallnum )
      ubz = bfz / (SQRT(Bsq) + smallnum )

! ----------------------------------------------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the x-direction

      IF ( idirection .EQ. 1) THEN

        aa = 2.0*rho
        bb = -(Bsq*ob4pi + 2.0*p_perp + ubx**2*(2.0*p_par-p_perp))
        bb = bb - (5.0/3.0)*p_elec
        cc = - ( p_perp**2 * ubx**2 * (1.0-ubx**2)
     1         - 3.0*p_par*p_perp*ubx**2 *( 2.0 - ubx**2)
     1         + 3.0*p_par**2*ubx**4 - 3.0*bfx**2*p_par*ob4pi
     1         ) / (2.0*rho)

        cc = cc + ( (5.0/3.0)*p_elec* ( bfx**2*ob4pi
     1               + 3.0*ubx**2*p_par
     1               + ubx**4*(p_perp - 4.0*p_par) ) ) / (2.0*rho)

! fast magneto
        m_f = SQRT( ( - bb + SQRT (bb**2 - 4.0*aa*cc)) / (2.0*rho) )

! alfven
        tempa = SQRT ( (bfx**2 * ob4pi / rho)
     1                  - (p_par - p_perp) * ubx**2 / rho )

        S_min = AMIN1 (vlx - tempa, vlx - m_f)

        S_max = AMAX1 (vlx + tempa, vlx + m_f)

      END IF

! ------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the y-direction

      IF ( idirection .EQ. 2) THEN

        aa = 2.0*rho
        bb = -(Bsq*ob4pi + 2.0*p_perp + uby**2*(2.0*p_par-p_perp))
        bb = bb - (5.0/3.0)*p_elec
        cc = - ( p_perp**2 * uby**2 * (1.0-uby**2)
     1         - 3.0*p_par*p_perp*uby**2 *( 2.0 - uby**2)
     1         + 3.0*p_par**2*uby**4 - 3.0*bfy**2*p_par*ob4pi
     1         ) / (2.0*rho)

        cc = cc + ( (5.0/3.0)*p_elec* ( bfy**2*ob4pi
     1               + 3.0*uby**2*p_par
     1               + uby**4*(p_perp - 4.0*p_par) ) ) / (2.0*rho)

! fast magneto
        m_f = SQRT( ( - bb + SQRT (bb**2 - 4.0*aa*cc)) / (2.0*rho) )

! alfven
        tempa = SQRT ( (bfy**2 * ob4pi / rho)
     1        - (p_par - p_perp) * uby**2 / rho)

        S_min = AMIN1 (vly - tempa, vly - m_f)

        S_max = AMAX1 (vly + tempa, vly + m_f)

      END IF

! ------------------------------
! Here, we find the extremal speeds "S_min, S_max" in the z-direction

      IF ( idirection .EQ. 3) THEN

        aa = 2.0*rho
        bb = -(Bsq*ob4pi + 2.0*p_perp + ubz**2*(2.0*p_par-p_perp))
        bb = bb - (5.0/3.0)*p_elec
        cc = - ( p_perp**2 * ubz**2 * (1.0-ubz**2)
     1         - 3.0*p_par*p_perp*ubz**2 *( 2.0 - ubz**2)
     1         + 3.0*p_par**2*ubz**4 - 3.0*bfz**2*p_par*ob4pi
     1         ) / (2.0*rho)

        cc = cc + ( (5.0/3.0)*p_elec* ( bfz**2*ob4pi
     1               + 3.0*ubz**2*p_par
     1               + ubz**4*(p_perp - 4.0*p_par) ) ) / (2.0*rho)

! fast magneto
        m_f = SQRT( ( - bb + SQRT (bb**2 - 4.0*aa*cc)) / (2.0*rho) )

! alfven
        tempa = SQRT( (bfz**2 * ob4pi / rho)
     1        - (p_par - p_perp) * ubz**2 / rho)

        S_min = AMIN1 (vlz - tempa, vlz - m_f)

        S_max = AMAX1 (vlz + tempa, vlz + m_f)

      END IF

#endif
! End of "#if ( IPDE_SYSTEM == 3)"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END SUBROUTINE PRIM_TO_WAVSPD




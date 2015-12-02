module hswf_mod

 use constants_mod,      only: grav, rdgas, rvgas, RADIAN, kappa, radius
 use mpp_domains_mod,    only: mpp_update_domains
 use time_manager_mod,   only: time_type, get_date, get_time
 use diag_manager_mod,   only: send_data
 use lin_cld_microphys_mod,         only: lin_cld_microphys_driver, sg_conv, qsmith

 use mpp_mod,            only: FATAL, mpp_error, stdlog, &
                               mpp_npes, mpp_get_current_pelist, get_unit
!epg: I added these for reading the namelist
 use fv_grid_tools_mod,  only: area
 use fv_grid_utils_mod,  only: g_sum
 use fv_mp_mod,          only: domain, gid
!epg: I added 'gid' for use in the Held_Suarez_init subroutine 
 use fv_diagnostics_mod, only: prt_maxmin
 use fv_timing_mod,      only: timing_on, timing_off

      implicit none
!-----------------------------------------------------------------------
      logical :: rf_initialized = .false.

      private
      public :: Held_Suarez_Strat, Held_Suarez_Tend, Sim_phys, age_of_air, &
           Held_Suarez_init
!epg: added Held_Suarez_init, so that we can use a namelist to control 
!     the forcing


!-----------------------------------------------------------------------

!epg: namelist parameters for the Held-Suarez troposphere
!     they are set to default values in the original H-S 1994 paper
      real :: ty = 60.0  ! equator-pole temp gradient (K)
      real :: tz = 10.0  ! vertical temp gradient (K)
      real :: t0 = 200.0 ! isothermal stratospheric temperature
                         ! (overwritten by t_tropopause if PK stratosphere invoked)
      real :: eps = 0.0  ! asymmetry in temperature forcing btw hemispheres (K)
      real :: tau_v = 1.0   ! surface Rayleigh damping time scale (days) 
      real :: tau_a = 40.0  ! free atmosphere Newtonian relaxation time scale (days)
      real :: tau_s = 4.0   ! tropical boundary layer Newt. relax. time scale (days)

!epg: namelist parameters for the Polvani-Kushner stratosphere
      logical :: pk_stratosphere = .false.  ! option to turn it on
      logical :: polar_vortex_only = .false. ! a flag that allows you to only
                                             ! add a cold vortex, but leave
                                             ! everything else isothermal at
                                             ! t0 K, as with the standard HS 
                                             ! forcing.
      real :: p_tropopause = 10000.0      ! lower boundary of the Kushner-Polvani
                                          !   stratospheric forcing (in Pa)
      real :: vortex_edge = 50.0        ! polar vortex edge location (in deg.)
      real :: vortex_edge_width = 10.0  ! polar vortex edge width
      real :: vortex_gamma = -2.e-3     ! polar vortex lapse rate (in K/m)

!epg: namelist parameters for tracers
      logical :: do_age_of_air = .false.   ! compute the age of air
      logical :: do_delta_tracer = .false. ! compute the age of air spectrum
      logical :: do_pulse_tracer = .false. ! a pulse of tracer, to test conservation
      real :: p_source = 70000.0     ! source/sinks of tracers are below this level

!epg: namelist parameters for global warming-like simulations (as explored by
!     Butler etal. 2010 and Wang et al. 2011
      logical :: do_ghg_heating = .false. ! flag to turn it on/off
      real :: ghg_heating_amplitude = 0.1 ! amplitude K/day (see Wang et al. 2011)
      real :: ghg_ypos = 0.0              ! meridional position, degrees
      real :: ghg_ywid = 22.91            ! meridional width, degrees
      real :: ghg_zpos = 0.300            ! vertical position, sigma levels
      real :: ghg_zwid = 0.11             ! vertical extent, sigma coords


contains


!epg: load in the namelist
 subroutine Held_Suarez_init

   ! for identifying the controlling processor
   integer, allocatable :: pelist(:)
   logical :: master
   integer :: commID

   ! for reading the namelist
   character(len=80) :: filename
   integer :: ios, f_unit
   logical :: exists
   integer :: unit

   !epg: for converting namelist variable units
   real :: deg_to_rad = 3.14159265358979/180.0

   !epg: a hswf namelists 
   namelist /hswf_nml/ty, tz, t0, eps, tau_v, tau_a, tau_s, pk_stratosphere, &
        polar_vortex_only, p_tropopause, vortex_edge, vortex_edge_width, &
        vortex_gamma, stratospheric_damping, strat_damping_pbottom, tau_strat, &
        do_age_of_air, do_delta_tracer, do_pulse_tracer, p_source, &
        do_ghg_heating, ghg_heating_amplitude, ghg_ypos, ghg_ywid, ghg_zpos, ghg_zwid

   allocate( pelist(mpp_npes()) )
   call mpp_get_current_pelist( pelist, commID=commID )
   master = gid==0

   !epg: read parameters from namelist
   filename = "input.nml"
   inquire(file=filename,exist=exists)

   if (.not. exists) then  ! This will be replaced with fv_error wrapper
      if(master) write(6,*) "file ",trim(filename)," doesn't exist"
      call mpp_error(FATAL,'FV core terminating')
   else
      f_unit = get_unit()
      open (f_unit,file=filename)
      ! Read hswf namelist
      rewind (f_unit)
      read (f_unit,hswf_nml,iostat=ios)
      if (ios .gt. 0) then
         if(master) write(6,*) 'hswf_nml ERROR: reading ',trim(filename),', iostat=',ios
         call mpp_error(FATAL,'FV core terminating')
      endif
      unit = stdlog()
      write(unit, nml=hswf_nml)
      close (f_unit)
   endif

   ! epg: for applying the global warming like perturbation
   if (do_ghg_heating) then
      ! convert degrees to radians
      ghg_ypos = ghg_ypos * deg_to_rad
      ghg_ywid = ghg_ywid * deg_to_rad
      ! convert K/day to K/second         
      ghg_heating_amplitude = ghg_heating_amplitude/86400.0
   end if


 end subroutine Held_Suarez_init

 
 subroutine Held_Suarez_Tend(npx, npy, npz, is, ie, js, je, ng, nq,   &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, u_dt, v_dt, t_dt, q_dt, agrid,  &
                              delz, hydrostatic, ak, bk, ks,   &
                              strat, rayf, master, Time, time_total)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq
      logical, intent(IN)    :: hydrostatic
      real   , INTENT(IN   ) ::  delz(is:ie,js:je,npz)

      real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
      real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
      real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
      real   , INTENT(INOUT) :: peln(is  :ie     ,1:npz+1,js   :je     )
      real   , INTENT(INOUT) ::  pkz(is  :ie     ,js   :je     ,1:npz)

      real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
      real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

! Tendencies:
      real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
      real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
      real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)


      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
      integer, INTENT(IN   ) :: ks

      real   , INTENT(IN   ) :: pdt
      logical, INTENT(IN   ) :: strat, rayf, master

      type(time_type), intent(in) :: Time
      real, INTENT(IN), optional:: time_total

! Local
      real pref(npz)
      integer  i,j,k,l
      !epg: these are namelist parameters now: t0, ty, tz
      real akap
      real  p0, sday, rkv, rka, rks, rkt, sigb, rsgb
      real  tmp
      real  ap0k, algpk
      real  tey, tez, fac, pw, sigl
      real  h0, dz
      real  dt_tropic
      real  rmr, rms
      real  relx, tau
      real  t_st, t_ms
      real  rdt, f1
      real  pc, c1

      real, allocatable, SAVE ::  rf(:)
      real :: frac(is-ng:ie+ng,js-ng:je+ng)
      real ::  pl(is:ie, js:je, 1:npz)
      real :: teq(is:ie, js:je, 1:npz)
      real rdg


! epg: for the Polvani-Kushner Stratosphere
      ! the following fields are used to construct a realistic stratosphere
      !----------- US Standard Atmospheric Temperature 1976 ------------------
      integer, parameter :: satm_levs = 8
      real, parameter, dimension(satm_levs) :: &
           satm_p = (/ 100000.00000000, &
                        22336.11050922, &
                         5403.29501078, &
                          856.66783593, &
                          109.45601338, &
                           66.06353133, &
                            3.90468337, &
                            0.00000001 /), &
           satm_t = (/ 288.150, &
                       216.650, &
                       216.650, &
                       228.650, &
                       270.650, &
                       270.650, &
                       214.650, &
                       186.946 /), &
            satm_g = (/ -6.5e-3, &
                        0.0e-3, &
                        1.0e-3, &
                        2.8e-3, &
                        0.0e-3, &
                        -2.8e-3, &
                        -2.0e-3, &
                        0.0e-3 /)
       real :: t_tropopause = 216.650
       real :: rad_to_deg = 180./3.14159265358979
       real :: lat_weight, t_satm, t_vortex

!epg: for stratospheric Rayleigh damping
       real :: rk_strat
       real :: pres
       real :: rayd




      rdg = -rdgas / grav

      !epg: these are set in the namelist
      !ty = 60.0
      !tz = 10.0             ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      !epg: set in namelist t0 = 200.         ! stratosphere temperature
      h0 = 7.
      sday = 24.*3600.
      rdt = 1. / pdt

      ! epg: these values are now set with the namelist
      rkv = pdt/(tau_v*sday)     ! Rayleigh friction near surface 
      rka = pdt/ (tau_a*sday)    ! free atmosphere Newtonian damping
      rks = pdt/ (tau_s*sday)    ! tropical boundary layer Newtonian damping
      !epg: also need time scale for stratospheric sponge
      if (tau_strat > 0) then
         rk_strat = pdt / (tau_strat*sday)
      else
         rk_strat = 0.
      endif

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)


      do k=1,npz
         pref(k) = ak(k) + bk(k)*1.E5
      enddo

!epg: w/ the Polvani-Kushner stratosphere, the temperature of the 
!     stratosphere is set to t_tropopause
      if ( pk_stratosphere .and. (.not. polar_vortex_only) ) then
         t0=t_tropopause
      endif

! Setup...
      if ( rayf .and. (.not. rf_initialized) ) then
          allocate( rf(npz) )
          c1 = 1. / (36.*3600)
          pc = 1.
          if(master) write(6,*) 'HSWF Forcing ...' 
          do k=1,ks
             tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
             rf(k) = c1*(1.+tanh(log10(pc/tmp)))
             if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
             rf(k) = 1./(1.+pdt*rf(k))
          enddo
          if(master) write(6,*) ' '
          rf_initialized = .true.
      endif

        do k=1,npz
           do j=js,je
              do i=is,ie
                 pl(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
              enddo
           enddo 
        enddo

        if ( .not. hydrostatic ) then
           do k=1,npz
              do j=js,je
                 do i=is,ie
                    pkz(i,j,k) = pl(i,j,k)**akap
                 enddo
              enddo 
           enddo
         endif
        
! Temperature forcing...
         do k=npz,1,-1
            do j=js,je
               do i=is,ie
                  !epg: we'll add the "eps" asymmetry factor
                  !tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) )
                  tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) &
                       - eps*SIN(agrid(i,j,2)) )
                  tez = tz*( ap0k/akap )*COS(agrid(i,j,2))*COS(agrid(i,j,2))

                  if (pk_stratosphere .and. pl(i,j,k) < p_tropopause ) then
! epg: the Polvani-Kushner temperature profile
                     lat_weight = 0.5 * ( 1. + &
                       tanh((agrid(i,j,2)*rad_to_deg-vortex_edge)/vortex_edge_width) )

                     if (polar_vortex_only) then
                        ! with this option, the stratosphere is isothermal except for 
                        ! the polar vortex
                        t_satm=t0
                     else
                        ! construct the US standard atmosphere temperature 
                        do l = 1, satm_levs
                           if ( pl(i,j,k) > satm_p(l+1)) then
                              t_satm = satm_t(l) * &
                                   (pl(i,j,k)/satm_p(l))**(-rdgas*satm_g(l)/grav)
                              !if (i == 1 .and. j == 1 .and. master) then
                              !   print *, k, t_satm
                              !endif
                              exit
                           endif
                        enddo
                     endif
                     ! compute vortex temp
                     t_vortex = t_tropopause * &
                          (pl(i,j,k)/p_tropopause)**(-rdgas*vortex_gamma/grav)
                     ! t_equilibrium is a weighted average of the US standard
                     ! atmosphere or the cold polar vortex
                     teq(i,j,k) = (1-lat_weight) * t_satm + lat_weight * &
                          t_vortex
                     ! lastly, computed the temperature tendency (implicit scheme)
                     t_dt(i,j,k) = rka*(teq(i,j,k)-pt(i,j,k))/(1.+rka) * rdt
 ! epg: end Polvani-Kushner temperature profile                    

                  elseif (strat .and. pl(i,j,k) < 10000.    &
                       .and. pl(i,j,k) > 100.  )  then
                     dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
!
! Lapse rate above tropic stratopause is 2.25 deg/km
! Relaxation time is t_st days at 100 mb (as H-S) and gradually
! decreases to t_ms Days at and above the stratopause
!
                     relx =  t_ms + tau*log(0.01*pl(i,j,k))
                     relx = pdt/(relx*sday)
                     dt_tropic = 2.25*COS(agrid(i,j,2)) * dz
                     teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
                     t_dt(i,j,k) = relx*(teq(i,j,k)-pt(i,j,k))/(1.+relx) * rdt
                  elseif (strat .and. pl(i,j,k) <= 100.)  then
!
! Mesosphere: defined as the region above 1 mb
!
                     dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
                     dt_tropic = -2.25*COS(agrid(i,j,2)) * dz
                     teq(i,j,k) = teq(i,j,k+1) + dt_tropic
                     t_dt(i,j,k) = ((pt(i,j,k)+rms*teq(i,j,k))*rmr - pt(i,j,k))*rdt
                  else

! Trop:  strictly Held-Suarez

                     sigl = pl(i,j,k)/pe(i,npz+1,j)
                     f1 = max(0., (sigl-sigb) * rsgb )
                     teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
                     teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
                     rkt = rka + (rks-rka)*f1*(COS(agrid(i,j,2))**4.0)
                     t_dt(i,j,k) = rkt*(teq(i,j,k)-pt(i,j,k))/(1.+rkt) * rdt
                  endif

                  !epg: this allows you to prescribe a CO2-like warming as in Wang et al 2011.
                  if (do_ghg_heating) then ! co2 warming
                     t_dt(i,j,k) = t_dt(i,j,k) + ghg_heating_amplitude &
                          *exp( -0.5*((agrid(i,j,2)-ghg_ypos)/ghg_ywid)**2 &
                          -0.5*((pl(i,j,k)/pe(i,npz+1,j)-ghg_zpos)/ghg_zwid)**2)
                  endif

               enddo     !i-loop
            enddo     !j-loop
         enddo     !k-loop

         if ( .not. hydrostatic ) then
!------------------------------------------------------------------------
! The nonhydrostatic pressure changes if there is heating (under constant
! volume and mass is locally conserved).
!------------------------------------------------------------------------
           do k=1,npz
              do j=js,je
                 do i=is,ie
                      pkz(i,j,k) = (rdg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k))**akap
                 enddo
              enddo 
           enddo
        endif


! Velocity dissipation damping
      do 2000 k=1,npz

        if (rayf .and. k <= ks) then
! Apply Rayleigh friction
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = ua(i,j,k)*(rf(k) - 1.) * rdt
                v_dt(i,j,k) = va(i,j,k)*(rf(k) - 1.) * rdt
             enddo
          enddo
        else
! Surface Rayleigh friction according to Held-Suarez
           do j=js,je
              do i=is,ie
                 sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
                 frac(i,j) = max(0., (sigl-sigb)*rsgb ) * rkv
                 if (frac(i,j) > 0.) then
                     u_dt(i,j,k) = -ua(i,j,k)*frac(i,j)/(1.+frac(i,j)) * rdt
                     v_dt(i,j,k) = -va(i,j,k)*frac(i,j)/(1.+frac(i,j)) * rdt
                 endif
               enddo
           enddo
        endif

2000  continue

!epg: tracers - age of air
        if ( do_age_of_air ) then
           ! make the tracer starts out at zero!
           if( time_total <= 3600 ) then
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       q(i,j,k,1) = 0.0
                    enddo
                 enddo
              enddo
           else
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       if( pe(i,k,j) >= p_source ) then
                          q(i,j,k,1) = time_total/86400.0
                       endif
                    enddo
                 enddo
              enddo
           endif

        endif

!epg: tracers - age spectrum
        if ( do_delta_tracer ) then

           if ( time_total/86400.0 > 330.0 ) then
              !remove any tracer below the source level
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       if( pe(i,k,j) >= p_source ) then
                          q(i,j,k,2) = 0.0
                       endif
                    enddo
                 enddo
              enddo

           elseif ( time_total/86400.0 > 300.0 ) then
              ! initialize the tracer
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       if( pe(i,k,j) >= p_source ) then
                          q(i,j,k,2) = 1.0
                       endif
                    enddo
                 enddo
              enddo
           else
              ! ensure that the tracer is initialy zero
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       q(i,j,k,2) = 0.0
                    enddo
                 enddo
              enddo
           endif

        endif

!epg: tracers - pulse tracer
        if ( do_pulse_tracer ) then

           if ( time_total/86400.0 > 330.0 ) then
              ! do nothing in this case

           elseif ( time_total/86400.0 > 300.0 ) then
              !initialize the pulse
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       if( pe(i,k,j) >= p_source ) then
                          q(i,j,k,3) = 1.0
                       endif
                    enddo
                 enddo
              enddo

           else
              !the tracer should start at zero
              do k=1,npz
                 do j=js,je
                    do i=is,ie
                       q(i,j,k,3) = 0.0
                    enddo
                 enddo
              enddo
           endif
        endif



#ifdef DO_AGE
      if( nq/=0 )     &
      call age_of_air(is, ie, js, je, npz, ng, time_total, pe, q(is-ng,js-ng,1,nq))
#endif


 end subroutine Held_Suarez_Tend

 subroutine Held_Suarez_Strat(npx, npy, npz, is, ie, js, je, ng, nq,  &
                              u, v, pt, q, pe, delp, peln, pkz, pdt,  &
                              ua, va, agrid, ak, bk, ks, strat,  &
                              rayf, master, Time, time_total)

      integer, INTENT(IN   ) :: npx, npy, npz
      integer, INTENT(IN   ) :: is, ie, js, je, ng, nq

	      real   , INTENT(INOUT) ::    u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
	      real   , INTENT(INOUT) ::    v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
	      real   , INTENT(INOUT) ::   pt(is-ng:ie+  ng,js-ng:je+  ng,npz)
	      real   , INTENT(INOUT) :: delp(is-ng:ie+  ng,js-ng:je+  ng,npz)
	      real   , INTENT(INOUT) ::    q(is-ng:ie+  ng,js-ng:je+  ng,npz, nq)
	      real   , INTENT(INOUT) ::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
	      real   , INTENT(INOUT) :: peln(is   :ie     ,1:npz+1,js   :je     )
	      real   , INTENT(INOUT) ::  pkz(is   :ie     ,js   :je     ,1:npz)

	      real   , INTENT(INOUT) ::   ua(is-ng:ie+ng,js-ng:je+ng,npz)
	      real   , INTENT(INOUT) ::   va(is-ng:ie+ng,js-ng:je+ng,npz)

	      real   , INTENT(IN   ) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
	      real   , INTENT(IN   ) :: ak(npz+1), bk(npz+1)
              integer, INTENT(IN   ) :: ks

	      real   , INTENT(IN   ) :: pdt
	      logical, INTENT(IN   ) :: strat, rayf, master

              type(time_type), intent(in) :: Time
	      real, INTENT(IN), optional:: time_total

! Local
	      real pref(npz)
              integer  i,j,k
	      real  ty, tz, akap 
	      real  p0, t0, sday, rkv, rka, rks, rkt, sigb, rsgb
	      real  tmp
	      real  ap0k, algpk
	      real  tey, tez, fac, pw, sigl
	      real  h0, dz
	      real  dt_tropic
	      real  rmr, rms
	      real  relx, tau
	      real  t_st, t_ms
	      real  f1
	      real  pc, c1

	      real, allocatable, SAVE ::  rf(:)
	      real :: frac(is-ng:ie+ng,js-ng:je+ng)
	      real ::  pl(is:ie, js:je, 1:npz)
	      real :: teq(is:ie, js:je, 1:npz)

      ty = 60.0
      tz = 10.0             ! Original value from H-S was 10.
      akap = 2./7.

      p0 = 100000.
      t0 = 200.
      h0 = 7.
      sday = 24.*3600.

      rkv = 0.5*pdt/sday
      rka = pdt/ (40.*sday)      ! was 40 days
      rks = pdt/ (4.*sday)       ! was 4 days

! For strat-mesosphere
      t_ms = 10.
      t_st = 40.
      tau = (t_st - t_ms) / log(100.)
      rms = pdt/(t_ms*sday)
      rmr =  1./(1.+rms)

      sigb = 0.7
      rsgb = 1./(1.-sigb)
      ap0k = 1./p0**akap
      algpk = log(ap0k)

      do k=1,npz
         pref(k) = ak(k) + bk(k)*1.E5
      enddo

! Setup...
      if ( rayf .and. (.not. rf_initialized) ) then
		  allocate( rf(npz) )
		  c1 = 1. / (36.*3600)
		  pc = 1.
		  if(master) write(6,*) 'HSWF Forcing ...' 
		  do k=1,ks
		     tmp = (ak(k+1)-ak(k))/log(ak(k+1)/ak(k))
		     rf(k) = c1*(1.+tanh(log10(pc/tmp)))
		     if(master) write(6,*) k, 0.01*tmp, 1./(rf(k)*sday)
		     rf(k) = 1./(1.+pdt*rf(k))
		  enddo
		  if(master) write(6,*) ' '
		  rf_initialized = .true.
		endif

		do k=1,npz
		   do j=js,je
		      do i=is,ie
			 pl(i,j,k) = delp(i,j,k) / ( peln(i,k+1,j)-peln(i,k,j))
		      enddo
		   enddo 
		enddo

	! Temperature forcing...
		do k=npz,1,-1
		   do j=js,je
		      do i=is,ie
			 tey = ap0k*( 315.0 - ty*SIN(agrid(i,j,2))*SIN(agrid(i,j,2)) )
			 tez = tz*( ap0k/akap )*COS(agrid(i,j,2))*COS(agrid(i,j,2)) 

			 if (strat .and. pl(i,j,k) < 10000.    &
				   .and. pl(i,j,k) > 100.  )  then
			   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
	!
	! Lapse rate above tropic stratopause is 2.25 deg/km
	! Relaxation time is t_st days at 100 mb (as H-S) and gradually
	! decreases to t_ms Days at and above the stratopause
	!
			   relx =  t_ms + tau*log(0.01*pl(i,j,k))
			   relx = pdt/(relx*sday)
			   dt_tropic = 2.25*COS(agrid(i,j,2)) * dz
			   teq(i,j,k)  =  teq(i,j,k+1) + dt_tropic
			   pt(i,j,k) = (pt(i,j,k)+relx*teq(i,j,k))/(1.+relx)
			 elseif (strat .and. pl(i,j,k) <= 100.)  then
	!
	! Mesosphere: defined as the region above 1 mb
	!
			   dz = h0 * log(pl(i,j,k+1)/pl(i,j,k))
			   dt_tropic = -2.25*COS(agrid(i,j,2)) * dz
			   teq(i,j,k) = teq(i,j,k+1) + dt_tropic
			   pt(i,j,k) = (pt(i,j,k)+rms*teq(i,j,k))*rmr
			 else

	! Trop:  strictly Held-Suarez

			   sigl = pl(i,j,k)/pe(i,npz+1,j)
			   f1 = max(0., (sigl-sigb) * rsgb )
			   teq(i,j,k) = tey - tez*(log(pkz(i,j,k))+algpk)
			   teq(i,j,k) = max(t0, teq(i,j,k)*pkz(i,j,k))
			   rkt = rka + (rks-rka)*f1*(COS(agrid(i,j,2))**4.0)
			   pt(i,j,k) = (pt(i,j,k)+rkt*teq(i,j,k))/(1.+rkt)
			 endif
		      enddo     !i-loop
		   enddo     !j-loop
		enddo     !k-loop

	! Velocity dissipation damping
	      do 2000 k=1,npz

		if (rayf .and. k <= ks) then

	! Apply Rayleigh friction
		  do j=js,je+1
		    do i=is,ie
		      u(i,j,k) = u(i,j,k)*rf(k)
		    enddo
		  enddo
		  do j=js,je
		    do i=is,ie+1
		      v(i,j,k) = v(i,j,k)*rf(k)
		    enddo
		  enddo

		else

	! Surface Rayleigh friction according to Held-Suarez
		  do j=js,je
		    do i=is,ie
		      sigl = 0.5*(pe(i,k,j)+pe(i,k+1,j)) / pe(i,npz+1,j)
		      frac(i,j) = max(0., (sigl-sigb)*rsgb )
		    enddo
		  enddo
#if defined(SPMD)
		  call mpp_update_domains( frac, domain )
#endif
	! Backward adjustment
		  do j=js,je+1
		    do i=is,ie+1
		      fac = frac(i,j)+frac(i,j-1)
		      if (fac > 0.) then
			u(i,j,k) = u(i,j,k)/(1.+rkv*fac)
		      endif
		    enddo
		  enddo
		  do j=js,je
		    do i=is,ie+1
		      fac = frac(i,j)+frac(i-1,j)
		      if (fac > 0.) then
			v(i,j,k) = v(i,j,k)/(1.+rkv*fac)
		      endif
		    enddo
		  enddo

		endif

	2000  continue

#ifdef DO_AGE
       if( nq/=0 )     &
           call age_of_air(is, ie, js, je, npz, ng, time_total, pe, q(is-ng,js-ng,1,nq))
#endif


	 end subroutine Held_Suarez_Strat


 subroutine Sim_phys(npx, npy, npz, is, ie, js, je, ng, nq,  &
                     u_dt, v_dt, t_dt, q_dt, u, v, w, ua, va,   &
                     pt, delz, q, pe, delp, peln, oro, hydrostatic, &
                     phys_hydrostatic, pdt, agrid, ak, bk, p_ref, fv_sg_adj,    &
                     master, Time, time_total)


 integer, INTENT(IN) :: npx, npy, npz
 integer, INTENT(IN) :: is, ie, js, je, ng, nq
 integer, INTENT(IN) :: fv_sg_adj
 real   , INTENT(IN) :: pdt
 real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
 real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
 logical, INTENT(IN) :: master
 real, INTENT(IN):: p_ref
 real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
 logical, INTENT(IN):: hydrostatic, phys_hydrostatic

 type(time_type), intent(in) :: Time
 real, INTENT(IN), optional:: time_total

 real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
 real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
 real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)
 real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
 real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )
 real, INTENT(INOUT):: delz(is   :ie   ,js   :je   ,npz)
 real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,npz)

! Tendencies:
 real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
 real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)
 real, INTENT(INOUT):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: va(is-ng:ie+ng,js-ng:je+ng,npz)
! Local
 real, parameter:: tice = 273.16
 real, dimension(is:ie,js:je,npz):: t3, p3, dz
 real, dimension(is:ie,js:je,npz+1):: zhalf, phalf
 real, dimension(is:ie,js:je,npz  ):: zfull
 real, dimension(is:ie,js:je,npz,nq):: q3
 real, dimension(is:ie,js:je):: rain, snow, ice, graupel, prec_mp, land
 real qs(is:ie)
 real pedge(npz+1)
 real pref(npz)
 real sst, sday, shour, rkv, rka, rks, rkt, zvir, rrg, tvm
 real rk1, tmp, rkq, pc, c1, cooling, frac
 real T0, pi, rdt, convt , tot_prec
 real fac_sm
 logical used
 logical :: do_lin_cld_microphys = .true.
 integer  i,j,k, iq, nq_con, k_mp
 integer  isd, jsd
 integer  tau_sg, seconds, days
 integer  nqv, nql, nqi


   isd = is - ng
   jsd = js - ng

! Factor for Small-Earth Appprox.
!  fac_sm = RADIUS / 6371.0e3
   fac_sm = 1.
   pi = 4.*atan(1.)

   rrg  = rdgas / grav

   sday  = 24.*3600.*fac_sm
   shour = 3600.*fac_sm

   rdt = 1. / pdt
   tau_sg = fv_sg_adj

   cooling = 1.5 * pdt / sday
!  cooling = 1.0 * pdt / sday
   rkt = pdt / (6.*sday)       ! temp relaxation
   rkq = pdt / (3.*shour)      ! sphum
   rk1 = pdt / (3.*shour)      ! sensible heat (from prescribed SST)
   rkv = pdt / (6.*shour)

   t_dt = 0.;  q_dt = 0.
   u_dt = 0.;  v_dt = 0.
      
   do k=1,npz+1
      pedge(k) = ak(k) + bk(k)*p_ref
   enddo
   do k=1,npz
      pref(k) = (pedge(k+1)-pedge(k)) / log(pedge(k+1)/pedge(k))
   enddo

! Determine k_mp: top layer for microphysics to operate
   k_mp = 1
!  do k=1,npz
!     if ( pref(k) > 1.E2 ) then
!          k_mp = k
!          exit
!     endif
!  enddo

   T0 = 200.
   do k=1,npz
      do j=js,je
         do i=is,ie
            tmp = pt(i,j,k) - cooling
            if ( tmp > T0 ) then
                 t3(i,j,k) = tmp
            else
                 t3(i,j,k) = (pt(i,j,k)+rkt*T0)/(1.+rkt)
            endif
         enddo 
      enddo 
   enddo      ! k-loop


!----------
! Setup SST:
!----------
   do j=js,je
      do i=is,ie
#ifdef UNIFORM_SST
!        sst = 302.
         sst = tice + 30.
#else
! Neale & Hoskins 2001
        if ( abs(agrid(i,j,2)) <  pi/3. ) then
             tmp = sin( 1.5*agrid(i,j,2) ) ** 2

! QOBS
!            sst = 273.16 + 27.*(1.-0.5*(tmp + tmp**2))
! FLAT_SST
!            sst = 273.16 + 27.*(1.-tmp**2)
! Default control case:
             sst = 273.16 + 27.*(1.-tmp)
         else
             sst = 273.16
         endif
#endif
         t3(i,j,npz) = (t3(i,j,npz)+rk1*sst)/(1.+rk1)
       enddo
   enddo

   do iq=1,nq
      do k=1,npz
         do j=js,je
            do i=is,ie
               q3(i,j,k,iq) = q(i,j,k,iq)
            enddo
         enddo
      enddo
   enddo
         
   do k=1,npz+1
      do j=js,je
         do i=is,ie
            phalf(i,j,k) = pe(i,k,j)
         enddo
      enddo
   enddo


! Compute zfull, zhalf
   do j=js,je
      do i=is,ie
         zhalf(i,j,npz+1) = 0.
      enddo
   enddo

   zvir = rvgas/rdgas - 1.
   do k=npz,1,-1
      do j=js,je
         do i=is,ie
#ifdef PHYS_NON_HYDRO
               dz(i,j,k) =  delz(i,j,k)
               p3(i,j,k) = -rrg*delp(i,j,k)/delz(i,j,k)*t3(i,j,k)*(1.+zvir*q3(i,j,k,1))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
            zfull(i,j,k) = 0.5*(zhalf(i,j,k)+zhalf(i,j,k+1))
#else
                     tvm = rrg*t3(i,j,k)*(1.+zvir*q3(i,j,k,1))
               dz(i,j,k) = -tvm*(peln(i,k+1,j)-peln(i,k,j))
               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
            zfull(i,j,k) = zhalf(i,j,k+1) + tvm*(1.-phalf(i,j,k)/p3(i,j,k))
            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
#endif
         enddo
      enddo
   enddo

!--------------------------
! Moist Physical processes:
!--------------------------

! Latent heat fluxes:
      do j=js,je
         call qsmith(ie-is+1, 1, 1, t3(is:ie,j,npz), p3(is:ie,j,npz), q3(is:ie,j,npz,1), qs)
         do i=is,ie
            rkt = rkq * max( 0.5, sqrt(ua(i,j,npz)**2+va(i,j,npz)**2) )
!                     * min(5.E2, delp(i,j,npz) ) / delp(i,j,npz)
           q3(i,j,npz,1) = (q3(i,j,npz,1) + rkt*qs(i)) / (1.+rkt)
         enddo
      enddo   ! j-loop

   land = 0.

   if ( tau_sg > 0 ) then
         nq_con = nq
         nqv = 1
         nql = 2
         nqi = 3

         call sg_conv(is, ie, js, je,  is, ie, js, je,    &
                      is, ie, js, je, npz, nq_con, pdt, tau_sg,  &
                      delp(is:ie,js:je,1:npz), phalf,    &
                      p3, zfull, zhalf, t3, q3,    &
                      ua(is:ie,js:je,1:npz),       &
                      va(is:ie,js:je,1:npz),       &
                       w(is:ie,js:je,1:npz),       &
                      u_dt(is:ie,js:je,1:npz),     &
                      v_dt(is:ie,js:je,1:npz),     &
                      t_dt, q_dt, 1, land, oro, nqv, nql, nqi, &
                      hydrostatic, phys_hydrostatic)
       do k=1,npz
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = (u_dt(i,j,k) - ua(i,j,k))*rdt                
                v_dt(i,j,k) = (v_dt(i,j,k) - va(i,j,k))*rdt                
                  t3(i,j,k) = t_dt(i,j,k)
                t_dt(i,j,k) = (t3(i,j,k) - pt(i,j,k)) * rdt
                  dz(i,j,k) = dz(i,j,k) * t3(i,j,k)/pt(i,j,k)
             enddo
          enddo
       enddo

       do iq=1,nq
       do k=1,npz
          do j=js,je
             do i=is,ie
                  q3(i,j,k,iq) = q_dt(i,j,k,iq)
                q_dt(i,j,k,iq) = (q3(i,j,k,iq) - q(i,j,k,iq))*rdt
             enddo
          enddo
       enddo
       enddo
     else
       do k=1,npz
          do j=js,je
             do i=is,ie
                u_dt(i,j,k) = 0.
                v_dt(i,j,k) = 0.
                t_dt(i,j,k) = (t3(i,j,k)-pt(i,j,k))*rdt
             enddo
          enddo
       enddo

       do iq=1,nq
       do k=1,npz
          do j=js,je
             do i=is,ie
                q_dt(i,j,k,iq) = (q3(i,j,k,iq)-q(i,j,k,iq))*rdt
             enddo
          enddo
       enddo
       enddo
     endif


   if ( do_lin_cld_microphys ) then
!---------------------------------------
! A 6-class cloud microphysics 
!---------------------------------------
      call lin_cld_microphys_driver(q3(:,:,:,1),   q3(:,:,:,2),   q3(:,:,:,3),  &
                     q3(:,:,:,4),   q3(:,:,:,5),   q3(:,:,:,6),  &
                     q3(:,:,:,7), q_dt(:,:,:,1), q_dt(:,:,:,2), q_dt(:,:,:,3), &
                     q_dt(:,:,:,4), q_dt(:,:,:,5), q_dt(:,:,:,6), &
                     q_dt(:,:,:,7), t_dt, t3, p3, dz, delp(is:ie,js:je,1:npz),  &
                     area, pdt, land, rain, snow, ice, graupel,  &
                     hydrostatic, phys_hydrostatic,  &
                     is,ie, js,je, 1,npz, k_mp,npz, Time)

#ifdef DEBUG_MP
      prec_mp(:,:) = rain(:,:) + snow(:,:) + ice(:,:) + graupel(:,:) 
      call prt_maxmin('prec_mp', prec_mp, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('Rain', rain, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('Snow', snow, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('ice',   ice, is, ie, js, je, 0, 1, 1., master)
      call prt_maxmin('Graupel', graupel, is, ie, js, je, 0, 1, 1., master)
#endif
    endif

!-------------------------
! one-layer (surface) drag 
!-------------------------
   k=npz
   do j=js,je
      do i=is,ie
                frac = rkv*max(0.5, sqrt(ua(i,j,npz)**2+va(i,j,npz)**2))
         u_dt(i,j,k) = u_dt(i,j,k) - ua(i,j,k)*frac/(1.+frac) * rdt
         v_dt(i,j,k) = v_dt(i,j,k) - va(i,j,k)*frac/(1.+frac) * rdt
      enddo
   enddo

 end subroutine Sim_phys

 subroutine age_of_air(is, ie, js, je, km, ng, time, pe, q)

      integer is, ie, js, je
      integer km
      integer ng

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

      real, intent(inout):: pe(is-1:ie+1, km+1, js-1:je+1)
      real, intent(in):: time        ! accumulated time since init
      real, intent(inout):: q(is-ng:ie+ng,js-ng:je+ng,km)

! Local
      integer i, j, k
      real p_source      ! source level (pa)
      real ascale
      real tiny
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( ascale = 5.e-6 / 60. )

!$omp parallel do default(shared) private(i, j, k)
      do k=1,km
        do j=js,je
            do i=is,ie
               if( time < tiny ) then
                   q(i,j,k) = 0.
               elseif( pe(i,k,j) >= p_source ) then
                   q(i,j,k) = ascale * time
               endif
            enddo
        enddo          ! j-loop
      enddo             ! k-loop

      end subroutine age_of_air

end module hswf_mod

c ======================================================================
c User Subroutine VUMAT for Johnson-Cook model.
c All rights of reproduction or distribution in any form are reserved.
c By Irfan Habeeb CN (Technion - IIT), cnirfan@gmail.com
c ======================================================================
      subroutine vumat(
C Read only -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew)
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname
      integer k, k1, k2, iter
      real*8 E, nu, A, B, C, m, n, rate0, Tr, Tm, rho, TQ, Cp, WH,
     1 sigeqv, f, epbar, ep, dep, eprateN, sighyd, sigyT, sigyN, dPwork,
     2 trInc, fT, fs, fN, fNs, sigdotp, epbarN, tempT, tempN, dWork,
     3 mu, alamda, tol,
     4 cmat(6,6), sTrial(6), sNew(6), devS(6), sTem(6), rd(6),
     5 L(6,6), sOld(6), sT(6), np(6)

c input parameters
      E = props(1)				! Young's modulus
      nu = props(2)				! Poisson's ratio
      A = props(3)				! JC parameters
      B = props(4)
      C = props(5)
      n = props(6)
      m = props(7)
      rate0 = props(8)
      Tr = props(9)				! reference temperature
      Tm = props(10)			! melting temperature
      rho = props(11)			! density
      TQ = props(12)			! Taylor-Q heat ener. from plastic deform.
      Cp = props(13)			! specific heat
      WH = props(14)			! work - heat conversion factor

c tolerance of the effe. stress, strain accuracy: ep ~ 1e-9
      tol = E*1e-9 									! strain accuracy or 1e-9
      Niter = 50                  	! max number of N-R iterations

c lame's parameters
      mu = E/(2.d0*(1.d0+nu))
      alamda = E*nu/((1.d0 + nu) * (1.d0 - 2.d0*nu))

c stiffness matrix
			cmat = 0.d0
			do k1 = 1, ndir
				do k2 = 1, ndir
					cmat(k1, k2) = alamda
				end do 
				cmat(k1, k1) = alamda + 2.d0*mu
			end do 
			do k1 = ndir+1, ndir+nshr
				cmat(k1, k1) = 2.d0*mu
			end do 

C -------------------- simulation first step & later -------------------
      do 30 k = 1, nblock
      if (stateOld(k, 1) .eq. 0.d0) then
        go to 10
      else 
        go to 20
      end if 
C----------------------------- initial state ---------------------------
 10   trInc = sum(strainInc(k, 1:3))
			do k1 = 1, ndir
				stressNew(k, k1) = stressOld(k, k1) + alamda*trInc + 
     1		2.d0*mu*strainInc(k, k1)
		 	end do 
		 	do k1 = ndir+1, ndir+nshr
		 		stressNew(k, k1) = stressOld(k, k1) + 
     1		2.d0*mu*strainInc(k, k1)
     	end do 

      stateNew(k, 1) = 1.d0           ! initiation check
      stateNew(k, 2) = 0.d0           ! plastic strain
      stateNew(k, 3) = Tr           	! initial temperature
      stateNew(k, 4) = 0.d0           ! yield stress
      stateNew(k, 5) = 0.d0						! plastic strain incre. last step
      stateNew(k, 6) = 1							! number of iterations
      tempNew(k) = Tr
      go to 30
C-------------------------- 2nd step and later -------------------------
 20   epbar = stateOld(k, 2)
      tempT = stateOld(k, 3)
      sigyT = stateOld(k, 4)
      sOld(1:6) = stressOld(k, 1:6)
      trInc = sum(strainInc(k, 1:3))

      do k1 = 1, ndir
        sT(k1) = stressOld(k,k1) +alamda*trInc +2.d0*mu*strainInc(k, k1)
      end do 
      do k1 = ndir+1, ndir+nshr
        sT(k1) = stressOld(k, k1) + 2.d0*mu*strainInc(k, k1)
      end do 

c to get the vector for radial reduction
      sighyd = sum(sT(1:ndir))/3.d0
      devS(1:ndir) = sT(1:ndir) - sighyd
      devS(ndir+1:ndir+nshr) = sT(ndir+1:ndir+nshr)

c equivalent stress
      sigeqv = sqrt(3.d0/2.d0 *(devS(1)**2.d0 + devS(2)**2.d0 +
     1 devS(3)**2.d0 + 2.d0*devS(4)**2.d0 + 2.d0*devS(5)**2.d0 +
     2 2.d0*devS(6)**2.d0 ))

c rd is constant over the iterations
			if (sigeqv .gt. 0.d0) then
      	rd = (3.d0*devS)/(2.d0*sigeqv)
      else
      	rd = 0.d0
      end if
			sTem = matmul(cmat, rd)

c initial value of the plastic strain increment and the bounds
      ep = 0.d0											! initial plast. strain incre.
      epmin = 0.d0
      epmax = sigeqv/(2.d0*mu) !max(1.d-5, 3.d0*1.d-5)

c NR - iteration starts
      iter = 0											! number of iterations
      do while (iter .lt. Niter) 
        iter = iter + 1
        if (iter .gt. Niter-1.) then
          print*, 'too many iterations, iter = ', iter
          call XPLB_EXIT
        end if
        
        epbarN = epbar + ep
        eprateN = ep/dt
        tempN = tempT 						! temperature is updated at the end

c updating stress
        sNew = sT - ep*sTem
        sighyd = sum(sNew(1:ndir))/3.d0

c deviatoric stress
        devS(1:ndir) = sNew(1:ndir) - sighyd
        devS(ndir+1:ndir+nshr) = sNew(ndir+1:ndir+nshr)

c equivalent stress
        sigeqv = sqrt(3.d0/2.d0 *(devS(1)**2.d0 + devS(2)**2.d0 +
     1  devS(3)**2.d0 + 2.d0*devS(4)**2.d0 + 2.d0*devS(5)**2.d0 +
     2  2.d0*devS(6)**2.d0 ))

c evaluating new yield stress for the increment 'ep'
	 			call func_syield(A, B, epbarN, n, Tr, Tm, tempN, m, C, eprateN,
     1			rate0, sigyN)
        if (sigyN .lt. sigyT) sigyN = sigyT

        f = sigeqv - sigyN
        if (abs(f) .lt. tol) exit 			! out of the iteration

c elastic criteria, ep = 0 & f < 0
        if ((ep .eq. 0.d0) .and. (f .lt. 0.d0)) go to 12

c rearranging the margins and new plastic strain increm.
        if ((f .ge. 0.d0) .and. (ep .ge. epmin)) epmin = ep
        if ((f .lt. 0.d0) .and. (ep .lt. epmax)) epmax = ep
        ep = 0.5d0 * (epmax + epmin)

c restoring the plastic strain increment from the last step
        if (iter .eq. 1) ep = stateOld(k, 5)
      end do 

 12   stressNew(k, 1:6) = sNew(1:6)

c work, plastic work and temp
      dWork = dot_product( 0.5d0*(sOld(1:6) + sNew(1:6)),
     1  strainInc(k, 1:6))
      dPwork = 0.5d0 * ep * sigeqv
      enerInternNew(k) = enerInternOld(k) + dWork/rho
      enerInelasNew(k) = enerInelasOld(k) + dPwork/rho

c updating state variables
 			stateNew(k, 1) = 1.d0							! to flag the initial step
			stateNew(k, 2) = epbarN						! plastic strain
			stateNew(k, 3) = tempT + WH*TQ*dPwork/rho/Cp	! temperature
			stateNew(k, 4) = sigyN						! yield stress
			stateNew(k, 5) = ep 							! plastic strain incre.
			stateNew(k, 6) = stateOld(k, 6) + iter ! number of iterations

 30   continue
      return
      end

c ----------------------------------------------------------------------
      subroutine func_syield(A, B, epbar, n, Tr, Tm, T, m, C, rate,
     1 rate0, sigbar)
		  include 'vaba_param.inc'
		  real*8 A, B, epbar, n, Tr, Tm, T, theta, m, C, rate, rate0, sigbar

		  if (T < Tr) then
		  	theta = 0.d0
		  else if (T > Tm) then
		  	theta = 1.d0
		  else 
		  	theta = (T - Tr)/(Tm - Tr)
		  end if 

		  if (rate == 0.d0) then
		  	rate = rate0
		  end if

		  sigbar = (A + B*epbar**n)*(1.d0 - theta**m)*
     1 (1.d0 + C*log(rate/rate0)) 
		  end subroutine
c ----------------------------------------------------------------------
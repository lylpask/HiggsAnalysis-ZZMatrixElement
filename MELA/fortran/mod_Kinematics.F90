MODULE ModKinematics
implicit none
save

public :: SetRunningScales,setPDFs,EvalAlphaS

CONTAINS



SUBROUTINE boost2Lab(x1,x2,NumPart,Mom)
implicit none
real(8) Mom(1:4,1:NumPart)
real(8) x1,x2
real(8) gamma,betagamma,MomTmp1,MomTmp4
integer :: i,NumPart

  gamma     = (x1+x2)/2d0/dsqrt(x1*x2)
  betagamma = (x2-x1)/2d0/dsqrt(x1*x2)

  do i=1,NumPart
      MomTmp1=Mom(1,i)
      MomTmp4=Mom(4,i)
      Mom(1,i)= gamma*MomTmp1 - betagamma*MomTmp4
      Mom(4,i)= gamma*MomTmp4 - betagamma*MomTmp1
  enddo

RETURN
END SUBROUTINE




!LORENTZ.F
!VERSION 20130123
!
!A subroutine that performs a general boost to a four vector
!(vector) based on another four vector (boost). The primed and
!unprimed frames have their axes in parallel to one another.
!Rotation is not performed by this subroutine.

      subroutine LORENTZ(vector, boost)

      implicit none

      double precision vector(4), boost(4) 
      double precision lambdaMtrx(4,4), vector_copy(4)
      double precision beta(2:4), beta_sq, gamma
      integer i,j
      double precision, parameter :: epsilon = 1d-13 !a small quantity slightly above machine precision

!      double precision KRONECKER_DELTA
!      external KRONECKER_DELTA

      do i=2,4
        beta(i) = boost(i)/boost(1)
      enddo

      beta_sq = beta(2)**2+beta(3)**2+beta(4)**2

  if(beta_sq.ge.epsilon)then

      gamma = 1d0/dsqrt(1d0-beta_sq)

      lambdaMtrx(1,1) = gamma

      do i=2,4
        lambdaMtrx(1,i) = gamma*beta(i)
        lambdaMtrx(i,1) = lambdaMtrx(1,i)
      enddo

      do i=2,4
      do j=2,4
        lambdaMtrx(i,j) = (gamma-1d0)*beta(i)*beta(j)/beta_sq + KRONECKER_DELTA(i,j)
      enddo
      enddo


!apply boost to vector1
      vector_copy = vector
      vector = 0d0
      do i=1,4
      do j=1,4
        vector(i) = vector(i) + lambdaMtrx(i,j)*vector_copy(j)
      enddo
      enddo
  endif

      return
      END subroutine LORENTZ



! KRONECKER_DELTA.F
!
! KRONECKER_DELTA(i,j)
! A function that returns 1 if i=j, and 0 otherwise.
      double precision function KRONECKER_DELTA(i,j)
      integer i,j
      if(i.eq.j)then
        KRONECKER_DELTA = 1d0
      else
        KRONECKER_DELTA = 0d0
      endif

      return
      end function KRONECKER_DELTA





subroutine SetRunningScales(p,id) ! p in JHU-GeV, id in JHUGen conventions
use ModParameters
use ModMisc
implicit none
real(dp), intent(in) :: p(1:4,4:6) ! No need to run the second index from 3 to 7: pH, pJ1, pJ2
integer, intent(in) :: id(4:7) ! id_JJH/id_JJVV, id_J1, id_J2, id_JJ (if applicable)  
real(8) :: polemass(3:7) ! mJJH, mH, mJ1, mJ2, mJJ (if applicable)
real(8) :: pJJHstar(4),pHstar(4),pJ(4,2),pJJ(4),pJHstar(4)
integer idx,ip

   pHstar(:) = 0d0
   pJJ(:) = 0d0
   polemass(3) = getMass(id(4)) ! Pole mass of the JJH system
   polemass(4) = M_Reso
   do idx=4,6
      if(idx.eq.4) then
         do ip=1,4
            pHstar(ip) = pHstar(ip) + p(ip,idx)
         enddo
      else
         polemass(idx) = getMass(id(idx))
         do ip=1,4
            pJJ(ip) = pJJ(ip) + p(ip,idx)
         enddo
      endif
   enddo
   polemass(7) = getMass(id(7)) ! Pole mass of the JJ system

   pJJHstar = pJJ + pHstar
   if(polemass(5).lt.polemass(6)) then
      pJ(:,1)=p(:,5)
      pJ(:,2)=p(:,6)
   else
      pJ(:,1)=p(:,6)
      pJ(:,2)=p(:,5)
      call swapr(polemass(5),polemass(6)) ! will use polemass(5) as the greater mass below
   endif
   pJHstar(1:4) = pJ(1:4,1) + pHstar(1:4)

   ! Determine the appropriate factorization scale for the chosen scheme from pole and invariant masses
   if(FacScheme .eq. kRenFacScheme_mhstar) then
      Mu_Fact = Get_MInv(pHstar(1:4))

   elseif(FacScheme .eq. -kRenFacScheme_mhstar) then
      Mu_Fact = polemass(4)

   elseif(FacScheme .eq. kRenFacScheme_mjjhstar) then
      Mu_Fact = Get_MInv(pJJHstar(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mjjhstar) then
      Mu_Fact = polemass(3)

   elseif(FacScheme .eq. kRenFacScheme_mjj_mhstar) then
      Mu_Fact = Get_MInv(pJJ(1:4))+Get_MInv(pHstar(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mjj_mhstar) then
      Mu_Fact = polemass(4)+polemass(7)
   elseif(FacScheme .eq. kRenFacScheme_mj_mj_mhstar) then
      Mu_Fact = Get_MInv(pJ(1:4,1))+Get_MInv(pJ(1:4,2))+Get_MInv(pHstar(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mj_mj_mhstar) then
      Mu_Fact = polemass(4)+polemass(5)+polemass(6)

   elseif(FacScheme .eq. kRenFacScheme_mjj) then
      Mu_Fact = Get_MInv(pJJ(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mjj) then
      Mu_Fact = polemass(7)
   elseif(FacScheme .eq. kRenFacScheme_mj_mj) then
      Mu_Fact = Get_MInv(pJJ(1:4))
   elseif(FacScheme .eq. -kRenFacScheme_mj_mj) then
      Mu_Fact = polemass(5)+polemass(6)

   elseif(FacScheme .eq. kRenFacScheme_mjhstar) then
      Mu_Fact = Get_MInv(pJHstar(1:4))
   elseif(FacScheme .eq. kRenFacScheme_mj_mhstar) then
      Mu_Fact = Get_MInv(pJ(1:4,1))+Get_MInv(pHstar(1:4))
   elseif((FacScheme .eq. -kRenFacScheme_mjhstar) .or. (FacScheme .eq. -kRenFacScheme_mj_mhstar)) then
      Mu_Fact = polemass(4)+polemass(5)
   elseif(FacScheme .eq. kRenFacScheme_mj) then
      Mu_Fact = Get_MInv(pJ(1:4,1))
   elseif(FacScheme .eq. -kRenFacScheme_mj) then
      Mu_Fact = polemass(5)
   endif

   ! Do the same for the renormalization scale
   if(RenScheme .eq. kRenFacScheme_mhstar) then
      Mu_Ren = Get_MInv(pHstar(1:4))

   elseif(RenScheme .eq. -kRenFacScheme_mhstar) then
      Mu_Ren = polemass(4)

   elseif(RenScheme .eq. kRenFacScheme_mjjhstar) then
      Mu_Ren = Get_MInv(pJJHstar(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mjjhstar) then
      Mu_Ren = polemass(3)

   elseif(RenScheme .eq. kRenFacScheme_mjj_mhstar) then
      Mu_Ren = Get_MInv(pJJ(1:4))+Get_MInv(pHstar(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mjj_mhstar) then
      Mu_Ren = polemass(4)+polemass(7)
   elseif(RenScheme .eq. kRenFacScheme_mj_mj_mhstar) then
      Mu_Ren = Get_MInv(pJ(1:4,1))+Get_MInv(pJ(1:4,2))+Get_MInv(pHstar(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mj_mj_mhstar) then
      Mu_Ren = polemass(4)+polemass(5)+polemass(6)

   elseif(RenScheme .eq. kRenFacScheme_mjj) then
      Mu_Ren = Get_MInv(pJJ(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mjj) then
      Mu_Ren = polemass(7)
   elseif(RenScheme .eq. kRenFacScheme_mj_mj) then
      Mu_Ren = Get_MInv(pJJ(1:4))
   elseif(RenScheme .eq. -kRenFacScheme_mj_mj) then
      Mu_Ren = polemass(5)+polemass(6)

   elseif(RenScheme .eq. kRenFacScheme_mjhstar) then
      Mu_Ren = Get_MInv(pJHstar(1:4))
   elseif(RenScheme .eq. kRenFacScheme_mj_mhstar) then
      Mu_Ren = Get_MInv(pJ(1:4,1))+Get_MInv(pHstar(1:4))
   elseif((RenScheme .eq. -kRenFacScheme_mjhstar) .or. (RenScheme .eq. -kRenFacScheme_mj_mhstar)) then
      Mu_Ren = polemass(4)+polemass(5)
   elseif(RenScheme .eq. kRenFacScheme_mj) then
      Mu_Ren = Get_MInv(pJ(1:4,1))
   elseif(RenScheme .eq. -kRenFacScheme_mj) then
      Mu_Ren = polemass(5)
   endif

   ! Never ever allow the scales to go negative
   Mu_Fact = abs(Mu_Fact) * MuFacMultiplier
   Mu_Ren = abs(Mu_Ren) * MuRenMultiplier

return
end subroutine SetRunningScales




SUBROUTINE setPDFs(x1,x2,pdf)
use ModParameters
implicit none
real(8) :: x1,x2,PDFScale
real(8) :: upv(1:2),dnv(1:2),usea(1:2),dsea(1:2),str(1:2),chm(1:2),bot(1:2),glu(1:2),phot(1:2),sbar(1:2),cbar(1:2),bbar(1:2)
integer,parameter :: swPDF_u=1, swPDF_d=1, swPDF_c=1, swPDF_s=1, swPDF_b=1, swPDF_g=1
real(8) :: pdf(-6:6,1:2),NNpdf(1:2,-6:7)

        PDFScale=Mu_Fact*100d0
        pdf(:,:) = 0d0
        
#if useLHAPDF==1
        call evolvePDF(x1,PDFScale,NNpdf(1,-6:7))
        call evolvePDF(x2,PDFScale,NNpdf(2,-6:7))
            NNpdf(1,-6:7) = NNpdf(1,-6:7)/x1
            NNpdf(2,-6:7) = NNpdf(2,-6:7)/x2
            
            pdf(Up_,1)   = NNpdf(1,+2)         * swPDF_u
            pdf(AUp_,1)  = NNpdf(1,-2)         * swPDF_u
            pdf(Dn_,1)   = NNpdf(1,+1)         * swPDF_d
            pdf(ADn_,1)  = NNpdf(1,-1)         * swPDF_d
            pdf(Chm_,1)  = NNpdf(1,+4)         * swPDF_c
            pdf(AChm_,1) = NNpdf(1,-4)         * swPDF_c
            pdf(Str_,1)  = NNpdf(1,+3)         * swPDF_s
            pdf(AStr_,1) = NNpdf(1,-3)         * swPDF_s
            pdf(Bot_,1)  = NNpdf(1,+5)         * swPDF_b
            pdf(ABot_,1) = NNpdf(1,-5)         * swPDF_b
            pdf(0,1)     = NNpdf(1,+0)         * swPDF_g            
            
            pdf(Up_,2)   = NNpdf(2,+2)         * swPDF_u
            pdf(AUp_,2)  = NNpdf(2,-2)         * swPDF_u
            pdf(Dn_,2)   = NNpdf(2,+1)         * swPDF_d
            pdf(ADn_,2)  = NNpdf(2,-1)         * swPDF_d
            pdf(Chm_,2)  = NNpdf(2,+4)         * swPDF_c
            pdf(AChm_,2) = NNpdf(2,-4)         * swPDF_c
            pdf(Str_,2)  = NNpdf(2,+3)         * swPDF_s
            pdf(AStr_,2) = NNpdf(2,-3)         * swPDF_s
            pdf(Bot_,2)  = NNpdf(2,+5)         * swPDF_b
            pdf(ABot_,2) = NNpdf(2,-5)         * swPDF_b
            pdf(0,2)     = NNpdf(2,+0)         * swPDF_g
            
#else

        if( PDFSet.eq.3 ) then

            call NNevolvePDF(x1,PDFScale,NNpdf(1,-6:7))
            call NNevolvePDF(x2,PDFScale,NNpdf(2,-6:7))
            NNpdf(1,-6:7) = NNpdf(1,-6:7)/x1
            NNpdf(2,-6:7) = NNpdf(2,-6:7)/x2
            
    !       PROTON CONTENT
            pdf(Up_,1)   = NNpdf(1,+2)         * swPDF_u
            pdf(AUp_,1)  = NNpdf(1,-2)         * swPDF_u
            pdf(Dn_,1)   = NNpdf(1,+1)         * swPDF_d
            pdf(ADn_,1)  = NNpdf(1,-1)         * swPDF_d
            pdf(Chm_,1)  = NNpdf(1,+4)         * swPDF_c
            pdf(AChm_,1) = NNpdf(1,-4)         * swPDF_c
            pdf(Str_,1)  = NNpdf(1,+3)         * swPDF_s
            pdf(AStr_,1) = NNpdf(1,-3)         * swPDF_s
            pdf(Bot_,1)  = NNpdf(1,+5)         * swPDF_b
            pdf(ABot_,1) = NNpdf(1,-5)         * swPDF_b
            pdf(0,1)     = NNpdf(1,+0)         * swPDF_g            
            
            pdf(Up_,2)   = NNpdf(2,+2)         * swPDF_u
            pdf(AUp_,2)  = NNpdf(2,-2)         * swPDF_u
            pdf(Dn_,2)   = NNpdf(2,+1)         * swPDF_d
            pdf(ADn_,2)  = NNpdf(2,-1)         * swPDF_d
            pdf(Chm_,2)  = NNpdf(2,+4)         * swPDF_c
            pdf(AChm_,2) = NNpdf(2,-4)         * swPDF_c
            pdf(Str_,2)  = NNpdf(2,+3)         * swPDF_s
            pdf(AStr_,2) = NNpdf(2,-3)         * swPDF_s
            pdf(Bot_,2)  = NNpdf(2,+5)         * swPDF_b
            pdf(ABot_,2) = NNpdf(2,-5)         * swPDF_b
            pdf(0,2)     = NNpdf(2,+0)         * swPDF_g            

        else
            print *, "PDFSet",PDFSet,"not available!"
            stop
        endif
#endif

        pdf(:,:) = dabs(pdf(:,:))
        RETURN

END SUBROUTINE



! QCD scale from MCFM
! Implementation into JHUGen by Ulascan Sarica, Dec. 2015
subroutine EvalAlphaS()
   use ModParameters
   IMPLICIT NONE
#if useLHAPDF==1
!--- This is simply a wrapper to the LHAPDF implementation of the running coupling alphas, in the style of the native MCFM routine
   DOUBLE PRECISION alphasPDF
   REAL(DP) :: Q
      Q = Mu_Ren/GeV
      alphas=alphasPDF(Q)
#else
!     Evaluation of strong coupling constant alphas
!     Original Author: R.K. Ellis
!     q -- Scale at which alpha_s is to be evaluated
!     alphas_mz -- ModParameters value of alpha_s at the mass of the Z-boson
!     nloops_pdf -- ModParameters value of the number of loops (1,2, or 3) at which the beta function is evaluated to determine running.
!     If you somehow need a more complete implementation, check everything at or before commit 28472c5bfee128dde458fd4929b4d3ece9519ab8
   INTEGER, PARAMETER :: NF6=6
   INTEGER, PARAMETER :: NF5=5
   INTEGER, PARAMETER :: NF4=4
   INTEGER, PARAMETER :: NF3=3
   INTEGER, PARAMETER :: NF2=2
   INTEGER, PARAMETER :: NF1=1
   
      IF (Mu_Ren .LE. 0d0) THEN 
         WRITE(6,*) 'ModKinematics::EvalAlphaS: Mu_Ren .le. 0, Mu_Ren (GeV) = ',(Mu_Ren*GeV)
         stop
      ENDIF
      IF (nQflavors_pdf .NE. NF5) THEN 
         WRITE(6,*) 'ModKinematics::EvalAlphaS: nQflavors_pdf invalid, nQflavors_pdf = ',nQflavors_pdf
         WRITE(6,*) 'ModKinematics::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
         stop
      ENDIF
      IF (nloops_pdf .NE. 1) THEN 
         WRITE(6,*) 'ModKinematics::EvalAlphaS: nloops_pdf invalid, nloops_pdf = ',nloops_pdf
         WRITE(6,*) 'ModKinematics::EvalAlphaS: Check 28472c5bfee128dde458fd4929b4d3ece9519ab8'
         stop
      ENDIF

      alphas=alphas_mz/(1_dp+alphas_mz*B0_PDF(NF5)*2_dp*log((Mu_Ren/zmass_pdf)))
#endif
      ! Temporary arrangement
      alphas_mz = 0.13229060d0
      alphas = 0.13229060d0
      ! Calculate the derived couplings
      call ComputeQCDVariables()
   RETURN
end subroutine EvalAlphaS




END MODULE

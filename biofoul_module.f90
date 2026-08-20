	MODULE BF_MOD

!  The Biofoul Module contains two procedures related to microplastic
!  biofouling:
!  Subroutine biofoul_subr computes the degree of biofouling.
!  Subroutine settling_vel calculates the particle settling velocity based on
!  plastic size and density and the degree of biofouling 
!
!  Created by:            Laura Sunberg
!  Created on:            2023
!  Last Modified on:      12 May 2026
! ROMSPath Version: 1.0.1

IMPLICIT NONE
PUBLIC

CONTAINS

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                     FUNCTION biofoul_subr                     ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE biofoul_subr(live_biofoul,dead_biofoul,total_biofoul,r_pl,settling_vel,dt,ng,Xpar,Ypar,Zpar,ix,ex,light,encounter,growth,grazing,mort,remineralization)
    USE INT_MOD,    ONLY: getInterp2D,getInterp3D,polintd,sinintd
    USE PARAM_MOD,  ONLY: adherence,v_A,cell_N,growth_factor,max_grazing_perday, & 
                         mortality_perday,remin_rate_perday,rho_bf,rho_pl, &
                         grazing_method,kp,K_NH4,K_NO3,Vp0,PhyIS,t_b,t_c,t_f
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: ng,dt
    DOUBLE PRECISION, INTENT(IN) :: r_pl,settling_vel,Xpar,Ypar,Zpar,light,ex(3),ix(3)
    DOUBLE PRECISION             :: PI, v_pl, SA_pl, r_A,phytopl,A_A, &
                                    v_bf,v_tot,r_tot,u_top,u_bot,v_top,v_bot,z_top,z_bot,z_dist,vel_dif,shear, & 
                                    temp,NH4,NO3,Vp,fac1,Epp,t_PPmax,cff1,cff2,inhNH4,N_Flux_NewProd,N_Flux_RegProd,PP,&
                                    beta_settling,beta_shear, beta_A,max_grazing, &
                                    zoopl,mp,rd,phytopl_local,depth,zetab,zetac,zetaf,ey(3)
                                    !beta_shear,encounter,mp,mortality,Q10,R20   
    DOUBLE PRECISION, INTENT(INOUT):: live_biofoul,dead_biofoul
    DOUBLE PRECISION, INTENT(INOUT) :: total_biofoul
    DOUBLE PRECISION, INTENT(OUT) :: encounter, growth, grazing, mort,remineralization
       
        !Load other variables 
        depth = DBLE(-1.0)* getInterp2D("depth",ng,Xpar,Ypar,t_c)
	zetab =  getInterp2D("zeta",ng,Xpar,Ypar,t_b)
	zetac =  getInterp2D("zeta",ng,Xpar,Ypar,t_c)
	zetaf =  getInterp2D("zeta",ng,Xpar,Ypar,t_f)
			  
	ey(1)=getInterp3d("phytoplankton",ng,Xpar,Ypar,Zpar,t_b,1,zetab,depth)
	ey(2)=getInterp3d("phytoplankton",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
	ey(3)=getInterp3d("phytoplankton",ng,Xpar,Ypar,Zpar,t_f,1,zetaf,depth)
	phytopl=polintd(ex,ey,3,ix(2))

	ey(1)=getInterp3d("temp",ng,Xpar,Ypar,Zpar,t_b,1,zetab,depth)
	ey(2)=getInterp3d("temp",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
	ey(3)=getInterp3d("temp",ng,Xpar,Ypar,Zpar,t_f,1,zetaf,depth)
	temp=polintd(ex,ey,3,ix(2))

	ey(1)=getInterp3d("NH4",ng,Xpar,Ypar,Zpar,t_b,1,zetab,depth)
	ey(2)=getInterp3d("NH4",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
	ey(3)=getInterp3d("NH4",ng,Xpar,Ypar,Zpar,t_f,1,zetaf,depth)
	NH4=polintd(ex,ey,3,ix(2))

	ey(1)=getInterp3d("NO3",ng,Xpar,Ypar,Zpar,t_b,1,zetab,depth)
	ey(2)=getInterp3d("NO3",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
	ey(3)=getInterp3d("NO3",ng,Xpar,Ypar,Zpar,t_f,1,zetaf,depth)
	NO3=polintd(ex,ey,3,ix(2))

	ey(1)=getInterp3d("zooplankton",ng,Xpar,Ypar,Zpar,t_b,1,zetab,depth)
	ey(2)=getInterp3d("zooplankton",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
	ey(3)=getInterp3d("zooplankton",ng,Xpar,Ypar,Zpar,t_f,1,zetaf,depth)
	zoopl=polintd(ex,ey,3,ix(2))

        ! encounters
        PI = 4*ATAN(1.d0)
        v_pl = (4./3.)*PI*r_pl**3
        SA_pl = 4.*PI*r_pl**2                                     ! plastic surface area
        r_A = (3.*v_A/(4.*PI))**(1./3.)                            ! radius algae
        A_A = phytopl * cell_N                                ! Ambient algae  convert mmol N/m^3 to No. Algae cells/m^3 with values from Marchetti & Harrison, 2007
        v_bf = total_biofoul * v_A*SA_pl                            ! volume biofouling. m^3
        v_tot = v_bf+v_pl                                       ! total plastic + biofouled volume. m^3
        r_tot = (3.*v_tot/(4.*PI))**(1./3.)                          ! total radius. m.
        shear = 1.D0 !vel_dif/ABS(z_dist)                             ! get shear!
        beta_shear = 1.3*shear*(r_tot+r_A)**3                      ! shear encounters
        if (settling_vel < 0.) then
                beta_settling = 0.5*PI*r_tot**2*ABS(settling_vel)        ! settling encounters
        else
                beta_settling = 0.
        end if
        beta_A = adherence*(beta_settling +beta_shear)                       ! encounter kernel
        encounter = beta_A*A_A/SA_pl                            ! change due to encounters
       
        !growth
        !temp = getInterp3D("temp",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
        !NH4 = getInterp3D("NH4",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)
        !NO3 = getInterp3D("NO3",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth) 
        Vp = Vp0*1.066**temp 
        fac1 = light*PhyIS
        Epp = Vp/SQRT(Vp**2+fac1**2)
        t_PPmax = Epp*fac1
        cff1 = NH4*K_NH4
        cff2 = NO3*K_NO3
        inhNH4 = 1/(1+cff1)
        N_Flux_NewProd = NO3*K_NO3*inhNH4/(1+cff2)*t_PPmax
        N_Flux_RegProd = NH4*K_NH4/(1+cff1)*t_PPmax
        PP = N_Flux_NewProd+N_Flux_RegProd
        growth = growth_factor*PP/24./60./60.*live_biofoul

        !grazing
        max_grazing = max_grazing_perday/24./60./60.
        !zoopl = getInterp3D("zooplankton",ng,Xpar,Ypar,Zpar,t_c,1,zetac,depth)     !get local zoplankton concentration
        SELECT CASE (grazing_method)
        CASE (1) !normalize by local phytoplankton to convert units
                grazing = max_grazing*zoopl/phytopl*live_biofoul 
        CASE (2) !use a constant value everywhere
                grazing = max_grazing*live_biofoul
        CASE (3) !normalize with plastic size
                grazing = max_grazing*zoopl*v_tot*cell_N/SA_pl 
        CASE (4) !normalize by local phytoplankton concentration, include [Phyt] factor
                grazing = max_grazing*(phytopl**2/(kp+phytopl**2))*zoopl/phytopl*live_biofoul 
        CASE (5) !normalize with plastic size, include [Phyt] factor
                phytopl_local = live_biofoul*SA_pl/(cell_N*V_tot)
                grazing = max_grazing*(phytopl_local**2/(kp+phytopl_local**2))*zoopl*v_tot*cell_N/SA_pl 
        CASE DEFAULT
                WRITE(*,*) "error in choosing grazing method"
        END SELECT

        !mortality
        mp = mortality_perday/24./60./60.
        mort = mp*live_biofoul


        live_biofoul = live_biofoul+dble(dt)*(encounter+growth-grazing-mort)                     ! update biofouling. no. algae cells/m^3

        if (live_biofoul < 0.) then
                live_biofoul = 0.D0
        end if

        !update dead biofouling

        !remineralization
        rd = remin_rate_perday/24./60./60.
        remineralization = rd*dead_biofoul
        
        dead_biofoul = dead_biofoul+dble(dt)*(mort-remineralization)
        
        if (dead_biofoul < 0.) then
                dead_biofoul = 0.D0
        end if

        total_biofoul = live_biofoul+dead_biofoul

        !WRITE(*,*) "remineralization: ", remineralization
        !WRITE(*,*) "live biofoul: ", live_biofoul
        !WRITE(*,*) "dead biofoul: ", dead_biofoul
        !WRITE(*,*) "total biofoul: ", total_biofoul

  END SUBROUTINE
  
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                     FUNCTION settling_vel_func                ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DOUBLE PRECISION FUNCTION settling_vel_func(r_pl,bf,ng,Ipar,Jpar,Zpar,t_c)
    USE INT_MOD,    ONLY: getInterp3D,getInterp2D
    USE PARAM_MOD,  ONLY: rho_bf, rho_pl, V_A

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: ng,t_c
    DOUBLE PRECISION, INTENT(IN) :: r_pl,bf,Ipar,Jpar,Zpar
    DOUBLE PRECISION             :: t_bf,rho_tot,rho_sw,g,wstar,nu_sw, &
                                    Dstar,V_bf,V_pl,V_tot,SA_pl,PI,r_tot,rho0, &
                                    rho_anom,test1,inner,temp,salt,mu_w,A,B,&
                                    mu_sw,depth,zeta
        
        depth = DBLE(-1.0)* getInterp2D("depth",ng,Ipar,Jpar,t_c) !new edit
	zeta =  getInterp2D("zeta",ng,Ipar,Jpar,t_c) !this is zetac!

        !convert biofouling cells/m^2 to biofouled layer thickness
        PI = 4*ATAN(1.d0)
        SA_pl = 4.*PI*r_pl**2
        V_bf = bf*V_A*SA_pl
        V_pl = (4./3.)*PI*r_pl**3
        V_tot = V_pl+V_bf
        t_bf = (V_tot*3./(4.*PI))**(1./3.)-r_pl

        !get total density
        r_tot = t_bf+r_pl
        rho_tot = (r_pl**3 * rho_pl+(r_tot**3-r_pl**3)*rho_bf)/r_tot**3

        !seawater density
        rho0 = 1000.
        rho_anom =  getInterp3D("rho",ng,Ipar,Jpar,Zpar,t_c,1,zeta,depth)           
        rho_sw = rho0+rho_anom

        !seawater viscosity (new!)
        temp =  getInterp3D("temp",ng,Ipar,Jpar,Zpar,t_c,1,zeta,depth)                   
        salt =  getInterp3D("salt",ng,Ipar,Jpar,Zpar,t_c,1,zeta,depth) !this may need conversion(?)
        salt = salt/1000.
        mu_w = 4.2844e-5 + 1./(0.157*(temp+64.993)**2 - 91.296)
        A = 1.541 + 1.998e-2*temp - 9.52e-5*temp**2
        B = 7.974 - 7.561e-2*temp + 4.724e-4*temp**2
        mu_sw = mu_w*(1.+A*salt+B*salt**2)
        nu_sw = mu_sw/rho_sw   

        !compute settling velocity
        r_tot = t_bf+r_pl
        g = 9.81
        !nu_sw = 1.06e-6
        Dstar = (rho_tot-rho_sw)*g*(2.*r_tot)**3/(rho_sw*nu_sw**2)
        if (Dstar < 0.05) then
                wstar = Dstar**2*1.71e-4
        else
                wstar = 10.**(-3.76715+1.92944*LOG10(Dstar) &
                -0.09815*LOG10(Dstar)**2 - 0.00575*LOG10(Dstar)**3 &
                + 0.00056*LOG10(Dstar)**4)
        endif
        inner = (rho_tot-rho_sw)*g*wstar*nu_sw/rho_sw
        settling_vel_func = -sign((abs(inner))**(1.0/3.0),inner)



 END FUNCTION settling_vel_func

END MODULE BF_MOD

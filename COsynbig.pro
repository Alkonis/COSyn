

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Script 1 of 4               ;;;;;;;;
;; load parameters             ;;;;;;;;  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;create an array of the parameters
;vary R_in, R_out, v_turb, inc, T_rot0, R_rot_alpha, H_den0, and H_den_alpha
best_diff=1d1

num_guesses=1501
matrix=FLTARR(7,num_guesses)
;diff_array=FLTARR(N_ELEMENTS(rbig),num_guesses)


qw=8L

matrix(0,*)=RANDOMU(qw+47L,num_guesses)*900.+100   ;range layers from 100-1000
matrix(1,*)=RANDOMU(qw+69L,num_guesses)*15.+6.     ;inner radius goes from 5-25
matrix(2,*)=RANDOMU(qw+10L,num_guesses)*76.+50     ;outer radius scales from 50-125
matrix(3,*)=RANDOMU(qw+42L,num_guesses)*5e5+1e5  ;v_turb goes from 0.5 - 4.5 km/s
matrix(4,*)=RANDOMU(qw-3L,num_guesses)*3500.+1000  ;fiducial rotational temperature
matrix(5,*)=RANDOMU(qw+37L,num_guesses)*.5+.2            ;power-law for temperature anging from 0-1
matrix(6,*)=RANDOMU(qw+38L,num_guesses)*99.+1.    ;Relative luminosity ranges from 1 - 100
;matrix(6,*)=2.5e10 ;10^(RANDOMU(x426L,num_guesses)*2+9.) ;fiducial density from 1e9 - 1e11
;matrix(7,*)=0.15 ;RANDOMU(x427L,num_guesses)*.5         ;power law for density ranging from 0.0 - 0.5

GOTO, skip_param
layers=300		;number of layers in disk (76 for b=2)
rel_lum=20.		;UV luminosity relative to HD141569
disk_in=13		;inner edge of disk
dist=1.496e13*disk_in	;convert inner edge of disk to cm
disk_out=100.0		;Outer edge of disk in AU
Mstar=2.4		;Stellar mass in solar units
inc=40.*!pi/180.	;Disk inclination
v_turb=3.e5		;Turbulent Velocity
T_rot0_fl=2.5e3		;Fiducial temp at 1AU
T_rot_alpha_fl=0.25	;Power law of rotational temp
T_rot0_cl=2.5e3		;Fiducial temp at 1AU
T_rot_alpha_cl=0.25
H_den0=2.5e10		;Fiducial density at 1AU
H_den_alpha=.15;0.05	;Power law of density
skip_param:

;FIXED PARAMETERS
f_i=1995.0		;Frequency range in wavenumbers
f_f=2179.0
rel_lum=20.
Mstar=2.4		;Stellar mass in solar units
H_den0=2.5e10		;Fiducial density at 1AU
H_den_alpha=.15;0.05	;Power law of density
inc=40.*!pi/180.	;Disk inclination
X12CO_13CO_fl =65./30.          ;C-12/13 ratio
X12CO_C18O_fl =550./16.25       ;O-16/18 ratio
X12CO_13CO_cl =65.		;C-12/13 ratio
X12CO_C18O_cl =560.
H_den0=2.5e10
H_den_alpha=0.1 

FOR bigi=0,num_guesses-1 DO BEGIN ;THE BIG LOOP

IF bigi EQ 0 THEN BEGIN ;start with best chi-by-eye fit....
layers=300
inc=40.*!pi/180.
disk_in=13.
dist=1.496e13*disk_in
disk_out=100.0
v_turb=3.e5
T_rot0_fl=2.5e3
T_rot_alpha_fl=0.25
T_rot0_cl=T_rot0_fl             ;Fiducial temp at 1AU
T_rot_alpha_cl=T_rot_alpha_fl
rel_lum=20.

ENDIF ELSE BEGIN
layers=matrix(0,bigi)
disk_in=matrix(1,bigi)
dist=1.496e13*disk_in
disk_out=matrix(2,bigi)
v_turb=matrix(3,bigi)
T_rot0_fl=matrix(4,bigi)
T_rot_alpha_fl=matrix(5,bigi)
T_rot0_cl=T_rot0_fl             ;Fiducial temp at 1AU
T_rot_alpha_cl=T_rot_alpha_fl
rel_lum=matrix(6,bigi)

	
ENDELSE

d=double(103.d*3.0856d18);Distance to star in cm
inst_res=6.0 		;Resolution of instrument

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Script 2 of 4
;this code reads in a set of einstein A's, 
;associated freqencies and wavelengths. It takes the UV luminosity and inner 
;disk edge from the user. 

; The code calculates the g-factors and solves the rate-balance equation
; for successive depths into a "slab" of gas. Each depth is N=2x10^12 cm-2
; which corresponds to tau=0.1 for the transition from v"=0 with the largest 
; oscillator strength.The relative population of each ground electronic state
; vibrational level is calculated. for each depth. The total column density of 
; each ring is then summed for each level

; This process is repeated for multiple slabs in steps of 1AU. 

;;
;;
;;
;Added a term to make general if UV luminosity is known relative to HD 141569, 
;and if ;distance is known

;;;;;;;;;;;;;;;;
;;;; STEP 1 ;;;;
;;;;;;;;;;;;;;;;
;define constants:

c=2.997924562e5			;speed of light in km/s
hc=6.626196e-27*2.997924562e10 	;erg*cm
hck=hc/1.380622e-16		;erg*K
cer=8.85282e-13			;pi*e^2/mc (cm)

;vib_einA=[30.84, 59.30, 85.40, 109.28, 131.01, 150.71, $
;	  169.65, 185.67, 199.93, 212.53] ; old values

;vib_einA=[32.87, 66.00, 97.25, 125.52, 152.16, 177.56, $
;	  202.07, 226.07, 249.79, 273.11] ; groundstate band A's at 2000K

vib_einA=[34.60, 67.68, 98.40, 126.99, 153.59, 178.31, 201.35, 223.10, 244.15, 265.21] ; groundstate band A's at 1000K

;vib_einA=[35.56, 68.91, 100.10, 129.21, 156.27, 181.36, $
;	  204.56, 225.92, 245.54, 263.44] ;ground state band A's at 200K



;;;;;;;;;;;;;;;;
;;;; STEP 4 ;;;;
;;;;;;;;;;;;;;;;
;Define annuli - eventually adjust to constant steps in velocity space (2km/s?)

;Create annuli in steps of 
;0.01 AU for r<.1 AU
;0.10 AU for 0.1<r<1
;1.00 AU for 1<r<disk_out

r_index_a	=0
r_index_b	=0
r_index_c	=0

r=disk_in
;note that ring width is redefined in COsynprof2! 
IF r LT .1 THEN BEGIN
REPEAT BEGIN
	r_index_a=r_index_a+1
	r=[r,disk_in+0.01*r_index_a]
ENDREP UNTIL MAX(r) GE 1.0 OR MAX(r) GE disk_out ;.1 OR MAX(r) GE disk_out
ENDIF

maxra=MAX(r)
;IF MAX(r) LT 10.0 AND MAX(r) GE 1.0 THEN BEGIN
IF MAX(r) LT 1 AND MAX(r) GE .1 THEN BEGIN
REPEAT BEGIN
	r_index_b=r_index_b+1
	r=[r,maxra+0.1*r_index_b]

;ENDREP UNTIL MAX(r) GE 10.0 OR MAX(r) GE disk_out
ENDREP UNTIL MAX(r) GE 1.0 OR MAX(r) GE disk_out
ENDIF

maxrb=MAX(r)

;IF MAX(r) LE disk_out AND MAX(r) GE 10.0 THEN BEGIN
IF MAX(r) LE disk_out AND MAX(r) GE 1.0 THEN BEGIN
REPEAT BEGIN
	r_index_c=r_index_c+1
	r=[r,maxrb+1.0*r_index_c]
ENDREP UNTIL MAX(r) GE disk_out
ENDIF

rdisk_index=WHERE(r LT 13.)

rdisk=r*1.496E13
steps=N_ELEMENTS(rdisk)


;ENTIRE CODE MUST BE RUN TWICE. ONCE WITH COLLISIONS AND ONCE WITHOUT.
;THIS IS A HACK TO PROPERLY TREAT RARE ISOTOPOLOGUES 

FOR coll_loop=0,1 DO BEGIN ; 0: include collisions, 1: neglect collisions

;;;;;;;;;;;;;;;;
;;;; STEP 3 ;;;;
;;;;;;;;;;;;;;;;

;Define arrays for rate equation calculation
;Account for lower and upper state as a function of depth into each annuli

EinA=fltarr(10,12)		;Einstein A in units of s^-1
lam_ang=fltarr(10,12)		;wavelength in units Angstrom
wavenum=fltarr(10,12)		;frequency in wavenumbers
tau_0=fltarr(10,12,layers)	;opacity of a given ring
dfdt_0=fltarr(10,12,layers)
dwdn=fltarr(10,12,layers)
Nv=fltarr(21,layers,steps)	;relative vibratinoal population for 
				;each depth into gas
g=fltarr(10,12,layers,steps)
rate_eqtn=DBLARR(21,21,layers,steps)
rate_eqtn(*,*,*,*)=double(0.00)
rate_eqtn(*,20,*,*)=double(1.00)


;;;;;;;;;;;;;;;;
;;;; STEP 4 ;;;;
;;;;;;;;;;;;;;;;
;read in tabulated data

openr,1,'inc/EinA.txt' 
readf,1,EinA 
close,1

openr,1,'inc/lambda.txt' 
readf,1,lam_ang 
close,1

openr,1,'inc/wavenum.txt'
readf,1,wavenum 
close,1

;Determine Luminosity of UV field from .13 to .23 um
;Assume polynomial with same shape as B9.5 star HD141569
;This needs to be modified to put in correct UV field for particular stars

;L=double(lam_ang)
;L=-double(3.8902018e+20)+double(1.2203011e+18*lam_ang)             $
;  -double(1.5831336e+15*lam_ang^2)+double(1.0886519e+12*lam_ang^3) $
;  -double(4.1867630e+8*lam_ang^4)+double(8.5411182e+4*lam_ang^5)   $
;  -double(7.2239110*lam_ang^6)
;
;L=double(L)*double(1d18) ; Avoid floating point overflow
;L=L*rel_lum
restore,filename='inc/HD100546_luminosity.dat' ;units are erg/s/Ang
L=L*rel_lum
;Calculate oscillator strengths from Einstein A's

fAX=EinA*(3.038/2.03)/wavenum^2	;Taken from Beegle et al. 1999
fXA=2*fAX			;Take into account stat weight

;;;;;;;;;;;;;;;;
;;;; STEP 5 ;;;;
;;;;;;;;;;;;;;;;
;Fill constants into rate equation outside of loop

for k=0,steps-1 do begin
	for j=0,layers-1 do begin
		rate_eqtn(11:20,0:10,j,k)=EinA(0:9,0:10)
	endfor
endfor

		FOR i=0,8 DO BEGIN
			rate_eqtn(i+1,i+1,*,*)=-vib_EinA(i)
			rate_eqtn(i+2,i+1,*,*)=vib_EinA(i+1)
		ENDFOR
		rate_eqtn(1,0,*,*)=vib_EinA(0)
		rate_eqtn(10,10,*,*)=-vib_EinA(9)

		FOR i=0,8 do begin
			rate_eqtn(i+11,i+11,*,*)=-TOTAL(EinA(i,*),2)
		ENDFOR






;;;;;;;;;;;;;;;;
;;;; STEP 6 ;;;;
;;;;;;;;;;;;;;;;
;Calculate temperature profile for disk from parameter file
;Calculate intrinsic line shape as function of temperature

;FOR k=1,steps-1 DO Rdisk(k)=Rdisk(k-1)+step_size
T_rot_fl=T_rot0_fl  *(1.5e13/rdisk)^T_rot_alpha_fl
;H_den=H_den0 *(1./rdisk)^H_den_alpha
T_rot_index=WHERE(T_rot_fl GE 3.5E3, T_rot_cnt)
IF T_rot_cnt NE 0 THEN T_rot_fl(T_rot_index)=3.5e3

;;;;;;col. gas.
T_rot_cl=T_rot0_cl  *(1.496e13/rdisk)^T_rot_alpha_cl
H_den=H_den0 *(1.496E13/rdisk)^H_den_alpha
;H_den(rdisk_index)=1e10
T_rot_index=WHERE(T_rot_cl GE 3.5E3, T_rot_cnt)
IF T_rot_cnt NE 0 THEN T_rot_cl(T_rot_index)=3.5e3



IF rel_lum LE 1E-3 THEN BEGIN	
	GOTO, skip_fluorcalc
ENDIF

b_tot=SQRT(2.*1.38e-16*6.02e23*T_rot_fl/28. + v_turb^2) ;UNITS=cm/s

;;;;;;;;;;;;;;;;
;;;; STEP 7 ;;;;
;;;;;;;;;;;;;;;;
;In general g=dW/dN *piF/hcnu
;to get dWdN we first need dFdt and thus F(Tau)
;For large values of tau, F(tau)~SQRT(ln(tau)). This is good to .1% for tau>20
;So get dFdtau numerically for tau<20

tau=findgen(2e3)/1e2
F_tau=fltarr(2000)
;Integrate function for given value of tau using trapezoid rule
sum=0
n=float(1E3)                   ; Number of steps for integral
a=float(0.0)                   ; Starting point
b=float(5.0)                  ; Make sure F(b) ~ 0   
delta=float((b-a)/n)           ; Step size of integral
FOR ii = 0., 1999. DO BEGIN
	sum=0.5*( (1 - exp(-tau(ii)*exp(-a^2))) + (1 - exp(-tau(ii)*exp(-b^2))) )
	FOR j = 1, n-1 DO BEGIN
	        x=a+delta*j
	        sum=sum+(1-exp(-tau(ii)*exp(-x^2)))
	ENDFOR
	sum=sum*delta
	F_tau(ii)=sum	
ENDFOR

;NOW WE HAVE F(tau) from tau=0-20
;WHEN TAU > 20, then F(tau)=SQRT(ln(tau))
dFdt=deriv(tau,F_tau)


;;;;;;;;;;;;;;;;
;;;; STEP 8 ;;;;
;;;;;;;;;;;;;;;;
;Create rate equation and solve for pops

FOR k=0,steps-1 DO BEGIN
	FOR i=0,7 DO BEGIN
           rate_eqtn(i+1,i+1,*,k)=-vib_EinA(i)
           rate_eqtn(i+2,i+1,*,k)=vib_EinA(i+1)
	ENDFOR
        rate_eqtn(1,0,*,k)=vib_EinA(0)
	rate_eqtn(9,9,*,k)=-vib_EinA(8)

IF coll_loop EQ 1 THEN goto, skip_coll        
        k_HCO_dn=(7.57e-15*T_rot_cl(k)/(1-exp(-3084./T_rot_cl(k))))*H_den(k)
        k_HCO_up=k_HCO_dn*exp(-3084./T_rot_cl(k))
	rate_eqtn(0,0,*,k)=rate_eqtn(0,0,*,k)-k_HCO_up
	rate_eqtn(1,0,*,k)=rate_eqtn(1,0,*,k)+k_HCO_dn
        
        FOR i=1,8 DO begin
		rate_eqtn(i-1,i,*,k)  =rate_eqtn(i-1,i,*,k)+k_HCO_up
		rate_eqtn(i  ,i,*,k)  =rate_eqtn(i,  i,*,k)-k_HCO_dn-k_HCO_up
		rate_eqtn(i+1,i,*,k)  =rate_eqtn(i+1,i,*,k)+k_HCO_dn
	ENDFOR
skip_coll:

;;;;;;;;;;NOW ADD UV PUMPING;;;;;;;;;;;;;;;
	Fuv=L/(4*!pi*rdisk(k)^2)
	FOR j=0,layers-1 DO BEGIN
		IF j GT 0 THEN BEGIN
			FOR i=0,9 DO BEGIN 
                          tau_0(i,*,j)=TOTAL(Nv(0:11,0:j,k),2)*7.55858e12*0.02654*fXA(i,*)*lam_ang(i,*)*1e-8/(SQRT(!pi)*b_tot(k))
                        ENDFOR
		ENDIF

		FOR ii=0,N_ELEMENTS(tau_0(0,*,j))-1 DO BEGIN	
			IF tau_0(0,ii,j) LT 20.0 THEN BEGIN 
				dFdt_0_index=WHERE(tau EQ ROUND(tau_0(0,ii,j)*10.)/10.,count)
				IF count NE 0 THEN dFdt_0(*,ii,j)=dFdt(dFdt_0_index) 
			ENDIF ELSE BEGIN
				dFdt_0(*,ii,j)=1./(2.*tau_0(0,ii,j)*SQRT(alog(tau_0(0,ii,j))))
                  
                ENDELSE	
		ENDFOR

		dWdN(*,*,j)=dfdt_0(*,*,j)*.02654*2.0*(lam_ang*1e-4)*(lam_ang*1e-8)*fXA/(SQRT(!pi)*c*1e5)
		g(*,*,j,k)=dWdN(*,*,j)*!pi*Fuv/(hc*wavenum)
;print,,j,k)

;if dfdt_0(0,0,j) EQ 0 then begin
;print,j
;read,x,prompt="?"
;endif
		;now add g-factors:
                rate_eqtn(0,0,j,k)=rate_eqtn(0,0,j,k)-TOTAL(g(*,0,j,k),1) ;sum of all
                                                        ;transitions from v=0 
							;ground state
		FOR i=1,10 DO rate_eqtn(i,i,j,k)=rate_eqtn(i,i,j,k) $
				-TOTAL(g(*,i,j,k),1)
		
		FOR i=0,8 DO rate_eqtn(0:10,11+i,j,k)=rate_eqtn(0:10,11+i,j,k)+g(i,0:10,j,k)
		rate_eqtn(1,10,j,k)=0.0

		;Now solve the system of equations to calculate the relative 
		;populations:

		SVDC, rate_eqtn(*,*,j,k), w, u, v, /DOUBLE
		z=fltarr(21) & z(*)=0.0 & z(20)=1.0 	;"solution" to system 
							;of equations
		Nv(*,j,k)=SVSOL(u, w, v, z, /DOUBLE) 	;these are the relative
							; populations
		Nv_index=WHERE(Nv(*,j,k) LT 0.)
		IF (Nv_index(0) NE -1) THEN Nv(Nv_index,j,k)=double(0)

	ENDFOR
ENDFOR
skip_fluorcalc:





IF coll_loop EQ 0 THEN BEGIN 
   Nv_coll=Nv
   R_ind=WHERE(rdisk LT 13.*1.496E13, count)
   IF count NE 0 THEN Nv(*,*,R_ind)=1.0*Nv(*,*,R_ind)
   tot_col_fluor        =total(Nv*7.55858e12,2)
   tot_col_fluor_back=tot_col_fluor
ENDIF

IF coll_loop EQ 1 THEN BEGIN
   Nv_nocoll=Nv
   R_ind=WHERE(rdisk LT 13.*1.496E13, count)
   IF count NE 0 THEN Nv(*,*,R_ind)=1.*Nv(*,*,R_ind)
   tot_col_fluor_nocoll =total(Nv*7.55858e12,2)
ENDIF

;total column density of each vibrational level in each ring
;for a flaring disk, the light enters the disk at a small angle.
;the column of excited gas is really tot_col_fluor*sin(angle)
;IF the disk is nearly isothermal (our CO is), then H goes as r^3/2. If the scale height at 10AU
; is 0.35 as claimed by Bouwman, then the angle is ~8deg. Since there is an inner wall, we 
;we assume the flux there is nearly normal

ENDFOR ; END OF COLL LOOP


;now correct for angle of incidence of light onto disk
;H(r)=SQRT(kTR^3/mum_H*GMstar)=5.59647E-10*SQRT(T/(Mstar/Msun))*R^3/2
;dH/dr=1.5*H/r (assuming T is nearly constant with respect to R). This
;is pretty close as T scales as R^-.3. If disk isn't really
;flared, then disk will be flatter and the light want illuminate outer
;disk. this is likely why UV fluorescence is relatively rare. This
;part of the code should be tweaked for each star...

m_disk=1.5*(5.59647E-10)*SQRT(T_rot_fl/Mstar)*(Rdisk)^(0.5) ;SLOPE of disk where UV impinges on disk dH(r)/dr
m_uv=(5.59647E-10)*SQRT(T_rot_fl/Mstar)*SQRT(Rdisk)       ;slope of uv rays H(r)/R   
phi=-atan((m_uv-m_disk)/(1+m_uv*m_disk))
phi(0)=!pi/2.d ;assume line hits inner annulus directly. 

xi=N_elements(tot_col_fluor(*,0))-1 ;vib level index
yi=N_ELEMENTS(tot_col_fluor(0,*))-1 ;annuli index

for j=1,yi do begin
 tot_col_fluor(0:xi,j)=sin(phi(j))*tot_col_fluor(0:xi,j)
 tot_col_fluor_nocoll(0:xi,j)=sin(phi(j))*tot_col_fluor_nocoll(0:xi,j)
ENDFOR
tot_col=tot_col_fluor(0:9,*)
tot_col_coll=tot_col



r=dist/1.496E13
Rdisk=Rdisk/1.496e13          ;convert to AU to keep consistent


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Script 3 of 4
;Compute intrinsic spectrum allowing for line opacity
;Accounts for line overlap along the lines of Carr code
;
;1) Calculate dispersion
;2) Calculate opacity at line center for each line in each ring
;	->tau0=N(CO)_j*(exp(E'/kT)/zCO)*0.026541*f_jk*
;		(1-exp(-1.44*nutild/T))/(nutild*SQRT(pi)*b)
;3) find all CO lines within +/-3b of line center for each wavelength step
;4) compute wavelength function
;5) calculate tau as function of wavelength over specified range

;define constants:
hc=6.626e-27*2.9979e10 	;erg*cm
hck=1.438		;erg*K
cer=8.85282e-13		;pi*e^2/mc (cm)
v	=double(2.5)		;velocity resolution in km/s
c	=double(2.997924562E5)	;speed of light in km/s
k	=1.380622E-16		;erg/K
B	=1.9225			; cm-1
;X12CO_13CO =60./3		;12/13 ratio
;X12CO_C18O =550./3		;12/18 ratio
;b_turb=3.0			;turbulent velocity in km/s
;b_tot=2.*SQRT(alog(2))*SQRT(T_rot*0.0046 + (line_rms/1e5)^2)*1e5 ;thermal+turb broadening for CO
;b_tot=SQRT(2.*k*6.02e23*T_rot_cl/28. + v_turb^2)
;b_tot(*)=line_rms

;rotational constants for each level
B10  =1.93124
B21  =1.91376
B32  =1.89628
B43  =1.87880
B54  =1.86132
B65  =1.84384
B76  =1.82636
B87  =1.80888
B98  =1.79140
B109 =1.77392
B1110=1.76644
B1310=1.847
B1321=1.830
B1332=1.813
B1810=1.840
Bv12=[B10,B21,B32,b43,b54,b65,b76,b87,b98,b109,b1110]
Bv13=[B1310,B1321,B1332]
Bv18=[1.840]

;Read in molecular data
;content is v", J", v', J', A(s-1), E" (cm-1), freq(cm-1) for each variable

;CHANGE PATH TO LOCATION OF MOLECULAR DATA
restore, filename='inc/CO_molecdat.dat'
;Use molecular data to calculate relative rotational populations 
; as a function of species, vibrational level and radius:

;1) loop over species
;  2) loop over ring
;    3) loop over vibrational level
;      
      
N_12CO_vJ=fltarr(120,n_elements(tot_col(*,0))-1,n_elements(tot_col(0,*))) ;we only want vibrational levels v=1-9
N_13CO_vJ=fltarr(120,3,n_elements(tot_col(0,*)))
N_C18O_vj=fltarr(120,1,n_elements(tot_col(0,*)))


FOR i=0,steps-1 DO BEGIN ;loop over rings
;first do 12CO

	FOR j=0,N_ELEMENTS(tot_col(*,0))-2 DO BEGIN ; loop over vib levels 1-9
	Jup =X12CO(j,3,*)
	Jdn =X12CO(j,1,*)
	wvn =X12CO(j,6,*)
	EinA=X12CO(j,4,*)	
;      print,tot_col_fluor(j+1,i)
;read,x,prompt='tot_col_fluor(j+1,i)'
		N_12CO_vJ(*,j,i)=tot_col_fluor(j+1,i)*(2*Jup+1.) $
				*exp(-hck*Bv12(j)*Jup*(Jup+1)/T_rot_fl(i)) $
				/(T_rot_fl(i)/(hck*Bv12(j))) $
                                + tot_col_coll(j+1,i)*(2*Jup+1.) $
				*exp(-hck*Bv12(j)*Jup*(Jup+1)/T_rot_cl(i)) $
				/(T_rot_cl(i)/(hck*Bv12(j)))
		
;print,tot_col_fluor(j+1,i)
;print,"---------------"
;print,Jup
;print,"--------"
;print,N_12CO_vJ(*,j,i)
;print,Jup
;print,Jdn
;print,wvn
;print,EinA
;read,x,prompt="tot_col_fluor(j+1,i)----Jup"
	ENDFOR

;Second do 13CO	
	FOR j=0,2 DO BEGIN
	Jup =X13CO(j,3,*)
	Jdn =X13CO(j,1,*)
	wvn =X13CO(j,6,*)
	EinA=X13CO(j,4,*)	
		N_13CO_vJ(*,j,i)=(tot_col_fluor_nocoll(j+1,i)/X12CO_13CO_fl)*(2*Jup+1.) $
				*exp(-hck*Bv13(j)*Jup*(Jup+1)/T_rot_fl(i)) $
				/(T_rot_fl(i)/(hck*Bv13(j))) ;$
;                                + (tot_col_coll(j+1,i)/X12CO_13CO_cl)*(2*Jup+1.) $
;				*exp(-hck*Bv13(j)*Jup*(Jup+1)/T_rot_cl(i)) $
;				/(T_rot_cl(i)/(hck*Bv13(j)))
;print,N_13CO_vJ(*,j,i)
;print,tot_col_fluor_nocoll(j+1,i)

;read,x,prompt="N_13COvj"
	ENDFOR

;Third do C18O
	FOR j=0,0 DO BEGIN
	Jup =XC18O(j,3,*)
	Jdn =XC18O(j,1,*)
	wvn =XC18O(j,6,*)
	EinA=XC18O(j,4,*)	
		N_C18O_vJ(*,j,i)=(tot_col_fluor_nocoll(j+1,i)/X12CO_C18O_fl)*(2*Jup+1.) $
				*exp(-hck*Bv18(j)*Jup*(Jup+1)/T_rot_fl(i)) $
				/(T_rot_fl(i)/(hck*Bv18(j))) ;$
;                               + (tot_col_coll(j+1,i)/X12CO_C18O_cl)*(2*Jup+1.) $
;				*exp(-hck*Bv18(j)*Jup*(Jup+1)/T_rot_cl(i)) $
;				/(T_rot_cl(i)/(hck*Bv18(j)))

;print,tot_col_fluor_nocoll(j+1,i)
;print,N_C18O_vJ(*,j,i)
;read,x,prompt="N_C18O"
	ENDFOR
ENDFOR

;Now have column density of CO as function of isotope, v, J, and R

; Construct frequency scale, set size of NIRSPEC is 4km/s per pixel
; Cover 2200cm-1 - 1800cm-1

freq_size	=alog10(f_f/f_i)/alog10(1+v/(3.*c))
freq		=dblarr(freq_size)		;frequency in wavenumbers
freq(0)		=double(f_i)
vfreq		=freq
vfreq(*)=0
FOR i=1,freq_size-1 DO freq(i)=freq(0)*(1+double(v/(3.*c)))^i
vfreq(*)=(freq(FIX(N_ELEMENTS(freq)/2.))-freq)*c/freq(FIX(N_ELEMENTS(freq)/2.))
freq_index=where(freq GE f_i AND freq LE f_f)


;Define annuli
annuli=dblarr(N_ELEMENTS(rdisk))
FOR i=0,steps-2 DO annuli(i)=!pi*(rdisk(i+1)^2-rdisk(i)^2)
annuli(steps-1)=!pi*((rdisk(steps-1)+1)^2-rdisk(steps-1)^2)
annuli=annuli*double(2.25e26)	;The area of each annulus in cm2.

;Now compute synthetic spectrum
stick_spec_12CO=dblarr(freq_size,steps)
stick_spec_12CO(*,*)=double(0.0)
stick_spec_13CO=stick_spec_12CO
stick_spec_C18O=stick_spec_12CO
stick_spec_tot=stick_spec_12CO
FOR i=0, steps-1 DO BEGIN					;Loop over annuli
	FOR j=0,N_ELEMENTS(tot_col(*,0))-2 DO BEGIN		;Loop over viblev
   	    Jup =X12CO(j,3,*)
	    Jdn =X12CO(j,1,*)
	    wvn =X12CO(j,6,*)
	    EinA=X12CO(j,4,*)
		FOR k=0,N_ELEMENTS(X12CO(0,0,*)) - 1 DO BEGIN 	;Loop over Rotlev
			IF wvn(k) GE f_i AND wvn(k) LE f_f THEN BEGIN
		  	   A0=N_12CO_vJ(k,j,i)*hc*wvn(k)*EinA(k)
			   A1=wvn(k)
			   A2=b_tot(i)*wvn(k)/(c*1e5)
;print,n_12co_vj(k,j,i)
;print,a0
;print,a1
;print,a2

;read,x,prompt="n_12co_vj(k.j,i),a0,a1,a2"
;print,-((a1-freq)/a2)^2
;read,x,prompt="inexp"
;print,exp(-((A1-freq)/A2)^2)
;read,x,prompt="exp"
			   stick_spec_12CO(*,i)=stick_spec_12CO(*,i)+(A0/(SQRT(!pi)*A2))*exp(-((A1-freq)/A2)^2)	
;print,exp(-((A1-freq)/A2)^2)
;read,x,prompt="exp factor:"
;print,stick_Spec_12CO(*,i)
;read,x,prompt="stick_Spec_12co(*,i)"
			ENDIF

		ENDFOR

	ENDFOR

	FOR j=0,2 DO BEGIN		;Loop over viblev
   	    Jup =X13CO(j,3,*)
	    Jdn =X13CO(j,1,*)
	    wvn =X13CO(j,6,*)
	    EinA=X13CO(j,4,*)
		FOR k=0,N_ELEMENTS(X13CO(0,0,*)) - 1 DO BEGIN 	;Loop over Rotlev
			IF wvn(k) GE f_i AND wvn(k) LE f_f THEN BEGIN
		  	   A0=N_13CO_vJ(k,j,i)*hc*wvn(k)*EinA(k)
			   A1=wvn(k)
			   A2=b_tot(i)*wvn(k)/(c*1e5)
			   stick_spec_13CO(*,i)=stick_spec_13CO(*,i)+(A0/(SQRT(!pi)*A2))*exp(-((A1-freq)/A2)^2)


;print,n_13co_vj(k,j,i)
;print,a0
;print,a1
;print,a2
;
;read,x,prompt="n_13co_vj(k.j,i),a0,a1,a2"
;print,-((a1-freq)/a2)^2
;read,x,prompt="inexp"
			ENDIF
		ENDFOR
	ENDFOR

	FOR j=0,0 DO BEGIN		;Loop over viblev
   	    Jup =XC18O(j,3,*)
	    Jdn =XC18O(j,1,*)
	    wvn =XC18O(j,6,*)
	    EinA=XC18O(j,4,*)
		FOR k=0,N_ELEMENTS(XC18O(0,0,*)) - 1 DO BEGIN 	;Loop over Rotlev
			IF wvn(k) GE f_i AND wvn(k) LE f_f THEN BEGIN
		  	   A0=N_C18O_vJ(k,j,i)*hc*wvn(k)*EinA(k)
			   A1=wvn(k)
			   A2=b_tot(i)*wvn(k)/(c*1e5)
			   stick_spec_C18O(*,i)=stick_spec_C18O(*,i)+(A0/(SQRT(!pi)*A2))*exp(-((A1-freq)/A2)^2)
			ENDIF
		ENDFOR
	ENDFOR

stick_spec_tot(*,i)=(stick_spec_12CO(*,i)+stick_spec_13CO(*,i)+stick_spec_C18O(*,i))
ENDFOR


print,stick_spec_tot
read,x,prompt="stick_Spec_tot"
;annuli(0)=annuli(0) ; Area of inner wall of H/R=.35 disk at 13AU
annuli(0)=1d28 ; Projected area of inner wall at 13AU
Iten_tot=stick_spec_tot
iten_tot(*,0)=iten_tot(*,0)*2.5
Lum_tot=Iten_tot
FOR i=0,steps-1 DO Lum_tot(*,i)=Iten_tot(*,i)*annuli(i)
Flux_tot=double(lum_tot)/(4.*!pi*double(d)^2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Script 4 of 4

;
; INPUTS:
; Temperature, 13CO/12CO ratio, C18O/12CO ratio, rms linewidth

;USER DEFINE PARAMETERS
;Mstar=0.5
;inc=45.*!pi/180.	;disk inclination
r_in=rdisk(0)		; Inner edge of disk
r_out=max(rdisk) 		; Outer edge of disk
;dv=v/2.			; step size in units of km/s (typically 1/2 rms of gas)
dv=1.

;Trot0=T_rot		;rotational temperature at inner edge
;inst_res=6.0		;resolution of inst

flux_tot_slit=flux_tot

;;;;Account for slit losses
;IF SEEING IS .6" THEN:
SLIT_LOSS=100.17*rdisk^(-1.260)
DISK_BEYOND_SLIT=63.

;IF SEEING IS .4" THEN:
;SLIT_loss=44.373*rdisk^(-1.179)
;DISK_BEYOND_SLIT=32.
FOR i=0,N_ELEMENTS(Rdisk) -1 DO BEGIN
   IF RDISK(i) GT DISK_BEYOND_SLIT THEN flux_tot_slit(*,i)=flux_tot_slit(*,i)*slit_loss(i)
ENDFOR

;n_rings=(ROUND(SQRT(887.2*Mstar/r_in)-SQRT(887.2*Mstar/r_out)))/dv + 1.
n_rings=N_ELEMENTS(rdisk)
;define width of rings
dr=rdisk 
FOR i=0,n_rings-1 DO BEGIN
	IF rdisk(i) LT 0.1 THEN dr(i)=0.01 ELSE $
		IF rdisk(i) LT 1.0 THEN dr(i)=0.1 $ 
		ELSE dr(i)=1.0
ENDFOR

vmax=fltarr(n_rings)
vmax=ROUND(SQRT(887.2*Mstar/rdisk))

; Construct frequency scale, set size of NIRSPEC is 4km/s per pixel
; Cover 2200cm-1 - 1800cm-1

vfreq		=freq
vfreq(*)=0
vfreq(*)=(freq(FIX(N_ELEMENTS(freq)/2.))-freq)*c/freq(FIX(N_ELEMENTS(freq)/2.))
freq_index=where(freq GE f_i AND freq LE f_f)
total_spec=Iten_tot
vel_spec=total_spec(*,0)


;;;;as a first pass, the user should select the line s/he
;;;;want's to generate the spectroastrometric signal for. 
num_lines=22
v_line5=(2032.3528-freq)*2.9979e5/2032.3528 ;start with P26 line! - nah! 
                                            ; add in the uv fluoresced lines 
                                            ; too for the fun of it!






v_line4=(2034.4083-freq)*2.9979e5/2034.4083
v_line1=(2033.4174-freq)*2.9979e5/2033.4174
v_line2=(2033.1423-freq)*2.9979e5/2033.1423
v_line3=(2032.8258-freq)*2.9979e5/2032.8258
v_line0=(2030.1586-freq)*2.9979e5/2030.1586

fco_list=[2034.9211, 2034.7209, 2034.1352, 2034.0473, 2032.8701, $
          2032.8096, 2032.5616, 2032.2011, 2030.5160, 2030.3076, $
          2029.6559, 2029.4427, 2029.4297, 2029.3679, 2029.2362, 2029.1276]

v_line_arr=FLTARR(N_ELEMENTS(freq),20)
FOR k=0,15 DO v_line_arr(*,k)=(fco_list(k)-freq)*2.9979e5/fco_list(k)




goto,skip_lowJ_lines
;Redefine for order near 2145cm-1
num_lines=6
v_line0=(2142.4729-freq)*2.9979e5/2142.4729 
v_line1=(2142.7200-freq)*2.9979e5/2142.7200
v_line2=(2142.9538-freq)*2.9979e5/2142.9538
;v_line3=(2149.4886-freq)*2.9979e5/2149.4886
v_line3=(2144.03-freq)*2.9979e5/2144.03
v_line4=(2145.9174-freq)*2.9979e5/2145.9174
v_line5=(2145.9988-freq)*2.9979e5/2145.9988
skip_lowJ_lines:

v_line=fltarr(N_ELEMENTS(v_line0),num_lines)
v_line(*,0)=v_line0
v_line(*,1)=v_line1
v_line(*,2)=v_line2
v_line(*,3)=v_line3
v_line(*,4)=v_line4
v_line(*,5)=v_line5
FOR k=0,15 DO v_line(*,6+k)=v_line_arr(*,k)

;2142.4729 ;v=2-1 P6
;2142.7200 ;v=3-2 R14
;2142.9538 ;v=4-3 P23
;2144.03   ;v=1-0 R13 (13CO)
;2145.9174 ;v=3-2 R15
;2145.9988 ;v=2-1 P7



grid=fltarr(26)

;f_fund		=[2027.6491, 2032.3528, 2037.0253]
;f_uv		=[2034.4083, 2033.4200, 2033.1425, 2032.8267,2030.1587 ]
v_line_index0=WHERE(v_line0 GT -15. AND v_line0 LT 15.)
v_line_index1=WHERE(v_line1 GT -15. AND v_line1 LT 15.)
v_line_index2=WHERE(v_line2 GT -15. AND v_line2 LT 15.)
v_line_index3=WHERE(v_line3 GT -15. AND v_line3 LT 15.)
v_line_index4=WHERE(v_line4 GT -15. AND v_line4 LT 15.)
v_line_index5=WHERE(v_line5 GT -15. AND v_line5 LT 15.)
v_line_index6=WHERE(v_line(*,6) GT -15. AND v_line(*,6) LT 15.)
v_line_index7=WHERE(v_line(*,7) GT -15. AND v_line(*,7) LT 15.)
v_line_index8=WHERE(v_line(*,8) GT -15. AND v_line(*,8) LT 15.)
v_line_index9=WHERE(v_line(*,9) GT -15. AND v_line(*,9) LT 15.)
v_line_index10=WHERE(v_line(*,10) GT -15. AND v_line(*,10) LT 15.)
v_line_index11=WHERE(v_line(*,11) GT -15. AND v_line(*,11) LT 15.)
v_line_index12=WHERE(v_line(*,12) GT -15. AND v_line(*,12) LT 15.)
v_line_index13=WHERE(v_line(*,13) GT -15. AND v_line(*,13) LT 15.)
v_line_index14=WHERE(v_line(*,14) GT -15. AND v_line(*,14) LT 15.)
v_line_index15=WHERE(v_line(*,15) GT -15. AND v_line(*,15) LT 15.)
v_line_index16=WHERE(v_line(*,16) GT -15. AND v_line(*,16) LT 15.)
v_line_index17=WHERE(v_line(*,17) GT -15. AND v_line(*,17) LT 15.)
v_line_index18=WHERE(v_line(*,18) GT -15. AND v_line(*,18) LT 15.)
v_line_index19=WHERE(v_line(*,19) GT -15. AND v_line(*,19) LT 15.)
v_line_index20=WHERE(v_line(*,20) GT -15. AND v_line(*,20) LT 15.)
v_line_index21=WHERE(v_line(*,21) GT -15. AND v_line(*,21) LT 15.)



iten_line0=TOTAL(Iten_tot(v_line_index0,*),1)*5.65e-3 ;stepsize in wavenumbers intensity is erg/s/cm2
iten_line1=TOTAL(Iten_tot(v_line_index1,*),1)*5.65e-3
iten_line2=TOTAL(Iten_tot(v_line_index2,*),1)*5.65e-3
iten_line3=TOTAL(Iten_tot(v_line_index3,*),1)*5.65e-3
iten_line4=TOTAL(Iten_tot(v_line_index4,*),1)*5.65e-3
iten_line5=TOTAL(Iten_tot(v_line_index5,*),1)*5.65e-3
iten_line6=TOTAL(Iten_tot(v_line_index6,*),1)*5.65e-3
iten_line7=TOTAL(Iten_tot(v_line_index7,*),1)*5.65e-3
iten_line8=TOTAL(Iten_tot(v_line_index8,*),1)*5.65e-3
iten_line9=TOTAL(Iten_tot(v_line_index9,*),1)*5.65e-3
iten_line10=TOTAL(Iten_tot(v_line_index10,*),1)*5.65e-3
iten_line11=TOTAL(Iten_tot(v_line_index11,*),1)*5.65e-3
iten_line12=TOTAL(Iten_tot(v_line_index12,*),1)*5.65e-3
iten_line13=TOTAL(Iten_tot(v_line_index13,*),1)*5.65e-3
iten_line14=TOTAL(Iten_tot(v_line_index14,*),1)*5.65e-3
iten_line15=TOTAL(Iten_tot(v_line_index15,*),1)*5.65e-3
iten_line16=TOTAL(Iten_tot(v_line_index16,*),1)*5.65e-3
iten_line17=TOTAL(Iten_tot(v_line_index17,*),1)*5.65e-3
iten_line18=TOTAL(Iten_tot(v_line_index18,*),1)*5.65e-3
iten_line19=TOTAL(Iten_tot(v_line_index19,*),1)*5.65e-3
iten_line20=TOTAL(Iten_tot(v_line_index20,*),1)*5.65e-3
iten_line21=TOTAL(Iten_tot(v_line_index21,*),1)*5.65e-3

;print,size(iten_line0)
;help,freq
;print,size(freq)
;read,x,prompt="?"
gs=FLTARR(2,11)
gs(0,*)=FINDGEN(11)-5.0 ; [-5,-4,-3,-2,-1,0,1,2,3,4,5]
gs(1,*)=exp(-gs(0,*)^2/(12./1.665)^2)/TOTAL(exp(-gs(0,*)^2/(6./1.665)^2))
;[0.000250768, 0.00321749, 0.0234145, 0.0966443, 0.226251, 0.300420, 0.226251, 0.0966443, 0.0234145, 0.00321749, 0.000250768]
FOR j=0,n_rings-1 DO BEGIN
;	n_seg=2.*ROUND(vmax(j))/dv +1.0
;        n_seg=2.*n_seg ;go all the way around!

        n_seg=4.*ROUND(vmax(j))/dv
        vseg=fltarr(n_seg)
        i=-1
        REPEAT BEGIN
           i=i+1
           vseg(i)=vmax(j)-i
        ENDREP UNTIL vseg(i) LE -vmax(j)
        REPEAT BEGIN
           i=i+1
           vseg(i)=vseg(i-1)+1
        ENDREP UNTIL vseg(i) EQ vmax(j)-1
;        vseg=-vseg
	phase=fltarr(n_seg)
	phase(0)=0.0
	phase(0:n_seg/2.)=acos(vseg(0:n_seg/2.)/vseg(0)) ;go from 0deg -> 180deg
        phase(n_seg/2.+1:n_seg-1)=2.*!pi-acos(vseg(n_seg/2.+1:n_seg-1)/vseg(0)) ;now finish the circle...
        dphase=fltarr(n_seg)
        
        FOR i=1,n_seg-2 DO dphase(i)=(phase(i+1)-phase(i))/2. + (phase(i)-phase(i-1))/2.
        dphase(0)=phase(1)
        dphase(n_seg-1)=dphase(1)

	area=fltarr(n_seg)
	area(*)=0
        


;        for k=0,n_seg - 2 do phase_mid(k)=phase(k)+(phase(k+1)-phase(k))/2.
;        phase_mid(n_seg-1)=phase(n_seg-1)+(!pi-phase(n_seg-1))/2.
	IF j NE n_rings-1 THEN BEGIN
		FOR i=0,n_seg-1 DO BEGIN
			area(i)=dphase(i)*(rdisk(j)+dr(j)/2.)*dr(j) ;units of square AU
;                        IF sin(phase(i))*rdisk(j) GT DISK_BEYOND_SLIT AND sin(phase(i+1))*rdisk(j) GT DISK_BEYOND_SLIT THEN area(i)=0 ;zero out segments not in slit

			vel_spec=vel_spec+interpol(total_spec(*,j)*area(i), $
                                                   freq+vseg(i)*sin(inc)*freq/c,freq)	
      
;print,vel_spec
;read,x,prompt="vel_spec  ^ "
                        grid=[[grid],[rdisk(j),vseg(i),phase(i),area(i)*2.25d26, $
                             iten_line0(j),iten_line1(j),iten_line2(j),iten_line3(j),iten_line4(j),iten_line5(j),iten_line6(j), $
                             iten_line7(j),iten_line8(j),iten_line9(j),iten_line10(j),iten_line11(j),iten_line12(j),iten_line13(j), $
                             iten_line14(j),iten_line15(j),iten_line16(j),iten_line17(j),iten_line18(j),iten_line19(j),iten_line20(j), $
                             iten_line21(j)]]
GOTO, skip_line_prof1
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,0),phase(i),area(i)*2.25d26, $
                             gs(1,0)*iten_line0(j),gs(1,0)*iten_line1(j), $
                             gs(1,0)*iten_line2(j),gs(1,0)*iten_line3(j), $
                             gs(1,0)*iten_line4(j),gs(1,0)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,1),phase(i),area(i)*2.25d26, $
                             gs(1,1)*iten_line0(j),gs(1,1)*iten_line1(j),$
                             gs(1,1)*iten_line2(j),gs(1,1)*iten_line3(j),gs(1,1)*iten_line4(j),gs(1,1)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,2),phase(i),area(i)*2.25d26,gs(1,2)*iten_line0(j), $
                             gs(1,2)*iten_line1(j),gs(1,2)*iten_line2(j),gs(1,2)*iten_line3(j),gs(1,2)*iten_line4(j),gs(1,2)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,3),phase(i),area(i)*2.25d26,gs(1,3)*iten_line0(j), $
                             gs(1,3)*iten_line1(j),gs(1,3)*iten_line2(j),gs(1,3)*iten_line3(j),gs(1,3)*iten_line4(j),gs(1,3)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,4),phase(i),area(i)*2.25d26,gs(1,4)*iten_line0(j),gs(1,4)*iten_line1(j), $
                             gs(1,4)*iten_line2(j),gs(1,4)*iten_line3(j),gs(1,4)*iten_line4(j),gs(1,4)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,5),phase(i),area(i)*2.25d26,gs(1,5)*iten_line0(j),gs(1,5)*iten_line1(j), $
                             gs(1,5)*iten_line2(j),gs(1,5)*iten_line3(j),gs(1,5)*iten_line4(j),gs(1,5)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,6),phase(i),area(i)*2.25d26,gs(1,6)*iten_line0(j),gs(1,6)*iten_line1(j),$
                            gs(1,6)*iten_line2(j),gs(1,6)*iten_line3(j),gs(1,6)*iten_line4(j),gs(1,6)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,7),phase(i),area(i)*2.25d26,gs(1,7)*iten_line0(j),gs(1,7)*iten_line1(j), $
                             gs(1,7)*iten_line2(j),gs(1,7)*iten_line3(j),gs(1,7)*iten_line4(j),gs(1,7)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,8),phase(i),area(i)*2.25d26,gs(1,8)*iten_line0(j),gs(1,8)*iten_line1(j), $
                             gs(1,8)*iten_line2(j),gs(1,8)*iten_line3(j),gs(1,8)*iten_line4(j),gs(1,8)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,9),phase(i),area(i)*2.25d26,gs(1,9)*iten_line0(j),gs(1,9)*iten_line1(j), $
                             gs(1,9)*iten_line2(j),gs(1,9)*iten_line3(j),gs(1,9)*iten_line4(j),gs(1,9)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,10),phase(i),area(i)*2.25d26,gs(1,10)*iten_line0(j),gs(1,10)*iten_line1(j), $
                             gs(1,10)*iten_line2(j),gs(1,10)*iten_line3(j),gs(1,10)*iten_line4(j),gs(1,10)*iten_line5(j)]]
skip_line_prof1:

		ENDFOR
	ENDIF ELSE BEGIN
		FOR i=0,n_seg-1 DO BEGIN
			area(i)=dphase(i)*(rdisk(j)+dr(j)/2)*dr(j)
;                        IF sin(phase(i))*rdisk(j) GT DISK_BEYOND_SLIT AND sin(phase(i+1))*rdisk(j) GT DISK_BEYOND_SLIT THEN area(i)=0 ;zero out segments not in slit
			vel_spec=vel_spec+interpol(total_spec(*,j)*area(i), $
				freq+vseg(i)*sin(inc)*freq/c,freq)	
                       grid=[[grid],[rdisk(j),vseg(i),phase(i),area(i)*2.25d26,iten_line0(j),iten_line1(j),iten_line2(j),iten_line3(j),iten_line4(j),iten_line5(j),iten_line6(j), $
                             iten_line7(j),iten_line8(j),iten_line9(j),iten_line10(j),iten_line11(j),iten_line12(j),iten_line13(j), $
                             iten_line14(j),iten_line15(j),iten_line16(j),iten_line17(j),iten_line18(j),iten_line19(j),iten_line20(j), $
                             iten_line21(j)]]
GOTO, skip_line_prof2
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,0),phase(i),area(i)*2.25d26,gs(1,0)*iten_line0(j), $
                                      gs(1,0)*iten_line1(j),gs(1,0)*iten_line2(j),gs(1,0)*iten_line3(j),gs(1,0)*iten_line4(j),gs(1,0)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,1),phase(i),area(i)*2.25d26,gs(1,1)*iten_line0(j), $
                                      gs(1,1)*iten_line1(j),gs(1,1)*iten_line2(j),gs(1,1)*iten_line3(j),gs(1,1)*iten_line4(j),gs(1,1)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,2),phase(i),area(i)*2.25d26,gs(1,2)*iten_line0(j), $
                                      gs(1,2)*iten_line1(j),gs(1,2)*iten_line2(j),gs(1,2)*iten_line3(j),gs(1,2)*iten_line4(j),gs(1,2)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,3),phase(i),area(i)*2.25d26,gs(1,3)*iten_line0(j), $
                                      gs(1,3)*iten_line1(j),gs(1,3)*iten_line2(j),gs(1,3)*iten_line3(j),gs(1,3)*iten_line4(j),gs(1,3)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,4),phase(i),area(i)*2.25d26,gs(1,4)*iten_line0(j), $
                                      gs(1,4)*iten_line1(j),gs(1,4)*iten_line2(j),gs(1,4)*iten_line3(j),gs(1,4)*iten_line4(j),gs(1,4)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,5),phase(i),area(i)*2.25d26,gs(1,5)*iten_line0(j), $
                                      gs(1,5)*iten_line1(j),gs(1,5)*iten_line2(j),gs(1,5)*iten_line3(j),gs(1,5)*iten_line4(j),gs(1,5)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,6),phase(i),area(i)*2.25d26,gs(1,6)*iten_line0(j), $
                                      gs(1,6)*iten_line1(j),gs(1,6)*iten_line2(j),gs(1,6)*iten_line3(j),gs(1,6)*iten_line4(j),gs(1,6)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,7),phase(i),area(i)*2.25d26,gs(1,7)*iten_line0(j), $
                                      gs(1,7)*iten_line1(j),gs(1,7)*iten_line2(j),gs(1,7)*iten_line3(j),gs(1,7)*iten_line4(j),gs(1,7)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,8),phase(i),area(i)*2.25d26,gs(1,8)*iten_line0(j), $
                                      gs(1,8)*iten_line1(j),gs(1,8)*iten_line2(j),gs(1,8)*iten_line3(j),gs(1,8)*iten_line4(j),gs(1,8)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,9),phase(i),area(i)*2.25d26,gs(1,9)*iten_line0(j), $
                                      gs(1,9)*iten_line1(j),gs(1,9)*iten_line2(j),gs(1,9)*iten_line3(j),gs(1,9)*iten_line4(j),gs(1,9)*iten_line5(j)]]
                        grid=[[grid],[rdisk(j),vseg(i)+gs(0,10),phase(i),area(i)*2.25d26,gs(1,10)*iten_line0(j), $
                                      gs(1,10)*iten_line1(j),gs(1,10)*iten_line2(j),gs(1,10)*iten_line3(j),gs(1,10)*iten_line4(j),gs(1,10)*iten_line5(j)]]
skip_line_prof2:
		ENDFOR
	ENDELSE

total_spec(*,j)=vel_spec

ENDFOR ;end of ring
index_grid=WHERE(grid(0,*) LE rdisk(0) AND grid(0,*) GT 0. AND grid(2,*) GT !pi) ;Elements of inner ring Not part of the wall... Will make all of these areas = 0.
;index_grid=WHERE(grid(0,*) EQ rdisk(0)) ;try without inner wall
grid(3,index_grid)=0.0;,WHERE(grid(0,*) LT 14. AND grid(0,*) GT 0. AND grid(2,*) GT !pi))

;Now plop a planet into the inner rim
;Let's assume the gas is gaussian and FWHM=3km/s and sigma=1.9km/s
;v0=6km/s = 30%
;v1=5/7   =22.6%
;v2=4/8   =9.7%
;v3=

;goto, skip_planet
;shift+6km/s at phase of 36deg (.62rad)
v_planet=5.;6.0 ;6.0
r_planet=13.
phase_planet=53.*!pi/180;acos(v_planet/13.);*!pi/180.
planet_inten=((2.*6.62d-27*2.9979d10^2*2030.^3/(exp(6.626e-27*2.9979e10*2030./(1.38e-16*1.e3))-1.d)))
planet_size=0. ;.1*2.25d26 ; area of disk in cm^2

;grid=[d,v,phi,area,I]
grid=[[grid],[r_planet,v_planet+gs(0,0),phase_planet,  $
            planet_size,0,0,0,0,0,planet_inten*gs(1,0),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,1),phase_planet,  $
            planet_size,0,0,0,0,0,planet_inten*gs(1,1),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,2),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,3),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,4),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,4),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,5),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,5),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,6),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,6),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,7),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,7),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,8),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,8),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,9),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,9),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
grid=[[grid],[r_planet,v_planet+gs(0,10),phase_planet, $
             planet_size,0,0,0,0,0,planet_inten*gs(1,10),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
;skip_planet:

centroid=freq
centroid(*)=0.0



;Properly counting the dilution of the signal is absolutely
;crucial. Need to ensure that the units are kept consistent do we get
;the right eqivalent width when we are done? probably worth checking
;on!


omega=55.*!pi/180.
Lc=5.13D-23*2.9979247e10*4.*!pi*(103.*3.08d18)^2*(.05/1.16) ;this is the luminosity of the continuum.
;print,"grid(1,*):"
;print,size(grid(1,*))
;read,x,prompt="grid(1,*)^"
FOR j=0,2.*MAX(grid(1,*)) DO BEGIN;vmax(0) DO BEGIN
   index1=WHERE(grid(1,*) LE MAX(grid(1,*))-j AND grid(1,*) GT MAX(grid(1,*))-(j+1),vel_count)
;print,index1
;print,"-----------------"
;print,2*max(grid(1,*))
;read,x,"index1, maxloop^"
   IF vel_count EQ 0 THEN GOTO, no_vel_elements
   i0=WHERE(v_line0 GT MEAN(grid(1,index1))-0.5 $
             AND v_line0 LE MEAN(grid(1,index1))+0.5,count0)
   i1=WHERE(v_line1 GT MEAN(grid(1,index1))-0.5 $
             AND v_line1 LE MEAN(grid(1,index1))+0.5,count1)
   i2=WHERE(v_line2 GT MEAN(grid(1,index1))-0.5 $
             AND v_line2 LE MEAN(grid(1,index1))+0.5,count2)
   i3=WHERE(v_line3 GT MEAN(grid(1,index1))-0.5 $
             AND v_line3 LE MEAN(grid(1,index1))+0.5,count3)
   i4=WHERE(v_line4 GT MEAN(grid(1,index1))-0.5 $
             AND v_line4 LE MEAN(grid(1,index1))+0.5,count4)
   i5=WHERE(v_line5 GT MEAN(grid(1,index1))-0.5 $
             AND v_line5 LE MEAN(grid(1,index1))+0.5,count5)
   i6=WHERE(v_line(*,6) GT MEAN(grid(1,index1))-0.5 $
            AND v_line(*,6) LE MEAN(grid(1,index1))+0.5,count6)
   i7=WHERE(v_line(*,7) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,7) LE MEAN(grid(1,index1))+0.5,count7)
   i8=WHERE(v_line(*,8) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,8) LE MEAN(grid(1,index1))+0.5,count8)
   i9=WHERE(v_line(*,9) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,9) LE MEAN(grid(1,index1))+0.5,count9)
   i10=WHERE(v_line(*,10) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,10) LE MEAN(grid(1,index1))+0.5,count10)
   i11=WHERE(v_line(*,11) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,11) LE MEAN(grid(1,index1))+0.5,count11)
   i12=WHERE(v_line(*,12) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,12) LE MEAN(grid(1,index1))+0.5,count12)
   i13=WHERE(v_line(*,13) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,13) LE MEAN(grid(1,index1))+0.5,count13)
   i14=WHERE(v_line(*,14) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,14) LE MEAN(grid(1,index1))+0.5,count14)
   i15=WHERE(v_line(*,15) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,15) LE MEAN(grid(1,index1))+0.5,count15)
   i16=WHERE(v_line(*,16) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,16) LE MEAN(grid(1,index1))+0.5,count16)
   i17=WHERE(v_line(*,17) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,17) LE MEAN(grid(1,index1))+0.5,count17)
   i18=WHERE(v_line(*,18) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,18) LE MEAN(grid(1,index1))+0.5,count18)
   i19=WHERE(v_line(*,19) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,19) LE MEAN(grid(1,index1))+0.5,count19)
   i20=WHERE(v_line(*,20) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,20) LE MEAN(grid(1,index1))+0.5,count20)
   i21=WHERE(v_line(*,21) GT MEAN(grid(1,index1))-0.5 $
             AND v_line(*,21) LE MEAN(grid(1,index1))+0.5,count21)
;print,count0
;print,count1
;print,count2
;print,count3
;print,count4
;print,count5
;print,count6
;print,count7
;print,count8
;print,count9
;read,x,prompt="counts^"

;Delta(y) is the projection of the velocity element along the axis of
;the slit. rp is the projected distance from the star to the disk on
;teh plane of the sky so that rp=ra*SQRT(cos(phase)^2+sin(phase)^2*cos(inc)^2)
;where ra is the distance from the star to the annulus, inc is the
;inclination of the disk, and phase is hte phase of the orbit.
;Delta(y) = rp*cos(theta+omega) where omega is the angle between the
;axis of the slit and the semimajor axis of the disk and theta is the
;projection of phi. 

   phi=grid(2,index1)   
   rp =grid(0,index1)*SQRT(cos(phi)^2+sin(phi)^2*cos(inc)^2) 
             ;projection of distance to annulus on plane of the sky
   theta=fltarr(N_ELEMENTS(index1))

   FOR i=0,N_ELEMENTS(phi)-1 DO BEGIN
      IF phi(i) LE !pi THEN theta(i)=ACOS( COS(phi(i)) /SQRT(COS(phi(i))^2 $
          + sin(phi(i))^2*COS(inc)^2) ) $
	ELSE theta(i)=2.*!pi - ACOS( COS(phi(i)) /SQRT(COS(phi(i))^2 $
                       + SIN(phi(i))^2*COS(inc)^2) )
   ENDFOR
   deltay=rp*cos(theta+omega)

   centroid(i0)=TOTAL(grid(4,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(4,index1)*grid(3,index1))+Lc)
   centroid(i1)=TOTAL(grid(5,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(5,index1)*grid(3,index1))+Lc)
   centroid(i2)=TOTAL(grid(6,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(6,index1)*grid(3,index1))+Lc)
   centroid(i3)=TOTAL(grid(7,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(7,index1)*grid(3,index1))+Lc)
   centroid(i4)=TOTAL(grid(8,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(8,index1)*grid(3,index1))+Lc)
   centroid(i5)=TOTAL(grid(9,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(9,index1)*grid(3,index1))+Lc)
   centroid(i6)=TOTAL(grid(10,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(10,index1)*grid(3,index1))+Lc)
   centroid(i7)=TOTAL(grid(11,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(11,index1)*grid(3,index1))+Lc)
   centroid(i8)=TOTAL(grid(12,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(12,index1)*grid(3,index1))+Lc)
   centroid(i9)=TOTAL(grid(13,index1)*grid(3,index1)*deltay) $
                /(TOTAL(grid(13,index1)*grid(3,index1))+Lc)
   centroid(i10)=TOTAL(grid(14,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(14,index1)*grid(3,index1))+Lc)
   centroid(i11)=TOTAL(grid(15,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(15,index1)*grid(3,index1))+Lc)
   centroid(i12)=TOTAL(grid(16,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(16,index1)*grid(3,index1))+Lc)
   centroid(i13)=TOTAL(grid(17,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(17,index1)*grid(3,index1))+Lc)
   centroid(i14)=TOTAL(grid(18,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(18,index1)*grid(3,index1))+Lc)
   centroid(i15)=TOTAL(grid(19,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(19,index1)*grid(3,index1))+Lc)
   centroid(i16)=TOTAL(grid(20,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(20,index1)*grid(3,index1))+Lc)
   centroid(i17)=TOTAL(grid(21,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(21,index1)*grid(3,index1))+Lc)
   centroid(i18)=TOTAL(grid(22,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(22,index1)*grid(3,index1))+Lc)
   centroid(i19)=TOTAL(grid(23,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(23,index1)*grid(3,index1))+Lc)
   centroid(i20)=TOTAL(grid(24,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(24,index1)*grid(3,index1))+Lc)
   centroid(i21)=TOTAL(grid(25,index1)*grid(3,index1)*deltay) $
            /(TOTAL(grid(25,index1)*grid(3,index1))+Lc)
no_vel_elements:
ENDFOR   

read,x,prompt="enter for centroid:"
print,centroid
read,x,prompt="centroid above"
FOR j=0,n_rings-1 DO total_spec(*,j)=total_spec(*,j)* $
            total(flux_tot_slit(*,j))/total(total_spec(*,j))

final_spec=total(total_spec,2)
inst_prof=final_spec
inst_prof(*)=0.
inst_prof=exp(-(vfreq)^2/(inst_res/1.665)^2)
inst_prof=SHIFT(inst_prof,FIX(N_ELEMENTS(vfreq)/2.))
inst_prof=inst_prof/total(inst_prof) 			;normalize profile
conv_spec=FFT(FFT(final_spec)*FFT(inst_prof),1)/2.	;account for reflection in FFT
cent_conv=FFT(FFT(centroid)*FFT(inst_prof),1)
cent_conv=cent_conv*total(abs(centroid))/total(abs(cent_conv))
conv_spec=conv_spec*total(flux_tot_slit)/total(conv_spec)

print,"cent_conv"
print,cent_conv
print,"conv_spec"
print,conv_spec
print,size(total_spec)
plot,conv_spec
plot,INTERPOL(conv_spec,freq-freq*5./2.9979e5,fbig)
diff_array(*,bigi)=((rbig_masked-1)*1.5e-12 - INTERPOL(conv_spec,freq-freq*5./2.9979e5,fbig))
index_diff=WHERE(FINITE(diff_array(*,bigi)) EQ 1)

IF TOTAL(ABS(diff_array(index_diff,bigi))) LT best_diff THEN BEGIN
   best_diff=TOTAL(ABS(diff_array(index_diff,bigi)))
   best_i=bigi
   cent_conv_best=cent_conv
   conv_spec_best=conv_spec
ENDIF

IF (bigi MOD 10) EQ 0 THEN $
   PRINT, FORMAT='("You are now at step",I5," out of",I5,".")',bigi,num_guesses-1

ENDFOR

diff_big=FLTARR(num_guesses)
chisq_big=diff_big

FOR i=0, num_guesses-1 DO BEGIN
   diff_big(i)=TOTAL(ABS(diff_array(index_diff,i)))
   chisq_big(i)=TOTAL(diff_array(index_diff,i)^2/1.4d-27)
ENDFOR

m8          =matrix
c8          =chisq_big
darr8       =diff_array
d8          =diff_big

save, m8,    $
      c8,    $
      darr8, $
      d8,    $
      filename='set8b.dat'


end_time=systime(/seconds)
hours=FIX((end_time-start_time)/3600.)
minutes=FIX((end_time-start_time)/60.-60*hours)
seconds=(end_time-start_time)-60*minutes-3600*hours
tot_time=STRCOMPRESS(STRING(hours)+' h '+STRING(minutes)+' m '+STRING(seconds)+' s ')
PRINT,'Time Elapsed: '+STRING(tot_time)
END

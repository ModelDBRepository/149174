TITLE ichan2.mod  
 
COMMENT
konduktivitas valtozas hatasa- somaban 

EAT 14Sep09 Seperated Na into a different file

ENDCOMMENT
 
UNITS {
        (mA) =(milliamp)
        (mV) =(millivolt)
        (uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
? interface 
NEURON { 
SUFFIX ichan2 
USEION kf READ ekf WRITE ikf  VALENCE 1
USEION ks READ eks WRITE iks  VALENCE 1
NONSPECIFIC_CURRENT il 
RANGE  gkf, gks
RANGE gkfbar, gksbar
RANGE gl, el
RANGE nfinf, nftau, ikf, nsinf, nstau, iks
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
        ekf  (mV)
	gkfbar (mho/cm2)
        eks  (mV)
	gksbar (mho/cm2)
	gl (mho/cm2)    
 	el (mV)
}
 
STATE {
	nf ns
}
 
ASSIGNED {
         
        gkf (mho/cm2)
        gks (mho/cm2)

        ikf (mA/cm2)
        iks (mA/cm2)


	il (mA/cm2)

	nfinf nsinf
 	nftau (ms) nstau (ms)
	nfexp nsexp
} 

? currents
BREAKPOINT {
	SOLVE states
        gkf = gkfbar*nf*nf*nf*nf
        ikf = gkf*(v-ekf)
        gks = gksbar*ns*ns*ns*ns
        iks = gks*(v-eks)

	il = gl*(v-el)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	
      nf = nfinf
      ns = nsinf

}

? states
PROCEDURE states() {	:Computes state variables m, h, and n 
        trates(v)	:      at the current v and dt.
        nf = nf + nfexp*(nfinf-nf)
        ns = ns + nsexp*(nsinf-ns)
}
 
LOCAL q10

? rates
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
       q10 = 3^((celsius - 6.3)/10)
                :"m" sodium activation system - act and inact cross at -40
             :"ns" sKDR activation system
        alpha = -0.028*vtrap((v+65-35),-6)
	beta = 0.1056/exp((v+65-10)/40)
	sum = alpha+beta        
	nstau = 1/sum      nsinf = alpha/sum
            :"nf" fKDR activation system
        alpha = -0.07*vtrap((v+65-47),-6)
	beta = 0.264/exp((v+65-22)/40)
	sum = alpha+beta        
	nftau = 1/sum      nfinf = alpha/sum
	
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
        TABLE nfinf, nfexp, nsinf, nsexp, nftau, nstau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

	       tinc = -dt * q10
	nfexp = 1 - exp(tinc/nftau)
	nsexp = 1 - exp(tinc/nstau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON


#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Global Fit 2>

Function calc_all (w)		//Calculate all of the FRET variables in the graphs (phase/modulation vs freq, FRET eff vs. r, dist histogram, phasor plot) using parameter wave w
	wave w
	wave phase = root:phase		//Declare all needed global variables
	wave modulation = root:modulation
	wave PhaseDonly = root:PhaseDonly
	wave ModulationDonly = root:ModulationDonly
	wave nf = root:nf
	wave df = root:df
	wave nfD = root:nfD
	wave dfD = root:dfD
	wave effr = root:effr
	wave forster = root:forster
	wave rR0 = root:rR0
	wave freq = root:freq
	wave GFit_eff2=root:GFit_eff2
	wave GFit_eff3=root:GFit_eff3

	k1=w[0]		//save current frac_D in parameters wave
	w[0]=1			//set frac_D to 100% donor-only calculations
	PhaseDonly=phase_func(w,x)
	ModulationDonly=mod_func(w,x)
	nfD=nf_func(w, freq)
	dfD=df_func(w, freq)
	w[0]=k1		//restore frac_D to parameters wave
	If (w[7]==0)	//If r1 is the only Gaussian
		GFit_eff2=eff_func(w)	//Set value for efficiency of 2 in plot
	Elseif (w[7]==1)	//If r2 is the only Gaussian
		GFit_eff3=eff_func(w)	//Set value for efficiency of 3 in plot
	Else
		GFit_eff2=NAN
		GFit_eff3=NAN
	Endif
	effr=effr_func(w,x)
	phase=phase_func(w,x)
	modulation=mod_func(w,x)
	nf=nf_func(w, freq)
	df=df_func(w, freq)
	forster=1/(1+(x/w[4])^6)
	rR0=x/w[4]/10
End

Function phase_func (w,f) : FitFunc		//Calculate the phase for frequency f using parameter wave w
	wave w; variable f
	wave avg_r1 = root:avg_r1
	wave avg_r2 = root:avg_r2
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = (atan(nf_func(w, f)/df_func(w, f))+2*pi*f*w[10])*180/pi
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ f
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = t0
	//CurveFitDialog/ w[11] = frac_back
	
	avg_r1=w[5]	//extract avg_r's for distance dependence plot
	avg_r2=w[8]
	return (atan2(nf_func(w, f),df_func(w, f))+2*pi*f*w[10])*180/pi		//Correct phase for time shift in the IRF, t0 (w[10])
End

Function mod_func (w,f) : FitFunc		//Calculate the modulation for frequency f using parameter wave w
	wave w; variable f
	wave avg_r1 = root:avg_r1
	wave avg_r2 = root:avg_r2
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = (nf_func(w, f)^2+df_func(w, f)^2)^0.5
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ f
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = t0
	//CurveFitDialog/ w[11] = frac_back
	
	avg_r1=w[5]	//extract avg_r's for distance dependence plot
	avg_r2=w[8]
	return (nf_func(w, f)^2+df_func(w, f)^2)^0.5
End

Function nf_func (w,f) : FitFunc		//Calculate the imaginary component (N or S) of the fluorescence response for frequency f using parameter wave w
	wave w; variable f
	wave pr = root:pr 			//Declare all needed global variables
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave nfr = root:nfr
	wave phaseB = root:phaseB
	wave modulationB = root:modulationB

	k0=(1-eff_func(w)+w[11]*eff_func(w))	//Correction factor to maintain frac_back the same for different FRET efficiencies

	pr=(1-w[7])*((1-w[0])*gauss(x, w[5], w[6])+w[0]*gauss(x,150,0.1))	//Apo Gaussian distance distribution with a component at 150 A which will not exhibit any FRET, i.e. donor-only
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	nfr=pr(x)*2*pi*f*1e6*tda1(x)^2/(1+(2*pi*f*1e6)^2*tda1(x)^2)	//distance dependence of N (N(r) to be integrated over r to get N
	
	pr=w[7]*((1-w[0])*(gauss(x, w[8], w[9]))+w[0]*gauss(x,150,0.1))	//Holo Gaussian distance distribution with a component at 150 A which will not exhibit any FRET, i.e. donor-only
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	nfr=nfr+pr(x)*2*pi*f*1e6*tda2(x)^2/(1+(2*pi*f*1e6)^2*tda2(x)^2)	//distance dependence of N (N(r) to be integrated over r to get N

	
	return (1-w[11]/k0)*area(nfr,0,inf)/jDA(w)+w[11]/k0*modulationB(f)*sin(phaseB(f)*pi/180)	//Correct N for the N of a background sample
End

Function df_func (w,f) : FitFunc		//Calculate the real component (D or G) of the fluorescence response for frequency f using parameter wave w
	wave w; variable f
	wave pr = root:pr		//Declare all needed global variables
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave dfr = root:dfr
	wave phaseB = root:phaseB
	wave modulationB = root:modulationB

	k0=(1-eff_func(w)+w[11]*eff_func(w))	//Correction factor to maintain frac_back the same for different FRET efficiencies
	
	pr=(1-w[7])*((1-w[0])*gauss(x, w[5], w[6])+w[0]*gauss(x,150,0.1))	//Apo Gaussian distance distribution with a component at 150 A which will not exhibit any FRET, i.e. donor-only
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	dfr=pr(x)*tda1(x)/(1+(2*pi*f*1e6)^2*tda1(x)^2)	//Distance dependence of D (D(r) to be integrated over r to get D
	
	pr=w[7]*((1-w[0])*(gauss(x, w[8], w[9]))+w[0]*gauss(x,150,0.1))	//Holo Gaussian distance distribution with a component at 150 A which will not exhibit any FRET, i.e. donor-only
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	dfr=dfr+pr(x)*tda2(x)/(1+(2*pi*f*1e6)^2*tda2(x)^2)	//Distance dependence of D (D(r) to be integrated over r to get D

	return (1-w[11]/k0)*area(dfr,0,inf)/jDA(w)+w[11]/k0*modulationB(f)*cos(phaseB(f)*pi/180) //Correct D for the D of a background sample
End

Function jDA (w) : FitFunc		//Calculate normalizaton factor J for the donor/acceptor sample using parameter wave w
	wave w
	wave pr = root:pr		//Declare all needed global variables
	wave tda1 = root:tda1
	wave tda2 = root:tda2
	wave dfr = root:dfr
	
	pr=(1-w[7])*((1-w[0])*gauss(x, w[5], w[6])+w[0]*gauss(x,150,0.1))	//Apo Gaussian distance distribution with a component at 150 A which will not exhibit any FRET, i.e. donor-only
	tda1=1/(1/w[1]+1/w[1]*(w[4]/x)^6)
	dfr=pr(x)*tda1(x)	//Distance dependence of D (D(r) to be integrated over r to get D
	
	pr=w[7]*((1-w[0])*(gauss(x, w[8], w[9]))+w[0]*gauss(x,150,0.1))	//Holo Gaussian distance distribution with a component at 150 A which will not exhibit any FRET, i.e. donor-only
	tda2=1/(1/w[3]+1/w[3]*(w[4]/x)^6)
	dfr=dfr+pr(x)*tda2(x)	//Distance dependence of D (D(r) to be integrated over r to get D

	return area(dfr,0,inf)	//J is just the unnormalized D at freq=0
End

Function jD (w) : FitFunc		//Calculate normalizaton factor J for the donor-only sample using parameter wave w
	wave w
	
	return (1-w[7])*w[1]+w[7]*w[3]
End

Function eff_func (w) : FitFunc	//Calculate the FRET efficiency as measured from the intensity using parameter wave w
	wave w
	
	return 1-jDA(w)/jD(w)
End

Function effnofracD_func (w,c) : FitFunc	//Calculate the FRET efficiency as measured from the intensity using parameter wave w without frac_D for fitting to SS FRET data
	wave w; variable c	//need vaiable c to use as a fitting function
	wave parameters = root:parameters
		
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = (nf_func(w, f)^2+df_func(w, f)^2)^0.5
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ f
	//CurveFitDialog/ Coefficients 12
	//CurveFitDialog/ w[0] = frac_D
	//CurveFitDialog/ w[1] = tauD1
	//CurveFitDialog/ w[2] = alphaD1
	//CurveFitDialog/ w[3] = tauD2
	//CurveFitDialog/ w[4] = r0
	//CurveFitDialog/ w[5] = avg_r1
	//CurveFitDialog/ w[6] = sig_r1
	//CurveFitDialog/ w[7] = frac_r2
	//CurveFitDialog/ w[8] = avg_r2
	//CurveFitDialog/ w[9] = sig_r2
	//CurveFitDialog/ w[10] = t0
	//CurveFitDialog/ w[11] = frac_back
	
	parameters=w
	parameters[0]=0
	return 1-jDA(parameters)/jD(parameters)
End

Function effr_func (w,c) : FitFunc	//Calculate the FRET efficiency at different values of avg_r1 (c) as measured from the intensity using parameter wave w
	wave w; variable c
	
	If (w[7]==0)	//If r1 is the only Gaussian
		k2=w[5]	//Save current avg_r1 in parameter wave
		w[5]=c		//Set avg_r1 to c
	Elseif (w[7]==1)	//If r2 is the only Gaussian
		k2=w[8]	//Save current avg_r2 in parameter wave
		w[8]=c		//Set avg_r2 to c
	Else
		return NAN
	Endif

	k1=1-jDA(w)/jD(w)
	
	If (w[7]==0)	//If r1 is the only Gaussian
		w[5]=k2	//Restore avg_r1 to parameter wave
	Elseif (w[7]==1)	//If r2 is the only Gaussian
		w[8]=k2	//Restore avg_r2 to parameter wave
	Endif
	
	return k1
End

Function background_correct (w, phase, modulation) //Correct the phase and modulation waves for the time shift (t0) and background in parameter wave w
	wave w, phase, modulation //Requires that phase and modulation waves have a scaled x value
	wave phaseB = root:phaseB		//Declare all needed global variables
	wave modulationB = root:modulationB
	duplicate/O phase, tempphase //Make a duplicate phase wave
	duplicate/O modulation,tempmod //Make a duplicate modulation wave

	k0=(1-eff_func(w)+w[11]*eff_func(w))	//Correction factor to maintain frac_back the same for different FRET efficiencies
	tempphase=(atan2(modulation(x)*sin(phase(x)*pi/180-2*pi*x*w[10])-w[11]/k0*modulationB(x)*sin(phaseB(x)*pi/180),modulation(x)*cos(phase(x)*pi/180-2*pi*x*w[10])-w[11]/k0*modulationB(x)*cos(phaseB(x)*pi/180)))*180/pi
	tempmod=1/(1-w[11]/k0)*((modulation(x)*sin(phase(x)*pi/180-2*pi*x*w[10])-w[11]/k0*modulationB(x)*sin(phaseB(x)*pi/180))^2+(modulation(x)*cos(phase(x)*pi/180-2*pi*x*w[10])-w[11]/k0*modulationB(x)*cos(phaseB(x)*pi/180))^2)^0.5
	phase=tempphase 		//Overwrite previous phase wave with corrected phase
	modulation=tempmod		//Overwrite previous modulation wave with the corrected modulation
End

Menu "Macros"		//Put macros in the Macro menu
	"Calc_all 1", calc_all(Coef_Modulation1)
	"Calc_all 2", calc_all(Coef_Modulation2)
	"Calc_all 3", calc_all(Coef_Modulation3)
	"Calc_all 4", calc_all(Coef_Modulation4)
	"Calc_all 5", calc_all(Coef_Modulation5)
	"Calc_all 6", calc_all(Coef_Modulation6)
	"FRET_efficiency 1", print effnofracD_func(coef_Phase1, 0)
	"FRET_efficiency 2", print effnofracD_func(coef_Phase2, 0)
	"FRET_efficiency 3", print effnofracD_func(coef_Phase3, 0)
	"FRET_efficiency 4", print effnofracD_func(coef_Phase4, 0)
	"FRET_efficiency 5", print effnofracD_func(coef_Phase5, 0)
	"FRET_efficiency 6", print effnofracD_func(coef_Phase6, 0)
	"Background_Correct",SetScale/P x X_Value1[0],X_Value1[1]-X_Value1[0],"", Modulation1,Modulation2,Modulation3,Modulation4,Modulation5,Modulation6,Phase1,Phase2,Phase3,Phase4,Phase5,Phase6;background_correct(coef_Modulation1,Phase1,Modulation1);background_correct(coef_Modulation1,GFit_Phase1,GFit_Modulation1);background_correct(coef_Modulation2,phase2,Modulation2);background_correct(coef_Modulation2,GFit_phase2,GFit_Modulation2);background_correct(coef_modulation3,phase3,modulation3);background_correct(coef_modulation3,GFit_phase3,GFit_modulation3);background_correct(coef_modulation4,phase4,modulation4);background_correct(coef_modulation4,GFit_phase4,GFit_modulation4);background_correct(coef_Modulation5,Phase5,Modulation5);background_correct(coef_Modulation5,GFit_Phase5,GFit_Modulation5);background_correct(coef_modulation6,phase6,modulation6);background_correct(coef_modulation6,GFit_phase6,GFit_modulation6)
	"Set scale BG",SetScale/P x X_ValueB[0],X_ValueB[1]-X_ValueB[0],"", PhaseB,ModulationB
End

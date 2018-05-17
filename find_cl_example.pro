;Example code to show how find_cl.py should be used to calculate confidence
;levels for data with a power-law power spectrum and uncertainties.


PRO find_cl_example

;Generate a synthetic red noise time series
	;Number of data points in time series
		npts = 200

	;Generate a white noise time series
		flux = RANDOMN(seed, npts) ;Normally distributed white noise
		time = findgen(npts)

	;Convert into a red noise time series
	for i=0, npts-2 do begin
  		flux[i+1] += 0.99*flux[i]
	endfor

;Calculate periodogram
	result = LNP_TEST(time, flux, wk1=freq, wk2=pow, ofac=1)
	npow = N_ELEMENTS(pow)

;Switch to log space, then power law is a straight line
	pow = ALOG10(pow)
	freq = ALOG10(freq)

;Fit power spectrum with a broken power law model
	err = REPLICATE(1., npow)
	p = [2.D, -1.D, 0.1D, 1.D] ;Fit initial guesses
	result = mpfitfun('broke_pow_law', freq, pow, err, p, yfit=pow_fit, /quiet, /NAN)

;Array containing some made-up uncertainties on the fitted power law model at each frequency index. For real data these could be estimated by performing Monte Carlo simulations
	pow_fit_err = REPLICATE(0.1, npow)

;Possible issues when uncertainties are too small, so limit how small they can be
	bad = WHERE(pow_fit_err LT 1.E-10, count)
	IF count NE 0 THEN pow_fit_err[bad] = 1.E-10

;Print to file so Python can work out confidence levels
	OPENW, 1, 'pow_fit_err.dat'
	PRINTF, 1, '# log Power fit error'
	FOR i=0, npow-1 DO BEGIN
 		PRINTF, 1, pow_fit_err[i], FORMAT='(E0)'
	ENDFOR
	CLOSE, 1

;Run Python script to calculate confidence levels
	SPAWN, 'python find_cl.py'

;Read result from Python:
	data = READ_ASCII('cl95.dat', COMMENT_SYMBOL='#')
	cl95 = data.field1
	data = READ_ASCII('cl99.dat', COMMENT_SYMBOL='#')
	cl99 = data.field1

;Normalise
	cl95 *= 0.5*MEAN(10.^(pow-pow_fit))
	cl99 *= 0.5*MEAN(10.^(pow-pow_fit))
	cl95 = ALOG10(cl95)+pow_fit
	cl99 = ALOG10(cl99)+pow_fit

	LOADCT,39
	WINDOW, 0
	PLOT, freq, pow, xstyle=1, ystyle=3, ytitle='log Power', xtitle='log Frequency'
	OPLOT, freq, pow_fit, color=250, thick=2
	OPLOT, freq, cl99, color=250, thick=2, linestyle=2
	OPLOT, freq, cl95, color=250, thick=2, linestyle=1

END

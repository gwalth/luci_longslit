;+
; NAME:
;   radec_greatcircle
;
; PURPOSE:
;   Solve for the great circle for a set of RA,DEC positions at a set of times.
;
; CALLING SEQUENCE:
;   fitval = radec_greatcircle(ralist, declist, xposlist, yposlist, $
;    [ timelist, start_parms, maxiter=, fixed=, muerr=, nuerr=, $
;    /debug, niter=, status=, muresid=, nuresid= ] )
;
; INPUTS:
;   ralist     - RA coordinates (degrees) [NPOINTS]
;   declist    - DEC coordinates (degrees) [NPOINTS]
;   xposlist   - X positions, perpendicular to scan direction [NPOINTS]
;   yposlist   - Y positions, along the scan direction [NPOINTS]
;
; OPTIONAL INPUTS:
;   timelist   - Time stamps.  If specified, then a time^2 term is
;                also fit along the scan direction, for a total of
;                6 fit parameters instead of 5 [NPOINTS].
;   start_parms- Initial guess for the great circle fit [NPARM].
;   maxiter    - Maximum number of iterations for fit; default to 1000
;   fixed      - If set, then fit for each parameter where this is 0,
;                and fit for each parameter where this is 1 [NPARM].
;   muerr      - Error of each point in MU coordinate; default to 1 arcsec;
;                either a scalar or of length [NPOINTS]; used for chi^2
;   nuerr      - Error of each point in NU coordinate; default to 1 arcsec;
;                either a scalar or of length [NPOINTS]; used for chi^2
;   debug      - If set, then plot the MU,NU deviations in units of arcsec
;                at each iteration of the fit.
;
; OUTPUTS:
;   fitval     - Best-fit parameters [NPARM].
;
; OPTIONAL OUTPUTS:
;   niter      - Number of iterations performed.
;   status     - Return status from MPFIT() function, set to a non-zero
;                value if an error occurred.
;   muresid    - MU coordinate residuals from best-fit
;   nuresid    - NU coordinate residuals from best-fit
;
; COMMENTS:
;   The fit parameters are as follows:
;     FITVAL[0] - Starting MU coordinate [degrees]
;     FITVAL[1] - Node of great circle [degrees]
;     FITVAL[2] - Inclination of great circle [degrees]
;     FITVAL[3] - Perpendicular offset of scan from the great circle [degrees];
;                 this is what is called XBORE for SDSS imaging scans.
;                 Fix this parameter to zero if you want to fit exactly along
;                 a great circle, and not offset from it.
;     FITVAL[4] - Tracking rate in degrees per unit of YPOSLIST;
;                 for example, if YPOSLIST is is units of pixels,
;                 this is in units of degrees/pixel
;     FITVAL[5] - Optional parameter for fitting a quadratic term in time
;
;   The set of parametric equations to convert from xposlist,yposlist,timelist
;   to ra,dec are as follows:
;     mu = FITVAL[0] + yposlist * FITVAL[4] + FITVAL[5] * timelist^2
;     nu = xposlist + FITVAL[3]
;     xx = cos(mu) * cos(nu)
;     yy = sin(mu) * cos(nu) * cos(FITVAL[2]) - sin(nu) * sin(FITVAL[2])
;     zz = sin(mu) * cos(nu) * sin(FITVAL[2]) + sin(nu) * cos(FITVAL[2])
;     ra = (180/!pi) * atan(yy,xx) + FITVAL[1]
;     dec = (180/!pi) * asin(zz)
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   cirrange
;   djs_plot
;   mpfit()
;   radec_to_munu
;
; INTERNAL SUPPPORT ROUTINES:
;   radec_gcfn()
;
; REVISION HISTORY:
;   21-Nov-2002  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function radec_gcfn, value, ralist=ralist, declist=declist, $
 xposlist=xposlist, yposlist=yposlist, timelist=timelist, $
 muerr=muerr, nuerr=nuerr, debug=debug, muresid=muresid, nuresid=nuresid

   mu_ref = value[0]
   node = value[1]
   incl = value[2]
   xbore = value[3]
   tracking = value[4]
   if (n_elements(value) GT 5) then qterm = value[5] $
    else qterm = 0

   ; Compute the mu,nu coordinates by projecting RA,DEC to the
   ; great circle specified by NODE,INCL.
   radec_to_munu, ralist, declist, mulist1, nulist1, node=node, incl=incl

   mulist2 = yposlist * tracking + mu_ref + qterm * timelist^2
   cirrange, mulist2
   nulist2 = xposlist + xbore

   muresid = mulist2 - mulist1
   muresid = muresid - 360. * (muresid GT 180) + 360. * (muresid LT -180)
   nuresid = nulist2 - nulist1

   muresid = muresid * 3600.d0
   nuresid = nuresid * 3600.d0

   muchi = muresid / muerr
   nuchi = nuresid / nuerr
   chivec = sqrt(muchi^2 + nuchi^2)

   if (keyword_set(debug)) then begin
      !p.multi = [0,1,3]
      xplot = yposlist * tracking
      djs_plot, xplot, muresid, psym=-4, charsize=2, $
       ytitle='\Delta\mu'
      djs_plot, xplot, nuresid, psym=-4, charsize=2, $
       ytitle='\Delta\nu'
      djs_plot, minmax(xplot), minmax([muchi,nuchi]), /nodata, charsize=2, $
       xtitle='(Y position) * Tracking', ytitle='\chi(\mu), \chi(\nu)'
      oplot, xplot, muchi, psym=-4
      oplot, xplot, nuchi, psym=-4
   endif

   return, chivec
end
;------------------------------------------------------------------------------
function radec_greatcircle, ralist, declist, xposlist, yposlist, timelist, $
 start_parms, maxiter=maxiter, fixed=fixed, debug=debug, $
 muerr=muerr, nuerr=nuerr, niter=niter, status=status, $
 muresid=muresid, nuresid=nuresid

   ;----------
   ; Set defaults

   if (NOT keyword_set(maxiter)) then maxiter = 1000
   if (NOT keyword_set(muerr)) then muerr = 1.0
   if (NOT keyword_set(nuerr)) then nuerr = 1.0

   npts = n_elements(ralist)
   if (n_elements(declist) NE npts) then $
    message, 'Number of elements in RALIST,DECLIST disagree'
   if (n_elements(xposlist) NE npts) then $
    message, 'Number of elements in RALIST,XPOSLIST disagree'
   if (n_elements(yposlist) NE npts) then $
    message, 'Number of elements in RALIST,YPOSLIST disagree'
   if (keyword_set(timelist)) then $
    if (n_elements(timelist) NE npts) then $
     message, 'Number of elements in RALIST,TIMELIST disagree'
   if (n_elements(muerr) GT 1 AND n_elements(muerr) NE npts) then $
    message, 'Number of elements in RALIST,MUERR disagree'
   if (n_elements(nuerr) GT 1 AND n_elements(nuerr) NE npts) then $
    message, 'Number of elements in RALIST,NUERR disagree'

   if (keyword_set(timelist)) then nparm = 6 $
    else nparm = 5
   parinfo = {value: 0.D, fixed: 0}
   parinfo = replicate(parinfo, nparm)

   if (keyword_set(start_parms)) then begin
      if (n_elements(start_parms) NE nparm) then $
       message, 'Wrong number of elements in START_PARMS'
      parinfo.value = start_parms
   endif
   if (keyword_set(fixed)) then begin
      if (n_elements(fixed) NE nparm) then $
       message, 'Wrong number of elements in FIXED'
      parinfo.fixed = fixed
   endif

   ;----------
   ; Do the fit

   functargs = { ralist: ralist, declist: declist, $
    xposlist: xposlist, yposlist: yposlist, timelist: timelist, $
    debug: keyword_set(debug), muerr: muerr, nuerr: nuerr }
   fitval = mpfit('radec_gcfn', parinfo=parinfo, functargs=functargs, $
    maxiter=maxiter, niter=niter, status=status)

   ;----------
   ; Make one last call to the function in order to return the residuals

   if (arg_present(muresid) OR arg_present(nuresid)) then begin
      junk = radec_gcfn(fitval, _EXTRA=functargs, $
      debug=0, muresid=muresid, nuresid=nuresid)
   endif

   return, fitval
end
;------------------------------------------------------------------------------

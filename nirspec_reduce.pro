;-------------------------------------------------------------------------

Function nr_makevariance,data,object=object,sky=sky,flatfield=flatfield,$
                         flatvar=flatvar,darkframe=darkframe,$
                         objcoadds=objcoadds,darkcoadds=darkcoadds,$
                         readnoise=readnoise,$
                         wave=wave

; Calculate variance array based on object counts, sky counts, flat
; field value, flat field variance, dark counts, and read noise.  All
; inputs are optional keywords.  Variables are as defined in comments
; below, except for 'data' which is the array for which the variance
; is to be computed (and is used only to count the number of pixels).
; At present, a gain of 5 e/ADU is assumed for calculating Poisson
; uncertainties.

gain = 5.
if keyword_set(readnoise) then rnoise = float(readnoise) else rnoise = 0.
rnoisevar = (rnoise^2)/gain ; Read noise variance in ADU

var  = 0.*data

if keyword_set(flatfield) then ffield = flatfield else ffield = 1. + 0.*data
; Object
if keyword_set(object) then begin
   objvar = (ffield^2)*abs(object)/gain
   print,'Object:',minmax(objvar)
   var = var + objvar
endif else object = 0.*data
; Sky
if keyword_set(sky) then begin
   skyvar = (ffield^2)*abs(sky)/gain
   print,'Sky:',minmax(skyvar)
   var = var + skyvar
endif else sky = 0.*data
; Flat
if keyword_set(flatvar) then begin
   flatvar = (object^2 + sky^2)*flatvar
   print,'Flat',minmax(flatvar)
   var = var + flatvar
endif
; Read noise - object frame
if keyword_set(objcoadds) then ocoadds = objcoadds else ocoadds = 1.
var = var + ocoadds*rnoisevar
; Dark - assume readnoise dominated
if keyword_set(darkcoadds) then dkcoadds = darkcoadds else dkcoadds = ocoadds
;if keyword_set(darkframe) then var = var + dkcoadds*rnoisevar

return,var
end

;-----------------------------------------------------------------------
;-----------------------------------------------------------------------

Pro nirspec_reduce,objframe,wavearr,slitarr,skymodel,objmodel,$
                   skyref=skyref,flatfield=flatfield,flatvar=flatvar,$
                   darkframe=darkframe,slitpos=slitpos,noobject=noobject,$
                   nomask=nomask,norotate=norotate,objcoadds=objcoadds,$
                   darkcoadds=darkcoadds,readnoise=readnoise,texp=texp,$
		   relgain=relgain,standard=standard,wavemodel=wavemodel,$
                   writewave=writewave,waveref=waveref,linelist=linelist,$
                   outfile=outfile,noshow=noshow

;+----------------------------------------------------------------------
;
; NIRSPEC_REDUCE     5/2004
;
; Iterative procedure to fit sky background and object spectra in
; NIRSPEC longslit data.  The program performs optimal sky
; subtraction, making use of input arrays that give proxies to the
; slit position and wavelength for each pixel.  After the initial sky
; subtraction is done, the object spectrum is extracted using an
; optimal subtraction routine.  A 2D model for the object spectrum is
; then subtracted from the original exposure and the sky subtraction
; repeated.  This process is iterated several times until
; self-consistent solutions for the sky background and object spectrum
; are found.
;
; The program produces 2D model of unbinned sky background and object
; spectrum, such that
;       raw object frame = sky model + obj model + dark + noise.
;
; Sky background modeling is done using a 2D b-spline fit, where the
; slit illumination function is fit with a low-order polynomial and a
; b-spline is fit in the wavelength direction.  Fitting is done using
; BSPLINE_ITERFIT, written by D. Schlegel & S. Burles, along with
; other supporting programs from the SDSS spectral reduction pipeline.
;
; The input proxy to the slit position for a given pixel is currently
; assumed to be the y-intercept of the hypothetical object trace that
; would pass through that pixel (assuming the dispersion direction is
; roughly the x-direction).  Similarly, the proxy to the wavelength is
; the x-intercept of the curve of constant wavelength which would pass
; through that pixel.
;
; Note: The program presently assumes that all input arrays except the
; object exposure have been rotated such that the spectral dispersion
; runs roughly in the x-direction, with bluer wavelengths on the
; left.  Output files will also be rotated.
;
; INPUTS
;        objframe  - Raw 2D NIRSPEC longslit object exposure (array or
;                    fits filename).
;        wavearr   - 2D array or fits filename, same size as exposure,
;                    whose pixel values are a proxy to the wavelength
;                    of the corresponding pixels in the object 
;                    exposure. 
;        slitarr   - 2D array or fits filename, same size as exposure,
;                    whose pixel values are a proxy to the slit position
;                    of the corresponding pixels in the object 
;                    exposure. 
;
; KEYWORDS   
;        skyref    - Raw 2D NIRSPEC longslit object exposure to be
;                    used as reference for primary sky subtraction
;                    (array or fits filename).  Should be an exposure 
;                    taken near in time to objframe, and in same
;                    setup.  Asuumed to have same orientation as objframe.
;        flatfield - 2D array or fits filename of flat field.  The
;                    program currently assumes that the spectral
;                    response, rough slit illumination function, and
;                    any defects due to dust or scratched on the dewar
;                    window have been removed, and that the mean value
;                    for the flatfield is 1.  Default is not to flat-field.
;        flatvar   - 2D array or fits filename of the expected variance
;                    in the flat field.
;        darkframe - 2D array or fits filename of dark frame with same
;                    exposure time as object frame, to be subtracted
;                    from object exposure.  Default is to not dark
;                    subtract. 
;        slitpos   - Approximate slit position of object, in units
;                    used in slitarr.  Useful for extracting multiple
;                    objects.
;        noobject  - Do not search for object or attempt to extract an object
;                    spectrum.  If this keyword is set, sky subtraction will be 
;                    performed across the entire slit, rather than in a narrow 
;                    region around an object trace.
;        nomask    - Do not explicitly mask pixels prior to fitting 
;                    b-spline to the sky background (which is done
;                    iteratively, with rejection) or performing
;                    optimal extraction (which also uses rejection).
;                    Default is to create an explicit mask based on
;                    outlying pixels in the flat and dark frames, if
;                    provided.  
;        norotate  - Do not rotate (already rotated) raw object array.      
;        objcoadds - Number of coadds for object exposure (used to
;                    determine expected read noise).
;        darkcoadds- Number of coadds for dark exposure (used to
;                    determine expected read noise).  Default is
;                    same number as for object.
;        readnoise - Read noise per coadd, in rms e/pixel.  Default
;                    is 10 (measured by hand).
;	 texp      - Total exposure time (seconds, int time x no. of coadds), 
;		     currently used to estimate dark counts.  Default is
;                    0 sec.
;        relgain   - Ratio of blue side-to-red side gains.  Default is
;                    1.0343, which was measured in Nov 2003 from 1.5
;                    sec x 7 coadd flats.   
;        standard  - If set, program looks for a standard star near the 
;                    "upper" (lefthand) edge of the slit.  The sky background 
;                    is then measured mainly to the inside of the star. 
;        wavemodel - Filname of wavelength fit to be used as an initial guess
;                    for automatic line identification.  Format is a three-
;                    column ASCII file, where the columns are wavelength-proxy 
;                    (as in the wavearr input), intensity, and physical 
;                    wavelength.  If this keyword is set, the program 
;                    'wave_calib' will attempt to automatically identify each 
;                    line in the linelist file (see linelist keyword).
;        writewave - Output an ASCII file giving the wavelength solution as
;                    a funciton of wavelength-proxy.  The format is that
;                    described above for the 'wavemodel' keyword, with
;                    entries in steps of 0.1 units used for the wavelength
;                    proxy (pixels).  This file can be used in future calls of 
;                    the 'wavemodel' keyword.  Filename will be 
;                    '<outfile>_wav.dat'.
;	 waveref   - Filename or 2D array of wavelength values
;                    to be used rather than fitting skylines (typically for
;                    a standard star taken with an exposure too short to
;                    accurately measure sky lines).
;        linelist  - Name of file containing line wavelengths.  Wavelengths
;                    should be in the first column.  Other columns ignored.
;                    Lines commented out with a '#' are ignored.
;        outfile   - Root name for output files, default is 'nirspec'.
;        noshow    - Do not display successive subtraction/extraction
;                    iterations. Default is to send output to windows
;                    10 and 11.
;
; OUTPUTS
;        skymodel  - 2D b-spline fit to sky background, array or fits
;                    filename to which array is to be written.
;        objmodel  - 2D model object spectrum, array or fits filename
;                    to which array is to be written.
;
; HISTORY
;        Written 5/19/2004 GDB
;        Keywords 'wavemodel', 'writewave', 'noobject', and 'linelist' 
;           added 8/2004 GDB
;----------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE:'
   print,'   nirspec_reduce,objframe,wavearr,slitarr,skymodel,objmodel,'
   print,'                  skyref=skyref,flatfield=flatfield,flatvar=flatvar,'
   print,'                  darkframe=darkframe,slitpos=slitpos,/noobject,'
   print,'                  /nomask,/norotate,objcoadds=objcoadds,'
   print,'                  darkcoadds=darkcoadds,readnoise=readnoise,'
   print,'                  texp=texp,relgain=relgain,/standard,'
   print,'                  wavemodel=wavemodel,/writewave,waveref=waveref,'
   print,'                  linelist=linelist,outfile=outfile,/noshow'
   return
endif

currentmulti = !p.multi

;;; Set gain (e/ADU)
gain = 5.

;;; 0.  Read in input arrys, if necessary

if (size(objframe,/tname) eq 'STRING') then begin
   checkexists = findfile(objframe,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Object file not found.'
      return
   endif else begin
       rawobj = readfits(objframe,hobj)
       fromfitsfile = 1 
       sxaddpar,hobj,'ORIGEXP',objframe,'Filename of original exposure.'
   endelse
endif else begin
   rawobj = objframe
   fromfitsfile = 0
   ; Create a basic FITS header
   fxhmake,hobj,rawobj,/date
endelse

;;; Reference sky (optional)
if keyword_set(skyref) then begin
  if (size(skyref,/tname) eq 'STRING') then begin
      checkexists = findfile(skyref,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Sky reference file not found.'
         return
      endif else begin
         skyrefarr = readfits(skyref)
         sxaddpar,hobj,'SKYREF',skyref,'Reference sky exposure'
      endelse
   endif else skyrefarr = skyref
endif else skyrefarr = 0.*rawobj

; Create subtracted sky frame
subobj = rawobj - skyrefarr

; Rotate object (optional)
if not(keyword_set(norotate)) then begin
   hrotate,rawobj,hobj,rawobj,hobj,3
   subobj = rotate(subobj,3)
endif

if (size(wavearr,/tname) eq 'STRING') then begin
   checkexists = findfile(wavearr,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Wavelength file not found.'
      return
   endif else begin
      wavegrid = readfits(wavearr)
      sxaddpar,hobj,'WAVEGRID',wavearr,'Wavelength-proxy grid'
   endelse
endif else wavegrid = wavearr

if (size(slitarr,/tname) eq 'STRING') then begin
   checkexists = findfile(slitarr,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Slit position file not found.'
      return
   endif else begin
      slitgrid = readfits(slitarr)
      sxaddpar,hobj,'SLITGRID',slitarr,'Slit-position-proxy grid'
   endelse
endif else slitgrid = slitarr

;;; 1. Dark-subtract (optional)

; Note: Due to apparent structure in the read-noise, dark-subtracting
; object frames has been found to add undue noise and interfere with
; sky-subtraction.  In addition, subtracting a refernce sky frame
; effectively gets rid of any dark levels.  Therefore, in order to
; measure the raw sky image (to get variance estimates), subtract a
; uniform dark level according to the exposure time.  The dark current
; is measured by GDB to be 0.55 e/sec/pixel, which is smaller than the
; documented 0.2 e/sec/pixel.  Guess is that dark current is
; overestimated in documentation due to row-to-row structure.

if keyword_set(darkframe) then begin
   if (size(darkframe,/tname) eq 'STRING') then begin
      checkexists = findfile(darkframe,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Dark exposure file not found.'
         return
      endif else begin
         darkarr = readfits(darkframe)
         sxaddpar,hobj,'DARKEXP',darkframe,$
            'Dark exposure (for bad pixel map only)'
      endelse
   endif else darkarr = darkframe
   ;rawobj = rawobj - darkarr
   ;;; Estimate dark level by exposure time
   if keyword_set(texp) then texposure = texp else texposure = 0.
   darkcurrent = 0.55  ; e/sec/pix, measured by GDB
   edarklevel = darkcurrent * texposure
   darklevel  = edarklevel / gain   
   print,'Dark counts (ADU) = ',darklevel
   sxaddpar,hobj,'DARKCNTS',darklevel,'Dark counts subtracted'
   rawobj = rawobj - darklevel
endif else darkarr = 0.

;;; 2. Flat-field (optional)

; Adjust relative gains between blue-side and red-side amplifiers
if keyword_set(relgain) then rel_gain = relgain else $
   rel_gain = 1.0343   ; Measured (blue gain / red gain) 
rawobj(0:511,*) = rawobj(0:511,*) / rel_gain
subobj(0:511,*) = subobj(0:511,*) / rel_gain

if keyword_set(flatfield) then begin
   if (size(flatfield,/tname) eq 'STRING') then begin
      checkexists = findfile(flatfield,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Flat field file not found.'
         return
      endif else begin
         flatarr = readfits(flatfield)
         sxaddpar,hobj,'FLATFLD',flatfield,'Flat field'
      endelse
   endif else flatarr = flatfield
   ; Fix zero-valued pixels by subbing in value small enough to be
   ; included in mask
   flatzerols = where(flatarr eq 0)
   if (total(flatzerols) ne -1) then flatarr(flatzerols) = 0.01
endif else begin
   ; If no input flat field, create a dummy array full of 1s
   flatarr = 1. + 0.*rawobj
endelse
rawobj = rawobj / flatarr
subobj = subobj / flatarr

if keyword_set(flatvar) then begin
   if (size(flatvar,/tname) eq 'STRING') then begin
      checkexists = findfile(flatvar,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Flat field variance file not found.'
         return
      endif else begin
         flatvararr = readfits(flatvar)
         sxaddpar,hobj,'FLATVAR',flatvar,'Flat field variance'
      endelse
   endif else flatvararr = flatvar
endif else flatvararr = 0.*rawobj

;;; 3. Create a bad pixel map by identifying outlying pixels in the
;;; dark and flat field (optional).  Note: The values for identifying
;;; outliers are currently hard-coded in and may need to be adjusted,
;;; although the ranges of allowable pixel values is already fairly
;;; wide.

badmap = 1. + 0.*rawobj     ; Default is no masking
if not(keyword_set(nomask)) then begin
   if keyword_set(darkframe) then begin
      darkbadls = where(darkarr lt 0 or darkarr gt 200)
      badmap(darkbadls) = 0.
   endif
   if keyword_set(flatfield) then begin
      flatbadls = where(flatarr lt 0.8 or flatarr gt 1.2)
      badmap(flatbadls) = 0.
   endif
endif

;;; 4. Locate the object by performing a rough fit to a narrow region
;;; region of the sky across the entire slit.  We assume that objects
;;; other than std stars have been kept away from the edges of the
;;; slit (we reserve the left/top edge of the slit for std stars).
;;; Use the reference sky-subtracetd frame (if one has been created).
;;;
;;; Sky-subtraction will be performed in regions of the slit to either
;;; side of the object.  The user may also elecct not to mask out an
;;; object, in which case the sky will be subtracted across the entire slit.

if keyword_set(noobject) then begin   ; No object - fit sky across entire slit.

   objheight   = 0.
   objposition = 0.
   objsig      = 0.
   refheight   = 0.
   refposition = 0.
   refsig      = 0.
   fit_entire_slit = 1

endif else begin    ; Locate object and mask when fitting sky

   fit_entire_slit = 0

   ; Fit parameters
   wavemin = 250
   wavemax = 400
   if keyword_set(standard) then begin  ; Assumes standards placed at edge
      slitmin = 400
      slitmax = 525
   endif else begin  ; Avoid standard star region, which may contain persistance
      slitmin = 245
      slitmax = 465
   endelse
   findls = where(wavegrid ge wavemin and wavegrid le wavemax and $
                  slitgrid ge slitmin and slitgrid le slitmax)
   wavelist    = wavegrid(findls)
   slitlist    = slitgrid(findls)
   subobjlist  = subobj(findls)
   badlist     = badmap(findls)
   order = sort(wavelist)
   findorder   = findls(order)
   waveorder   = wavelist(order)
   slitorder   = slitlist(order)
   subobjorder = subobjlist(order)
   badorder    = badlist(order)
   invvar = 0.*subobjorder
   invvar(*) = 0.025        ; Reasonable value for H-band sky counts in 10 min
   maskls = where(badorder eq 0)      ; Obvious bad points
   invvar(maskls) = 0.
   bk_space = 0.7
   rej = 5
   n_poly = 2
   xtwo = slitorder
   max_iter = 20.
   ; Fit sky in narrow region
   if keyword_set(standard) then skyfitorder = 0 else begin
      sset=bspline_iterfit(waveorder,subobjorder,x2=xtwo,npoly=n_poly,$
                        invvar=invvar,bkspace=bk_space,yfit=skyfitorder,$
                        upper=rej,lower=rej,maxiter=max_iter)
   endelse
   noskyorder = subobjorder - skyfitorder
   
   ; Divide sky-subtracted region of the slit into several chunks and
   ; compute the median value in each chunk.
   slitlen  = slitmax - slitmin
   nchunks  = 30.
   chunklen = slitlen/nchunks
   poslist  = fltarr(nchunks)
   medlist  = fltarr(nchunks)
   for i=0,nchunks-1 do begin
      minpos = slitmin + i*chunklen
      maxpos = slitmin + (i+1)*chunklen
      poslist(i) = 0.5*(minpos+maxpos)
      thischunk  = where(slitorder ge minpos and slitorder lt maxpos)
      medlist(i) = median(noskyorder(thischunk))
   endfor

   ; Perform a Gaussian fit to the object profile using the location of
   ; the maximum in the above median list or a user input value as a
   ; guess at the location of the object.
   heightguess = max(medlist)
   if keyword_set(slitpos) then centguess = slitpos else $
      centguess = poslist(where(medlist eq max(medlist)))
   goodls = where(noskyorder ge -200 and noskyorder le (5*heightguess > 200))
   goodslitorder  = slitorder(goodls)
   goodnoskyorder = noskyorder(goodls)
   gfit = gaussfit(goodslitorder,goodnoskyorder,gterms,nterms=3,yerror=gerr,$
                   estimates=[heightguess,centguess,2])
   ; Reject outliers and refit
   if not(keyword_set(standard)) then begin
      resid = goodnoskyorder - gfit
      goodls = where(abs(resid) le 3*gerr)
      goodslitorder = goodslitorder(goodls)
      goodnoskyorder = goodnoskyorder(goodls)
      gfit = gaussfit(goodslitorder,goodnoskyorder,gterms,nterms=3,yerror=gerr,$
                      estimates=gterms)
   endif
   ; Record temporary trace parameters
   objheight   = gterms(0)
   objposition = gterms(1)
   objsig      = abs(gterms(2))
   
   ; Return if no reasonable object profile is found.
   if (objposition lt slitmin or objposition gt slitmax or $
       objsig lt 0.5 or objsig gt 20.) then begin
      print,'No object profile found.  Try entering parameters by hand.'
      return
   endif

   ; Display profile (optional).
   if not(keyword_set(noshow)) then begin
      chan,10
      plot,slitorder,noskyorder,psym=3,yrange=[-3*objheight,3*objheight]
      oplot,goodslitorder,gfit,psym=3
   endif
   print,'Initial fit terms: ',gterms

   ; If using a reference sky exposure, locate and fit object trace in
   ; reference exposure.
   if keyword_set(skyref) then begin
      refheightguess = min(medlist)
      refcentguess   = poslist(where(medlist eq min(medlist)))
      goodls = where(noskyorder le 200 and $
                     noskyorder ge (5*refheightguess < (-200)))
      goodslitorder  = slitorder(goodls)
      goodnoskyorder = noskyorder(goodls)
      refgfit = gaussfit(goodslitorder,goodnoskyorder,refgterms,nterms=3,$
                         yerror=gerr,estimates=[refheightguess,refcentguess,2])
      ; Reject outliers and refit
      if not(keyword_set(standard)) then begin
         resid = goodnoskyorder - refgfit
         goodls = where(abs(resid) le 3*gerr)
         goodslitorder = goodslitorder(goodls)
         goodnoskyorder = goodnoskyorder(goodls)
         refgfit = gaussfit(goodslitorder,goodnoskyorder,refgterms,nterms=3,$
                            yerror=gerr,estimates=refgterms)
      endif
      ; Record temporary trace parameters
      refheight   = refgterms(0)
      refposition = refgterms(1)
      refsig      = abs(refgterms(2))
      ; Return if no reasonable profile is found.
      if (refposition lt slitmin or refposition gt slitmax or $
          refsig lt 0.5 or refsig gt 20.) then begin
         print,'No reference object profile found.'
         return
      endif
      ; Testing
      if not(keyword_set(noshow)) then $
         oplot,goodslitorder,refgfit,psym=3,color=100
      print,'Reference object fit terms: ',refgterms
   endif else begin
      refheight   = 0.
      refposition = 0.
      refsig      = 0.
   endelse

endelse

;;;
;;; 5. Sky subtraction loop.
;;;

; Sky-subtraction / object extraction loop: First, iterate fit in
; dispersion direction only until fit converges.  Then, fit b-spline
; in dispersion direction and a low-order polynomial along slit.

; Begin by fitting non reference sky-subtracted image in order to
; determine the expected variance in the counts.  When that fit has
; converged, fit reference sky-subtracetd image and extract object
; spectrum.

; Determine parameters for computing expected noise.  

; Number of coadds for object (could also read from header)     
if keyword_set(objcoadds) then ncoadds = objcoadds else ncoadds=1.
; Number of coadds for dark (could also read from header)     
if keyword_set(darkcoadds) then dkcoadds = darkcoadds else dkcoadds=objcoadds
; Read noise, in ADU
if keyword_set(readnoise) then erdnoise = readnoise else erdnoise=10. ;25.

; Determine which pixels to fit.
slitedge1   = 230
if keyword_set(standard) then slitedge2 = 525 else slitedge2 = 511
if (fit_entire_slit) then begin
   slitedge2 = 518
   slitmin  = slitedge1
   slitmax  = slitedge2
endif else begin
   subwidth = 100.  ; Width of slit section around object in which to fit sky
   slitmin  = objposition - subwidth/2. > slitedge1
   slitmax  = (objposition + subwidth/2.) < slitedge2
endelse
exposed     = where(slitgrid ge slitmin and slitgrid le slitmax)
nexposed    = n_elements(exposed)
wavelist    = wavegrid(exposed)
slitlist    = slitgrid(exposed)
rawobjlist  = rawobj(exposed)
subobjlist  = subobj(exposed)
darklist    = darkarr(exposed)
flatlist    = flatarr(exposed)
flatvarlist = flatvararr(exposed)
badlist     = badmap(exposed)
; Sort arrays according to wavelength
order        = sort(wavelist)
exposedorder = exposed(order)
waveorder    = wavelist(order)
slitorder    = slitlist(order)
rawobjorder  = rawobjlist(order)
subobjorder  = subobjlist(order)
darkorder    = darklist(order)
flatorder    = flatlist(order)
flatvarorder = flatvarlist(order)
badorder     = badlist(order)

; Testing
nosky = rawobj
skymodel = 0.*rawobj
noskynoobj = rawobj

niter = 0         ; Counter for number of sky-subtraction iterations
wavebin = 1.      ; Initial binning for extracted object spectrum
firstsubfit = 0   ; Flag for the first fit to the ref-sky-subtracted image
fitconverged = 0  ; Flag that current sky fit has converged
objfitmade  = 0   ; Flag that object spectrum has been extracted
extractobj   = 0  ; Flag to extract object
pgrid = 0.*rawobj ; Probability grid
skyfactor = 1.    ; Multiplicative constant by which to adjust sky value
coaddfactor = 1.  ; Multiplicative constant by which to adjust # coadds
done  = 0         ; Flag that that sky-subtraction/object-extraction loop is done.

; Set initial sky and object fits to zero such that the variance will
; first be computed assuming that the array is everywhere readnoise
; dominated.
objfitorder    = 0.*rawobjorder
skyfitorder    = 0.*rawobjorder
rawskyfitorder = 0.*rawobjorder

; Identify bad pixels
badls = where(badorder eq 0)

; Begin by fitting raw (un-ref-sky-subtracted) frame
fitraw = 1

; If fitting standard, only fit raw sky once
if keyword_set(standard) then niter=5

while not(done) do begin

   ; Increment count of iterations
   niter = niter + 1

   print,'Skyfit iteration '+strtrim(string(niter),2)
   if (fitraw) then begin
      print,'Fitting unsubtracted image...'
      skyorder = rawobjorder
      n_poly   = 0.
   endif else begin
      print,'Fitting subtracted image...'
      skyorder = subobjorder
      n_poly = floor(niter/2.)
      print,'n_poly = '+strtrim(string(n_poly))
      ; If a sky reference is being used, increase the expected photon
      ; noise under the assumption that the two frames have roughly
      ; equal counts.  Also, double the number of exposures to
      ; account for readnoise.
      if keyword_set(skyref) then begin
         skyfactor = 2.
         coaddfactor = 2.
      endif
   endelse   

   ; Create inverse variance array based on updated raw sky and object fits
   varorder = nr_makevariance(rawobjorder,object=objfitorder,$
                              sky=skyfactor*rawskyfitorder,$
                              flatvar=flatvarorder,$
                              darkframe=darkarr,readnoise=erdnoise,$
                              objcoadds=coaddfactor*ncoadds,$
                              darkcoadds=dkcoadds,wave=waveorder)
   invvarorder = 1./varorder

   ; Mask bad pixels and those with suspicious values
   if not(keyword_set(standard)) then begin
      if (niter gt 1 and not(firstsubfit)) then begin
         outlierls = where(abs(skyorder-skyfitorder) gt 3*sqrt(varorder))
      endif else begin
         outlierls = where(skyorder ge 60000 or skyorder le -60000)
      endelse
   endif else outlierls = -1

   if (total(outlierls) ne -1) then invvarorder(outlierls) = 0.
   if (total(badls)     ne -1) then invvarorder(badls) = 0.

   ; Mask out the object and object in reference sky exposure (if any)
   if not(keyword_set(noobject)) then begin
      print,'Masking object.'
      if keyword_set(standard) then nsig = 5. else nsig = 3.
      objls = where(slitorder ge objposition - nsig*objsig and $
                    slitorder le objposition + nsig*objsig)
      refls = where(slitorder ge refposition - nsig*refsig and $
                    slitorder le refposition + nsig*refsig)
      invvarorder(objls) = 0.
      if (total(refls) ne -1) then invvarorder(refls) = 0.
   endif 

   ; Perform 2D b-spline fit to sky
   every_n  = 25.
   rej      = 5.
   max_iter = 20.
   if (n_poly gt 0) then x_2 = slitorder else x_2 = 0
   if (n_poly eq 2) then done = 1

   sset=bspline_iterfit(waveorder,skyorder,x2=x_2,npoly=n_poly,$
                        invvar=invvarorder,everyn=every_n,yfit=skyfitorder,$
                        upper=rej,lower=rej)

   ; If fitting the raw (un-reference-subtracted) frame, record sky
   ; fit to be used in computing variance.
   if (fitraw) then rawskyfitorder = skyfitorder

   ; For standards, pretend sky is zero
   if keyword_set(standard) then skyfitorder(*)=0.

   ; Testing - display sky-subtracted frame
   if (fitraw) then begin
      nosky = rawobj
      nosky(exposedorder) = rawobj(exposedorder) - skyfitorder
      noskyorder = rawobjorder - skyfitorder
   endif else begin
      nosky = subobj
      nosky(exposedorder) = subobj(exposedorder) - skyfitorder
      noskyorder = subobjorder - skyfitorder
   endelse
   skymodel(exposedorder) = skyfitorder
   if not(keyword_set(noshow)) then begin
      chan,11
      doit,nosky(*,130:511),-80,80
   endif

   ; Create a new variance array based on the sky fit.
   varorder = nr_makevariance(rawobjorder,object=objfitorder,$
                              sky=skyfactor*rawskyfitorder,$
                              flatvar=flatvarorder,$
                              darkframe=darkarr,readnoise=erdnoise,$
                              objcoadds=coaddfactor*ncoadds,$
                              darkcoadds=dkcoadds,wave=waveorder)
   invvarorder = 1./varorder
   if (total(badls)     ne -1) then invvarorder(badls) = 0.

   ; Test: always extract object, unless user specifies no object present
   if keyword_set(noobject) then extractobj = 0 else extractobj = 1

   ; Extract object spectrum
   if (extractobj) then begin
      ; Subtract model sky
      if (fitraw) then noskyorder = rawobjorder - skyfitorder $
         else noskyorder = subobjorder - skyfitorder
      ; Refine object profile using Gaussian fit
      ; If fitting raw sky, or if reducing a standard, 
      ; fit a narrow section in wavelength.  Otherwise, use
      ; the the extracted spectrum to normalize the profile
      ; at all wavelengths.  Assumes that if fitting the subtracted
      ; sky, a rough object spectrum has already been extracted.
      if (fitraw or keyword_set(standard)) then begin
         section = where(waveorder ge 350 and waveorder le 400 and $
                         abs(noskyorder) le 5*objheight)
         slitsection  = slitorder(section)
         noskysection = noskyorder(section)
      endif else begin 
         specorder    = interpol(spec,speclam,waveorder)
         normalized   = noskyorder / specorder
         goodnormls   = where(abs(normalized) le 1)
         slitsection  = slitorder(goodnormls)
         noskysection = normalized(goodnormls)
         ; "Normalize" expected peak intensity, if necessary.
         if (objheight gt 1) then objheight = 0.5
      endelse
      gfit = gaussfit(slitsection,noskysection,gterms,nterms=3,$
                      estimates = [objheight,objposition,objsig],$
                      yerror=gfiterr)
      if not(keyword_set(standard)) then begin
         resid = noskysection - gfit
         goodls = where(abs(resid) le 4*gfiterr)
         goodslitsection = slitsection(goodls)
         goodnoskysection = noskysection(goodls)
         gfit = gaussfit(goodslitsection,goodnoskysection,gterms,nterms=3,$
                         estimates=gterms,yerror=gfiterr)
      endif
      ; Update trace parameters
      objheight   = gterms(0)
      objposition = gterms(1)
      objsig      = abs(gterms(2))
      ; Construct probability grid
      porder = 0*slitorder
      const = 1./(objsig*sqrt(2*!pi))
      for i=long(0),long(nexposed-1) do porder(i) = $
         const * exp( -1*((slitorder(i)-objposition)^2)/(2*objsig^2))
      ; Standard star - bspline fit profile
      if keyword_set(standard) then begin
         invvarsection = invvarorder(section)
         sset_profile = bspline_iterfit(slitsection,noskysection,$
                            bkspace=1,invvar=invvarsection,yfit=bfit)
         profilemax   = max(bfit(where(slitsection le slitedge2)))
         profileorder = bspline_valu(slitorder,sset_profile)
         porder       = profileorder / profilemax
         ; Suppress spurious structure outside of main peak
         wheremax = where(bfit eq profilemax)
         profilecenter = slitsection(wheremax)
         profilecenter = profilecenter(0)  ; Force conversion to scalar.
         rightpixels = where(slitsection gt profilecenter and $
                             bfit lt 0.01*profilemax)
         if (total(rightpixels) ne -1) then $
            rightedge = min(slitsection(rightpixels)) else $ 
            rightedge = objposition + 5*objsig
         rightls = where(slitorder ge rightedge)
         if (total(rightls) ne -1) then porder(rightls) = 0.
         leftpixels = where(slitsection lt profilecenter and $
                            bfit lt 0.01*profilemax)
         if (total(leftpixels) ne -1) then $
            leftedge = max(slitsection(leftpixels)) else $
            leftedge = objposition - 5*objsig
         leftls = where(slitorder le leftedge)
         if (total(leftls) ne -1) then porder(leftls) = 0.
         pordersection = porder(section)
      endif
      nearzerols = where(porder le 1e-6)
      if (total(nearzerols) ne -1) then begin
         porder(nearzerols) = 0.
      endif else begin
         print,'Error - Object profile extends over entire slit region.'
         return
      endelse
      ; Extract
      spec = opt_extract(noskyorder,waveorder,invvarorder,porder,$
                         ;speclam=speclam,bin=wavebin,maxcount=10*objheight,$
                         speclam=speclam,bin=wavebin,maxcount=60000.,$
                         modelcounts=objfitorder,specvar=specvar)

      ; Display extracted spectrum (optional)
      if not(keyword_set(noshow)) then begin
         chan,10
         y_min = -0.5*median(spec)
         y_max = 2.5*median(spec)
         plot,speclam,spec,yrange=[y_min,y_max],xstyle=1
         oplot,speclam,sqrt(specvar),color=100
         xyouts,800,0.9*y_max,'niter  = '+strtrim(string(niter),2),/data
         xyouts,800,0.8*y_max,'n_poly = '+strtrim(string(n_poly),2),/data
      endif

      objfitmade = 1

   endif

   ; Has current fit converged?  For now, specify number of 
   ; iterations.
   if (niter eq 6) then begin
      fitconverged = 1
      niter = 0
   endif

   ; Continue to fit raw sky, ref-subtracted sky, or done?
   if (fitraw and fitconverged) then begin
      fitraw = 0
      fitconverged = 0
   endif else begin
      if (not(fitraw) and not(fitconverged)) then begin
         if (firstsubfit) then firstsubfit = 0  ; First sub fit done.
      endif
      if (not(fitraw) and fitconverged) then done = 1
   endelse

endwhile

;;;
;;; Wavelength calibrate
;;;

if keyword_set(waveref) then begin  ; pre-made wavelength fit (optional)

   if (size(waveref,/tname) eq 'STRING') then begin
      checkexists = findfile(waveref,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Wave reference file not found.'
         return
      endif else begin
          wavereffile = 1
          waverefarr = readfits(waveref,hwaveref)
          sxaddpar,hobj,'WAVEREF',waveref,'Filename of wavelength reference.'
      endelse
   endif else begin
      wavereffile = 0
      waverefarr = waveref
   endelse
   ; Possible that this is a subarray.  If smaller than 1024 x 1024,
   ; patch into a full array.
   if (n_elements(waverefarr) lt 1024L^2) then begin
      sz = size(waverefarr)
      if (sz(2) ne 450) then begin
         print,'Error - Unexpected size for wavelength reference array.'
         return
      endif
      ; Fill remainder of wavelength array with dummy (high) values
      tmp_waverefarr = 1e9 + 0*rawobj
      tmp_waverefarr(*,100:549) = waverefarr
      waverefarr = tmp_waverefarr
   endif
   wavefitlist  = waverefarr(exposed)
   wavefitorder = wavefitlist(order)
   ; Record in header the heliocentric velocity correction applied to 
   ; reference wavelength array.  In general, standards will have
   ; coordinates and will be taken at similar times as the science
   ; objects serving as wavelength references, therefore heliocentric
   ; velocities should be very similar.
   if (wavereffile) then begin
      vhelio = sxpar(hwaveref,'HELIOVEL',count=heliovelcount)
      if (heliovelcount eq 1) then begin
         sxaddpar,hobj,'HELIOCNT','Heliocentric correction from waveref.'
         sxaddpar,hobj,'HELIOVEL',vhelio,'lam_hcnt = lam_obs*(1+heliovel/c)'
         print,'Heliocentric velocity borrowed from reference array.'
      endif
   endif

endif else begin   

   ; Perform skyine fit.  Note - tabulated OH line wavelengths are vacuum.
   ; Use model wavelength fit for H-band.  The model wavelength fit
   ; is different from the reference wavelength array above.  It is a 
   ; 3-column ASCII table that gives the intensity and wavelength over
   ; a range in wavelength-proxies (the values so far used as 
   ; "wavelengths").  This fir is used as an estimate to compute a 
   ; wavelength fit for the current array.

   if keyword_set(linelist) then line_file = linelist else $
      line_file = '/home/tallis/gdb/nirspec/oh_lists/'+$
	          'lowd_ir_ohlines_hband_best.lst' 

   ;; Automatic line identification?
   if keyword_set(wavemodel) then begin
      wave_model = wavemodel
      ;Wavelength fit for NIRSPEC H-band, on tallis:
      ;wave_model  = '/scr2/gdb/keck/nirspec/nov03/raw/calibrations/spec/'+$
      ;              'skyline_calib_ref.dat' 
      wave_calib,waveorder,rawskyfitorder,wavefitorder,reference=wave_model,$
              linefile=line_file,coeffs=wavecoeffs,/nopersist
   endif else begin
      ; Manual line identification
      wave_calib,waveorder,rawskyfitorder,wavefitorder,$
              linefile=line_file,coeffs=wavecoeffs
   endelse

   ; Record in header that vacuum wavelengths used.
   if (fromfitsfile) then sxaddpar,hobj,'VACUUM',$
      'Vacuum wavelengths used for OH line wave calibration.'

   ;;; Apply heliocentric velocity correction
   if (fromfitsfile) then begin
      ; Coordinates for Keck
      kecklong = 155.478  ; degrees
      kecklat  = 19.8283  ; degrees
      keckalt  = 4160.    ; meters
      ; Get UT modified Julian date and convert to UT Julian date.
      ; Note: could also read UT date and time and convert to Julian
      ; date using JULDAY.
      mjday = sxpar(hobj,'MJD-OBS')
      jday  = mjday + 2400000.5
      ; Get telescope coordinates and epoch
      objra  = sxpar(hobj,'RA')  ; Degrees, for NIRSPEC
      objdec = sxpar(hobj,'DEC') ; Degrees, for NIRPSEC
      epoch  = sxpar(hobj,'EQUINOX')
      ; Get heliocentric velocity correction
      vcorr = heliocentric(objra,objdec,epoch,jd=jday,long=kecklong,$
                           lat=kecklat,altitude=kecklat)
      ; Change sign such that the velocity is that towards the object
      ; with respect to the sun
      vhelio = -1*vcorr
      ; Adjust wavelengths
      c = 299792d0
      wavefitorder = wavefitorder * (1. + vhelio/c)
      ; Add keywords similar to those added by MAKEE.
      sxaddpar,hobj,'HELIOCNT','Heliocentric correction applied.'
      sxaddpar,hobj,'HELIOVEL',vhelio,'lam_hcnt = lam_obs*(1+heliovel/c)'
      print,'Heliocentric velocity correction applied.'
   endif else vhelio = 0.

endelse


; Extract object - binning in angstroms.  For R ~ 2000, lam ~ 15000,
; res element ~ 8 ang
if not(keyword_set(noobject)) then begin
   wavebinang=3.
   spec = opt_extract(noskyorder,wavefitorder,invvarorder,porder,$
                      speclam=speclam,bin=wavebinang,maxcount=60000.,$ 
                      modelcounts=objfitorder,specvar=specvar)
   if not(keyword_set(noshow)) then begin
      chan,10
      plot,speclam,spec,yrange=[-0.5*median(spec),2.5*median(spec)],xstyle=1
      oplot,speclam,sqrt(specvar),color=100
   endif
endif

;;;
;;; Write outputs - For 2D Outputs, extract 1023 x 550 arrays
;;;

if keyword_set(outfile) then rootname = outfile else rootname = 'nirspec'

; 2D Subtracted frame
print,'Writing subtracetd frame...'
nosky = subobj
nosky(exposedorder) = subobj(exposedorder) - skyfitorder
hextract,nosky,hobj,smnosky,smhobj,0,1023,100,549
writefits,rootname+'_sub.fits',smnosky,smhobj
; 2D Sky model
print,'Writing sky model frame...'
skymodel = 0.*rawobj
skymodel(exposedorder) = skyfitorder
hextract,skymodel,hobj,smskymodel,smhobj,0,1023,100,549
writefits,rootname+'_sky.fits',smskymodel,smhobj
; 2D Variance array
print,'Writing variance frame...'
variance2d = 0.*rawobj
variance2d(exposedorder) = varorder
hextract,variance2d,hobj,smvariance2d,smhobj,0,1023,100,549
writefits,rootname+'_var.fits',smvariance2d,smhobj
; 2D Wavelength array - reconstruct over entire array
print,'Writing wavelength frame...'
if keyword_set(waveref) then begin
   wave2d = waverefarr
endif else begin 
   nwavecoeffs = n_elements(wavecoeffs)
   wave2d = 0d0*rawobj
   for i=0,nwavecoeffs-1 do wave2d = wave2d + wavecoeffs(i)*(wavegrid^i)
   ; Heliocentric correction
   wave2d = wave2d * (1. + vhelio/c)
endelse
hextract,wave2d,hobj,smwave2d,smhobj,0,1023,100,549
writefits,rootname+'_wav.fits',smwave2d,smhobj
; ASCII file containing wavelength fit (optional)
if keyword_set(writewave) then begin
   ; Create a list of wavelenth-proxy values, in steps of 0.1 (pixels)
   wavegridstep       = 0.1
   nwavesteps         = ceil((max(waveorder) - min(waveorder)) / wavegridstep)
   tmp_waveorder      = min(waveorder) + wavegridstep*findgen(nwavesteps)
   tmp_wavefitorder   = interpol(wavefitorder,waveorder,tmp_waveorder)
   tmp_rawskyfitorder = interpol(rawskyfitorder,waveorder,tmp_waveorder)
   print,'Writing wavelength fit to ASCII file...'
   openw,1,rootname+'_wav.dat'
   for i=0,nwavesteps-1 do printf,1,tmp_waveorder(i),tmp_rawskyfitorder(i),$
                                    tmp_wavefitorder(i)
   close,1
endif
; Outputs for extracted spectrum, if any
if not(keyword_set(noobject)) then begin
   ; 1D Extracted spectrum - NAXIS keywords updated automatically.
   hspec = hobj
   sxaddpar,hspec,'CRVAL1',min(speclam)
   sxaddpar,hspec,'CDELT1',wavebinang
   sxaddpar,hspec,'CD1_1',wavebinang
   sxaddpar,hspec,'CRPIX1',1
   sxaddpar,hspec,'CTYPE1','LAMBDA'
   sxaddpar,hspec,'CTYPE2','PIXEL'
   writefits,rootname+'_spec.fits',spec,hspec
   ; 1D Error array
   specerr = sqrt(specvar)
   writefits,rootname+'_err.fits',specerr,hspec
   ; PS Plot - 1D Spectrum + Error
   ps_openfile,rootname+'_spec'
   plot,speclam/10000.,spec,yrange=[0,2.5*median(spec)],xstyle=1,$
        title=rootname,xtitle='wavelength (microns)',ytitle='flux'
   oplot,speclam/10000.,sqrt(specvar),color=100
   ps_closefile
   ; PS Plot - Object Profile
   ps_openfile,rootname+'_profile',/color
   !p.multi = [0,1,1]
   x_min   = min(slitsection)
   x_max   = max(slitsection)
   y_min   = -0.5*objheight
   y_max   = 1.5*objheight
   y_range = [y_min,y_max]
   if keyword_set(standard) then begin
      plot,slitsection,noskysection,psym=3,xstyle=1,$
         yrange=y_range,ystyle=2,title='Object Profile - '+rootname,$
         xtitle = 'slit position',ytitle ='ADU'
      oplot,slitsection,pordersection*profilemax,psym=3,color=100
   endif else begin
      ; Only plot every other point
      npoints = n_elements(slitsection)
      plotls  = 2*findgen(floor(npoints/2.))
      plot,slitsection(plotls),noskysection(plotls),psym=3,xstyle=1,$
         yrange=y_range,ystyle=2,title='Object Profile - '+rootname,$
         xtitle = 'slit position',ytitle ='ADU'
      ngoodpoints = n_elements(goodslitsection)
      goodplotls  =  2*findgen(floor(ngoodpoints/2.))
      oplot,goodslitsection(goodplotls),gfit(goodplotls),psym=3,color=100
   endelse
   objfwhm = 2.*sqrt(2.*alog(2))*objsig
   xyouts,x_min+0.1*(x_max-x_min),0.9*y_max,$
      "FWHM = "+strtrim(string(objfwhm),2),/data
   ps_closefile
endif

!p.multi=currentmulti

return
end

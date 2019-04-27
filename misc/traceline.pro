Pro traceline,arr,x0,y0,xtrace,ytrace,ctrace,xdir=xdir,ydir=ydir,range=range,$
        bound=bound,nsteps=nsteps,nside=nside,subpix=subpix,bymax=bymax,$
        bygauss=bygauss,nterms=nterms,median=median,coord=coord,cbound=cbound,$
        merit=merit,noplot=noplot,modelprofile=modelprofile,$
        modelindex=modelindex

;+---------------------------------------------------------------
;
; TRACELINE        2/2004
;
; Trace lines in a 2D array.  A typical use would be to trace an
; object or an arc line in a 2D spectrum..
;
; INPUTS
;       arr       - 2D spectrum
;       x0,y0     - Either scalar pixel coordinates of a point on the line
;                   or 1D arrays giving x and y positions of several
;                   points (generally on different lines).
;
; KEYWORDS
;       xdir      - Assume that the lines runs in the x-direction;
;                   step along in x and find maximum points in y
;       ydir      - Assume that the lines runs in the y-direction;
;                   step along in y and find maximum points in x
;                   *NOTE: Either /xdir or /ydir must be set.*
;       range     - 2-element vector giving the minimum and maximum
;                   values in the trace direction (x if /xdir is set,
;                   y if /ydir is set) over which to trace line
;       bound     - 2-element vector giving the minimum and maximum
;                   values perpendicular to the trace direction (y if
;                   /xdir is set, x if /ydir is set) over which to trace
;                   line.
;       nsteps    - Number of steps to make away from startingt
;                   position.  If scalar, will step this many pixels
;                   in both directions (2*nsteps+1 positions total).
;                   If two-element array [negnstep,posnstep], will go 
;                   negnstep pixels in the negative direction and
;                   posnstep pixels in the positive direction
;                   (negnstep+posnstep+1 pixels total.  If both range 
;                   and nsteps keywords are set, program will stop at 
;                   whichever limit is hit first.
;       nside     - Number of pixels to either side of the nominal
;                   position to use in crosscorrelating with the 
;                   adjacent position.  If tracing by finding maxima, 
;                   number of pixels to either side of the nominal
;                   position in which to look for the line
;                   maxiumum at any position.  Default is 5.  Minimum is 2.
;       subpix    - Number of subpixels per pixel over which to
;                   interpolate when computing the crosscorrelation
;                   between slices.  Default is 1 (no interpolation).
;                   Will not affect results if /bymax keyword is set.
;                   Should be set to match slope of line, i.e. if
;                   slope is about 0.1 (as for NIRSPEC arc/sky lines)
;                   then use subpix=10
;       bymax     - Trace lines by finding a maximum at every step.
;                   Default is to crosscorrelate short slices
;                   of the array at each step.  Note: For tracing
;                   strong lines, /bymax seems to work as well
;                   or better than computing crosscorrelations
;                   and is faster.  However, when tracing an image at
;                   a "random" position it will tend to
;                   jump over to nearby strong lines.
;       bygauss   - Trace lines by fitting a Guassian at every step.
;                   Gaussian will have four terms only (no linear term).
;                   This seems to work very well for bright standard stars.
;       nterms    - Scalar, number of terms for Gaussian fit (see GAUSSFIT).
;                   Must be between 3 and 5.  Default is 5, which includes
;                   a scalar offset and a linear slope.
;       median    - Take a median slice of this many pixels thickness when
;                   computing the central value at a given position.
;                   When using crosscorrelations, the reference slice
;                   will also be computed from the median of a slice of
;                   this thickness.  Must be an odd number.
;       coord     - 2D array of auxiliary coordinate values (e.g., a 2D
;                   wavelength or slit position solution).  If set, the
;                   program will compute the value of this coordinate at
;                   each trace (x,y) position.
;       cbound    - 2-element vector giving the minimum and maximum
;                   auxilliary coordinate values over which to trace
;                   line.
;       merit     - Array of figures of merit for each trace point.
;                   If the /bygauss keyword is set, then this is the
;                   inverse of the squared standard deviation from the
;                   fit.  If the /bymax keyword is set, then this is 1
;                   when a maximum is found.  If the default
;                   crosscorrelation is used, then this will be the
;                   value of the normalized crosscorrelation at each
;                   point.  In all cases, if a trace point is not
;                   found within some reasonable range, this value
;                   will be 0.  The array will have the same
;                   dimensions as xtrace and ytrace.
;       noplot    - Suppresses plotting of trace points.
;       modelprofile - A 1-D model of the expected line intensity, to
;                   be used as a reference for crosscorreIating with
;                   cuts along the line.  Typically, the model should
;                   be more finely sampled than the data to be traced.
;                   If this keyword is set, the modelindex keyord must 
;                   also be set.  No affect if the /bymax or /bygauss 
;                   keyword is set.  Must also set the modelindex keyword.
;       modelindex - 1-D array of x- or y-values values (not necessarily 
;                   integer) corresponding to the model.  If the /ydir
;                   keyword is set, then this is a list of x-values at which
;                   the model is given.  If the /xdir keyword is set, then
;                   this is a list of y-values.  This keyword must be set if
;                   the /modelprofile keyword is set.
;
; OUPUTS
;       xtrace,ytrace - If x0 and y0 are scalars, array of pixel coordinates
;                       of the trace along the line.  If x0 and y0 are
;                       N-dimanesional arrays, nxN arrays of pixel x
;                       and y values where n is the number of points
;                       along a trace.
;       ctrace        - Optional - array of trace auxiliary
;                       coordinates (keyword coord must be set).  In
;                       the case where the x and y trace values are
;                       not integers, the coordinate value will be
;                       linearly interpolated from neighboring pixels
;                       (along the x-direction if the keyword /ydir is
;                       set, or along the y-direction if the keyword
;                       /xdir is set).
;
; HISTORY
;       Written 2/15/2004 GDB
;       BYGAUSS keyword added 8/19/2004 GDB
;       COORD and CTRACE added 8/30/2004 GDB
;       MODELPROFILE, MODELINDEX keywords added 12/9/2004 GDB
;       Bounding keywords added 1/30/2005 GDB
;       NTERMS keyword added 4/25/2007 GDB
;---------------------------------------------------------

if (n_params() lt 3) then begin
   print,'CALLING SEQUENCE: traceline,arr,x0,y0,xtrace,ytrace,ctrace,/xdir,'
   print,'    /ydir,range=range,bound=bound,nsteps=nsteps,nside=nside,'
   print,'    subpix=subpix,/bymax,/bygauss,nterms=nterms,median=median,'
   print,'    coord=coord,cbound=cbound,merit=merit,/noplot,'
   print,'    modelprofile=modelprofile,model=modelindex'
   return
endif

; Determine size of array.
sz = size(arr)
if (sz(0) ne 2) then begin
   print,'TRACELINE: Error - Input array should be two-dimensional.'
   xtrace = -1
   ytrace = -1
   ctrace = -1
   return
endif
xarrsz = sz(1)
yarrsz = sz(2)

; If using an auxiliary coordinate, array must be same size as input
; array
if keyword_set(coord) then begin
   coordsz = size(coord)
   if (coordsz(1) ne sz(1) or coordsz(2) ne sz(2)) then begin
      print,'TRACELINE: Error - Coordinate array must be same size as '+$
            'input array.'
      xtrace = -1
      ytrace = -1
      ctrace = -1
      return
   endif
   ; Upper and lower boundaries in the auxilliary coordinate
   if keyword_set(cbound) then begin
      cmin = min(cbound)
      cmax = max(cbound)
   endif else begin
      cmin = min(coord)
      cmax = max(coord)
   endelse
endif

; Determine how many lines to trace
x0=[x0]
y0=[y0]
if (n_elements(x0) ne n_elements(y0)) then begin
   print,'TRACELINE: x0 and y0 must have the same number of elements!'
   xtrace = -1
   ytrace = -1
   ctrace = -1
   return
endif
nlines = n_elements(x0)
; Tile plots of trace points
ntiles2 = floor(sqrt(nlines))
ntiles1 = ceil((1.*nlines)/ntiles2)
pmultisave=!p.multi
!p.multi=[0,ntiles1,ntiles2]

; If tracing a specific number of pixels, make sure nsteps is either a
; scalar or a 2-elements array
if keyword_set(nsteps) then begin
   if (n_elements([nsteps]) gt 2) then begin
      print,'TRACELINE: nsteps must either be a scalar or a two-element array.'
      xtrace=-1
      ytrace=-1
      ctrace=-1
      return
   endif
endif

; At each step, search at most +/- this many pixels
; around the nominal position
if keyword_set(nside) then side=nside else side=5
if (nside lt 2) then begin
   print,'TRACELINE: Minimum value for nside is 2.'
   return
endif

; Interpolate over this many subpixels per pixel in each slice
if keyword_set(subpix) then nsubpix=float(subpix) else nsubpix = 1.

; Create median slices of this thickness
if keyword_set(median) then begin
   if (median mod 2 eq 0) then begin
      print,'TRACELINE: Keyword "median" must be odd!'
      return
   endif else nmed = abs(median)
endif else nmed = 1

; Location with respect to nominal position
intlocation = findgen(2*side+1)-side ; Integer number of pixels
location    = (findgen(2*side*nsubpix+1)-side*nsubpix)/nsubpix ; Fractional pixels

; Will compute crosscorr up to +/- 2 integer pixels in lag
nlag        = 2*(2 < side)*nsubpix+1
subpixlag   = findgen(nlag)-2*nsubpix ; Lag in number of sub-pixels
laglocation = (findgen(nlag)-2*nsubpix)/nsubpix ; Lag in fractional pix 

; Provide boost to crosscorr to compensate for the number of (sub)pixels
; used to compute the crosscorr at each lag
ccboost = abs(subpixlag) + 2*side*nsubpix+1 
ccboost = ccboost/min(ccboost)

; Set direction in which to step through array
if not(keyword_set(xdir) or keyword_set(ydir)) then begin
   print,'TRACELINE: Error - Must set either /xdir or /ydir keyword.'
   xtrace = -1
   ytrace = -1
   ctrace = -1
   return
endif
if keyword_set(xdir) then begin ; Trace line in x-direction
   xline=1
   yline=0
   xstep=nmed
   ystep=0.
   xside=(nmed-1)/2.
   yside=side
   meddim=1
   if keyword_set(range) then begin
      xmin = min(range)
      xmax = max(range)
   endif else begin
      xmin = 0
      xmax = xarrsz-1
   endelse
   if keyword_set(bound) then begin
      ymin = min(bound)
      ymax = max(bound)
   endif else begin
      ymin = 0
      ymax = yarrsz-1
   endelse
endif else begin ; Trace line in y-direction
   xline=0
   yline=1
   xstep=0
   ystep=nmed
   xside=side
   yside=(nmed-1)/2.
   meddim=2
   if keyword_set(range) then begin
      ymin = min(range)
      ymax = max(range)
   endif else begin
      ymin = 0
      ymax = yarrsz-1
   endelse
   if keyword_set(bound) then begin
      xmin = min(bound)
      xmax = max(bound)
   endif else begin
      xmin = 0
      xmax = xarrsz-1
   endelse
endelse

; Store trace coordinates in arrays long enough to
; accommodate the longest possible lines.
if keyword_set(xdir) then length_bound = xarrsz else length_bound = yarrsz
xtrace = fltarr(length_bound,nlines)
ytrace = fltarr(length_bound,nlines)
ctrace = fltarr(length_bound,nlines)
merit  = fltarr(length_bound,nlines)
; Keep track of the maximum actual line length.
max_length = 0

;; Flag that the first line is being traced
;firstline = 1

for i=0,nlines-1 do begin
   
   ;;; Start at the user-specified coordinate and proceed in the
   ;;; positive direction until the end of the line (not yet
   ;;; implemented), edge of the array, or user-specified maximum is
   ;;; reached.  At that point return to the nominal position and
   ;;; proceed in the negative direction.  .  At each point, extract
   ;;; a small 1-D slice of the array, locate the maximum (or find the
   ;;; center/offset by some other means) and record its position.

   ; Set initial nominal position
   xnom = x0(i)
   ynom = y0(i)

   donext = 1 ; Flag to proceed to the next position
   first  = 1 ; Flag to indicate that we are finding the first trace coordinates

   stepstaken = 0 ; Counter to record how many steps taken in this direction

   ; If the initial nominal position is at the low edge of the array,
   ; then do not proceed in the negative direction.
   donegative = 1
   if (xline and xnom-xside le 0) then donegative = 0
   if (yline and ynom-yside le 0) then donegative = 0
   
   ; If the initial nominal position is at the high edge of the array,
   ; begin by stepping in the negative direction
   if (xline and xnom+xside ge xarrsz-1) then xstep = -xstep
   if (yline and ynom+yside ge yarrsz-1) then ystep = -ystep

   ; If tracing by crosscorrelating adjacent slices, intialize the
   ; trace coordinates with the current nominal positions.  Set the
   ; reference slice to be the slice at the user-defined coordinates.
   if not(keyword_set(bymax) or keyword_set(bygauss)) then begin
      xtrace0 = [xnom]
      ytrace0 = [ynom]
      merit0   = [1]
      ; Make sure not starting too close to the edge of the arra
      if (xtrace0-xside lt 0 or xtrace0+xside gt xarrsz-1 or $
          ytrace0-yside lt 0 or ytrace0+yside gt yarrsz-1) then begin
         print,'TRACELINE: When using crosscorrelation, starting point'
         print,'           must be away from the edge of the array.'
         xtrace = -1
         ytrace = -1
         ctrace = -1
         return
      endif
      ; Use a model for the refence slice?
      if keyword_set(modelprofile) then begin
         if keyword_set(ydir) then begin
            reflo = xnom-xside
            refhi = xnom+xside
            refls = where(modelindex ge reflo and modelindex le refhi)
            modellocation = modelindex(refls) - xnom
            modelslice    = modelprofile(refls)
            refslice      = interpol(modelslice,modellocation,location)
         endif else begin
            reflo = ynom-yside
            refhi = ynom+yside
            refls = where(modelindex ge reflo and modelindex le refhi)
            modellocation = modelindex(refls) - ynom
            modelslice    = modelprofile(refls)
            refslice      = interpol(modelslice,modellocation,location)
         endelse
      endif else begin
         ; Use the first slice for the reference slice.
         if (nmed gt 1) then begin
            refslice = median(arr(xnom-xside:xnom+xside,ynom-yside:ynom+yside),$
                              dimension=meddim)
         endif else begin
            refslice = arr(xnom-xside:xnom+xside,ynom-yside:ynom+yside)
         endelse
         refslice = interpol(refslice,intlocation,location)
      endelse
      ; Make sure there are no spurious (NaN) points in the referece slice
      temp_refslice = refslice
      goodls = where(refslice ge -1e29 and refslice le 1e29)
      refslice(*) = 0.
      if (total(goodls) ne -1) then refslice(goodls) = temp_refslice(goodls)
      ; Iterate count.
      first = 0
      xnom = xnom+xstep
      ynom = ynom+ystep
      stepstaken = stepstaken+1
   endif

   while donext do begin

      ; Force conversion to scalar
      xnom = round(xnom(0))
      ynom = round(ynom(0))

      ; Slice will ultimately be 1-D since either xside or yside is 0,
      ; or a median will be taken along one of the directions.
      xlo   = (xnom - xside) > xmin
      xhi   = (xnom + xside) < xmax
      ylo   = (ynom - yside) > ymin
      yhi   = (ynom + yside) < ymax
      if (nmed gt 1 and xlo ne xhi and ylo ne yhi) then do_median = 1 $
         else do_median = 0
      if do_median then begin
         slice = median(arr(xlo:xhi,ylo:yhi),dimension=meddim)
      endif else begin
         slice = arr(xlo:xhi,ylo:yhi)
      endelse

      ; If part of the slice falls outside the 2D array, shorten the
      ; location arrays.
      if xline then begin
         min_location = ylo - ynom
         max_location = yhi - ynom
      endif else begin
         min_location = xlo - xnom
         max_location = xhi - xnom
      endelse
      int_ls = where(intlocation ge min_location and $
                     intlocation le max_location)
      sub_ls = where(location ge min_location and $
                     location le max_location)
      lag_ls = where(laglocation ge min_location and $
                     laglocation le max_location)
      sm_intlocation = intlocation(int_ls)
      sm_location    = location(sub_ls)
      sm_laglocation = laglocation(lag_ls)
   
      ;; Optional: If part of the slice falls outside the allowed
      ;; range, move the slice over.
      ;if keyword_set(cbound) then begin
      ;   if do_median then begin
      ;      coord_slice = median(coord(xlo:xhi,ylo:yhi),dimension=meddim)
      ;   endif else begin
      ;      coord_slice = coord(xlo:xhi,ylo:yhi)
      ;   endelse
      ;   below_ls = where(coord_slice lt cmin or coord_slice gt cmax,$
      ;                       n_notbound)
      ;   if xline then begin
      ;
      ;   endif else begin
      ;      
      ;   endelse
      ;endif

      ;; For now, pad points that fall outside the array with zeros.
      ;; This should only affect the initial position.
      ;if xline then begin
      ;   if (ynom-yside lt 0) then begin
      ;      diff = ynom-yside
      ;      dummyvalue = slice(0)
      ;      for d=0,abs(diff)-1 do slice = [dummyvalue,slice]
      ;   endif
      ;   if (ynom+yside gt yarrsz-1) then begin
      ;      diff = ynom+yside - (yarrsz-1)
      ;      dummyvalue = slice(n_elements(slice)-1)
      ;      for d=0,abs(diff)-1 do slice = [slice,dummyvalue]
      ;   endif
      ;endif else begin
      ;   if (xnom-xside lt 0) then begin
      ;      diff = xnom-xside
      ;      dummyvalue = slice(0)
      ;      for d=0,abs(diff)-1 do slice = [dummyvalue,slice]
      ;   endif
      ;   if (xnom+xside gt xarrsz-1) then begin
      ;      diff = xnom+xside - (xarrsz-1)
      ;      dummyvalue = slice(n_elements(slice)-1)
      ;      for d=0,abs(diff)-1 do slice = [slice,dummyvalue]
      ;   endif
      ;endelse

      ; Interpolate onto sub-sampled mesh.  If finding the maximum,
      ; use spline interpolation.
      if keyword_set(bymax) then begin
         yinit = spl_init(sm_intlocation,slice)
         slice = spl_interp(sm_intlocation,slice,yinit,sm_location)
      endif else begin
         slice = interpol(slice,sm_intlocation,sm_location)
      endelse          

      ; Optional: Allow only points that fall within the specified
      ; range of auxilliary coordinates.
      if keyword_set(cbound) then begin
         if do_median then $
            coord_slice = median(coord(xlo:xhi,ylo:yhi),dimension=meddim) $
            else coord_slice = coord(xlo:xhi,ylo:yhi)
         coord_slice = interpol(coord_slice,sm_intlocation,sm_location)
         bound_ls    = where(coord_slice ge cmin and coord_slice le cmax,$
                             n_bound)
         if (n_bound gt 0) then begin
            tmp_slice       = slice(bound_ls)
            if keyword_set(refslice) then tmp_refslice = refslice(bound_ls)
            tmp_location    = sm_location(bound_ls)
            tmp_laglocation = sm_laglocation
         endif else begin
            print,'TRACELINE: Array slice falls outside of allowed bounds.'
            
         endelse
      endif else begin
         tmp_slice       = slice
         if keyword_set(refslice) then tmp_refslice = refslice
         tmp_location    = sm_location
         tmp_laglocation = sm_laglocation
      endelse

      ; Locate the peak point.
      if keyword_set(bymax) then begin
         offset    = tmp_location(where(tmp_slice eq max(tmp_slice)))
         conf      = 1
      endif else begin
         if keyword_set(bygauss) then begin
            ; Require enough pixels to fit
            if (n_elements(tmp_slice) le 5) then begin
               offset = 0
               conf   = 0
            endif else begin
               ; Guess at fit parameters.  If this is the first point,
               ; use the input nominal position.
               maxpt  = where(tmp_slice eq max(tmp_slice))
               t1     = max(tmp_slice)
               t2     = tmp_location(maxpt)
               if (n_elements(t2) gt 1) then t2=median(t2,/even)
               ;If this is the first point, use the input nominal position.
               if first then t2=0.
               t3     = 2.
               t4     = 0.
               t5     = 0.
               if keyword_set(nterms) then n_terms = nterms $
                  else n_terms = 5
               all_estimates = [t1,t2,t3,t4,t5]
               estimates     = all_estimates(0:n_terms-1)
               gfit   = gaussfit(tmp_location,tmp_slice,terms,nterms=n_terms,$
                                 sigma=sig,estimates=estimates)
               offset = terms(1)
               conf   = 1. / (sig(1)^2)
               ; Reject anything greater than 500-sigma in offset, with
               ; FWHM < 1 pixel, or with a central value outside the
               ; selected range  as spurious.
               if (abs(terms(1)/sig(1)) ge 500 or $
                   2*sqrt(2*alog(2))*abs(terms(2)) lt 1 or $
                   offset lt min(tmp_location) or $
                   offset gt max(tmp_location)) then begin
                  offset = 0
                  conf   = 0
               endif 
           endelse
        endif else begin
            ; Require enough pixels to perform reasonable cross-correlation
            if ((n_elements(tmp_slice) le 5) or $
                 n_elements(tmp_refslice) lt n_elements(subpixlag)) then begin
               offset = 0
               conf   = 0
            endif else begin
               cc = c_correlate(tmp_refslice,tmp_slice,subpixlag)
               ccadjust = ccboost*cc
               offset = tmp_laglocation(where(ccadjust eq max(ccadjust)))
               conf   = cc
            endelse
         endelse
      endelse
      ; If more than one offset point is found, take the median.
      if (n_elements(offset) gt 1) then offset = median(offset,/even)

      ;; Expect the offset from the previous row to be with in several
      ;; pixels.  If it is not, reject it as a spurious result and sub
      ;; in the previous value.
      ;if not(first) then begin
      ;   if (abs(offset) gt 3) then begin
      ;      offset = 0.
      ;      conf   = 0.
      ;   endif
      ;endif

      ; Record array coordinates of max point and define next
      ; nominal position
      if xline then begin
         if first then begin
            xtrace0=[xnom]
            ytrace0=[ynom+offset]
            merit0  =[conf]
            first=0
         endif else begin
            xtrace0=[xtrace0,xnom]
            ytrace0=[ytrace0,ynom+offset]
            merit0  =[merit0,conf]
         endelse
         xnom=round(xnom+xstep)
         ynom=round(ynom+offset)
         stepstaken = stepstaken+1
      endif else begin
         if first then begin
            xtrace0=[xnom+offset]
            ytrace0=[ynom]
            merit0  =[conf]
            first=0
         endif else begin
            xtrace0=[xtrace0,xnom+offset]
            ytrace0=[ytrace0,ynom]
            merit0  =[merit0,conf]
         endelse
         xnom=round(xnom+offset)
         ynom=round(ynom+ystep)    
         stepstaken = stepstaken+1
      endelse

      ; Check to see if we've hit a limit.  If so, return to the starting
      ; position and continue in the negative direction, or terminate if
      ; we've already gone in the negative direction.  Default is to proceed.
      limit = 0
      ; Edge of array or a user defined limit?
      ;if (xnom-xside lt xmin or xnom+xside gt xmax or $
      ;    ynom-yside lt ymin or ynom+yside gt ymax) then begin
      if (xnom lt xmin or xnom gt xmax or $
          ynom lt ymin or ynom gt ymax) then begin
         limit = 1
      endif else begin
         if ((xline and (xnom-xside lt xmin or xnom+xside gt xmax)) or $
             (yline and (ynom-yside lt ymin or ynom+yside gt ymax))) then begin
            limit = 1
         endif else begin
            if (keyword_set(coord) and keyword_set(cbound)) then begin
               next_xlo = (xnom - xside) > xmin
               next_xhi = (xnom + xside) < xmax
               next_ylo = (ynom - yside) > ymin
               next_yhi = (ynom + yside) < ymax
               if (nmed gt 1 and next_xlo ne next_xhi and $
                   next_ylo ne next_yhi) then begin
                  next_coord_slice = median(coord(next_xlo:next_xhi,$
                                             next_ylo:next_yhi),dimension=meddim)
               endif else begin
                  next_coord_slice = coord(next_xlo:next_xhi,next_ylo:next_yhi)
               endelse
               next_bound_ls = where(next_coord_slice ge cmin and $
                                     next_coord_slice le cmax,n_next_bound)
               if (n_next_bound eq 0) then limit = 1
            endif
         endelse
      endelse
      ; Maximum number of steps?
      if keyword_set(nsteps) then begin
         if (n_elements([nsteps]) eq 2) then begin
            ; Going in positive direction?
            if ((xstep gt 0) or (ystep gt 0)) then begin
               if (stepstaken gt nsteps(1)) then limit=1
            endif else begin ; Negative direction
               if (stepstaken gt nsteps(0)) then limit=1
            endelse
         endif else begin
            if (stepstaken gt nsteps) then limit=1
         endelse
      endif
      ; If fitting Gaussian profiles, within 1.5-sigma of the edge?
      if keyword_set(bygauss) then begin
         if xline then begin
            if (ynom - 1.5*terms(2) lt ymin or $
                ynom + 1.5*terms(2) gt ymax) then limit = 1
         endif else begin
            if (xnom - 1.5*terms(2) lt xmin or $
                xnom + 1.5*terms(2) gt xmax) then limit=1
         endelse
      endif

      ; End of line?  

      if limit then begin
         ; Have we been going in the positive direction?
         if ((xstep gt 0 or ystep gt 0) and donegative) then begin
            if xline then begin
               xstep = -xstep
               xnom  = round(xtrace0(0)+xstep)
               ynom  = round(ytrace0(0))
               stepstaken=1
            endif else begin
               ystep = -ystep
               xnom  = round(xtrace0(0))
               ynom  = round(ytrace0(0)+ystep)
               stepstaken=1
            endelse
            ; If going in the negative direction will take us beyond
            ; the allowed range, then go on to the next line.
            if (xnom lt xmin or xnom gt xmax or $
                ynom lt ymin or ynom gt ymax) then begin
               donext = 0
               xstep  = abs(xstep)
               ystep  = abs(ystep)
            endif
         endif else begin
               donext = 0 ; Already gone in both pos and neg directions.
               xstep  = abs(xstep) ; Restore original direction
               ystep  = abs(ystep) ; for the next line.
         endelse
      endif

   endwhile

   ; Sort the positions according to the orientation of the line
   if xline then sortorder=sort(xtrace0) else sortorder=sort(ytrace0)
   xtrace0 = xtrace0(sortorder)
   ytrace0 = ytrace0(sortorder)
   merit0  = merit0(sortorder)

   ; Record trace values for this line in the big arrays
   ;if (firstline) then begin
   ;   ntracepoints = n_elements(xtrace0)
   ;   xtrace = fltarr(ntracepoints,nlines)
   ;   ytrace = fltarr(ntracepoints,nlines)
   ;   ctrace = fltarr(ntracepoints,nlines) ; Auxiliary coordinate (optional)
   ;   merit  = fltarr(ntracepoints,nlines)
   ;   firstline = 0
   ;endif
   line_length = n_elements(xtrace0)
   xtrace(0:line_length-1,i) = xtrace0
   ytrace(0:line_length-1,i) = ytrace0
   merit(0:line_length-1,i)  = merit0

   ; Keep track of the longlest line traced
   if (line_length gt max_length) then max_length = line_length

   ; Get rid of negative x and y values
   xtrace0 = xtrace0 > 0
   ytrace0 = ytrace0 > 0

   ; Optional - Find the auxiliary coordinate for each trace point
   ; Use linear interpolation for non-integer x,y coords.
   if keyword_set(coord) then begin
      ctrace0 = 0.*xtrace0
      for k=0,line_length-1 do begin
         xround = round(xtrace0(k))
         yround = round(ytrace0(k))
         xfloor = floor(xtrace0(k))
         yfloor = floor(ytrace0(k))
         xceil  = ceil(xtrace0(k))
         yceil  = ceil(ytrace0(k))
         if (xfloor eq xceil and yfloor eq yceil) then begin
            ctrace0(k) = coord(xtrace0(k),ytrace0(k))
         endif else begin
            if xline then begin
               ;coord_slice = coord(xround,yround-1:yround+1)
               ;y_slice     = yround + [-1,0,1]
               coord_slice = coord(xround,yfloor:yceil)
               y_slice     = [yfloor,yceil]
               ctrace0(k)  = interpol(coord_slice,y_slice,ytrace0(k))
            endif else begin
               ;coord_slice = coord(xround-1:xround+1,yround)
               ;x_slice     = xround + [-1,0,1]
               coord_slice = coord(xfloor:xceil,yround)
               x_slice     = [xfloor,xceil]
               if (n_elements(coord_slice) ne n_elements(x_slice)) then begin
                  ctrace0(k) = -9999
                  merit(k)   = 0
               endif else begin
                  ctrace0(k)  = interpol(coord_slice,x_slice,xtrace0(k))
               endelse
            endelse
         endelse
      endfor
      ctrace(0:line_length-1,i) = ctrace0
   endif

   ; Plot trace coordinates
   x_plot = xtrace0
   y_plot = ytrace0
   if keyword_set(coord) then begin
      if xline then y_plot = ctrace0 else x_plot = ctrace0
   endif

   if not(keyword_set(noplot)) then $
      plot,x_plot,y_plot,xrange=[min(x_plot),max(x_plot)],$
         yrange=[min(y_plot),max(y_plot)],psym=3

endfor

; Trim arrays to the longest line actually traced.
xtrace = xtrace(0:max_length-1,*)
ytrace = ytrace(0:max_length-1,*)
ctrace = ctrace(0:max_length-1,*)
merit  = merit(0:max_length-1,*)

!p.multi=pmultisave

return
end

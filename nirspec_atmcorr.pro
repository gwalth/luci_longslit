Pro nirspec_atmcorr,object,star,corr_object,error=error,a0=a0,$
                    response=response,interactive=interactive

;+-----------------------------------------------------------------
;
; NIRSPEC_ATMCORR    6/2004 GDB
;
; Correct for atmospheric absorption in NIRSPEC (H band)
; low-resolution spectra.  
;
; Atmospheric absorption as a function of wavelength is taken from
; the the normalized spectrum of a reference star (typically A0, F,
; or G) which has been divided by tabulated standard values for
; the intrinsic flux of that spectral type (e.g., Meyer et al. 1998 -
; see NIRSPEC_DIVSTD.pro).  
;
; The normalized reference star spectrum (corrected for intrinsic
; features) is converted to a flux decrement via
;                 decrement = 1 - (corrected star spectrum).
; The decrement is then scaled to the airmass of the science object
; as        decr(object) = decr(star) * air(object) / air(star).
; Note that this is an approximation for small optical depths.  The
; object spectrum is then divided by 1-decr(object) to recover the
; unabsorbed spectrum.
;
; NOTE: A model transmission spectrum may be used in place of the 
; standard star.  It should also be in a FITS file with wavelength
; informaiton in the header.  If no airmass information is avilable
; in the header, the model is assumed to be at airmass = 1.  The 
; user may then adjust the airmass, as well as the resolution and
; velocity offset, by setting the /interactive keyword.
;
; Meyer et al. (1998) tabulated standards only cover 15154.8 - 
; 17860.0 AA.  To make up for this at bluest wavelengths, the 
; absorption spectrum as determined from an A0 star can be substituted
; for lambda < 15400 AA (where there have been some issues in getting
; a reasonable absorption spectrum from other spectral types).  A0
; stars cab be well approximated as featureless in this range.  The 
; flux decrement observed in the A) star will be similarly scaled to
; the airmass at which the object spectrum was taken.
;
; NOTE: Using an A0 star from a different night (and likely from a
; different part of the night) to cover the blue end doesn't seem to
; work very well.  The amount of absorption at those wavelengths will
; be sensitive to the amount of water vaopr in the air.
;
; The program does not attempt to correct for any heliocentric
; velocity shifts that may have been applied.  The largest possible 
; difference between two exposures is ~60 km/s, which at R~1600 is
; ~1/3 of a resolution element.  However, in general the stellar
; spectrum is taken at roughly the same place in the sky and at the
; same time as the object spectrum, so the heliocentric velocity
; corrections will be very similar.
;
; Optionally, the spectrum can also be corrected for a user-supplied
; response function.
;
; INPUTS
;       object  - Filename of object spectrum
;       star    - Filename of normalized stellar spectrum from which
;                 intrinsic stellar features have been removed.
;
; KEYWORDS
;       error   - Filename of object error spectrum.  Will be
;                 corrected for atmospheric absorption in the same
;                 manner as the object spectrum.  The resulting
;                 error spectrum will be written to a file with the
;                 same name as the object output, with ".fits" 
;                 replaced by "e.fits".
;       a0      - Not recommended!  Filename of normalized A0 star, to
;                 be used to compute absorption at lambda < 15400.
;                 The default is to use the input stellar spectrum at all
;                 wavelengths. 
;       response - Filename of a response function (fits file), giving
;                 instrumental response as a function of wavelength, by
;                 which output spectrum will be divided.  Response
;                 function should cover minimum and maximum
;                 wavelengths covered by object spectrum.  Values for
;                 the reponse function will be interpolated onto the
;                 wavelength array of the object.
;       interactive - If set, user may manually adjust the velocity
;                 offset and intensity of the stellar spectrum 
;                 (atmospheric trnamission).  The use may also smooth
;                 the stellar spectrum to achieve a lower resolution.
;
; OUTPUTS
;       corr_object - Filename of corrected object spectrum.  Program
;                 will automatically append '.fits'.
;
; HISTORY
;       Written 6/28/2004 GDB
;       INTERACTIVE keyword added 1/17/2005 GDB
;-----------------------------------------------------------------------

if (n_params() lt 3) then begin
   print,'CALLING SEQUENCE: nirspec_atmcorr,object,star,corr_object,'
   print,'                        error=error,a0=a0,response=response,'
   print,'                        interactive=interactive'
   return
endif

c = 2.998e5

;;; Read in files
obj = rflam(object,hobj,wobj)
str = rflam(star,hstr,wstr)
if keyword_set(error) then err = rflam(error,herr,werr) $
   else err = 1 + 0*obj
if keyword_set(a0) then bluestr = rflam(a0,hbluestr,wbluestr) $
 else begin
   bluestr  = str
   hbluestr = hstr
   wbluestr = wstr
endelse

;;; Obtain airmasses - Default is 1.
am_obj = sxpar(hobj,'AIRMASS',count=n_am_obj)
am_str = sxpar(hstr,'AIRMASS',count=n_am_str)
am_bluestr = sxpar(hbluestr,'AIRMASS',count=n_am_bluestr)
if (n_am_obj eq 0) then am_obj = 1.
if (n_am_str eq 0) then am_str = 1.
if (n_am_bluestr eq 0) then am_bluestr = 1.

;;; Convert stellar spectra to flux dectrements
strdec     = 1. - str
bluestrdec = 1. - bluestr

;;; Scale decrements to airmass of object.  If airmasses are
;;; available, then ignore airmass corrections.
if (n_am_obj eq 0 or n_am_str eq 0) then scale = 1. else $
   scale = am_obj/am_str
if (n_am_obj eq 0 or n_am_str eq 0) then bluescale = 1. else $
   bluescale = am_obj/am_bluestr
dec     = strdec * scale
bluedec = bluestrdec * bluescale

;;; Default is no velocity shift and no smoothing (can be adjusted in
;;; interactive mode)
shift  = 0.0  ; Velocity shift (km/s)
smooth = 0.0  ; Smoothing (Gaussian kernal FWHM, km/s)

;;; Interpolate the decrememt arrays onto the object wavelengths
interp_dec     = interpol(dec,wstr,wobj)
edge_ls        = where(wobj gt max(wstr) or wobj lt min(wstr),n_edge)
if (n_edge ne 0) then interp_dec(edge_ls) = 0. 
interp_bluedec = interpol(bluedec,wbluestr,wobj)
blueedge_ls    = where(wobj gt max(wbluestr) or wobj lt min(wbluestr),n_blueedge)
if (n_blueedge ne 0) then interp_bluedec(blueedge_ls) = 0. 

;;; Create a master flux decrements array. using the flux decrement
;;; from the A0 star (if any) for lambda < 15400.
bigdec = 0.*obj
redls  = where(wobj ge 15400.,complement=bluels)
bigdec(redls)  = interp_dec(redls)
bigdec(bluels) = interp_bluedec(bluels)

;bigdec = 0.*obj
;wmax   = max(wobj) < max(wstr)
;strls  = where(wstr ge 15400. and wstr le wmax)
;objls  = where(wobj ge 15400. and wobj le wmax)
;bigdec(objls) = dec(strls)
;wmin   = min(wobj) > min(wbluestr)
;bluestrls = where(wbluestr ge wmin and wbluestr lt 15400.)
;blueobjls = where(wobj ge wmin and wobj lt 15400.)
;bigdec(blueobjls) = bluedec(bluestrls)

;;; Divide object and error by expected throughput
objcorr = obj / (1.-bigdec)
errcorr = err / (1.-bigdec)

;;; Interactively adjust velocity offset, intensity, and resolution
;;; of atmospheric transmission function (optional).
if keyword_set(interactive) then begin
   key = ''
   obj_pixsz     = wobj(1) - wobj(0)  ; Pixel size for object
   str_pixsz     = wstr(1) - wstr(0)  ; Pixel size for atmospheric transmission
   bluestr_pixsz = wbluestr(1) - wbluestr(0)
   old_p  = !p
   !p.multi = [0,1,4]
   xr = [min(wobj),max(wobj)]
   delta_y = max(obj) - min(obj)
   yr = [min(obj),max(obj)+0.2*delta_y]
   while (key ne 'q') do begin
      ; Plot atmospheric trnamission, uncorrected object spectrum, 
      ; and corrected object spectrum
      ; Ignoring A0 star in the blue for now!
      plot,wobj,objcorr,xrange=xr,yrange=yr,title='Corrected',charsize=1.5
      plot,wobj,obj,xrange=xr,yrange=yr,title='Uncorrected',charsize=1.5     
      plot,wobj,1-scale*bigdec,xrange=xr,$
         title='Interpolated Atmospheric Transmission',charsize=1.5 
      plot,wstr,1-dec,xrange=xr,title='Raw Atmospheric Transmission',charsize=1.5 
      print,'Scale  = ',strtrim(scale,2)
      print,'Shift  = ',strtrim(shift,2)+' km/s'
      print,'Smooth = ',strtrim(smooth,2)+ 'km/s'
      ; Get keyboard input
      print,'Enter key to adjust atmospheric transmission (? for help)'
      key = get_kbrd(1)
      case key of
         '?': begin
                 print,'j : Increase velocity shift by 5 km/s'
                 print,'l : Decrease velocity shift by 5 km/s'
                 print,'x : Perform cross-correlation to find best velocity shift'
                 print,'i : Increase scale factor by 0.05'
                 print,'k : Decrease scale factor by 0.05'
                 print,'r : Increase smoothing by 10 km/s'
                 print,'e : Decrease smoothing by 10 km/s'
                 print,'s : Enter smoothing FWHM'
                 print,'? : Print help'
              end
         'j': shift = shift - 5.
         'l': shift = shift + 5.
         'i': scale = scale + 0.05
         'k': scale = scale - 0.05
         'e': smooth = (smooth - 10.) > 0.
         'r': smooth = smooth + 10.
         's': read,smooth,prompt='Enter smoothing kernel FWHM (km/s): '
         'x': begin
                 ; Cross-correlate smoothed atmospheric transmission
                 ; with object spectrum to obtain best velocity shift.
                 xcorr_ls = where(wobj ge 14610. and wobj le 15130,n_xcorr)
                 if (n_xcorr eq 0) then begin
                    print,'Cross-correlation wavelength range not covered.' 
                 endif else begin
                    obj_xcorr  = obj(xcorr_ls)
                    dec_xcorr  = bigdec(xcorr_ls)
                    ; Sub-sample
                    x    = findgen(n_xcorr)
                    xsub = findgen(10.*n_xcorr)/10.
                    sub_obj_xcorr = interpol(obj_xcorr,x,xsub)
                    sub_dec_xcorr = interpol(dec_xcorr,x,xsub)
                    lag = findgen(401.) - 20  ; Lag, in 0.1 pixels
                    xcorr = c_correlate(sub_obj_xcorr,1-sub_dec_xcorr,lag)
                    opt_lag = lag(where(xcorr eq max(xcorr))) / 10.
                    opt_lag = mean(opt_lag)
                    ; Translate optimum lag into a velocity shift
                    temp_shift = c*opt_lag*obj_pixsz / 14870.
                    shift = shift - temp_shift
                 endelse
              end
         'q': print,'Done.'
         else: begin
                  print,'Unknown key'
                  key = 'unknown'
               end
      endcase
      if (key ne 'q' and key ne 'unknown' and key ne '?') then begin
         ; Smooth atmospheric transmission
         if (smooth gt 0.) then begin
            ker_fwhm = (smooth/c)*mean(wobj)/str_pixsz
            sm_dec   = convol(dec,gaussker(fwhm=ker_fwhm))
         endif else sm_dec = dec
         ; Shift wavelengths
         w_shift  = wstr*(1+shift/c)
         ; Interpolatet onto object wavelengths
         bigdec = interpol(sm_dec,w_shift,wobj)
         ; Repeat above steps for blue side,  if necessary
         if keyword_set(a0) then begin
            if (smooth gt 0.) then begin
               blueker_fwhm = (smooth/300000.)*mean(wobj)/bluestr_pixsz
               sm_bluedec   = convol(bluedec,gaussker(fwhm=blueker_fwhm))
            endif else sm_bluedec = bluedec
            wblue_shift  = wbluestr*(1+shift/c)
            interp_bluedec = interpol(bluedec,wblue_shift,wobj)
            bigdec(bluels) = interp_bluedec(bluels)
         endif
         ; Divide object and error by scaled expected throughput
         objcorr = obj / (1.-scale*bigdec)
         errcorr = err / (1.-scale*bigdec)
      endif
   endwhile
endif

;;; Optional - Divide by reponse function
if keyword_set(response) then begin
   resp = rflam(response,hresp,wresp)
   respmatch = interpol(resp,wresp,wobj)
   objcorr   = objcorr / respmatch
   errcorr   = errcorr / respmatch
endif

;;; Update headers and write fits files
sxaddpar,hobj,'ATMSTAR',star,'Stellar spectrum used for atmosph abs correction.' 
sxaddpar,hobj,'AIR_STAR',am_str,'Airmass at which stellar spectrum taken.'
if keyword_set(a0) then begin
   sxaddpar,hobj,'A0_STAR',star,'A0 spectrum used for atmosph abs below 15400 AA.'
   sxaddpar,hobj,'AO_AIR',am_bluestr,'Airmass at which A0 spectrum taken.'
endif
if keyword_set(response) then $
   sxaddpar,hobj,'RESPONS',response,'Reponse function by which spectrum divided.'
sxaddpar,hobj,'ATMSCALE',scale,'Scale factor applied to atmospheric correction'
sxaddpar,hobj,'ATMSHIFT',shift,'Velocity shift applied to atmospheric correction (km/s)'
sxaddpar,hobj,'ATMSMTH',smooth,'Smoothing applied to atmospheric corretion (km/s)'

writefits,corr_object+'.fits',float(objcorr),hobj

if keyword_set(error) then begin
   sxaddpar,herr,'ATMSTAR',star,'Stellar spectrum used for atmosph abs correction.' 
   sxaddpar,herr,'AIR_STAR',am_str,'Airmass at which stellar spectrum taken.'
   if keyword_set(a0) then begin
      sxaddpar,herr,'A0_STAR',star,'A0 spectrum used for atmosph abs below 15400 AA.'
      sxaddpar,herr,'AO_AIR',am_bluestr,'Airmass at which A0 spectrum taken.'
   endif
   if keyword_set(response) then $
      sxaddpar,herr,'RESPONS',response,'Reponse function by which spectrum divided.'
   writefits,corr_object+'e.fits',float(errcorr),herr
   sxaddpar,herr,'ATMSCALE',scale,'Scale factor applied to atmospheric correction'
   sxaddpar,herr,'ATMSHIFT',shift,'Velocity shift applied to atmospheric correction (km/s)'
   sxaddpar,herr,'ATMSMTH',smooth,'Smoothing applied to atmospheric corretion (km/s)'
endif

return
end

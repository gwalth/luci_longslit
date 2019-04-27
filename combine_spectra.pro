Pro combine_spectra,speclist,outspec,outerr,outwave,fluxtag=fluxtag,$
                    errtag=errtag,wavetag=wavetag,bin=bin,loglin=loglin,$
                    medrange=medrange,scale=scale,noscale=noscale,rej=rej,$
                    altrej=altrej,chop=chop,flux_arr=flux_arr,err_arr=err_arr,$
                    wave_arr=wave_arr

;+------------------------------------------------------------------------
;
; COMBINE_SPECTRA     6/2004
;
; Generic proceedure for combining 1D spectra.  Input FITS files
; should contain wavelength information in their headers, unless
; separate wavelength files are provided.  Input spectra need not
; cover the same wavelength interval, however it is assumed that the
; binning is the same for all spectra and that the wavelength bins are
; equally centered and spaced.
;
; May also directcly input arrays of flux, errors, and wavelength
; using the flux_arr, err_arr, and wave_arr keyords.
;
; By default, program will multiply each input spectrum by a scalar
; such that its median is equal to that of the spectrum with the
; highest median.  Scaling may be turned off (to preserve absolute
; flux calibration) by setting the 'noscale' keyword.
; 
;
; INPUTS
;        speclist   - Array of filenames of spectra to be combined;
;                     may either be full filenames or root names
;                     if the fluxtag keyword is set.  Headers
;                     should contain wavelength information, unless the
;                     wavetag keyword is set.
;
; KEYWORDS
;        fluxtag    - Suffix for flux filenames, to be appended to
;                     entries in speclist, e.g., '_flux.fits'.
;        errtag     - Suffix for error filenames, to be appended to
;                     entries in speclist, e.g., '_err.fits'.
;                     Error arrays are used to compute weighted means.
;                     Non-positive errors are treated as missing data.
;        wavetag    - Suffix for wavelength filenames, to be appended to
;                     entries in speclist, e.g., '_wave.fits'.  If this
;                     keyword is supplied, wavelengths are read from
;                     a file rather than determined from information
;                     in the FITS headers of the flux files.  Pixels
;                     with wavelengths <= 0 are treated as missing data.
;        bin        - Wavelength bin size.  Default is to use the wavelength
;                     bin size contained in the header or the bin size
;                     measured from the input wavelength arrays.  Currently
;                     assumes infinitely narrow input pixels (no sub-pixel
;                     interpolation).
;        loglin     - If set, binning is loglinear, and bin is assumed
;                     to be in km/s.  Must also set the bin keyword.
;        medrange   - 2-elements array giving the minimum and maximum
;                     wavelengths within which to compute the median
;                     flux value in each individual spectrum before
;                     combining.
;        scale      - Array of scale values, one per input flux array,
;                     by which to multiply flux and error values prior to
;                     combining.  Default is to comput the median in
;                     in each array and scale the arrays such that the
;                     medians are equal.
;        noscale    - Do not scale input flux and error values before 
;                     combining.
;        rej        - Scalar giving the number of standard deviations
;                     for iteratively rejecting outliers.  If flux
;                     errors are provided, a pixel will be rejected it
;                     its values is rej * the error in flux value for
;                     that pixel value away from the mean.  If errors are not
;                     not provided, then a pixel will be rejected if it
;                     is rej standard deviations away from the mean.
;                     Default is no rejection.
;        altrej     - Scalar, alternate rejection value.  A pixel
;                     will be rejected if it is altrej * the error
;                     in the mean away from the mean.  Default is
;                     no rejection.  Must also set the rej keyword
;                     (may be a high value) and input flux errors.
;        chop       - Two-element array giving the fraction of the
;                     input wavelength range to be masked at the
;                     beginning and end of each spectrum.  If the flux
;                     data is in a multi-row array (e.g., multiple
;                     orders), the beginning and end of each row will
;                     be chopped.  Default is no explicit masking.
;        flux_arr, err_arr, wave_arr - 2D or 3D arrays of flux, error, and
;                     wavelength (1 row per spectrum).  If provided,
;                     program will use these values instead of reading
;                     arrays from FITS files.  Input speclist should
;                     be set to a dummy value.
;
; OUTPUTS
;        outspec    - Variable or filename of combined spectrum
;        outerr     - (Optional) Variable or filename of 1-sigma
;                     combined error array.  If errtag keyword is set,
;                     then errors in combined spectrum are computed
;                     from input errors.  Otherwise, the standard
;                     deviation in the input spectra at each
;                     wavelength bin are returned. 
;        outwave    - (Optional) Variable or filename of wavelength
;                     array for combined spectrum.  Note: If writing to
;                     files, headers for outspec and outerr will
;                     already contain wavelength information.
;
; Side effects: Output flux and error arrays are converted to single-
;               precision floating point.
;
; HISTORY
;        Written 6/18/2004 GDB
;        WAVETAG,MEDRANGE,REJ,CHOP keywords added, modified to allow
;           2D input spectra 3/2/2005 GDB
;        Modified method for rejecting outliers 6/13/2005 GDB
;        LOGLIN, FLUX_ARR, ERR_ARR, and WAVE_ARR keywords added 5/18/2007 GDB
;----------------------------------------------------------------------- 

if (n_params() lt 2) then begin
   print,'CALLING SEQUENCE: '
   print,'      combine_specrtra,speclist,outspec,outerr,fluxtag=fluxtag,'
   print,'                       errtag=errtag,wavetag=wavetag,bin=bin,'
   print,'                       /loglin,medrange=medrange,scale=scale,'
   print,'                       /noscale,rej=rej,altrej=altrej,chop=chop,'
   print,'                       flux_arr=flux_arr,err_arr=err_arr,'
   print,'                       wave_arr=wave_arr'
   print,''  
   print,'   Note: If flux_arr, err_arr, and wave_arr keywords are set,'
   print,'         then speclist should be set to a dummy value.'
   return
endif

;;; Get flux, error, and wavelength arrays

if keyword_set(flux_arr) then begin

   ; Direct input from keywords (optional)
   fluxarr = flux_arr
   if not(keyword_set(wave_arr)) then begin
      print,'If using flux_arr keyword, must also set wave_arr keyword.'
      return
   endif

   ; If these are 2D arrays, embed them in 3D arrays with second
   ; dimension length = 1.
   sz = size(fluxarr)
   case sz(0) of
      1:    begin
               nspec   = 1
               xsz_max = sz(1)
               ysz_max = 1
            end
      2:    begin
               nspec   = sz(2)
               xsz_max = sz(1)
               ysz_max = 1
            end
      3:    begin
               nspec   = sz(4)
               xsz_max = sz(1)
               ysz_max = sz(3)
            end
      else: begin
               print,'Input arrays may have at most 3 dimensions.'
               return
            end
   endcase
   fluxarr = fltarr(xsz_max,ysz_max,nspec)
   errarr  = fltarr(xsz_max,ysz_max,nspec) - 1  ; Initialize with dummy values
   wavearr = dblarr(xsz_max,ysz_max,nspec)
   goodarr = intarr(xsz_max,ysz_max,nspec) + 1
   fluxarr(*) = flux_arr
   if keyword_set(err_arr) then begin
      have_errs = 1
      errarr(*) = err_arr
   endif else have_errs = 0
   wavearr(*) = wave_arr
   ; Default bin size
   binsz = abs(wavearr(1,0,0) - wavearr(0,0,0))
       
endif else begin

   ;;; Read in flux, error, and wavelength information from files
   ;;; (default)

   ;;; 0. Determine full filenames
   
   nspec = n_elements(speclist)
   
   if keyword_set(fluxtag) then fluxlist = speclist + fluxtag $
      else fluxlist = speclist
   ;Read in error files?
   if keyword_set(errtag) then begin
      errlist   = speclist + errtag
      have_errs = 1
   endif else begin
      have_errs = 0
   endelse
   ; Read in wavelength files?
   if keyword_set(wavetag) then begin
      wavelist  = speclist + wavetag
      have_wave = 1
   endif else begin
      have_wave = 0
   endelse
   
   ;print,'*** Rounding wavelengths to 0.01 ***'
   
   ;;; 1. Determine minimum and maximum wavelengths.  Double check that
   ;;;    the wavelength bin size is the same for all spectra.  Also get
   ;;;    the maximum x (and possibly y) dimensions of the input arrays.
   
   ; Get starting values from first input.
   if have_wave then begin
      wav = readfits(wavelist(0),hwav)
      ; ****
      if (min(wav) ge (2L^15)) then begin
         print,'*** Subtracting 2^15 from wavelength ***'
         wav=wav-2L^15
      endif
      ; ****
      ;wav = round(wav*100)/1d2
      pos_ls = where(wav gt 200,n_pos)
      minwav = min(wav(pos_ls))
      maxwav = max(wav(pos_ls))
      binsz   = (wav(pos_ls))(1) - (wav(pos_ls))(0)
   endif else begin
      flux   = rflam(fluxlist(0),hflux,wav,/double)
      minwav = min(wav)
      maxwav = max(wav)
      binsz   = sxpar(hflux,'CDELT1')
   endelse
   sz = size(wav)
   xsz_max = sz(1)
   if (sz(0) eq 2) then ysz_max = sz(2) else ysz_max = 1
   
   ; Loop through other spectra
   for i=1,nspec-1 do begin
      if have_wave then begin
         wav = readfits(wavelist(i),hwav)
         ; ****
         if (min(wav) ge 2L^15) then begin
            print,'*** Subtracting 2^15 from wavelength ***'
            wav = wav - 2L^15
         endif
         ; ****
         ;wav = round(wav*100)/1d2
         pos_ls = where(wav gt 200,n_pos)
         if (n_pos gt 0) then begin
            if (min(wav(pos_ls)) lt minwav) then minwav = min(wav(pos_ls))
            if (max(wav(pos_ls)) gt maxwav) then maxwav = max(wav(pos_ls))
            temp_binsz = (wav(pos_ls))(1) - (wav(pos_ls))(0)
            if (abs(temp_binsz-binsz) gt (1e-2)*binsz and $
                not(keyword_set(loglin))) then begin
               print,'Warning: Input spectra have unequal size wavelength bins.'
            endif
         endif
      endif else begin
         flux = rflam(fluxlist(i),hdr,wav,/double)
         if (min(wav) lt minwav) then minwav = min(wav)
         if (max(wav) gt maxwav) then maxwav = max(wav)
         if (sxpar(hdr,'CDELT1') ne binsz) then begin
            print,'Warning: Input spectra have unequal size wavelength bins.'
         endif
      endelse
      sz = size(wav)
      if (sz(1) gt xsz_max) then xsz_max = sz(1)
      if (sz(0) eq 2) then begin
         if (sz(2) gt ysz_max) then ysz_max = sz(2)
      endif
   endfor
   
   ;;; 2. Create master arrays of flux and error values
   
   fluxarr = fltarr(xsz_max,ysz_max,nspec)
   ; Initialize error array with dummy values
   errarr  = fltarr(xsz_max,ysz_max,nspec) - 1
   wavearr = dblarr(xsz_max,ysz_max,nspec)
   ; Keep track of good pixels for each spectrum
   goodarr = intarr(xsz_max,ysz_max,nspec) + 1
   
   ; Read in arrays
   for i=0,nspec-1 do begin
      if have_wave then begin
         temp_flux = readfits(fluxlist(i),out_hdr)
         temp_wave = readfits(wavelist(i),htemp_wave)
         ; ****
         if (min(temp_wave) ge 2L^15) then begin
            print,'*** Subtracting 2^15 from wavelength ***'
            temp_wave = temp_wave - 2L^15
         endif
         ; ****
         ;temp_wave = round(temp_wave*100)/1d2
      endif else begin
         temp_flux = rflam(fluxlist(i),out_hdr,temp_wave)
      endelse
      temp_sz = size(temp_flux)
      temp_xsz = temp_sz(1)
      if (temp_sz(0) eq 2) then temp_ysz = temp_sz(2) else temp_ysz = 1
      fluxarr(0:temp_xsz-1,0:temp_ysz-1,i) = temp_flux
      wavearr(0:temp_xsz-1,0:temp_ysz-1,i) = temp_wave
      if have_errs then begin
         temp_err = readfits(errlist(i),htemp_err)
         errarr(0:temp_xsz-1,0:temp_ysz-1,i) = temp_err
      endif
   endfor

endelse

; Mask pixels with err <=0 or wavelength <= 200
neg_wave_ls = where(wavearr le 200,n_neg_wave)
if (n_neg_wave ne 0) then goodarr(neg_wave_ls) = 0
if have_errs then begin
   neg_err_ls = where(errarr le 0,n_neg_err)
   if (n_neg_err ne 0) then goodarr(neg_err_ls) = 0
endif

; Optional - Mask the ends of each spectrum or order
if keyword_set(chop) then begin
   if (n_elements(chop) ne 2) then begin
      print,'COMBINE_SPECTRA: Keyword chop must be a 2-element array.'
      return
   endif
   lo_chop = chop(0)
   hi_chop = chop(1)
   for i=0,nspec-1 do begin
      for j=0,ysz_max-1 do begin
         temp_wavearr = wavearr(*,j,i)
         temp_goodarr = goodarr(*,j,i)
         ls = where(temp_goodarr eq 1,n_ls)
         if (n_ls ne 0) then begin
            temp_wmin = min(temp_wavearr(ls))
            temp_wmax = max(temp_wavearr(ls))
            delta_w   = temp_wmax - temp_wmin
            lo_bound  = temp_wmin + lo_chop*delta_w
            hi_bound  = temp_wmax - hi_chop*delta_w
            chop_ls   = where(temp_wavearr le lo_bound or $
                              temp_wavearr ge hi_bound,n_chop)
            if (n_chop ne 0) then goodarr(chop_ls,j,i) = 0
         endif
      endfor
   endfor
endif


;;; 3. Scale spectra such that the median values are equal (optional).

if keyword_set(scale) then scalelist = scale else begin
   if keyword_set(medrange) then begin
      med_wlo = min(medrange)
      med_whi = max(medrange)
   endif else begin
      med_wlo = min(wavearr)
      med_whi = max(wavearr)
   endelse
   medlist  = fltarr(nspec)
   for i=0,nspec-1 do begin
      temp_fluxarr = fluxarr(*,*,i)
      temp_wavearr = wavearr(*,*,i)
      temp_goodarr = goodarr(*,*,i)
      temp_ls = where(temp_wavearr ge med_wlo and temp_wavearr le med_whi and $
                     temp_goodarr eq 1)
      medlist(i) = median(temp_fluxarr(temp_ls),/even)
   endfor
   maxmed    = max(medlist)
   scalelist = maxmed / medlist
endelse

if not(keyword_set(noscale)) then begin
   for i=0,nspec-1 do begin
      fluxarr(*,*,i) = fluxarr(*,*,i) * scalelist(i)
      if (have_errs) then errarr(*,*,i) = errarr(*,*,i) * scalelist(i)
   endfor
endif   

;;; 4. Computed weighted (optional) mean in each wavelength bin 

; Weight flux values by their inverse variance (optional)
if (have_errs) then weightarr = 1./(double(errarr)^2) else $
   weightarr = 1. + 0.*fluxarr

; Keep track of pixels that have been rejected
maskarr = 0*goodarr
initial_mask_ls = where(goodarr eq 0,n_initial_mask,$
                        complement=initial_include_ls)
if (n_initial_mask ne 0) then maskarr(initial_mask_ls) = 1

; Compute output wavelengths
maxwav_comb = max(wavearr(initial_include_ls))
minwav_comb = min(wavearr(initial_include_ls))
if keyword_set(loglin) then begin
   if keyword_set(bin) then binsz = double(bin) else begin
      print,'Must set bin keyword if /loglin keyword is set.'
      return
   endelse
   c = 2.99792458d5
   logbinsz    = alog10(binsz/c + 1)
   lamcent     = minwav_comb
   loglamcent  = alog10(lamcent)
   j_min       = floor((alog10(minwav_comb/lamcent)/logbinsz) + 0.5)
   j_max       =  ceil((alog10(maxwav_comb/lamcent)/logbinsz) - 0.5)
   nwavbins    = j_max - j_min + 1
   j_list      = j_min + dindgen(nwavbins)
   out_wave    = 10^(loglamcent + j_list*logbinsz)
   out_wave_lo = 10^(loglamcent + (j_list-0.5)*logbinsz)
   out_wave_hi = 10^(loglamcent + (j_list+0.5)*logbinsz)
endif else begin
   if keyword_set(bin) then binsz = double(bin)
   nwavbins    = round((maxwav_comb - minwav_comb) / binsz) + 1
   out_wave    = minwav_comb + binsz*dindgen(nwavbins)
   out_wave_lo = minwav_comb + binsz*(dindgen(nwavbins) - 0.5)
   out_wave_hi = minwav_comb + binsz*(dindgen(nwavbins) + 0.5)
endelse
out_spec = 0.*out_wave
out_err  = 0.*out_wave

; Speed things up by dividing the spectrum into wavelength regions
; *** Should replace this with HISTOGRAM method. ***
bins_per_region = round(sqrt(nwavbins))
n_regions = ceil((1.0*nwavbins)/bins_per_region)
for r=0,n_regions do begin
   i_lo = r*bins_per_region < (nwavbins - 1)
   i_hi = ((r+1)*bins_per_region) < (nwavbins - 1)
   ;region_wlo = minwav_comb + (i_lo - 0.5)*binsz
   ;region_whi = minwav_comb + (i_hi + 0.5)*binsz
   region_wlo = out_wave_lo(i_lo)
   region_whi = out_wave_hi(i_hi)
   region = where(wavearr ge region_wlo and wavearr le region_whi and $
                  goodarr eq 1,n_region)
   if (n_region gt 0) then begin
      region_flux   = fluxarr(region)
      region_wave   = wavearr(region)
      region_weight = weightarr(region)
      for i=i_lo,i_hi do begin
         ; Identify pixels in this wavelength bin
         w = out_wave(i)
         ;covered_ls = where(region_wave ge w-0.5*binsz and $
         ;                   region_wave lt w+0.5*binsz,n_covered)
         covered_ls = where(region_wave ge out_wave_lo(i) and $
                            region_wave lt out_wave_hi(i),n_covered)
         ; Compute mean and iteratively reject outliers (optional)
         if (n_covered ne 0) then begin
            ; One-at-a-time rejection scheme
            done = 0
            while not(done) do begin
               fluxvals  = region_flux(covered_ls)
               weights   = region_weight(covered_ls)
               errvals   = 1./sqrt(weights)
               temp_mean = total(fluxvals * weights) / total(weights)
               if (have_errs) then begin                
                  temp_err = sqrt(1./total(weights))    
               endif else begin                         
                  if (n_covered gt 1) then $            
                     temp_err = stddev(fluxvals) else $ 
                     temp_err = -1.                     
               endelse 
               if (keyword_set(rej) and n_covered gt 1) then begin
                  ; Identify the largest outlier
                  resids    = abs(fluxvals - temp_mean)
                  max_resid = max(resids)
                  ; Compute the mean without this pixel
                  keep     = where(resids ne max_resid,n_keep,$
                                   complement=discard,ncomplement=n_discard)
                  ; Reject only one pixel at a time
                  if (n_discard gt 1) then begin
                     keep      = [keep,discard(1:n_discard-1)]
                     n_keep    = n_elements(keep)
                     discard   = discard(0)
                     n_discard = n_elements(discard)
                  endif
                  if (n_discard ne n_covered and n_discard ne 0) then begin
                     new_mean = total(fluxvals(keep) * weights(keep)) / $
                                   total(weights(keep))
                     if (have_errs) then begin
                        new_err = sqrt(1./total(weights(keep)))
                     endif else begin
                        if (n_keep gt 1) then $
                           new_err = stddev(fluxvals(keep)) else $
                           new_err = -1
                     endelse
                     new_resids = abs(fluxvals - new_mean)
                     ; If the pixel is more than an n-sigma outlier
                     ; from the new mean, then discard it.
                     if (have_errs) then pixel_err = errvals(discard) else $
                         pixel_err = new_err
                     discard_pixel = 0
                     if keyword_set(altrej) then begin
                        if (new_resids(discard) gt altrej*new_err) then $
                           discard_pixel = 1
                     endif
                     if (new_resids(discard) gt rej*pixel_err and $
                         pixel_err ne -1) then discard_pixel = 1
                     if (discard_pixel) then begin
                        ; Flag this pixel as discarded
                        maskarr(region(covered_ls(discard))) = 1
                        ; Update list of pixels in this wavelength bin
                        covered_ls = covered_ls(keep)
                        n_covered  = n_keep
                     endif else done = 1
                  endif else done = 1
               endif else done = 1
            endwhile
            out_spec(i) = temp_mean
            out_err(i)  = temp_err
         endif else begin
            out_spec(i) = 0
            out_err(i)  = -1
         endelse
      endfor
   endif else begin
      out_spec(i_lo:i_hi) = 0
      out_err(i_lo:i_hi)  = -1
   endelse
endfor

;;; 5. Write to files or return variables

; Convert output to single-precision floating point
out_spec = float(out_spec)
out_err  = float(out_err)

; Modify header for output

if keyword_set(loglin) then begin
   sxaddpar,out_hdr,'CRVAL1',alog10(minwav_comb)
   sxaddpar,out_hdr,'CDELT1',alog10(binsz/c + 1)
   sxaddpar,out_hdr,'CD1_1',alog10(binsz/c + 1)
   sxaddpar,out_hdr,'CRPIX1',1
   sxaddpar,out_hdr,'CTYPE1','LINEAR'
   sxaddpar,out_hdr,'CTYPE2','PIXEL'
endif else begin
   sxaddpar,out_hdr,'CRVAL1',minwav_comb
   sxaddpar,out_hdr,'CDELT1',binsz
   sxaddpar,out_hdr,'CD1_1',binsz
   sxaddpar,out_hdr,'CRPIX1',1
   sxaddpar,out_hdr,'CTYPE1','LAMBDA'
   sxaddpar,out_hdr,'CTYPE2','PIXEL'
endelse
if not(keyword_set(flux_arr)) then begin
   sxaddhist,'Combined spectrum created from the following files:',out_hdr
   for i=0,nspec-1 do sxaddhist,'  '+fluxlist(i),out_hdr
endif

if (size(outspec,/tname) eq 'STRING') then begin
   writefits,outspec,out_spec,out_hdr
endif else outspec = out_spec

if (n_params() ge 2) then begin
   if (size(outerr,/tname) eq 'STRING') then begin
      writefits,outerr,out_err,out_hdr
   endif else outerr = out_err
endif

if (n_params() ge 3) then begin
   if (size(outwave,/tname) eq 'STRING') then begin
      writefits,outwave,out_wave,out_hdr
   endif else outwave = out_wave
endif

return
end

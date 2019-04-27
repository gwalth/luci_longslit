;-------------------------------------------------------------------------

Function lr_makevariance,data,object=object,sky=sky,flatfield=flatfield,$
                         flatvar=flatvar,darkframe=darkframe,$
                         objcoadds=objcoadds,darkcoadds=darkcoadds,$
                         gain=gain,readnoise=readnoise,biasnoise=biasnoise

; Calculate variance array based on object counts, sky counts, flat
; field value, flat field variance, dark counts, and read noise.  All
; inputs are optional keywords.  Variables are as defined in comments
; below, except for 'data' which is the array for which the variance
; is to be computed (and is used only to count the number of pixels).
; Gain is in e/ADU for calculating Poisson uncertainties.  Sky and object
; counts (in electrons) are assumed to obey Poisson stats.

print,'Contributions to variance...'

; Set gain
if keyword_set(gain) then gn = gain else gn = 1.

; Set flat-field (arbitrary units)
if keyword_set(flatfield) then ffield = flatfield else ffield = 1. + 0.*data

; Convert read noise from electrons to flat-fielded ADU
if keyword_set(readnoise) then rnoise = float(readnoise) / (gn * ffield) $
   else rnoise = 0.*data  ; Read noise rms in ADU
rnoisevar = rnoise^2      ; Read noise variance in ADU
print,'Read',minmax(rnoisevar)

var  = 0.*data

; Calculate contributions to variance in pixel values from uncertainties
; in the following.

; NOTE: Object, sky, etc. counts already normalized!

; Object
if keyword_set(object) then begin
   objvar = abs(object)/(gn*ffield)
   print,'Object:',minmax(objvar)
   var = var + objvar
endif else object = 0.*data
; Sky
if keyword_set(sky) then begin
   skyvar = abs(sky)/(gn*ffield)
   print,'Sky:',minmax(skyvar)
   var = var + skyvar
endif else sky = 0.*data
; Flat
if keyword_set(flatvar) then begin
   ffieldvar = ((object + sky)^2)*flatvar/(ffield^2)
   print,'Flat',minmax(ffieldvar)
   var = var + ffieldvar
endif
; Read noise - object frame
if keyword_set(objcoadds) then ocoadds = objcoadds else ocoadds = 1.
var = var + ocoadds*rnoisevar
; Noise in 2D bias
if keyword_set(biasnoise) then begin
   biasnoisevar = (biasnoise/ffield)^2
   var = var + biasnoisevar
   print,'Bias',minmax(biasnoisevar)
endif
; Dark - Poisson noise
if keyword_set(darkframe) then begin
   darkvar = abs(darkframe)/(gn*ffield)
   print,'Dark:',minmax(darkvar)
   var = var + darkvar
endif
; Read noise - dark frame.  Detect zero-valued keyword.
if (n_elements(darkcoadds) ne 0) then begin
  var = var + darkcoadds*rnoisevar
endif

; Edit out NaN pixels
num_ls = where(var ge min(var) and var le max(var),n_num,$
               complement=nan_ls,ncomplement=n_nan)
if (n_nan ne 0) then var(nan_ls) = max(var)

return,var
end

;-----------------------------------------------------------------------

Function lr_set_bkpoints,x,y,dispersion,bk_space,min_bk_space

; Place points break points.  Attempt to identify and place extra
; break points around skylines.  x=wavelength, y=sky counts

; Compute median sky counts in bins of wavelength
bin   = dispersion
h     = histogram(x,binsize=bin,omin=xmin,rev=x_ri)
n_bin = n_elements(h)
xbin  = xmin + bin*(dindgen(n_bin)+0.5)
ybin  = fltarr(n_bin)
for i=0,n_bin-1 do if (x_ri(i+1) gt x_ri(i)) then $
   ybin(i) = median(y(x_ri(x_ri(i):x_ri(i+1)-1)))

; Make sure bk_space is greater than min_bk_space
bk_space = (bk_space > min_bk_space)

; Buffer pixels at the edges of lines and in slowly varying regions
; with break points at an intermediate spacing
intermed_bk_space = min_bk_space + 0.3*(bk_space-min_bk_space)

; Compute local quantities for groups pixels
local_sz     = n_bin/10.     ; Number of bins per local region
nside        = 2
binpergroup  = 2*nside+1     ; Ideal number of bins per group
group_nbin   = fltarr(n_bin) ; Actual number of bins in each group
group_mean   = fltarr(n_bin)
group_median = fltarr(n_bin)
local_median = fltarr(n_bin)
local_sdev   = fltarr(n_bin)
for i=0,n_bin-1 do begin                                  
   ; Quantities in the local region (long stretch)
   j_lo = 0 > (i-local_sz/2)
   j_hi = (j_lo+local_sz) < (n_bin-1)
   local_median(i) = median(ybin(j_lo:j_hi))
   local_sdev(i)   = robust_sigma(ybin(j_lo:j_hi))
   ; Quantities for this group of pixels
   ilo = (i-nside) > 0
   ihi = (i+nside) < (n_bin-1)
   group_nbin(i)   = ihi-ilo+1
   group_mean(i)   = mean(ybin(ilo:ihi))
   group_median(i) = median(ybin(ilo:ihi))
endfor

; Default is the max break point spacing
bkspace_list = replicate(bk_space,n_bin)

; Groups that look like emission lines
emline_thresh = local_median+3*local_sdev/sqrt(group_nbin)
emline_ls = where(group_mean ge emline_thresh and $
                  group_median ge local_median+local_sdev,n_emline)
if (n_emline ne 0) then begin
   for j=0,n_emline-1 do begin
      ; Note: Adding extra pixels to the side of groups may cause
      ; ringing in the sky fit at line edges
      i   = emline_ls(j)
      ilo = (i-nside) > 0
      ihi = (i+nside) < (n_bin-1)
      bkspace_list(i)       = min_bk_space
      bkspace_list(ilo:ihi) = intermed_bk_space < bkspace_list(ilo:ihi)
   endfor
endif

; Groups that look like absorption lines
absline_thresh = local_median-3*local_sdev/sqrt(group_nbin)
absline_ls = where(group_mean le absline_thresh and $
                   group_median le local_median-local_sdev,n_absline)
if (n_absline ne 0) then begin
   for j=0,n_absline-1 do begin
      i   = absline_ls(j)
      ilo = (i-nside) > 0
      ihi = (i+nside) < (n_bin-1)
      bkspace_list(i)       = min_bk_space
      bkspace_list(ilo:ihi) = intermed_bk_space < bkspace_list(ilo:ihi)
   endfor
endif

; Additional regions that do not appear constant
region_nside = 10
for i=0,n_bin-1 do begin
   ilo = (i-region_nside) > 0
   ihi = (i+region_nside) < (n_bin-1)
   region_ybin       = ybin(ilo:ihi)
   region_group_mean = group_mean(ilo:ihi)
   region_bkspace    = bkspace_list(ilo:ihi)
   ; Bins where spacing is still the maximum
   max_ls = where(region_bkspace eq bk_space,n_max)
   if (n_max ge region_nside) then begin
      ; Determine whether scatter goes down with binning
      sdev_bin   = stddev(region_ybin(max_ls))
      sdev_group = stddev(region_group_mean(max_ls))
      if (sdev_bin/sdev_group le 0.7*sqrt(binpergroup)) then begin
         ; Scatter not decreasing with binning as expected for a
         ; constant sky bg
         bkspace_list(ilo:ihi) = region_bkspace < intermed_bk_space
      endif
   endif
endfor

; Create break spaces in each stretch of equal spacing
done_bkspace  = 0
first_stretch = 1
i_hi          = -1
while not(done_bkspace) do begin					  
   i_lo = i_hi + 1	   						  
   current_bkspace = bkspace_list(i_lo)					  
   i_hi = min(where(xbin ge xbin(i_lo) and $
              bkspace_list ne current_bkspace)) - 1			  
   if (i_hi eq -2) then begin 						  
      i_hi = n_bin-1							  
      done_bkspace = 1							  
   endif								  
   temp_minwave  = xbin(i_lo) - 0.5*bin					  
   temp_maxwave  = xbin(i_hi) + 0.5*bin					  
   temp_delta_wave  = temp_maxwave - temp_minwave			  
   temp_n_bkpts  = ceil(temp_delta_wave/current_bkspace)		  
   temp_bkspace  = temp_delta_wave/temp_n_bkpts				  
   temp_bkpoints = temp_minwave + $
                        temp_bkspace*(0.5+findgen(temp_n_bkpts))	  
   if first_stretch then begin						  
      bkpoints      = temp_bkpoints					  
      first_stretch = 0							  
   endif else begin							  
      bkpoints = [bkpoints,temp_bkpoints]				  
   endelse 	      							  
endwhile

return,bkpoints
end

;-----------------------------------------------------------------------

Function lr_refine_bkpoints,old_bkpoints,redchisq,$
                            min_bkspace=min_bkspace,status=status

; Identify regions where the fit to the sky is poor.  In each region,
; replace the current break points with twice as many, equally
; spaced between the bounding break points where the fit is good.

bkpoints = old_bkpoints

; Set minimum spacing between break points
if keyword_set(min_bkspace) then min_bk_space = min_bkspace else $
   min_bk_space = 0.1

; Identify regions where the fit is poor.
refine_ls = where(redchisq gt 2.0,n_refine)

; Flag that the minimum break point spacing has been reached
min_bk_space_reached = 0

; Need to refine break point grid?
if (n_refine eq 0) then begin
   print,'Fit looks okay.  No changes to break points.'
   status = 0
endif else begin
   status = 1
   print,'Bad fit detected.  Refining break points...'
   ; Flag points to be refined    
   refine_flag = 0.*bkpoints
   refine_flag(refine_ls) = 1
   ; Identify blocks of points to refine
   refine_done = 0
   hi_bound    = min(bkpoints) - 1
   ; In each block, double the number of break points, spacing
   ; them out evenly.
   while not(refine_done) do begin
      n_bkpoints = n_elements(bkpoints)
      temp_refine_ls = where(refine_flag eq 1 and bkpoints gt hi_bound,$
                             n_temp_refine)
      if (n_temp_refine eq 0) then begin
         refine_done = 1 
      endif else begin
         lo_ref_index = min(temp_refine_ls)
         if (lo_ref_index eq 0) then begin
            lo_bound_index = -1
            lo_bound       = 2*bkpoints(0) - bkpoints(1)
         endif else begin
            lo_bound_index = lo_ref_index - 1
            lo_bound       = bkpoints(lo_bound_index)
         endelse
         temp_no_refine_ls = where(refine_flag eq 0 and $
                                   bkpoints gt lo_bound,n_temp_no_refine)
         if (n_temp_no_refine eq 0) then begin
            hi_bound_index = -1
            hi_bound       = 2*bkpoints(n_bkpoints-1) $
                                - bkpoints(n_bkpoints-2)
            hi_ref_index   = n_elements(bkpoints) - 1   
         endif else begin
            hi_bound_index = min(temp_no_refine_ls)
            hi_bound       = bkpoints(hi_bound_index)
            hi_ref_index   = hi_bound_index - 1
         endelse
         n_to_replace   = hi_ref_index - lo_ref_index + 1
         total_space    = hi_bound - lo_bound
         n_new_bkpoints = 2*n_to_replace
         new_bk_space   = total_space / (1 + n_new_bkpoints)
         ; If the new spacing is smaller than the minimum break space,
         ; fit in as many points as are allowed.
         if (new_bk_space lt min_bk_space) then begin
            min_bk_space_reached = 1
            n_new_bkpoints = (floor(total_space/min_bk_space) - 1) > 0
            new_bk_space   = total_space / (1 + n_new_bkpoints)
         endif
         if (n_new_bkpoints gt 0) then begin
            new_bkpoints   = lo_bound + new_bk_space*(1+findgen(n_new_bkpoints))
            new_flag       = 1 + 0.*new_bkpoints
            if (lo_ref_index eq 0) then begin
               if (hi_ref_index eq n_bkpoints-1) then begin
                  bkpoints    = new_bkpoints
                  refine_flag = new_flag
               endif else begin
                  upper_bkpoints = bkpoints(hi_bound_index:n_bkpoints-1)
                  upper_flag     = refine_flag(hi_bound_index:n_bkpoints-1)
                  bkpoints       = [new_bkpoints,upper_bkpoints]
                  refine_flag    = [new_flag,upper_flag]
               endelse
            endif else begin
               lower_bkpoints = bkpoints(0:lo_bound_index)
               lower_flag     = refine_flag(0:lo_bound_index)
               if (hi_ref_index eq n_bkpoints-1) then begin
                  bkpoints    = [lower_bkpoints,new_bkpoints]
                  refine_flag = [lower_flag,new_flag]
               endif else begin
                  upper_bkpoints = bkpoints(hi_bound_index:n_bkpoints-1)
                  upper_flag     = refine_flag(hi_bound_index:n_bkpoints-1)
                  bkpoints = [lower_bkpoints,new_bkpoints,upper_bkpoints]
                  refine_flag = [lower_flag,new_flag,upper_flag]
               endelse
            endelse
         endif
      endelse
   endwhile
   ; Was the minimum break point spacing reached?
   if min_bk_space_reached then print,'Minimum break point spacing reached.'
   ; Were any new break points actually added?
   n_old = n_elements(old_bkpoints)
   n_new = n_elements(bkpoints)
   if (n_old eq n_new) then begin
      ; The old and new break point arrays have the same number of
      ; elements.  Have any changed?
      differences = where(n_new ne n_old,n_differences)
      if (n_differences eq 0) then begin
         print,'Could not add any new break points.'
         status = 0
      endif else status = 1
   endif else status = 1
endelse

return,bkpoints

end

;-----------------------------------------------------------------------

Function lr_traceobject,array,slitgrid,position,direction,$
                        range=range,bound=bound,tracerange=tracerange,$
                        slitrange=slitrange,$
                        objtype=objtype,spectrum=spectrum,speclam=speclam,$
                        specvar=specvar,wavegrid=wavegrid,pixtr=pixtr,$
                        slittr=slittr,slitfit=slitfit,waverange=waverange,$
                        nside=nside,noshow=noshow,$
                        fitterms=fitterms

; NOTE: range      - min/max x or y in dispersion direction over which to
;                    trace object
;       bound      - min/max x or y in spatial direction over which to 
;                    trace object
;       tracerange - Also min/max x or y in dispersion direction over which
;                    to trace object, set in call to longslit_reduce
;       slitrange  - min/max slit position over which to trace object

print,'Tracing object...'

; Default is no good trace
fitterms = -1

; What type of object
if keyword_set(objtype) then obj_type=objtype else obj_type='ordinary'

; Locate pixels near the object position
trace_ls = where(slitgrid ge position-15 and slitgrid le position+15,n_trace)
if (n_trace eq 0) then begin
   print,'LR_TRACEOBJECT: No pixels near object position.'
   return,-1
endif
indices = array_indices(array,trace_ls)
slit    = slitgrid(trace_ls)
counts  = array(trace_ls)
x_ls    = indices(0,*)
y_ls    = indices(1,*)

; Default is to trace over the entire array
sz = size(array)
if (direction eq 'x') then begin
   coord1       = x_ls  ; Coordinate parallel to the trace
   coord2       = y_ls  ; Coordinate perpendicular to the trace
   coord1_limit = [0,sz(1)]
   coord2_limit = [0,sz(2)]
   slit_corr_sz = sz(1)
endif else begin
   coord1       = y_ls
   coord2       = x_ls
   coord1_limit = [0,sz(2)]
   coord2_limit = [0,sz(1)]
   slit_corr_sz = sz(2)
endelse
slit_limit = [-1e6,1e6]

; Limit region of the array over which to trace object.  Default is to
; trace over the entire array.  (Naming conventions for 'bound' and
; 'range' are from traceline.pro.)
if keyword_set(bound)     then coord2_limit = [min(bound),max(bound)]
if keyword_set(range)     then coord1_limit = [min(range),max(range)]
if keyword_set(slitrange) then slit_limit   = [min(slitrange),max(slitrange)]
good_ls = where(coord1 ge coord1_limit(0) and coord1 le coord1_limit(1) and $
                coord2 ge coord2_limit(0) and coord2 le coord2_limit(1) and $
                slit   ge slit_limit(0)   and slit   le slit_limit(1),n_good)
if (n_good eq 0) then begin
   print,'LR_TRACEOBJECT: No pixels within specified bounds.'
   return,-1
endif
trace_ls = trace_ls(good_ls)
slit     = slit(good_ls)
counts   = counts(good_ls)
coord1   = coord1(good_ls)
coord2   = coord2(good_ls)

; Limit wavelength interval over which to trace object (optional)
if (keyword_set(waverange) and keyword_set(wavegrid)) then begin
   trace_wave = wavegrid(trace_ls)
   good_ls    = where(trace_wave ge min(waverange) and $
                      trace_wave le max(waverange),n_good)
   if (n_good eq 0) then begin
      print,'LR_TRACEOBJECT: No pixels within wavelength range.'
      return,-1
   endif
   trace_ls = trace_ls(good_ls)
   slit     = slit(good_ls)
   counts   = counts(good_ls)
   coord1   = coord1(good_ls)
   coord2   = coord2(good_ls)
endif

; Trace only over regions where the spectrum has significant flux 
; (optional).
if keyword_set(spectrum) then begin
   smspec       = median(spectrum,5)
   specerr      = sqrt(specvar)
   nsigma       = smspec / specerr
   trace_wave   = wavegrid(trace_ls)
   trace_nsigma = interpol(nsigma,speclam,trace_wave)
   good_ls      = where(trace_nsigma ge 3 and $
                        trace_wave ge min(speclam) and $
                        trace_wave le max(speclam),n_good)
   if (n_good eq 0) then begin
      print,'LR_TRACEOBJECT: No pixels with significant flux.'
      return,-1
   endif
   trace_ls = trace_ls(good_ls)
   slit     = slit(good_ls)
   counts   = counts(good_ls)
   coord1   = coord1(good_ls)
   coord2   = coord2(good_ls)
endif

; NOTE: If spectrum is noisy or has too many cosmic rays, median
; filtering of the 2D array parallel to the trace direction may
; be helpful.

; Fit Gaussians to small slices of the object profile along the
; slit
slice_sz = 3
h        = histogram(coord1,omin=min_coord1,binsize=slice_sz,$
                     locations=slice_start,rev=trace_ri)
pixtr    = slice_start + 0.5*slice_sz
n_slice  = n_elements(h)
slittr   = fltarr(n_slice)
fit_fwhm = fltarr(n_slice)
fit_nsig = fltarr(n_slice)
for i=0,n_slice-1 do begin
   if (h(i) lt 10) then continue
   slice_ls     = trace_ri(trace_ri(i):trace_ri(i+1)-1)
   slice_slit   = slit(slice_ls)
   slice_counts = counts(slice_ls)
   ; Omit strongest pixels when estimating fit parameters.  If these
   ; are actual object counts, they will be at the center anyway.
   order        = sort(slice_counts)
   tmp_slit     = (slice_slit(order))(0:0.95*h(i))
   tmp_counts   = (slice_counts(order))(0:0.95*h(i))
   cent_guess   = total(tmp_slit*tmp_counts)/total(tmp_counts)
   height_guess = max(tmp_counts)
   sig_guess    = 1.5
   estimates    = [height_guess,cent_guess,sig_guess]
   gfit         = gaussfit(slice_slit,slice_counts,terms,nterms=3,$
                           est=estimates,sigma=terms_sig)
   ; Reject if getting near the edge
   if (terms(1) le min(slice_slit)+3 or $
       terms(1) ge max(slice_slit)-3) then continue
   slittr(i)    = terms(1)
   fit_fwhm(i)  = 2.*sqrt(2*alog(2)) * terms(2)
   fit_nsig(i)  = terms(0) / terms_sig(0)
endfor

; Evaluate fits - Note: Restriction on fit_fwhm may omit
; good trace points from very narrow traces.
accept_ls   = where(fit_fwhm gt 1.5 and fit_fwhm lt 10 and $
                    fit_nsig ge 5 and $
                    abs(slittr-position) le 15,n_accept)
if (n_accept ne 0) then begin
   trace_merit = replicate(1,n_accept)
   pixtr       = pixtr(accept_ls)
   slittr      = slittr(accept_ls)
endif else begin
   ; To be compatible with following lines
   trace_merit  = 0
endelse

;;;;; Old and problematic method
;; Limit region of the array over which to trace the object
;sz = size(array)
;if keyword_set(bound) then array_bound = [min(bound),max(bound)] else begin
;   if (direction eq 'x') then array_bound = [0,sz(2)] else $
;      array_bound = [0,sz(1)]
;endelse
;if keyword_set(range) then $
;   chip_range = [min(range),max(range)] else begin
;   if (direction eq 'x') then chip_range = [0,sz(1)] else $
;      chip_range = [0,sz(2)]
;endelse
;
;; Choose wavelength regions over which to trace object (optional)
;if keyword_set(spectrum) then begin
;   ; Locate regions of the spectrum with significant flux
;   smspec            = median(spectrum,5)
;   specerr           = sqrt(specvar)
;   w_lo              = -1
;   w_hi              = -1
;   looking_for_start = 1
;   looking_for_end   = 0
;   ;if keyword_set(waverange) then begin
;   ;   i_min = min(where(speclam ge min(waverange)))
;   ;   i_max = max(where(speclam le max(waverange)))
;   ;endif else begin
;      i_min = 0
;      i_max = n_elements(speclam)-1
;   ;endelse
;   for i=i_min,i_max do begin
;      if looking_for_start then begin    ; Locate starting point
;         if (smspec(i) ge 3*specerr(i)) then begin
;            w_lo = [w_lo,speclam(i)]
;            looking_for_start = 0
;            looking_for_end   = 1
;         endif
;      endif else begin                ; Locate ending point
;         if (smspec(i) lt 3*specerr(i)) then begin
;            w_hi = [w_hi,speclam(i-1)]
;            looking_for_start = 1
;            looking_for_end   = 0
;         endif
;      endelse
;   endfor
;   if looking_for_end then w_hi = [w_hi,speclam(i_max)]
;   if (total(w_lo) eq -1) then begin
;      print,'LR_TRACEOBJECT: Error - No regions of significant flux.'
;      return,-1
;   endif else begin
;      w_lo    = w_lo[1:n_elements(w_lo)-1]
;      w_hi    = w_hi[1:n_elements(w_hi)-1]
;      n_chunk = n_elements(w_lo)
;   endelse
;endif else begin
;   ; Treat entire spectrum as a single chunk
;   w_lo = [-1e6]
;   w_hi = [1e6]
;   n_chunk = 1
;endelse
;
;; Edges of the array
;if (direction eq 'x') then begin
;   x_limit = chip_range
;   y_limit = array_bound
;endif else begin
;   x_limit = array_bound
;   y_limit = chip_range
;endelse
;
;; Offset in slitgrid from object position
;slitdiff = abs(slitgrid-position)
;      
;; Trace object over each wavelength chunk
;first_chunk = 1
;; Default is no good trace
;trace_merit = 0
;for i=0,n_chunk-1 do begin
;   ; Choose a starting point for trace along expected object 
;   ; position
;   position_tol = 0.5
;   wave_tol     = 0.5*(w_hi(i) - w_lo(i))
;   mean_wave    = 0.5*(w_hi(i) + w_lo(i))
;   slit_ls      = where(slitdiff le position_tol,n_slit)
;   if (n_slit gt 0) then begin
;      wavediff = abs(wavegrid(slit_ls)-mean_wave)
;      wave_ls  = where(wavediff le wave_tol,n_wave)
;      if (n_wave gt 0) then begin
;         center_ls = slit_ls(wave_ls)
;         n_center  = n_wave
;      endif else begin
;         center_ls = -1
;         n_center  = 0
;      endelse
;   endif else begin
;      center_ls = -1
;      n_center  = 0
;   endelse
;   if (n_center eq 0) then begin
;      ;print,'LR_TRACEOBJECT: No pixels within specified range.'
;   endif else begin
;      ; Stay away from the edges of the array
;      indices = array_indices(array,center_ls)
;      x_ls    = indices(0,*)
;      y_ls    = indices(1,*)
;      buffer  = 8
;      good_ls = where(x_ls ge min(x_limit)+buffer and $
;                      x_ls le max(x_limit)-buffer and $
;                      y_ls ge min(y_limit)+buffer and $
;                      y_ls le max(y_limit)-buffer,n_good)
;      if (n_good eq 0) then begin
;         ;print,'LR_TRACEOBJECT: Object too close to specified edge.'
;      endif else begin
;         n_center  = n_good
;         center_ls = center_ls(good_ls)
;         x_ls      = x_ls(good_ls)
;         y_ls      = y_ls(good_ls)
;         peakvals  = array(center_ls)
;         ; Eliminate the brightest 20% of the selected pixels to avoid
;         ; cosmic rays and then select the remaining brightest pixel
;         ; as a starting trace point.
;         cutoff      = cutvalue(peakvals,80,/incl)
;         peakvalue   = max(peakvals(where(peakvals le cutoff)))
;         peak_pix    = where(peakvals eq peakvalue)
;         peak_pix    = peak_pix(0)
;         xstart      = x_ls(peak_pix)
;         ystart      = y_ls(peak_pix)
;         ; Trace object with respect to slit position.  Assume that 
;         ; the objetc trace will deviate by no more than 10 pixels along
;         ; the slit.
;         if keyword_set(nside) then n_side=floor(nside) else n_side=20
;         sz = size(array)
;         if (direction eq 'x') then begin
;            traceline,array,xstart,ystart,xtr,ytr,temp_slittr,/xdir,/bygauss,$
;               nside=n_side,coord=slitgrid,range=[min(x_ls),max(x_ls)],/noplot,$
;               merit=temp_trace_merit,median=5,cbound=position+[-15,15],$
;               bound=array_bound
;            temp_pixtr = xtr
;            slit_corr_sz = sz(1)
;         endif
;         if (direction eq 'y') then begin
;            traceline,array,xstart,ystart,xtr,ytr,temp_slittr,/ydir,/bygauss,$
;               nside=n_side,coord=slitgrid,range=[min(y_ls),max(y_ls)],/noplot,$
;               merit=temp_trace_merit,median=5,cbound=position+[-15,15],$
;               bound=array_bound
;            temp_pixtr = ytr
;            slit_corr_sz = sz(2)
;         endif
;         if first_chunk then begin
;            pixtr       = temp_pixtr
;            slittr      = temp_slittr
;            trace_merit = temp_trace_merit
;            first_chunk = 0
;         endif else begin
;            pixtr       = [pixtr,temp_pixtr]
;            slittr      = [slittr,temp_slittr]      
;            trace_merit = [trace_merit,temp_trace_merit]
;         endelse
;      endelse
;   endelse
;endfor

; Check that at least some ok points were found
ok_merit_ls = where(trace_merit ne 0,n_ok_merit)
if (n_ok_merit lt 10) then begin
   if (n_ok_merit eq 0) then $
      print,'LR_TRACEOBJECT: Error - No good trace points found.' $
      else print,'LR_TRACEOBJECT: Error - Not enough good trace points found.'
   return,-1
endif
; Restrict range over which to fit trace?
if keyword_set(tracerange) then begin
   if (n_elements(tracerange) ne 2) then begin
      print,'LR_TRACEOBJECT: Error - Keyword tracerange must be a '+$
            '2-element array.'
      return,-1
   endif
   trace_lo = min(tracerange)
   trace_hi = max(tracerange)
endif else begin
   trace_lo = -1e9
   trace_hi =  1e9
endelse
; Do a rough fit first.  Include only those trace points 
; within 5 pixels of the median (unless the object is
; bright, in which cases all trace points are 
; likely to be good).
if (obj_type eq 'bright') then begin
   tmp_ls  = where(pixtr ge trace_lo and pixtr le trace_hi and $
                   trace_merit ne 0,n_tmp)
endif else begin
   tmp_ls  = where(abs(slittr-position) le 5 and $
                   pixtr ge trace_lo and pixtr le trace_hi and $
                   trace_merit ne 0,n_tmp)
endelse

if (n_tmp eq 0) then begin
   print,'LR_TRACEOBJECT: Error - No trace points found wthin range.'
   return,-1
endif

;; This is now redundant.
;; If an object spectrum has been extracted, fit trace
;; at wavelengths where the spectrum shows significant counts.
;if keyword_set(spectrum) then begin
;   ; Filter spectrum first
;   smspec = median(spectrum,5)
;   ; Get wavelengths of trace points
;   wavetr = 0*xtr
;   for t=0,n_elements(xtr)-1 do wavetr(t) = wavegrid(xtr(t),ytr(t))
;   ; Interpolate the extracted spectrum over trace wavelengths
;   spectr    = interpol(smspec,speclam,wavetr)
;   specvartr = interpol(specvar,speclam,wavetr)
;   ; Find 5-sigma points, if any
;   signal = where(spectr ge 5*sqrt(specvartr),n_signal)
;   if (n_signal ne 0) then tmp_ls = signal else begin
;      print,'No significant counts found.'
;      return,-1
;   endelse
;endif

; Set order of the polynomial fit
merit_ls = where(trace_merit ne 0,n_merit)
pixrange = max(pixtr(merit_ls)) - min(pixtr(merit_ls))
case obj_type of
   'faint'   : pfit_order = 1
   'bright'  : pfit_order = 4 < ceil(pixrange/1000.)
   else      : pfit_order = 4 < ceil(pixrange/1000.) ;2 < ceil(pixrange/1000.)
endcase
pfit = poly_fit(pixtr(tmp_ls),slittr(tmp_ls),pfit_order,$
                yfit=tmp_slitfit,yerror=sliterr) 
; Reject outliers and refit until fit converges
old_ls        = 0
slit_done     = 0
max_slit_iter = 20
slit_iter     = 1
while not(slit_done) do begin
   slitfit = 0.
   for h=0,pfit_order do slitfit = slitfit + (pixtr^h)*pfit(h)
   ; Try not to exlude too many points
   sliterr = sliterr > 0.1
   new_ls = where(abs(slittr-slitfit) le 2*sliterr and $
                  pixtr ge trace_lo and pixtr le trace_hi and $
                  trace_merit ne 0,complement=omit_ls,ncomplement=n_omit)
   if (n_elements(new_ls) eq n_elements(old_ls) and $
       max(abs(new_ls-old_ls)) eq 0) then begin
      slit_done = 1
   endif else begin
      pfit = poly_fit(pixtr(new_ls),slittr(new_ls),pfit_order,$
                      yfit=tmp_slitfit,yerror=sliterr)
      old_ls = new_ls
   endelse
   slit_iter = slit_iter + 1
   if (slit_iter ge max_slit_iter) then slit_done=1
endwhile
; Compute and apply correction (offset along the slit as a function of
; pixel index in the dispersion direction) to give an object-centric
; slit position.
slit_correction  = fltarr(slit_corr_sz)
dispersion_index = findgen(slit_corr_sz)
for jj=0,pfit_order do slit_correction = slit_correction + $
                                     (dispersion_index^jj)*pfit(jj)
; Define correction with respect to one edge of the array
slit_correction = slit_correction - slit_correction(0)
obj_slitgrid    = slitgrid
if (direction eq 'x') then $
   for m=0,slit_corr_sz-1 do obj_slitgrid(m,*) = slitgrid(m,*) - $
                                                 slit_correction(m)
if (direction eq 'y') then $
   for m=0,slit_corr_sz-1 do obj_slitgrid(*,m) = slitgrid(*,m) - $
                                                 slit_correction(m)
; Display result of trace
if not(keyword_set(noshow)) then begin
   chan,12
   loadct,39,/silent
   plot,pixtr,slittr,psym=3,/xstyle,/ystyle,xtitle='pixel',$
      ytitle='slit position',yrange=[median(slittr)-3*stddev(slittr),$
      median(slittr)+3*stddev(slittr)]
   plotsym,0
   if (n_omit ne 0) then $
      oplot,pixtr(omit_ls),slittr(omit_ls),psym=8,symsize=0.5,color=150
   oplot,pixtr,slitfit,linestyle=1,color=250
endif

; Return polynomial fit coefficients
fitterms = pfit
print,'Trace fit coeffients: ',strtrim(fitterms,2)+' '

return,obj_slitgrid
end

;-----------------------------------------------------------------------

Function lr_probability_fn,slit,terms

; Calculate the expected flux probability function given the slit
; positions and the object profile parameters.

; Gaussian parameters
t0 = terms(0) ; Height
t1 = terms(1) ; Position
t2 = terms(2) ; Width (sigma)

const = 1./(t2*sqrt(2*!pi))     ; Normalized

prob = const * exp(-1*((slit-t1)^2)/(2*t2^2))

; Check - Does the object profile go to zero?
nearzero_ls = where(prob le (1e-6)*max(prob),n_nearzero)
if (n_nearzero eq 0) then begin
   print,'WARNING - Object profile extends over entire slit region.'
endif

return,prob
end

;-----------------------------------------------------------------------

Function lr_gaussprofile,slit,counts,wave,invvar,terms,waverange=waverange,$
                         bestprofile=bestprofile,spectrum=spectrum,$
                         specwave=specwave,specerr=specerr,$
                         dispersion=dispersion,estimates=estimates,$
                         scaledcounts=scaledcounts,fitlist=fitlist,$
                         result=result,noshow=noshow

if keyword_set(bestprofile) then insert=' best ' else insert=' '

print,'Getting'+insert+'Gaussian object profile...'

; Perform a Gaussian fit to the object profile along the slit.
; Optionally use the extracted spectrum to normalize the profile at
; all wavelengths and make a fit over the entire wavelength range.
; Otherwise, fit over a specified wavelength range.

if keyword_set(estimates) then begin
   est_height = estimates(0)
   est_terms  = estimates
   n_terms    = n_elements(estimates)
endif else begin
   est_height = 1e9
   est_terms  = [0,0,0]
   n_terms    = 3
endelse

; This is a kludge to reduce faint orders in MIKE.  Using
; only three terms can force the Gaussian fit to be too
; wide.
n_terms = 4
est_terms = [est_terms,0]

; Default keyword outputs
terms        = est_terms
scaledcounts = counts
result       = 0

; Choose pixels for profile fitting
if keyword_set(bestprofile) then begin
   interp_spec  = interpol(spectrum,specwave,wave)
   interp_err   = interpol(specerr,specwave,wave)
   norm_counts  = counts / (interp_spec*dispersion)
   norm_invvar  = invvar * ((interp_spec*dispersion)^2)
   ; Fit only over regions with detected flux and avoid cosmic rays.
   section      = where(interp_spec ge 3*interp_err and $
                        abs(norm_counts) le 2       and $
                        wave ge min(waverange)      and $
                        wave le max(waverange),n_section)
   if (n_section eq 0) then begin
      print,'LR_GAUSSPROFILE: No points to fit.'
      result = 0
      return,-1
   endif
   temp_counts  = norm_counts(section)
   temp_invvar  = norm_invvar(section)
   ; "Normalize" estimate of peak intensity.
   if (keyword_set(estimates) and estimates(2) ne 0) then $
      est_height = 1./sqrt(2*!pi)/estimates(2) $
      else est_height = 0.5
   est_terms(0) = est_height
   ; Return points for plotting
   scaledcounts = norm_counts
endif else begin
   section = where(wave ge min(waverange)      and $
                   wave le max(waverange)      and $
                   abs(counts) le 5*est_height,n_section)
   if (n_section eq 0) then begin
      print,'LR_GAUSSPROFILE: No points to fit.'
      result = 0
      return,-1
   endif
   temp_counts  = counts(section)        
   temp_invvar  = invvar(section)
   ; Return points for plotting
   scaledcounts = counts
endelse
temp_slit = slit(section)

terms = est_terms

; Make sure there are a reasonable number of points to fit.
if (n_section lt 20) then begin
   print,'LR_GAUSSPROFILE: Too few points to fit.'
   result = 0
   return,-1
endif

; Fit a Gaussian.  Reject outliers and refit.
max_gfit_iter = 10
n_gfit_iter   = 0
gfit_done     = 0
old_gfit      = 0
while not(gfit_done) do begin
   if (n_gfit_iter eq 0) then begin
      good_slit   = temp_slit
      good_counts = temp_counts
   endif else begin
      resid       = temp_counts - gfit
      good_ls     = where(abs(resid) le 7*gfiterr,n_good)
      if (n_good lt 10) then begin
         print,'LR_GAUSSPROFILE: Profile fit failed.  Too few points.'
         return,-1
      endif
      good_slit   = temp_slit(good_ls)
      good_counts = temp_counts(good_ls)
   endelse
terms(2) = 1.5
   temp_gfit = gaussfit(good_slit,good_counts,terms,nterms=n_terms,$
                        estimates=terms,yerror=gfiterr,sigma=sig)
   temp_z    = (temp_slit-terms(1))/terms(2)
   gfit      = terms(0) * exp(-(temp_z^2)/2.)
   gfit_diff = gfit - old_gfit
   if (max(abs(gfit_diff)) le 0.01*max(gfit)) then gfit_done = 1
   n_gfit_iter = n_gfit_iter + 1
   if (n_gfit_iter gt max_gfit_iter) then begin
      print,'LR_GAUSSPROFILE: Too many iterations.'
      gfit_done = 1
   endif
   old_gfit = gfit
endwhile

; Force sigma to be positive
terms(2) = abs(terms(2))

; Flag suspicious fits.  Default is an acceptable fit.
result = 1
if (sig(0) ge 0.2*terms(0)     or $  ; Large uncertainties
    sig(2) ge 0.2*terms(2)     or $
    sig(1) ge 0.5              or $
    terms(1) lt min(good_slit) or $  ; Center off edge of slit
    terms(1) gt max(good_slit) or $
    terms(2) lt 1/2.35         or $  ; Unreasonable FWHM
    terms(2) gt 10/2.35)       then result = 0
if (keyword_set(bestprofile)   and $
    (terms(0) lt 0.01          or $  ; Unreasonable amplitude
     terms(0) gt 1))           then result = 0
if (result eq 0) then print,'LR_GAUSSFIT: Supicious profile'

; Display (optional)
if not(keyword_set(noshow)) then begin
   chan,13
   loadct,39,/silent
   gmax = max(gfit)
   plot,temp_slit,temp_counts,psym=3,/xstyle,/ystyle
   oplot,temp_slit,gfit,psym=3,color=250
endif

; Compute profile over the entire input wave range
g_profile = lr_probability_fn(slit,terms)

; Return list of fitted points
fitlist = section

; Return profile
return,g_profile

end

;-----------------------------------------------------------------------

Function lr_splineprofile,slit,counts,wave,invvar,waverange=waverange,$
                          bestprofile=bestprofile,spectrum=spectrum,$
                          specwave=specwave,specerr=specerr,$
                          dispersion=dispersion,bad=bad,$
                          scaledcounts=scaledcounts,fitlist=fitlist,$
                          variable=variable,noshow=noshow

; Fit a b-spline to the object profile.

; NOTE: The spline profile will automatically be normalized if the
;       bestprofile keyword is set and the input "guess" spectrum
;       was made using boxcar extraction.  A normalized profile
;       will give the same count levels as boxcar extraction.

if keyword_set(bestprofile) then insert=' best ' else insert=' '

print,'Getting'+insert+'spline profile...'

; Default keyword outputs
scaledcounts = counts
result       = 0

; Default is no bad pixels
if keyword_set(bad) then badpixels = bad else badpixels = 1 + 0*slit

; Flag pixels with suspicious dispersion
if keyword_set(dispersion) then begin
   if (n_elements(dispersion) gt 1) then begin
      meddisp     = median(dispersion(where(dispersion ne 0)))
      bad_disp_ls = where(abs(dispersion) gt 10*abs(meddisp),n_bad_disp)
      if (n_bad_disp ne 0) then badpixels(bad_disp_ls) = 0
   endif
endif else dispersion = 1

; Choose pixels to fit
if keyword_set(bestprofile) then begin
   interp_spec  = interpol(spectrum,specwave,wave)
   interp_err   = interpol(specerr,specwave,wave)
   norm_counts  = counts / (interp_spec*dispersion)
   norm_invvar  = invvar * ((interp_spec*dispersion)^2)
   ; Fit only over regions with detected flux and avoid cosmic rays.
   section      = where(interp_spec ge 3*interp_err and $
                        abs(norm_counts) le 2       and $
                        badpixels gt 0              and $
                        wave ge min(waverange)      and $
                        wave le max(waverange),n_section)
   if (n_section eq 0) then begin
      print,'LR_SPLINEPROFILE: No points to fit.'
      return,-1
   endif
   ; To speed things up, limit the number of fitted pixels.
   max_n_per_unit_slit = 500.
   delta_slit    = max(slit) - min(slit)
   max_n_section = max_n_per_unit_slit * delta_slit
   if (n_section ge max_n_section) then begin
      temp_ls = round((n_section/max_n_section)*findgen(max_n_section))
      section = section(temp_ls)
   endif
   temp_counts  = norm_counts(section)
   temp_invvar  = norm_invvar(section)
   ; Return points for plotting
   scaledcounts = norm_counts
endif else begin
   section = where(wave ge min(waverange) and $
                   wave le max(waverange) and $
                   badpixels ne 0,n_section)
   if (n_section eq 0) then begin
      print,'LR_SPLINEPROFILE: No points to fit.'
      return,-1
   endif
   temp_counts  = counts(section)
   temp_invvar  = invvar(section)
   ; Return points for plotting
   scaledcounts = counts
endelse
temp_slit = slit(section)
temp_wave = wave(section)

; Place breakpoints at every half integer slit position
bk_space = 0.5

; Compute fit.  Reject outliers and refit.
max_sfit_iter = 3
n_sfit_iter   = 0
sfit_done     = 0
old_profile   = 0
while not(sfit_done) do begin
   if (n_sfit_iter eq 0) then begin
      good_slit   = temp_slit
      good_counts = temp_counts
      good_invvar = temp_invvar
      good_wave   = temp_wave
      ; Step through with a median filter and reject outliers.
      ; Also compute the local std deviation.
      h = histogram(good_slit,binsize=bk_space,location=step_start,$
                    rev=filt_ri)
      n_step = n_elements(h)
      local_sdev = fltarr(n_step)
      local_slit = step_start + 0.5*bk_space
      for i=0,n_step-1 do begin
         if (h(i) le 2) then continue
         temp_ls       = filt_ri(filt_ri(i):filt_ri(i+1)-1)
         temp_sdev     = robust_sigma(good_counts(temp_ls))
         if (temp_sdev eq -1) then $
            temp_sdev  = stddev(good_counts(temp_ls))
         temp_med      = median(good_counts(temp_ls))
         temp_resid    = good_counts(temp_ls) - temp_med
         rej_ls        = where(abs(temp_resid) gt 5*temp_sdev,n_rej)
         local_sdev(i) = temp_sdev
         if (n_rej ne 0) then good_invvar(temp_ls(rej_ls)) = 0
      endfor
      no_sdev_ls = where(local_sdev eq 0,n_no_sdev,complement=sdev_ls)
      if (n_no_sdev ne 0) then $
         local_sdev(no_sdev_ls) = median(local_sdev(sdev_ls))
      local_sdev = median(local_sdev,5)
      temp_sdev  = interpol(local_sdev,local_slit,temp_slit)
   endif else begin
      ; Reject points more aggressively when getting best
      ; profile
      if keyword_set(bestprofile) then nsig_rej = 2 else nsig_rej = 5
      resid     = temp_counts - new_profile
      good_ls   = where(abs(resid) le nsig_rej*temp_sdev,n_good)
      if (n_good eq 0) then begin
         print,'LR_SPLINEPROFILE: Profile fit failed.'
         return,-1
      endif
      good_slit   = temp_slit(good_ls)
      good_counts = temp_counts(good_ls)
      good_invvar = temp_invvar(good_ls)
      good_wave   = temp_wave(good_ls)
   endelse
   if keyword_set(bestprofile) then begin
      good_x2 = good_wave
      temp_x2 = temp_wave
      ; Allow profile to vary in the dispersion direction?
      if keyword_set(variable) then n_poly = 3 else n_poly = 1
   endif else begin
      good_x2 = 0
      temp_x2 = 0
      n_poly  = 0
   endelse
   profile_sset = bspline_iterfit(good_slit,good_counts,invvar=good_invvar,$
                                  bkspace=bk_space,yfit=temp_profile,$
                                  x2=good_x2,npoly=n_poly)
   new_profile  = bspline_valu(temp_slit,profile_sset,x2=temp_x2)
   profile_diff = new_profile - old_profile
   if (max(abs(profile_diff)) le 0.01*max(new_profile)) then sfit_done = 1
   n_sfit_iter  = n_sfit_iter + 1
   if (n_sfit_iter gt max_sfit_iter) then begin
      print,'LR_SPLINEPROFILE: Too many iterations.'
      sfit_done = 1
   endif
   old_profile = new_profile
endwhile

; Compute the profile over the entire input wave range
if keyword_set(bestprofile) then x_2 = wave else x_2 = 0
s_profile = bspline_valu(slit,profile_sset,x2=x_2)

; Display (optional)
if not(keyword_set(noshow)) then begin
   chan,13
   loadct,39,/silent
   smax = max(s_profile)
   ymin = -0.5*smax > min(temp_counts)
   ymax = 1.5*smax < max(temp_counts)
   plot,temp_slit,temp_counts,psym=3,/xstyle,yrange=[ymin,ymax]
   oplot,good_slit,temp_profile,psym=3,color=200
endif

; Return list of fitted points
fitlist = section

return,s_profile

end

;-----------------------------------------------------------------------
;-----------------------------------------------------------------------

Pro longslit_reduce,objframe,wavearr,slitarr,nosky,noskyvar,wavefit,$
                    skymodel,nirspec=nirspec,n_echelle=n_echelle,$
                    luci=luci,lrisr=lrisr,lrisb=lrisb,$
                    fors2_1=fors2_1,fors2_2=fors2_2,gmos=gmos,hires=hires,$
                    ldss3=ldss3,dbsp=dbsp,mike=mike,mage=mage,$
                    ccd=ccd,hdr=hdr,skyref=skyref,$
                    bias=bias,flatfield=flatfield,flatvar=flatvar,$
                    darkframe=darkframe,nodarksub=nodarksub,nobias=nobias,$
                    gn=gn,rdnoise=rdnoise,bkspace=bkspace,$
                    min_bkspace=min_bkspace,set_bkpts=set_bkpts,$
                    illumpoly=illumpoly,slitpos=slitpos,fwhm=fwhm,$
                    subwidth=subwidth,subrange=subrange,noobject=noobject,$
                    maskall=maskall,objcoadds=objcoadds,darkcoadds=darkcoadds,$
                    readnoise=readnoise,texp=texp,fixgain=fixgain,$
                    bright=bright,faint=faint,boxcar=boxcar,$
                    sprofile=sprofile,findrange=findrange,trace_obj=trace_obj,$
                    tracerange=tracerange,stdtrace=stdtrace,stdpos=stdpos,$
                    stdoffset=stdoffset,autoid=autoid,wavemodel=wavemodel,$
                    writewave=writewave,waveref=waveref,waveshift=waveshift,$
                    meanshift=meanshift,linelist=linelist,vacuum=vacuum,$
                    waveinfo=waveinfo,nohelio=nohelio,bin=bin,loglin=loglin,$
                    outfile=outfile,noshow=noshow,maxiter=maxiter,force=force,$
                    noskysub=noskysub

;+----------------------------------------------------------------------
;
; LONGSLIT_REDUCE     5/2004
;
; General purpose longslit spectrum reduction package.
;
; Instrument-specific parameters are set in the file 'longslit_reduce.inc'.
;
; The program performs optimal sky subtraction using two-dimensional
; maps of slit position and wavelength.  Object spectra are optionally
; extracted using an optimal subtraction routine.
;
; Sky background modeling is done using a 2D b-spline fit, where the
; slit illumination function is fit with a low-order polynomial and a
; b-spline is fit in the wavelength direction.  Fitting is done using
; BSPLINE_ITERFIT, written by D. Schlegel & S. Burles, which is part
; of the SDSS spectral reduction pipeline.
;
; INPUTS
;        objframe  - 2D array or FITS filename, raw object exposure.
;        wavearr   - 2D array or FITS filename, same size as 'objframe',
;                    containing proxies to the wavelengths of pixels in the 
;                    input 'objframe'. 
;        slitarr   - 2D array or FITS filename, same size as 'objframe',
;                    containing proxies to the slit position of pixels in 
;                    the input 'objframe'. 
;
; KEYWORDS   
;        nirspec   - If set, data are assumed to be NIRSPEC (default)
;        n_echelle - If set, data are assumed to be NIRSPEC Echelle
;        lrisr     - If set, data are assumed to be LRIS-Red
;        lrisb     - If set, data are assumed to be LRIS-Blue
;        fors2_1   - If set, data are assumed to be FORS2, chip 1
;        fors2_2   - If set, data are assumed to be FORS2, chip 2
;        gmos      - If set, data are assumed to be GMOS
;        hires     - If set, data are assumed to be HIRES (upgraded detector)
;        ldss3     - If set, data are assumed to be LDSS3
;        dbsp      - If set, data are assumed to be DBSP
;        mike      - If set, data are assumed to be MIKE
;        mage      - If set, data are assumed to be MagE
;        luci      - If set, data are assumed to be LUCIFER
;        ccd       - Scalar, number of the CCD (1,2,3, etc.) on which to
;                    perform sky subtraction.  Note:  The first chip should be
;                    numbered 1, not 0.  The default is to model the
;                    sky across all ccds (those included in the input
;                    objframe).  However, this should be avoided when
;                    there is a significant difference in QE between
;                    chips (e.g., LRIS-B chips at lambda < 4000
;                    angstroms).  If an object spectrum is to be
;                    extracted, the program will automatically select
;                    which ccd.
;        hdr       - String, FITS header for input objframe.  Useful
;                    if objframe is an array (rather than a FITS
;                    filename) and heliocentric velocity correction is
;                    desired.  If objframe is a FITS file, hobj will
;                    be overwritten with the header from the file. On
;                    output, hobj will contain information provided by
;                    the program.
;        skyref    - 2D array of FITS filename, raw 2D longslit object
;                    exposure to be used as reference for primary sky
;                    subtraction.  Should be an exposure taken near in
;                    time to objframe, and in same setup.  Asuumed to
;                    have same orientation as objframe.
;        bias      - 2D array or FITS filename, bias.  For LRIS-R,
;                    LRIS-B, LDSS3, and FORS2, the default is to
;                    subtract the bias using the pre/overscan region.
;                    Otherwise, the default is no bias subtraction.
;        flatfield - 2D array or FITS filename, flat field.  If the
;                    flat field has been normalized it will be used to
;                    create a bad pixel map.  Default is no
;                    flat-fielding.
;        flatvar   - 2D array or FITS filename, expected
;                    variance in the flat field, in ADU.  Default is
;                    to assume that the flat contain zero error.
;        darkframe - 2D array or FITS filename, dark frame to be used
;                    to create bad pixel map ONLY.  For NIRSPEC data,
;                    subtracting 2D darkframes has been found to give
;                    unsatisfactory results in the sky subtraction.
;                    Therefore, dark counts to be subtracted are
;                    calculated from the 'texp' keyword and the
;                    'darkcurrent' parameter set in the included file.
;        nodarksub - If set, program will not subtract dark counts
;                    prior to doing sky subtraction.  However, the
;                    expected nonise in the dark counts will still
;                    be calculated based on the 'texp' keyword and the
;                    'darkcurrent' parameter set in the included file.  
;                    This is useful if the object frame has been
;                    dark-subtracted prior to calling this program.
;        nobias    - If set, do not perform bias subtraction (optical
;                    data only).
;        gn        - Scalar, or array with one value per amplifier,
;                    gain (e-/ADU).  Default is the 'gain' value(s)
;                    set in the included file.
;        rdnoise   - Scalar, or array with value per amplifier,
;                    readnoise (e-).  Default is the 'erdnoise'
;                    value(s) set in the included file.
;        bkspace   - Scalar, initial spacing of b-spline break points for
;                    fitting sky, in units of the input 'wavearr'.
;                    Default is the 'bk_space' value set in the
;                    included file.  A typically value is ~ 1 pixel.
;        min_bkspace - Scalar, minimum spacing of b-spline break
;                    points for fitting sky, in units of the input
;                    'wavearr'.  Default is the 'min_bk_space' value
;                    set in the included file.  If min_bkspace <
;                    bkspace, then program will attempt to add
;                    breakpoints where it finds a bad fit, unless
;                    the set_bkpts keyword is set (recommended).
;        set_bkpts - Recommended.  If set, program will attempt to
;                    position b-spline break points such that the
;                    minimum break point spacing ('min_bk_space' in
;                    the included file, or min_bkspsace keyword) is
;                    used around skylines, and the max break point
;                    space ('bk_space' in the included file, or
;                    bkspace keywoprd) is used elsewhere.  Default is
;                    to start by using the maximum break point spacing
;                    everywhere.  In the default mode, if min_bk_space
;                    is less than bk_space, then the program will
;                    attempt to add break points where it thinks the
;                    fit to the sky is poor.  However, this generally
;                    produces mixed results and is slower.  Using the
;                    /set_bkpts keyword is especially recommended in
;                    cases where the skylines are well resolved from
;                    the general sky background.
;        illumpoly - Scalar, order + 1 of the polynomial to be used to
;                    fit the slit illumination function.  Default is 2
;                    (linear).  For no slit illumination fit, set
;                    illumpoly = -1.  The numbering scheme follows the
;                    convention of the bspline fitting package (1 is
;                    constant, 2 linear, 3 quad, etc.).
;        slitpos   - Scalar, approximate slit position of object, in
;                    units of the input 'slitarr'.  Useful for
;                    extracting multiple objects.  If the 'maskall'
;                    keyword is set, then the spectrum of the object
;                    at this position will be extracted.  If zero on
;                    input, will contain the slit position at the center
;                    of the object profile on output.
;        fwhm      - Scalar, estimated FWHM of a Gaussian object
;                    profile in units of the input 'slitarr'.  Unless
;                    the force keyword is set, the program will
;                    automatically attemp to fit a Gaussian profile to
;                    the object in order to mask the object out, and
;                    possibly to use for optimal extraction.  If the
;                    fit fails, then this value (if non-zero) will be
;                    used instead.  On output, the keyword will
;                    contain the value used.
;        subwidth  - Scalar, range along the slit over which to fit
;                    the sky, in units of the input 'slitarr'.  The
;                    sky will be fitted to pixels within +/-
;                    subwidth/2 of the object position.  Default is
;                    the value of 'sub_width' in the included file.
;        subrange  - Two-element array, minimum and maximum slit
;                    values over which to fit the sky, in units of the
;                    input 'slitarr'.  
;                    *If both the 'subwidth' and 'subrange' keywords
;                    are set, then the sky will be fit to pixels that
;                    fall within both ranges.
;        noobject  - If set, program will not search for an object to
;                    mask (unless the maskall keyword is set) or
;                    attempt to extract an object spectrum.  Sky
;                    subtraction will be performed across the entire
;                    slit unless the 'subrange' keyword is set.
;        maskall   - If set, program will attempt to mask all objects
;                    found along the slit.  The brightest object will
;                    be extracted unless the 'slitpos' keyword or
;                    /noobject keyword is set.
;        objcoadds - Scalar, number of coadds for the input 'objframe'
;                    (used to determine expected read noise).  Default
;                    is 1.
;        darkcoadds- Scalar, number of coadds for dark exposure (used
;                    to determine expected read noise).  Default is to
;                    set this value equal to the 'objcoadds' keyword
;                    value.  NOTE: This keyword currently has no
;                    effect since the dark count level is calculated
;                    from the exposure time and the expected amount of
;                    dark current.  Therefore, additional read noise
;                    from a dark exposure is not included.
;        readnoise - Scalar, read noise per exposure, or per coadd, if
;                    applicable, in rms e/pixel.  Default is the value
;                    of 'erdnoise' set in the included file.
;        texp      - Scalar, exposure time in seconds int time.  If
;                    multiple coadds were taken then the total time
;                    for all coadds should be used.  The exposure time
;                    is currently used only to estimate the expected
;                    number of dark counts.  Default is 0 sec.
;        fixgain   - If set, program will adjust counts in the object
;                    frame for gain differences prior to
;                    flat-fielding.  Counts in the region covered by
;                    amplifier i will be divided by
;                    gain(i)/mean(gain).  This keyword should be set
;                    only if the flatfield has had the gain
;                    differences taken out.
;        bright    - If set, program will use bright object settings:
;                    Mask out an extra-wide region aound the object
;                    when fitting the sky background, use a 3rd order
;                    fit of the object trace with respect to slit
;                    position (default is 2nd order if tracing is
;                    done), and, if performing optimal extraction, use
;                    a b-spline fit to the object profile.  Profile
;                    may also vary in the dispersion direction, depending
;                    on what has been commented out below.
;        faint     - If set, use faint object settings: Trace the
;                    object with respect to slit position to 1st order
;                    only (default is one order per 1000 pixels, up to
;                    4th order, if tracing is done).
;        boxcar    - If set, spectral extraction is done using boxcar
;                    extraction.  Default is optimal extraction.
;        sprofile  - If set, optimal extaction is performed using a
;                    b-spline fit to the object profile.  It is
;                    recommended to set bestprofile = 1 in the
;                    longslit_redeuce.inc when using this option.
;                    Default is to use a Gaussian profile, unless the
;                    'bright' keyword is set.
;        findrange - 2-element array, minimum and maximum x-value (if
;                    the dispersion is in the x-direction) or y-value
;                    (if the dispersion is in the y-direction) over
;                    which to locate objects.  If set, the wavelengths
;                    bound by this range overide the 'wavemin' and
;                    'wavemax' values in the included file.
;        trace_obj - If set, program will attempt to compute
;                    'object-centric' slit coordinates by tracing the
;                    object.  If trace_obj is < 0, then object will
;                    not be traced.  Default is the behavior set by
;                    the traceobj parameter in the included file.
;        tracerange- 2-element array, minimum and maximmum pixel
;                    indices in the dispersion direction over which to fit
;                    object trace with respect to slit position.
;        stdtrace  - 2D array or FITS filename, exposure containing a
;                    bright source to be used to determine the object
;                    trace.  If set, the science object will not be
;                    traced separately, but will be assumed to have a
;                    parallel trace to the standard object.  The
;                    standard object will be traced as a bright object
;                    (i.e., to higher order than a typical object).
;        stdpos    - Scalar, approximate slit position of the standard
;                    trace object, in units of the input 'slitarr'.
;        stdoffset - Scalar, estimate of the offset from the standard
;                    trace to the science object trace, in units of
;                    the input 'slitarr'.  Unless the force keyword is
;                    set, the program will automatically attempt to
;                    locate the science object.  If that fails, the
;                    input offset value (if non-zero) will be used
;                    instead.  On output, the keyword will contain the
;                    value of the offset used.
;        autoid    - If set, wavelength fitting program will
;                    automatically identify skylines.  This keyword
;                    should be set only when the input wavearr
;                    contains values that are close to the expected
;                    final calibrated values.  The default is to
;                    identify skylines by hand, unless the 'wavemodel'
;                    keyword is set.
;        wavemodel - String, filname of tabulated wavelength fit to be
;                    used as an initial guess for automatic line
;                    identification.  Format is a three column ASCII
;                    file, where the columns are: uncalibrated
;                    wavelengths (as in the input 'wavearr'),
;                    intensity, and calibrated wavelengths.  If this
;                    keyword is set, the program 'wave_calib' will
;                    attempt to automatically identify each line in
;                    the linelist file (see 'linelist' keyword).
;        writewave - If set, program will output an ASCII file giving
;                    a tabulated wavelength solution as a funciton of
;                    the uncalibrated input wavelength values.  The
;                    format is the same as that described above for
;                    the 'wavemodel' keyword.  Sampling rate is at
;                    least 2 steps per pixel.  This file can be used
;                    in future calls of the 'wavemodel' keyword.
;                    Filename is '<outfile>_wav.dat'.
;        waveref   - 2D array or FITS filename, calibrated wavelength
;                    values to be used rather than calibrating
;                    wavelengths to identified skylines.  Useful for a
;                    standard star taken with an exposure too short to
;                    accurately measure sky lines.  If a heliocentric
;                    velocity value is present in the FITS header,
;                    that velocity shift will be undone before
;                    applying a heliocentric correction for the new
;                    array, unless the /nohelio keyword is set.  If
;                    the /nohelio keyword is set, no heliocentric
;                    velocity correction will be applied or undone.
;        waveshift - If set, program will use identified skylines to
;                    compute a linear shift in wavelength only.  This
;                    assumes that the input 'wavearr' is in physical
;                    units (e.g., angstroms or microns).  This option
;                    should be used for data where wavelength
;                    calibration has been performed using an arc lamp
;                    and a small offset may be required for each
;                    science exposure.  The fitted wavelengths will be
;                    a linear function of the input wavelengths if
;                    multiple lines are found, or a constant offset
;                    from the input wavelengths if only one line is
;                    found.  Default is to compute a third-order
;                    polynomial fit of the input wavelengths to the
;                    wavelengths of the identified skylines.
;        meanshift - Scalar, mean wavelength shift.  If skylines are
;                    identified, keyword returns the mean change in
;                    wavelength before applying a heliocentric
;                    velocity correction.  Otherwise, this shift is
;                    applied to the input wavelengths prior to
;                    correcting for heliocentric velocity.  Useful if
;                    the rough wavelength shift is known but cannot be
;                    measured from skylines.
;        linelist  - String, name of file containing line wavelengths.
;                    Wavelengths should be in the first column.  Other
;                    columns are ignored.  Lines commented out with a
;                    '#' are ignored.  
;                    NOTE: No air to vacuum conversion is performed
;                    unless the 'vacuum' keyword is set.  Default is a
;                    list of near-IR OH sky lines in vacuum
;                    wavelengths.
;        vacuum    - If set, wavelengths are converted from air to
;                    vacuum following skyline calibration.  To use
;                    this option, ANY and ALL wavelength inputs
;                    ('waveref', 'wavemodel', 'linelist', etc.) should
;                    be in air wavelengths.  The final 2D wavelength
;                    array and 1D spec wavelengths will be written in
;                    vacuum wavelengths.  
;                    NOTE: Unless either the nehlio keyword is set, a
;                    helio-centric velocity correction is always
;                    applied when the object input is a FITS file with
;                    the necessary header information, whether the
;                    output is in air or vacuum wavelengths.
;        nohelio   - If set, no heliocentric velocity correction is
;                    applied to the calibrated wavelengths.  Default
;                    is to apply a coorection if the necessary
;                    information is present in the input FITS headers.
;                    If the 'reference' keyword is set, the
;                    heliocentric velocity correctin applied to the
;                    reference wavelength array will not be undone,
;                    nor will a new correction for the current array
;                    be applied.
;        bin       - Scalar, wavelength binning for the final
;                    extracted spectrum, in calibrated wavelength
;                    units, or, if the /loglin keyword is set, in
;                    km/s.  Default is the 'value of the 'wavebin'
;                    parameter set in the included file.
;        loglin    - If set, extraced spectrum uses bins of constant
;                    velocity width, and the value of bin is assumed
;                    to be in km/s.  Default is constant wavelength units.
;        outfile   - String, root name for output files.  Default is
;                    'longslit'.
;        noshow    - If set, program will not display successive
;                    subtraction/extraction iterations. Default is to
;                    display results to the screen.
;        maxiter   - Scalar, maximum number of sky-fitting iterations.
;                    Default is 20.
;        force     - Force the program to use the input fwhm and
;                    slitpos or stdoffset values (if non-zero), rather
;                    than attempting to fit a Gaussian to the object
;                    profile.
;        noskysub  - Do not subtract sky.  Useful for short exposures
;                    of bright objects that fill the entire slit.
;                    Program will still process the raw image, trace
;                    the object, extract a spectrum etc., as desired.
;                    NOTE: Since skylines will not be available for
;                    fitting, the waveref keyword must also be set.

;
; OUTPUT FILES
;        outfile_sub.fits   - 2D processed sky-subtraced frame
;        outfile_var.fits   - 2D formal variance in the sky-subtracted
;                             frame
;        outfile_sky.fits   - 2D b-spline sky model
;        outfile_wav.fits   - 2D wavelength solution
;        outfile_wav.dat    - (Optional) Tabulated wavelength solution,
;                             written if the 'writewave' keyword
;                             is set.
;        outfile_spec.fits  - 1D extracted spectrum (header contains
;                             wavelength information)
;        outfile_err.fits   - 1D formal 1-sigma uncertainty in the
;                             extracted spectrum (header contains
;                             wavelength information)
;        outfile_spec.ps    - Postscript plot of the extracted spectrum
;                             and error
;        outfile_profile.ps - Postscript plot of a cut along the
;                             slit and the fitted object profile.
;        outfile_trace.ps   - Postscript plot of the fit to the object
;                             trace 
;
; OPTIONAL OUTPUTS
;        nosky     - 2D array, processed, sky-subtracted frame
;        noskyvar  - 2D array, formal variance in the 2D processed, 
;                    sky-subtracetd frame
;        skymodel  - 2D array, b-spline fit to sky background
;        wavefit   - 2D array, final wavelength solution
;
; HISTORY
;        Written 5/19/2004 GDB
;        Keywords 'wavemodel', 'writewave', 'noobject', and 'linelist' 
;           added 8/2004 GDB
;        Converted from 'NIRSPEC_REDUCE' 8/2004 GDB
;        Updated 4/2005 GDB
;        SET_BKPTS keyword added 3/9/2007 GDB
;        NOSKYSUB keyword added 5/16/2007 GDB
;----------------------------------------------------------------------

if (n_params() lt 3) then begin
   print,'CALLING SEQUENCE:'
   print,'longslit_reduce,objframe,wavearr,slitarr,nosky,noskyvar,wavefit,'
   print,'          skymodel,/nirspec,/n_echelle,/lrisr,/lrisb,/fors2_1,'
   print,'          /fors2_2,/gmos,/hires,/ldss3,/dbsp,/mike,/mage,'
   print,'          ccd=ccd,hdr=hdr,skyref=skyref,bias=bias,'
   print,'          flatfield=flatfield,flatvar=flatvar,darkframe=darkframe,'
   print,'          /nodarksub,/nobias,gn=gn,rdnoise=rdnoise,bkspace=bkspace,'
   print,'          min_bkspace=min_bkspace,/set_bkpts,illumpoly=illumpoly,'
   print,'          slitpos=slitpos,fwhm=fwhm,'
   print,'          subwidth=subwidth,subrange=subrange,/noobject,/maskall,'
   print,'          objcoadds=objcoadds,darkcoadds=darkcoadds,'
   print,'          readnoise=readnoise,texp=texp,/fixgain,/bright,/faint,'
   print,'          /boxcar,/sprofile,findrange=findrange,/trace_obj,'
   print,'          tracerange=tracerange,stdtrace=stdtrace,stdpos=stdpos,'
   print,'          stdoffset=stdoffset,/autoid,wavemodel=wavemodel,'
   print,'          /writewave,waveref=waveref,waveshift=waveshift,'
   print,'          meanshift=meanshift,linelist=linelist,/vacuum,'
   print,'          waveinfo=waveinfo,/nohelio,bin=bin,/loglin,'
   print,'          outfile=outfile,/noshow,maxiter=maxiter,/force,/noskysub'
   return
endif

;;; System variables
currentmulti = !p.multi
!p.multi     = [0,1,1]

;;; Plot in color
device,decomp=0

;;; Display progress to the screen?
if keyword_set(noshow) then no_show = 1 else no_show = 0

;;; Constants
c = 2.99792458d5     ; Speed of light (km/s)

;;; Default outputs
nosky    = -1
noskyvar = -1
wavefit  = -1
skymodel = -1
speclam  = [0]
spec     = [0]
specvar  = [0]

;;; Return slit position of the extracted object?
if not(keyword_set(slitpos)) then return_slitpos = 1 else return_slitpos = 0

;;; Determine which instrument.  Default is unknown.

instr = 'unknown'
if keyword_set(nirspec)   then instr = 'nirspec'
if keyword_set(n_echelle) then instr = 'n_echelle'
if keyword_set(lrisr)     then instr = 'lrisr'
if keyword_set(lrisb)     then instr = 'lrisb'
if keyword_set(fors2_1)   then instr = 'fors2_1'
if keyword_set(fors2_2)   then instr = 'fors2_2'
if keyword_set(gmos)      then instr = 'gmos'
if keyword_set(hires)     then instr = 'hires'
if keyword_set(ldss3)     then instr = 'ldss3'
if keyword_set(dbsp)      then instr = 'dbsp'
if keyword_set(mike)      then instr = 'mike'
if keyword_set(mage)      then instr = 'mage'
if keyword_set(luci)      then instr = 'lucifer'
;;; Check that if the sky is not being fit, that a wavelength
;;; solution has been included.

if (keyword_set(noskysub) and not(keyword_set(waveref))) then begin
   print,'Must set waveref keyword when using /noskysub option.'
   return
endif

;;;
;;; Include instrument-specific parameters
;;;

; Include file should be aliased to an instrument-specific setup
; file.
@longslit_reduce.inc


;;;
;;; Read in input arrys (if necessary)
;;;

if (size(objframe,/tname) eq 'STRING') then begin
   checkexists = findfile(objframe,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Object file not found.'
      return
   endif else begin
       print,' '
       print,'Reducing '+objframe
       print,' '
       rawobj = readfits(objframe,hobj)
       sxaddpar,hobj,'ORIGEXP',filename(objframe),$
          'Filename of original exposure.'
   endelse
endif else begin
   ; If type is unsigned integer, convert to long
   if (size(objframe,/tname) eq 'UINT') then rawobj = long(objframe) $
      else rawobj = objframe
   ; Header information provided?
   if keyword_set(hdr) then begin
      hobj = hdr
   endif else begin
     ; Create a basic FITS header
      fxhmake,hobj,rawobj,/date
   endelse
endelse

; Check BZERO keyword
if (sxpar(hobj,'BZERO') ne 0) then begin
   print,'Setting BZERO keyword to 0.0'
   sxaddpar,hobj,'BZERO',0.0
endif

;;; Reference sky (optional)
if keyword_set(skyref) then begin
  if (size(skyref,/tname) eq 'STRING') then begin
      checkexists = findfile(skyref,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Sky reference file not found.'
         return
      endif else begin
         skyrefarr = readfits(skyref)
         sxaddpar,hobj,'SKYREF',filename(skyref),'Reference sky exposure'
      endelse
   endif else begin
      ; If type is unsigned integer, convert to long
      if (size(skyref,/tname) eq 'UINT') then skyrefarr = long(skyref) $
         else skyrefarr = skyref
   endelse
   ; Mask for reference objects (default is no masking)
   ref_mask = 1 + 0*skyrefarr
endif else begin
   skyrefarr = 0.*rawobj
   ref_mask  = 1 + 0*rawobj
endelse

;;; Wavelength-proxy array
if (size(wavearr,/tname) eq 'STRING') then begin
   checkexists = findfile(wavearr,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Wavelength file not found.'
      return
   endif else begin
      wavegrid = readfits(wavearr)
      sxaddpar,hobj,'WAVEGRID',filename(wavearr),'Wavelength-proxy grid'
   endelse
endif else wavegrid = wavearr

;;; Get size of arrays
sz  = size(wavegrid)
xsz = sz(1)
ysz = sz(2)

;;; Slit-position-proxy array
if (size(slitarr,/tname) eq 'STRING') then begin
   checkexists = findfile(slitarr,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Slit position file not found.'
      return
   endif else begin
      slitgrid = readfits(slitarr)
      sxaddpar,hobj,'SLITGRID',filename(slitarr),'Slit-position-proxy grid'
   endelse
endif else slitgrid = slitarr

;;; Create array of gain values
if keyword_set(gn) then begin
   if (n_elements(gn) eq 1) then gain(*) = gn(0) else gain = gn
endif
gainarr    = 0.*rawobj
gainarr(*) = 1.
for i=0,namp-1 do $
   gainarr(amp_xmin(i):amp_xmax(i),amp_ymin(i):amp_ymax(i)) = gain(i)

;;;
;;; Bias Subtract (optional)
;;;

;;; Bias-subtract.
biasnoise_adu = 0.   ; Default is no noise from the bias
if not(keyword_set(nobias)) then begin
   if keyword_set(bias) then begin
      if (size(bias,/tname) eq 'STRING') then begin
         checkexists = findfile(bias,count=nfiles)
         if (nfiles eq 0) then begin
            print,'Bias file not found.'
            return
         endif else begin
            bias2darr = readfits(bias)
            sxaddpar,hobj,'BIAS',filename(bias),'2D bias array'
         endelse
       endif else bias2darr = bias
       print,'Subtracting 2D bias...'
       rawobj    = rawobj - bias2darr
       skyrefarr = skyrefarr - bias2darr
       ; Calculate the noise in the bias frame.  If more than one
       ; chip or amp is used, take the noise measurement from the
       ; first one (this will be an approximation for the other
       ; amps/chips, but shouldn't be too far off).
       if trim then begin
          xlo = amp_xmin(0) - chip_xmin(0)
          ylo = amp_ymin(0) - chip_ymin(0)
       endif else begin
          xlo = amp_xmin(0)
          ylo = amp_ymin(0)
       endelse
       xhi = xlo + amp_xmax(0) - amp_xmin(0)
       yhi = ylo + amp_ymax(0) - amp_ymin(0) 
       biasnoise_adu = robust_sigma(bias2darr(xlo:xhi,ylo:yhi))
   endif else begin
      case instr of
         'luci' : begin
                       print, 'No bias subtraction for LUCIFER'
                  end
         'nirspec': begin
                       print,'No bias subtraction for NIRSPEC'
                    end
         'lrisr':   begin
                       print,'Subtracting LRIS-R bias...'
                       rawobj = lris_subbias(rawobj,/red)
                       if keyword_set(skyref) then $
                          skyrefarr = lris_subbias(skyrefarr,/red,npoly=1)
                    end
         'lrisb':   begin
                       print,'Subtracting LRIS-B bias...'
                       rawobj = lris_subbias(rawobj,/blue)
                       if keyword_set(skyref) then $
                          skyrefarr = lris_subbias(skyrefarr,/blue,npoly=4)
                    end
         'ldss3':   begin
                       print,'Subtracting LDSS3 bias...'
                       rawobj = ldss3_subbias(rawobj)
                       if keyword_set(skyref) then $
                          skyrefarr = ldss3_subbias(skyrefarr)
                    end
         'mike':    begin
                       ; Bias subtracted in mike_reduce.pro
                    end
         'mage':    begin
                       ; Bias subtracted in mage_reduce.pro
                    end
         'dbsp':    begin
                       print,'Subtracting DBSP bias...'
                       rawobj = dbsp_subbias(rawobj)
                       if keyword_set(skyref) then $
                          skyrefarr = dbsp_subbias(skyrefarr)
                    end
         'fors2_1': begin
                       print,'Subtracting FORS2 (CCD1) bias...'
                       rawobj = fors2_subbias(rawobj,ccd=1)
                       if keyword_set(skyref) then $
                          skyrefarr = fors2_subbias(skyrefarr,ccd=1)
                    end
         'fors2_2': begin
                       print,'Subtracting FORS2 (CCD2) bias...'
                       rawobj = fors2_subbias(rawobj,ccd=2)
                       if keyword_set(skyref) then $
                          skyrefarr = fors2_subbias(skyrefarr,ccd=2)
                    end
         else:      print,'Bias subtraction: Instrument name not recognized.'
      endcase
   endelse
endif

;;; Trim pre/overscan regions?
if trim then begin
   rawobj    = rawobj(min(chip_xmin):max(chip_xmax),$
                      min(chip_ymin):max(chip_ymax))
   skyrefarr = skyrefarr(min(chip_xmin):max(chip_xmax),$
                         min(chip_ymin):max(chip_ymax))
   gainarr   = gainarr(min(chip_xmin):max(chip_xmax),$
                       min(chip_ymin):max(chip_ymax))
endif


; Create subtracted sky frame
subobj = rawobj - skyrefarr

; Rotate object (optional).
if rotate then begin
   hrotate,rawobj,hobj,rawobj,hobj,3
   subobj = rotate(subobj,3)
   gainarr = rotate(gainarr,3)
endif

;;;
;;; Dark-subtract (optional)
;;;

; Note for NIRSPEC: Due to apparent structure in the read-noise,
; dark-subtracting object frames has been found to add undue noise and
; interfere with sky-subtraction.  In addition, subtracting a refernce
; sky frame effectively gets rid of any dark levels.  Therefore, in
; order to measure the raw sky image (to get variance estimates),
; subtract a uniform dark level according to the exposure time.  Use
; 2D dark frame for bad pixel map only.  If the /nodarksub keyword is set,
; then calculate but do not subtract dark counts.

if keyword_set(darkframe) then begin
   if (size(darkframe,/tname) eq 'STRING') then begin
      checkexists = findfile(darkframe,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Dark exposure file not found.'
         return
      endif else begin
         darkarr = readfits(darkframe)
         sxaddpar,hobj,'DARKEXP',filename(darkframe),$
            'Dark exposure (for bad pixel map only)'
      endelse
   endif else darkarr = darkframe
endif else darkarr = 0.

;;; Estimate dark level by exposure time.  Default is no dark counts.
darkmodel  = 0.
edarklevel = 0.
if keyword_set(texp) then texposure = texp else texposure = 0.
if (total(darkcurrent) ne 0) then begin
   if (n_elements(darkcurrent)) eq 1 then begin
      dc_arr = fltarr(n_elements(namp))
      dc_arr(*) = darkcurrent(0)
   endif else dc_arr = darkcurrent
   edarklevel = dc_arr * texposure
   ; Compute expected dark counts for each amplifier
   darkmodel = 0.*rawobj
   for i=0,namp-1 do  begin
      x0 = amp_xmin(i)
      x1 = amp_xmax(i)
      y0 = amp_ymin(i)
      y1 = amp_ymax(i)  
      darkmodel(x0:x1,y0:y1) = edarklevel(i)/gain(i)
   endfor
   ; Subtract dark current calculated for each amplifier(optional)
   if keyword_set(nodarksub) then begin
      print,'Calculated dark counts not subtracted.  Will be used for'
      print,'   error estimates only.'
      print,' '
   endif else begin
      print,'Subtracting calculated dark counts'
      print,' '
      rawobj = rawobj - darkmodel
   endelse
endif
darkcounts = edarklevel / gain
print,'Dark counts (ADU) = ',darkcounts
darkprint = ''
for i=0,n_elements(darkcounts)-1 do darkprint = darkprint + $
                                            strtrim(darkcounts(i),2) + ' '
sxaddpar,hobj,'DARKCNTS',darkprint,'Calculated dark counts'

;;; Read noise: If read noise is differt for each amplifier, then 
;;; construct an array of read noise values.
; User-specified read noise?
if keyword_set(rdnoise) then erdnoise = rdnoise
if (n_elements(erdnoise) gt 1) then begin
   rdnoise_arr    = 0.*rawobj
   for i=0,namp-1 do $
      rdnoise_arr(amp_xmin(i):amp_xmax(i),amp_ymin(i):amp_ymax(i)) = erdnoise(i)
   multi_rdnoise = 1
endif else multi_rdnoise = 0

;;; Set number of coadds
; Number of coadds for object
if keyword_set(objcoadds) then ncoadds = objcoadds else ncoadds=1.
;; Number of coadds for dark.
;if keyword_set(darkcoadds) then dkcoadds = darkcoadds else dkcoadds=ncoadds
; Currently using model dark counts only - no read noise from dark.
dkcoadds = 0.

;;;
;;; Flat-field (optional)
;;;

; Optional - Adjust counts to account for differences in amplifier gains.
if keyword_set(fixgain) then begin
   gainratio = gain / mean(gain)
   for i=0,namp-1 do  begin
      x0 = amp_xmin(i)
      x1 = amp_xmax(i)
      y0 = amp_ymin(i)
      y1 = amp_ymax(i)
      rawobj(x0:x1,y0:y1) = rawobj(x0:x1,y0:y1) / gainratio(i)
      subobj(x0:x1,y0:y1) = subobj(x0:x1,y0:y1) / gainratio(i)
   endfor      
endif

if keyword_set(flatfield) then begin
   if (size(flatfield,/tname) eq 'STRING') then begin
      checkexists = findfile(flatfield,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Flat field file not found.'
         return
      endif else begin
         flatarr = readfits(flatfield)
         sxaddpar,hobj,'FLATFLD',filename(flatfield),'Flat field'
      endelse
   endif else flatarr = flatfield
   ;; Fix zero-valued pixels by subbing in value small enough to be
   ;; included in mask
   ;flatzerols = where(flatarr eq 0,n_flatzero)
   ;if (n_flatzero ne 0) then flatarr(flatzerols) = 0.01
endif else begin
   ; If no input flat field, create a dummy array full of 1s
   flatarr = 1. + 0.*subobj
endelse
nonzeroflat_ls = where(flatarr gt 0)
rawobj(nonzeroflat_ls) = rawobj(nonzeroflat_ls) / flatarr(nonzeroflat_ls)
subobj(nonzeroflat_ls) = subobj(nonzeroflat_ls) / flatarr(nonzeroflat_ls)

if keyword_set(flatvar) then begin
   if (size(flatvar,/tname) eq 'STRING') then begin
      checkexists = findfile(flatvar,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Flat field variance file not found.'
         return
      endif else begin
         flatvararr = readfits(flatvar)
         sxaddpar,hobj,'FLATVAR',filename(flatvar),'Flat field variance'
      endelse
   endif else begin
      if (n_elements(flatvar) eq 1) then begin
         if (flatvar eq -1) then flatvararr = 0.*rawobj
      endif else begin
         flatvararr = flatvar
      endelse
   endelse
endif else begin
   flatvararr = 0.*rawobj
endelse

;;;
;;; Create a bad pixel map by identifying outlying pixels in the
;;; dark (optional) and flat field (optional).  Also flag 
;;; pre/overscan regions and known bad pixels.
;;;

; Default is no masking
badmap = 1 + 0*rawobj

;;; Bad pixels are flagged with a 0.  Non-illuminated regions are flagged
;;; with a -1.  Masking should be done in that order.

;;; Flag bad pixels.

if keyword_set(darkframe) then begin
   print,'Creating bad pixel map from dark frame...'
   for i=0,namp-1 do  begin
      x0 = amp_xmin(i)
      x1 = amp_xmax(i)
      y0 = amp_ymin(i)
      y1 = amp_ymax(i)
      dark_section = darkarr(x0:x1,y0:y1)
      med_dark  = median(dark_section)
      sdev_dark = robust_sigma(dark_section)
      darkbadls = where(dark_section ge med_dark + 5*sdev_dark or $
                        dark_section le med_dark - 5*sdev_dark,n_darkbad)
      if (n_darkbad ne 0) then begin
         badmap_section = badmap(x0:x1,y0:y1)
         badmap_section(darkbadls) = 0
         badmap(x0:x1,y0:y1) = badmap_section
      endif
   endfor
endif

if keyword_set(flatfield) then begin
   ; Rough test of whether flat has been normlized: 
   one_ls = where(flatarr ge 0.9 and flatarr le 1.1,n_one)
   if n_one ge 0.1*n_elements(flatarr) then flatnormal=1 else flatnormal=0
   ; If the flat has been normalized, use it to locate bad pixels.
   if flatnormal then begin
      print,'Normalized flat field detected.'
      print,'Creating bad pixel map from flat field...'
      flatbadls = where((flatarr gt 0 and flatarr lt 0.7) or $
                        flatarr gt 1.3,n_flatbad)
      if (n_flatbad ne 0) then badmap(flatbadls) = 0.
   endif
   ; Flag pixels where the flat is <= 0 or the flatvar is < 0 
   ; as non-illuminated
   flatzero_ls = where(flatarr le 0 or flatvararr lt 0,n_flatzero)
   if (n_flatzero ne 0) then badmap(flatzero_ls) = -1
endif else flatnormal=0

;;; Special masking for specific instruments
if (instr eq 'lrisr') then begin
   ; Omit first and last columns due to edge effects
   case xsz of
      2248: begin
               ; First and last columns
               badmap(40,*)   = 0.
               badmap(2087,*) = 0.
               ; Bad columns and pixels
               if (ysz eq 1000) then begin   ; Longslit windowing
                  badmap(996,0:444) = 0.
                  for i=0,99 do badmap(996,10*i+4) = 0.
                  badmap(1000,0:446) = 0.
               endif else begin              ; Full-frame windowing
                  badmap(996,0:1493)  = 0.   
                  badmap(1000,0:1495) = 0.
                  badmap(862,0:722)   = 0.
               endelse
            end
      2048: begin
               ; First and last columns
               badmap(0,*)   = 0.
               badmap(2047,*) = 0.
               ; Bad columns and pixels
               if (ysz eq 1000) then begin   ; Longslit windowing
                  badmap(956,0:444) = 0.
                  for i=0,99 do badmap(946,10*i+4) = 0.
                  badmap(960,0:446) = 0.
               endif else begin
                  badmap(956,0:1493) = 0.
                  badmap(960,0:1495) = 0.
                  badmap(842,0:722)  = 0.
               endelse
            end
      else: print,'LRIS-R size not recognized.  Skipping bad column masking.'
   endcase
endif

;;; Flag non-illuminated areas.

;;; Flag prescan and overscan regions
if not(trim) then begin
   if (total(prescan_x) ne -1) then begin
      badmap(prescan_x(0):prescan_x(1),prescan_y(0):prescan_y(1)) = -1
   endif
   if (total(oscan_x) ne -1) then begin
      badmap(oscan_x(0):oscan_x(1),oscan_y(0):oscan_y(1)) = -1
   endif
endif

;;; Flag non-finite values in the calibration arrays
finite_ls = where(finite(wavegrid) eq 1 and $
                  finite(slitgrid) eq 1,n_finite,$
                  complement=nonfinite_ls,$
                  ncomplement=n_nonfinite)
if (n_nonfinite ne 0) then badmap(nonfinite_ls) = -1

;;; Optional - constrain region of slit over which to locate
;;; objects.
if keyword_set(subrange) then begin
   if (n_elements(subrange) ne 2) then begin
      print,'Keyword subrange must be a 2-element array!'
      return
   endif else begin
      slitmin = subrange(0)
      slitmax = subrange(1)
   endelse
endif

;;; Optional - specify the area of the detector over which to
;;; locate objects.
if keyword_set(findrange) then begin
   if (dispdir eq 'x') then begin
      find_wavegrid = wavegrid(findrange(0):findrange(1),*)
      find_slitgrid = slitgrid(findrange(0):findrange(1),*)
   endif  else begin 
      find_wavegrid = wavegrid(*,findrange(0):findrange(1))
      find_slitgrid = slitgrid(*,findrange(0):findrange(1))
   endelse
   ; Constrain slit position over which to bound the search
   ; wavelength?
   if keyword_set(stdpos) then begin
      slit_lobound = stdpos-1
      slit_hibound = stdpos+1
   endif else begin
      if keyword_set(slitpos) then begin
         slit_lobound = slitpos-1
         slit_hibound = slitpos+1
      endif else begin
         slit_lobound = slitmin
         slit_hibound = slitmax
      endelse
   endelse
   find_region = where(find_slitgrid ge slit_lobound and $
                       find_slitgrid le slit_hibound,n_find_region)
   if (n_find_region ne 0) then begin
      wavemin = min(find_wavegrid(find_region))
      wavemax = max(find_wavegrid(find_region))
   endif else begin
      print,'ERROR - Boundaries set by findrange do not cover any of the'
      print,'        expected slit range.'
      return
   endelse
endif

;;;
;;; Optional - Define object-centric slit coordinates by tracing a
;;; standard star with respect to the input slit position grid.
;;;

; Default is successful (or no) object trace
trace_ok    = 1
trace_terms = -1

; Flag that the object the standard star has been traced
object_traced = 0

; Limits for tracing objects.  Avoid running off the side of the chip
; or onto margins of the chip that are fully masked.
xcheck_badmap = fltarr(xsz)
for i=0,xsz-1 do xcheck_badmap(i) = max(badmap(i,*))
ycheck_badmap = fltarr(ysz)
for i=0,ysz-1 do ycheck_badmap(i) = max(badmap(*,i))
unmasked_xmin = min(where(xcheck_badmap gt 0))
unmasked_xmax = max(where(xcheck_badmap gt 0))
unmasked_ymin = min(where(ycheck_badmap gt 0))
unmasked_ymax = max(where(ycheck_badmap gt 0))
if (dispdir eq 'x') then begin
   chip_range = [chip_xmin > unmasked_xmin,chip_xmax < unmasked_xmax]
   chip_bound = [chip_ymin > unmasked_ymin,chip_ymax < unmasked_ymax]
endif else begin
   chip_range = [chip_ymin > unmasked_ymin,chip_ymax < unmasked_ymax]
   chip_bound = [chip_xmin > unmasked_xmin,chip_xmax < unmasked_xmax]
endelse
if keyword_set(tracerange) then $
   trace_range = tracerange else $
   trace_range = 0
if keyword_set(stdtrace) then begin
   print,' '
   print,'Getting trace from standard...'
   if (size(stdtrace,/tname) eq 'STRING') then begin
      checkexists = findfile(stdtrace,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Trace file not found.'
         return
      endif else begin
          stdtracearr = readfits(stdtrace,htrace)
          sxaddpar,hobj,'TRACE',filename(stdtrace),'Standard star trace.'
      endelse
   endif else stdtracearr = long(stdtrace)
   ; Look for trace object in the same region as the science object.
   findls = where(wavegrid ge wavemin and wavegrid le wavemax and $
                  slitgrid ge slitmin and slitgrid le slitmax)
   wavelist      = wavegrid(findls)
   order         = sort(wavelist)
   waveorder     = wavelist(order)
   findorder     = findls(order)
   slitorder     = slitgrid(findorder)
   stdtraceorder = stdtracearr(findorder)
   ; Divide region into chunks and compute the median value in each
   ; chunk.
   ; *** Should optimize this using HISTOGRAM ***
   slitlen  = slitmax - slitmin
   chunklen = 5.
   nchunks  = floor(slitlen/chunklen)
   poslist  = fltarr(nchunks)
   medlist  = fltarr(nchunks)
   for i=0,nchunks-1 do begin
      minpos = slitmin + i*chunklen
      maxpos = slitmin + (i+1)*chunklen
      poslist(i) = 0.5*(minpos+maxpos)
      thischunk  = where(slitorder ge minpos and slitorder lt maxpos,$
                         n_thischunk)
      ; Some points may fall in chip gaps
      if (n_thischunk ne 0) then begin
         medlist(i) = median(stdtraceorder(thischunk))
      endif else medlist(i) = 0. 
   endfor
   ; Pick out the maximum point
   heightguess = (max(medlist))
   centguess   = poslist(where(medlist eq max(medlist)))
   centguess   = centguess(0)
   ; Optional - user-specified position for the standard
   if keyword_set(stdpos) then centguess = stdpos
   ; Fit a Gaussian at that point and record trace object profile
   ; parameters.
   gfit = gaussfit(slitorder,stdtraceorder,gterms,nterms=4,$
                   yerror=gerr,estimates=[heightguess,centguess,2,0])
   stdheight   = abs(gterms(0))
   stdposition = gterms(1)
   stdsig      = abs(gterms(2))
   ; Trace the input trace object.
   obj_type='bright'
   ; Use a narrow window for tracing a standard.
   n_side = 15 < (max(slitorder) - min(slitorder))/3.
   obj_slitgrid = lr_traceobject(stdtracearr,slitgrid,stdposition,$
                                 dispdir,range=chip_range,$
                                 slitrange=[slitmin,slitmax],$
                                 tracerange=trace_range,$
                                 objtype=obj_type,wavegrid=wavegrid,$
                                 pixtr=pixtr,slittr=slittr,slitfit=slitfit,$
                                 ;waverange=[wavemin,wavemax],$
                                 nside=n_side,noshow=no_show,$
                                 bound=chip_bound,fitterms=trace_terms)
   trace_ok = 1
   if (total(obj_slitgrid) eq -1) then begin
      print,'WARNING: Object trace failed.  Using input slit positions.'
      obj_slitgrid = slitgrid
      trace_ok = 0
   endif
   object_traced = 1
   ; Use these object-centric coordinates for all subsequent
   ; masking and extraction.  Do not trace the science object.
   traceobj = 0
endif else begin
   ; Default is to start by assuming that objects follow lines
   ; of constant slit position (no atmospheric diffraction).
   obj_slitgrid = slitgrid
endelse

; Optional - specify via keyword whether to trace object
if keyword_set(trace_obj) then begin
   if (trace_obj lt 0) then traceobj = 0 else traceobj = 1
endif

;;;
;;; Locate objects by performing a rough fit to the sky over a narrow
;;; wavelength range.  Use the reference sky-subtracetd frame if one
;;; has been created.
;;;

; Optional - Set spacing of break points via keywords (defaults are 
; values in included file).
if keyword_set(bkspace)     then bk_space     = bkspace
if keyword_set(min_bkspace) then min_bk_space = min_bkspace

; Optional - Use keyword to specify maximum polynomial order used to fit 
; slit illumination function.  Default is value in included file.
; 1 = constant, 2 = linear, etc.
; Note: 'cholesky_band' (part of bspline fitting package) seems
; to sometimes crash if you attempt polynomials in xtwo > 2.  
if keyword_set(illumpoly) then illum_poly = illumpoly  > 1
if keyword_set(noskysub) then illum_poly = 1

; If an individual CCD is speficied, set the min/max slit positions to
; be the min/max that are covered by that chip.
if (keyword_set(ccd)) then begin
   if (ccd gt nchip) then begin
      print,'Warning: CCD number too large!'
      return
   endif 
   chipindex  = ccd - 1
   if trim then begin
      chip_xlo  = 0.
      chip_xhi  = chip_xmax(chipindex) - chip_xmin(chipindex)
      chip_ylo  = 0.
      chip_yhi  = chip_ymax(chipindex) - chip_ymin(chipindex)
   endif else begin
      chip_xlo  = chip_xmin(chipindex)
      chip_xhi  = chip_xmax(chipindex)
      chip_ylo  = chip_ymin(chipindex)
      chip_yhi  = chip_ymax(chipindex)
   endelse
   slitedge1 = min(slitgrid(chip_xlo:chip_xhi,chip_ylo:chip_yhi)) > slitedge1
   slitedge2 = max(slitgrid(chip_xlo:chip_xhi,chip_ylo:chip_yhi)) < slitedge2
   ; Mask out other chips, in case the slit positions covered by the 
   ; specified chip extend over other chips.
   for i=0,nchip-1 do begin
      if (i ne chipindex) then begin
         chip_xlo  = chip_xmin(i)
         chip_xhi  = chip_xmax(i)
         chip_ylo  = chip_ymin(i)
         chip_yhi  = chip_ymin(i)
         badmap(chip_xlo:chip_xhi,chip_ylo:chip_yhi) = -1 ; Non-illuminated
      endif
   endfor
endif

; Default is no masking
objheight       = 0.
objposition     = 0.
objsig          = 0.
fit_entire_slit = 0.
if keyword_set(noobject) then begin   ; No object
   if not(keyword_set(subrange)) then fit_entire_slit = 1
endif

; Parameters for estimating variance
if keyword_set(skyref) then begin
   skyfactor  = 2. 
   biasfactor = 0. ; Both frames have had the same bias subtracted
   darkfactor = 2.
endif else begin
   skyfactor  = 1.  
   biasfactor = 1.
   darkfactor = 2.
endelse

; Number of sigma around object position to mask
if keyword_set(bright) then begin
   case instr of
      'mike':      nsig = 5.
      'hires':     nsig = 5.
      'n_echelle': nsig = 5.
      'mage'     : nsig = 5.
      else:        nsig = 7.
   endcase
endif else begin
   case instr of
      'mike':      nsig = 2. ;3.5   ; Narrow slit
      'hires':     nsig = 3.5 
      'n_echelle': nsig = 3.
      'mage':      nsig = 2.5 ;3.5
      else:        nsig = 5.
   endcase
endelse

; Default is not to locate and mask any objects
domasking = 0
if keyword_set(maskall) then domasking = 1
if not(keyword_set(noobject)) then begin
   ; Object to extract on slit.  Fit profile?
   if (keyword_set(slitpos) and keyword_set(fwhm) and $
       keyword_set(force)) then begin
      print,'Using user-specified object profile parameters'
      objposition = slitpos
      objsig      = fwhm / (2.*sqrt(2.*alog(2)))
   endif else begin
      domasking = 1
   endelse
endif

if domasking then begin

   ; Locate object(s) and mask when fitting sky

   ; Plot in color
   loadct,39,/silent

   ; Identify region in which to look for object
   findls = where(wavegrid ge wavemin and wavegrid le wavemax and $
                  slitgrid ge slitmin and slitgrid le slitmax and $
                  badmap ne -1)

   ; Look for object in reference sky-subtracted frame (if any)
   wavelist      = wavegrid(findls)
   order         = sort(wavelist)
   waveorder     = wavelist(order)
   findorder     = findls(order)
   slitorder     = slitgrid(findorder)
   obj_slitorder = obj_slitgrid(findorder)
   rawobjorder   = rawobj(findorder)
   subobjorder   = subobj(findorder)
   gainorder     = gainarr(findorder)
   flatorder     = flatarr(findorder)
   flatvarorder  = flatvararr(findorder)
   badorder      = badmap(findorder)
   if multi_rdnoise then rdnoiseorder = rdnoise_arr(findorder) $
      else rdnoiseorder = erdnoise

   ; Estimate variance in the difference array, or the raw frame, if
   ; no differencing done.
   varorder = lr_makevariance(rawobjorder,$
                              flatfield=flatorder,$
                              sky=skyfactor*rawobjorder,$  ; Estimate
                              flatvar=flatvarorder,$
                              darkframe=darkfactor*darkmodel,$
                              readnoise=rdnoiseorder,$
                              objcoadds=skyfactor*ncoadds,$
                              darkcoadds=darkfactor*dkcoadds,$
                              gain=gainorder,$
                              biasnoise=biasfactor*biasnoise_adu)
   invvarorder = 1./varorder
   
   ; Mask bad pixels
   badls = where(badorder le 0,n_bad)
   if (n_bad ne 0) then invvarorder(badls) = 0.
   ; Fit sky in narrow region
   rej      = 5.        ; Number of sigma for rejecting pixels
   ; Use at most a linear illumination function
   n_poly   = 2 < illum_poly
   xtwo     = slitorder
   max_iter = 20.
   if keyword_set(noskysub) then begin
      skyfitorder = 0.*subobjorder
   endif else begin
      sset=bspline_iterfit(waveorder,subobjorder,x2=xtwo,npoly=n_poly,$
                           invvar=invvarorder,bkspace=bk_space,$
                           yfit=skyfitorder,upper=rej,lower=rej,$
                           maxiter=max_iter)
   endelse
   noskyorder = subobjorder - skyfitorder
   ; Divide sky-subtracted region of the slit into small chunks and
   ; compute the median value in each chunk.
   obj_slitmin = min(obj_slitorder)
   obj_slitmax = max(obj_slitorder)
   slitlen  = obj_slitmax - obj_slitmin
   chunklen = 2. ;4.
   nchunks  = floor(slitlen/chunklen)
   poslist  = fltarr(nchunks)
   medlist  = fltarr(nchunks)
   badlist  = fltarr(nchunks)
   for i=0,nchunks-1 do begin
      minpos = obj_slitmin + i*chunklen
      maxpos = obj_slitmin + (i+1)*chunklen
      poslist(i) = 0.5*(minpos+maxpos)
      thischunk  = where(obj_slitorder ge minpos and $
                         obj_slitorder lt maxpos,n_thischunk)
      ; Some points may fall in chip gaps
      if (n_thischunk ne 0) then begin
         medlist(i) = median(noskyorder(thischunk))
         badlist(i) = median(badorder(thischunk))
      endif else medlist(i) = 0. 
   endfor

   ; Comment out to be able to fit orders near edge of array
   ;; Flag bad regions
   ;temp_ls = where(badlist eq 0,n_temp)
   ;if (n_temp ne 0) then medlist(temp_ls) = -99999

   ; Display cross-section (opional)
   if not(no_show) then begin
      chan,10
      ymax = 3*max(medlist) > 3*robust_sigma(noskyorder)
      plot,obj_slitorder,noskyorder,psym=3,/xstyle,yrange=[-ymax,ymax]
   endif

   ; Flags for multiple object masking
   n_objects_masked    = 0
   masking_ref_objects = 0
   done_masking        = 0

   ; If the object position has been set with respect to a standard
   ; star trace and the width of the object profile has also been set,
   ; then skip masking the first object.
   if (keyword_set(stdtrace) and keyword_set(stdoffset) and $
       keyword_set(fwhm)) then begin
      objposition      = stdposition + stdoffset
      objsig           = fwhm / (2.*sqrt(2*alog(2)))
      n_objects_masked = 1
      temp_ls = where(poslist ge objposition-7*objsig and $
                      poslist le objposition+7*objsig,n_temp)
      if (n_temp ne 0) then medlist(temp_ls) = -99999
      if not(keyword_set(maskall)) then done_masking = 1
   endif

   ; Optional - force the program to use user-provided values
   if (keyword_set(force) and keyword_set(slitpos) and keyword_set(fwhm)) $
      then begin
      objposition = slitpos
      objsig      = fwhm / (2.*sqrt(2*alog(2)))
      n_objects_masked = 1
      temp_ls = where(poslist ge objposition-7*objsig and $
                      poslist le objposition+7*objsig,n_temp)
      if (n_temp ne 0) then medlist(temp_ls) = -99999
      if not(keyword_set(maskall)) then done_masking = 1
   endif

   ; Enter loop to locate and fit Gaussian profiles to objects

   while not(done_masking) do begin

      if masking_ref_objects then print,'Masking reference object' else $
         print,'Masking object'

      ; If a slit position is specified, mask that object first.
      ; Otherwise, the brightest object present will be the one
      ; ultimately extracted. 
      ; ****Note: The specified slit position may be somewhat
      ;           offset from the "object-centric" slit position.
      ;           If the difference is too great the program may not
      ;           find the object, but it should be obvious if this
      ;           happens.
      if (keyword_set(slitpos) and n_objects_masked eq 0 and $
          not(masking_ref_objects)) then begin
         centguess   = slitpos
         heightguess = medlist(where(abs(poslist-centguess) eq $
                                     min(abs(poslist-centguess))))
         heightguess = heightguess(0)
         goodls      = where(noskyorder ge -2*heightguess and $
                             noskyorder le  10*heightguess and $
                             badorder gt 0 and $
                             abs(obj_slitorder-centguess) lt 100) 
         user_input  = 1
      endif else begin
         if not(masking_ref_objects) then begin
            ; Look for positive profiles
            heightguess = max(medlist(where(medlist ne -99999)))
            centguess   = poslist(where(medlist eq heightguess))
            centguess   = centguess(0)
            goodls      = where(noskyorder ge -2*heightguess and $
                                noskyorder le  10*heightguess and $
                                badorder gt 0 and $
                                abs(obj_slitorder-centguess) lt 100)
            user_input  = 0
         endif else begin
            ; Look for negative profiles
            heightguess = min(medlist(where(medlist ne -99999)))
            centguess   = poslist(where(medlist eq heightguess))
            centguess   = centguess(0)
            goodls      = where(noskyorder le -2*heightguess and $
                                noskyorder ge  10*heightguess and $
                                badorder gt 0 and $
                                abs(obj_slitorder-centguess) lt 100)
            user_input  = 0
         endelse
      endelse
      ; If masking multiple objects, then make sure that this at 
      ; least a 4-sigma outlier from the noise in the median list.
      sdev_median = stddev(medlist(where(medlist ne -99999 and $
                                         medlist ne heightguess)))
      if (not(keyword_set(maskall) and abs(heightguess) lt 4*sdev_median) or $
          user_input) then begin
         goodobj_slitorder  = obj_slitorder(goodls)
         goodnoskyorder     = noskyorder(goodls)
         gfit = gaussfit(goodobj_slitorder,goodnoskyorder,gterms,nterms=4,$
                         yerror=gerr,estimates=[heightguess,centguess,2,0])
         ; Reject outliers and refit
         if not(keyword_set(bright)) then begin
            resid = goodnoskyorder - gfit
            goodls = where(abs(resid) le 3*gerr)
            goodobj_slitorder = goodobj_slitorder(goodls)
            goodnoskyorder    = goodnoskyorder(goodls)
            gfit = gaussfit(goodobj_slitorder,goodnoskyorder,gterms,nterms=4,$
                            yerror=gerr,estimates=gterms,sigma=sig_gterms)
         endif
         print,'Object profile fit terms: ',gterms
         ; If the object profile is not reasonable, discard it and
         ; quit masking objects.  Otherwise, continue.
         if ((gterms(1) lt slitmin or gterms(1) gt slitmax or $
              abs(gterms(2)) lt 0.5 or abs(gterms(2)) gt 20. or $
              (masking_ref_objects and gterms(0) gt 0) or $
              (not(masking_ref_objects) and gterms(0) lt 0)) $
             and not(user_input)) then begin
            print,'Unreasonable object profile found - discarding and '$
                  +'exiting masking loop.'
            print,'Try entering parameters by hand.'
            print,'Note: This may be an artifact of the program if the object'
            print,'      profile is very wide.'
            done_masking=1
         endif else begin
            ; If this is the first object being masked and extraction
            ; is being done, record temporary profile parameters for
            ; later refinement and extraction.  Otherwise, flag pixels
            ; to be permanently masked.  Optionally use user-specified
            ; values.
            if (n_objects_masked eq 0 and not(masking_ref_objects) and $
                not(keyword_set(noobject))) then begin
               objheight   = abs(gterms(0))
               if (keyword_set(stdtrace) and keyword_set(stdoffset)) then $
                  objposition = stdposition + stdoffset else $
                  objposition = gterms(1)
               if keyword_set(fwhm) then $
                  objsig = fwhm / (2.*sqrt(2*alog(2))) else $
                  objsig = abs(gterms(2))
            endif else begin
               temp_objls = where($
                            obj_slitgrid ge gterms(1)-nsig*abs(gterms(2)) and $
                            obj_slitgrid le gterms(1)+nsig*abs(gterms(2)) and $
                            badmap ne -1,temp_n_objpix)
               if (temp_n_objpix ne 0) then begin
                  if masking_ref_objects then ref_mask(temp_objls) = 0 else $
                     badmap(temp_objls) = 0
               endif
            endelse
            ; Suppress those points in the median list covered by the 
            ; masked object.
            temp_ls = where(poslist ge gterms(1)-nsig*abs(gterms(2)) and $
                            poslist le gterms(1)+nsig*abs(gterms(2)),n_temp)
            if (n_temp ne 0) then medlist(temp_ls) = -99999
            ; Suppress those points in the fitted region covered by the
            ; masked object, in order to avoid masking objects more
            ; than once.
            temp_bad_ls = $
                where(obj_slitorder ge gterms(1)-nsig*abs(gterms(2)) and $
                      obj_slitorder le gterms(1)+nsig*abs(gterms(2)),$
                      n_temp_bad)
            if (n_temp_bad ne 0) then badorder(temp_bad_ls) = 0.
            ; Display profile (optional).
            if not(no_show) then begin
               chan,10
               plot_order = sort(goodobj_slitorder)
               oplot,goodobj_slitorder(plot_order),gfit(plot_order),color=250
            endif
            ; Keep count of the total number of objects masked.
            n_objects_masked = n_objects_masked + 1
         endelse
         ; Are there any unmasked pixels left?
         unmasked_ls = where(medlist ne -99999,n_unmasked)
         if (n_unmasked eq 0) then done_masking = 1
      endif else begin
         ; No more significant object candidates.
          if not(masking_ref_objects) then begin
             print,'Number of objects masked = '+strtrim(n_objects_masked,2)
             if keyword_set(skyref) then begin
                masking_ref_objects = 1
                n_objects_masked    = 0
             endif else done_masking = 1
          endif else begin
             print,'Number of reference objects masked = ' + $
                strtrim(n_objects_masked,2)
             done_masking = 1
          endelse
      endelse

      ; Mask more than one object?
      if not(keyword_set(maskall)) then begin
         if (not(masking_ref_objects) and keyword_set(skyref)) then begin
            masking_ref_objects = 1
            n_objects_masked    = 0
         endif else begin
            done_masking = 1 
         endelse
      endif 

      ; Too many objects masked?
      if (n_objects_masked ge 20) then begin
         print,'Masking too many objects!'
         done_masking = 1
      endif

   endwhile

   ; If multiple CCDs are included in the input array, limit sky
   ; subtration to only those CCDs on which the object was found
   ; (preferably just one).  ***There may be an offset between 
   ; the object-centric slit position and the normal slit position
   ; but it shouldn't matter here.***
   if (nchip gt 1 and not(keyword_set(noobject))) then begin
      firstgoodchip = 1
      for i=0,nchip-1 do begin
         chip_xlo  = chip_xmin(i)
         chip_xhi  = chip_xmax(i)
         chip_ylo  = chip_ymin(i)
         chip_yhi  = chip_ymin(i)
         chip_slitmin = min(slitgrid(chip_xlo:chip_xhi,chip_ylo:chip_yhi))
         chip_slitmax = max(slitgrid(chip_xlo:chip_xhi,chip_ylo:chip_yhi))
         if (objposition ge chip_slitmin and $
             objposition le chip_slitmax) then begin
            if firstgoodchip then begin
               temp_slitmin = chip_slitmin
               temp_slitmax = chip_slitmax
               firstgoodchip = 0
            endif else begin
               temp_slitmin = temp_slitmin < chip_slitmin
               temp_slitmax = temp_slitmax > chip_slitmax
            endelse
         endif
      endfor
      slitedge1 = slitedge1 > temp_slitmin
      slitedge2 = slitedge2 < temp_slitmax
   endif

endif

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

; Refine which pixels to fit based on object position (if any)

; Default is to fit the entire slit
slitmin  = slitedge1
slitmax  = slitedge2

; Allow keywords to override slit min and max
if not(fit_entire_slit) then begin
   ; Take the maximum allowable range as determined by all
   ; keywords.
   if keyword_set(subrange) then begin
      if (n_elements(subrange) ne 2) then begin
         print,'Keyword subrange must be a 2-element array!'
         return
      endif else begin
         slitmin = subrange(0) ;> slitmin
         slitmax = subrange(1) ;< slitmax
      endelse
   endif
   if not(keyword_set(noobject)) then begin
      if keyword_set(subwidth) then sub_width = subwidth
      slitmin  = (objposition - sub_width/2.) > slitmin
      slitmax  = (objposition + sub_width/2.) < slitmax
   endif
endif

; Only fit wavelength ranges over which at least 50% of the designated
; slit range is covered by usable pixels.
temp_exposed  = where(slitgrid ge slitmin and slitgrid le slitmax and $
                      badmap gt 0,n_exposed)
if (n_exposed eq 0) then begin
   print,'LONGSLIT_REDUCE: Error - No exposed pixels within range.'
   return
endif
temp_wavelist = wavegrid(temp_exposed)
temp_slitlist = slitgrid(temp_exposed)
; Calculate the rough dispersion in units of input wavelength/pixel.
; Ignores the fact that the dispersion may not run exactly along x or
; y.
index = array_indices(rawobj,temp_exposed)
if (dispdir eq 'x') then temp_dispindex = double(index(0,*))
if (dispdir eq 'y') then temp_dispindex = double(index(1,*))
temp_waveterms = poly_fit(temp_dispindex,temp_wavelist,3)
temp_displist = dblarr(n_exposed)
for i=1,2 do $
   temp_displist = temp_displist + i*temp_waveterms(i)*(temp_dispindex^(i-1))
temp_displist  = abs(temp_displist)
temp_disp      = median(temp_displist)
temp_wavelist  = wavegrid(temp_exposed)
temp_slitlist  = slitgrid(temp_exposed)
lo_wavebound   = min(temp_wavelist)
hi_wavebound   = max(temp_wavelist)
required_frac  = 0.5
lo_limit_found = 0
hi_limit_found = 0
i=0L
while not(lo_limit_found) do begin
   w_lo = min(temp_wavelist) + i*temp_disp
   w_hi = min(temp_wavelist) + (i+1)*temp_disp
   if (w_lo ge max(temp_wavelist)) then begin
      print,'Error: Slit range is nowhere more than 50% on the chip'
      return
   endif
   region = where(temp_wavelist ge w_lo and temp_wavelist le w_hi,n_region)
   if (n_region gt 0) then begin
      slit_region = temp_slitlist(region)
      if (max(slit_region)-min(slit_region) ge $
          required_frac*(slitmax-slitmin)) then begin
         lo_limit_found = 1
         lo_wavebound   = w_lo
      endif
   endif
   i = i+1
endwhile
i=0
while not(hi_limit_found) do begin
   w_hi   = max(temp_wavelist) - i*temp_disp
   w_lo   = max(temp_wavelist) - (i+1)*temp_disp
   if (w_hi le min(temp_wavelist)) then begin
      print,'Error: Slit range is nowhere more than 50% on the chip'
      return
   endif
   region = where(temp_wavelist ge w_lo and temp_wavelist le w_hi,n_region)
   if (n_region ne 0) then begin
      slit_region = temp_slitlist(region)
      if (max(slit_region)-min(slit_region) ge $
          required_frac*(slitmax-slitmin)) then begin
         hi_limit_found = 1
         hi_wavebound   = w_hi
      endif
   endif
   i = i+1
endwhile

exposed = where(slitgrid ge slitmin and slitgrid le slitmax and $
                wavegrid ge lo_wavebound and wavegrid le hi_wavebound and $
                badmap ne -1,n_exposed)
wavelist      = wavegrid(exposed)
order         = sort(wavelist)
exposedorder  = exposed(order)
waveorder     = wavelist(order)
slitorder     = slitgrid(exposedorder)
obj_slitorder = obj_slitgrid(exposedorder)
rawobjorder   = rawobj(exposedorder)
subobjorder   = subobj(exposedorder)
darkorder     = darkarr(exposedorder)
flatorder     = flatarr(exposedorder)
flatvarorder  = flatvararr(exposedorder)
gainorder     = gainarr(exposedorder)
obj_badorder  = badmap(exposedorder)
ref_maskorder = ref_mask(exposedorder)

if multi_rdnoise then rdnoiseorder = rdnoise_arr(exposedorder) $
   else rdnoiseorder = erdnoise

; Calculate the rough dispersion in units of input wavelength/pixel.
; Ignores the fact that the dispersion may not run exactly along x or
; y.
indexorder = array_indices(rawobj,exposedorder)
if (dispdir eq 'x') then dispindex = double(indexorder(0,*))
if (dispdir eq 'y') then dispindex = double(indexorder(1,*))
disp_porder = 3
waveterms   = poly_fit(dispindex,waveorder,disp_porder)
disporder   = dblarr(n_exposed)
for i=1,disp_porder do disporder = disporder + i*waveterms(i)*(dispindex^(i-1))
disporder   = abs(disporder)
disp        = median(disporder)

; Start by assuming that a standard star has been traced to give
; object-centric slit positions, or that the object trace follows a
; line of constant slit position.  If traceobj=1, the object will
; later be traced in order to produce an object-centric slit position.
slitfit = objposition

; Initialize sky fit results
nosky    = rawobj
skymodel = 0.*rawobj

; Initialize extracted spectrum results
spec    = [0]
specvar = [0]
speclam = [0]

firstsubfit  = 0  ; Flag for the first fit to the ref-sky-subtracted image
fitconverged = 0  ; Flag that current sky fit has converged
skyfitmade   = 0  ; Flag that a sky fit has been made
objfitmade   = 0  ; Flag that object spectrum has been extracted
extractobj   = 0  ; Flag to extract object
skyfactor    = 1. ; Multiplicative constant for computing variance
done         = 0  ; Flag that sky-subtraction/object-extraction loop is done.

; Set initial sky and object fits to zero.
objfitorder     = 0.*rawobjorder
skyfitorder     = 0.*rawobjorder
rawskyfitorder  = 0.*rawobjorder

; Identify bad pixels
obj_badls = where(obj_badorder le 0,n_bad_obj)
ref_badls = where(obj_badorder le 0 or ref_maskorder le 0,n_bad_ref)

; Begin by not masking objects in the sky reference frame (if any)
badorder = obj_badorder
badls    = obj_badls
n_bad    = n_bad_obj

; Begin by fitting raw (un-ref-sky-subtracted) frame
fitraw = 1

; Begin by using a uniform array of evenly spaced break points
refined_bkpoints = 0

; Counter for number of sky-subtraction iterations
niter = 0

; Remember previous sky fit in order to gauge whether fit has converged.
old_skyfitorder = 0.

; Extract object spectrum?
if keyword_set(noobject) then extractobj = 0 else extractobj = 1

; Optional - Do boxcar extraction
if keyword_set(boxcar) then doboxcar=1 else doboxcar = 0

; Create an initial formal variance array
varorder = lr_makevariance(rawobjorder,object=objfitorder,$
                           flatfield=flatorder,$
                           sky=skyfactor*rawobjorder,$  ; Estimate***
                           flatvar=flatvarorder,$
                           darkframe=darkfactor*darkmodel,$
                           readnoise=rdnoiseorder,$
                           objcoadds=skyfactor*ncoadds,$
                           darkcoadds=darkfactor*dkcoadds,$
                           gain=gainorder,$
                           biasnoise=biasfactor*biasnoise_adu)
; *** Note: In previous version, was using sky=skyfactor*rawskyfitorder,
;     which would be zero at this point.
; Boost variance to include as many pixels as possible in the initial fit
; (optional). Needed?
if (flatnormal) then varorder = varorder + satlevel / flatorder

; Identify potential outlying pixels before fitting sky.  These
; pixels will be masked only on the first sky fit iteration.
; Do two filtering passes, offset by half a filter bin length.
; Only points that are rejected in both passes should be filtered.
; This is to avoid cases where the edge of a skyline would be
; filtered out if it occured at the edge of a bin.
print,' '
print,'Pre-filtering sky...'
; Default is not to filter any points
filter_ls   = -1
filter_ls1  = -1
filter_ls2  = -1
; *** Does this filtering scheme work in cases where a reference 
; *** sky frame is being subtracted?  Or if there is no object?
if not(keyword_set(noskysub)) then begin
   filter_sz   = 5*min_bk_space
   minwave     = min(waveorder)
   temp_sky_ls = where(obj_slitorder lt objposition-3*objsig or $
                       obj_slitorder gt objposition+3*objsig,n_sky)
   temp_wave   = waveorder(temp_sky_ls)
   temp_sky    = rawobjorder(temp_sky_ls)
   if (n_sky ne 0) then for pass=1,2 do begin
      if (pass eq 1) then filter_offset = 0. else filter_offset = 0.5*filter_sz
      filt_hist = histogram(temp_wave,min=minwave+filter_offset,$
                            binsize=filter_sz,rev=filt_ri)
      first_filter = 1
      for i=0,n_elements(filt_hist)-1 do begin
         if (filt_ri(i+1) gt filt_ri(i)) then begin
            bin_indices   = filt_ri(filt_ri(i):filt_ri(i+1)-1)
            bin_sky       = temp_sky(bin_indices)
            bin_median    = median(bin_sky)
            bin_sdev      = robust_sigma(bin_sky)
            bin_filter_ls = where(abs(bin_sky-bin_median) ge 3*bin_sdev,$
                                         n_bin_filter)
            if (n_bin_filter gt 0) then begin
               if first_filter then begin
                  temp_filter_ls = temp_sky_ls(bin_indices(bin_filter_ls))
                  first_filter   = 0
               endif else begin
                  temp_filter_ls = [temp_filter_ls,$
                                    temp_sky_ls(bin_indices(bin_filter_ls))]
               endelse
            endif
         endif
      endfor
      case pass of
         1: filter_ls1 = temp_filter_ls
         2: filter_ls2 = temp_filter_ls
      endcase
   endfor
endif

; Filter only points that were flagged on both passes
filterorder = 0*waveorder
if (total(filter_ls1) ne -1 and total(filter_ls2) ne -1) then begin
   ; Use the fact that array indices are integers
   big_filter_ls = [filter_ls1,filter_ls2]
   big_filt_hist = histogram(big_filter_ls,omin=filt_min)
   repeat_ls     = where(big_filt_hist ge 2)
   filter_ls     = repeat_ls + filt_min
   filterorder(filter_ls) = 1
endif else begin
   filter_ls = -1
endelse

; Set maximum number of sky-fitting iterations.
if keyword_set(maxiter) then niter_max = maxiter else niter_max = 20

; Make sure there are sky pixels away from the object that can be fit
if not(keyword_set(noobject)) then begin
   sky_ls = where((obj_slitorder lt objposition - nsig*objsig or $
                   obj_slitorder gt objposition + nsig*objsig) and $
                  obj_badorder  eq 1 and $
                  filterorder   eq 0,n_sky)
endif else begin
   sky_ls = where(obj_badorder  eq 1 and $
                  filterorder   eq 0,n_sky)
endelse
if (n_sky eq 0) then begin
   print,'Error - Object extends across entire slit.  Exiting.'
   return
endif

; Add extra break points around sky lines (optional)
if keyword_set(set_bkpts) then begin
   print,'Setting break points...'
   bkpoints = lr_set_bkpoints(waveorder(sky_ls),rawobjorder(sky_ls),$
                              disp,bk_space,min_bk_space)
endif

; Check object profile and trace after sky fit converges?  Skip check
; if no object is being extracted, or if the object profile and trace
; are completely specified by the user.
check_profile = 1
if keyword_set(noobject) then check_profile = 0
if (keyword_set(force) and $
    keyword_set(fwhm)  and $
    ((keyword_set(slitpos) and not(traceobj)) or $
     (keyword_set(std) and keyword_set(stdoffset)))) then check_profile = 0


while not(done) do begin

   ; Increment count of iterations
   niter = niter + 1
   print,' '
   print,'Skyfit iteration '+strtrim(niter,2)
   print,'Slit polynomial order = '+strtrim(illum_poly,2)

   ; Determine whether fitting the raw or ref-sky-subtracted
   ; frame, and specify the maximum number of iterations
   if fitraw then begin
      print,'Fitting raw frame.'
      skyorder = rawobjorder
      skyfactor   = 1.
      biasfactor  = 1.
      darkfactor  = 1.
   endif else begin
      print,'Fitting differenced frame.'
      skyorder = subobjorder
      ; Increase the expected photon noise under the assumption that
      ; the two frames have roughly equal counts.  Also, double the
      ; number of exposures to account for readnoise.
      skyfactor  = 2.
      ; Ignore errors in the bias, since both frames will have had
      ; the same bias level subtracted.
      biasfactor = 0.
      ; Increase the expected noise in the dark counts.
      darkfactor = 2.
      ; Mask objects in reference array
      badorder   = obj_badorder * ref_maskorder
      badls      = ref_badls
      n_bad      = n_bad_ref
   endelse

   ; Calculate the inverse variance array
   invvarorder = 1./varorder

   ; Mask bad pixels and those with suspicious values.  Do not mask
   ; pixels that have been flagged to be "unmasked" unless they exceed
   ; the saturation level.
   if (niter gt 1) then begin
      ; Include more points on the 2nd pass than on subsequent passes
      if (i eq 2) then nsig_out=5 else nsig_out=3
      outlierls = where(abs(skyorder-skyfitorder) gt nsig_out*sqrt(varorder) $
                        or abs(skyorder) ge satlevel,n_outlier)
   endif else begin
      outlierls = where(abs(skyorder) ge satlevel,n_outlier)
   endelse

   if (n_outlier ne 0) then invvarorder(outlierls) = 0.
   if (n_bad     ne 0) then invvarorder(badls) = 0.

   if (niter eq 1 and total(filter_ls) ne -1) then begin
     print,'Masking filtered points.'
     invvarorder(filter_ls) = 0.
   endif

   ; Mask out the object.  Atmospheric refraction has been taken into
   ; account if a standard star has been traced, or will be when the
   ; object is traced relative to the physical slit position.

   n_objpix = 0
   n_refpix = 0
   if not(keyword_set(noobject)) then begin
      print,'Masking object.'
      ; Use 'object' slit positions
      objls = where(obj_slitorder ge objposition - nsig*objsig and $
                    obj_slitorder le objposition + nsig*objsig,n_objpix)
      if (n_objpix ne 0) then invvarorder(objls) = 0.
      ; Note: Additional objects, including objects in the reference
      ; sky frame, if any, have been masked above.
   endif 

   ; Perform 2D b-spline fit to sky
   ; Spacing of break points set in included file.
   rej      = 5.   ; Number of sigma for pixel rejection
   max_iter = 20.
   ; When fitting sky, use true slit positions
   x_2 = slitorder

   if keyword_set(noskysub) then begin
      print,'Skipping sky subtraction.'
      skyfitordder = 0.*skyorder
      sset         = 0
      bkpoints     = 0
   endif else begin
      if (not(keyword_set(set_bkpts)) and refined_bkpoints eq 0) then begin
         sset=bspline_iterfit(waveorder,skyorder,x2=x_2,npoly=illum_poly,$
                              invvar=invvarorder,bkspace=bk_space,$
                              yfit=skyfitorder,upper=rej,lower=rej,$
                              maxiter=max_iter)
         ; Record initial break points.
         bkpoints  = sset.fullbkpt
      endif else begin
         sset=bspline_iterfit(waveorder,skyorder,x2=x_2,npoly=illum_poly,$
                              invvar=invvarorder,bkpt=bkpoints,$
                              yfit=skyfitorder,upper=rej,lower=rej,$
                              maxiter=max_iter)
      endelse
      skyfitmade = 1
   endelse

   ; If fitting the raw (un-reference-subtracted) frame, record sky
   ; fit to be used in computing the variance and for skyline
   ; identification
   if fitraw then begin
      rawskyfitorder  = skyfitorder
      sset_wave = sset
   endif

   ; Display the sky-subtracted frame
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
   if not(no_show) then begin
      chan,11
      if keyword_set(noskysub) then begin
         disprange = 2*median(noskyorder)
      endif else begin
         disprange = 5*sqrt(median(abs(rawskyfitorder)))
      endelse
      if (instr eq 'nirspec' or instr eq 'n_echelle') then disprange = 80
      if (disprange lt 1) then disprange = median(rawskyfitorder)
      ; Crop display range to fit in screen
      device,get_screen_size=ssz
      xdisp_min   = x_display(0)
      delta_xdisp = (x_display(1) - x_display(0)) < (ssz(0) - 50)
      xdisp_max   = xdisp_min + delta_xdisp
      ydisp_min   = y_display(0)
      delta_ydisp = (y_display(1) - y_display(0)) < (ssz(1) - 50)
      ydisp_max   = ydisp_min + delta_ydisp
      loadct,0,/silent
      doit,nosky(xdisp_min:xdisp_max,ydisp_min:ydisp_max),$
           -disprange,disprange
   endif

   ; Create a new variance array based on the sky fit.
   print,'Recalculating variance based on sky fit.'
   varorder = lr_makevariance(rawobjorder,object=objfitorder,$
                              flatfield=flatorder,$
                              sky=skyfactor*rawskyfitorder,$
                              flatvar=flatvarorder,$
                              darkframe=darkfactor*darkmodel,$
                              readnoise=rdnoiseorder,$
                              objcoadds=skyfactor*ncoadds,$
                              darkcoadds=darkfactor*dkcoadds,$
                              gain=gainorder,$
                              biasnoise=biasfactor*biasnoise_adu)

   ; Re-calculate inverse variance array and mask object (if any)
   invvarorder = 1./varorder
   print, varorder
   if (n_bad    ne 0) then invvarorder(badls) = 0.
   ;if (n_objpix ne 0) then invvarorder(objls) = 0.
   ;if (n_refpix ne 0) then invvarorder(refls) = 0.
   

   ; Has the current fit converged?  
   
   ; For preliminary fits to the raw sky, require that the fit
   ; converge to 5%.  For the final fit, require the fit to converge
   ; to 1%.
   if (fitraw and keyword_set(skyref)) then frac=0.05 else frac=0.01
   diff_skyfitorder = abs(skyfitorder - old_skyfitorder)
   notconverged_ls  = where(diff_skyfitorder ge frac*abs(skyfitorder) and $
                            diff_skyfitorder ge sqrt(varorder),$
                            n_notconverged)
   ;; For omitting bright skyline when reducing FORS2 data
   ;notconverged_ls  = where(diff_skyfitorder ge frac*abs(skyfitorder) and $
   ;                         diff_skyfitorder ge 0.5*sqrt(varorder) and $
   ;                         (waveorder lt 5573 or waveorder gt 5582),$
   ;                         n_notconverged)

   print,'Number of sky fit pixels not yet converged: ',$
      strtrim(n_notconverged,2)
   if (n_notconverged eq 0) then fitconverged = 1 else fitconverged = 0

   ; Kludge for some instruments to keep program from hammering away
   ; at the last few pixels (This intended for full-frame reductions.)
   if (instr eq 'fors2_1' or instr eq 'fors2_2') then begin
      if (n_notconverged lt 1000) then fitconverged = 1 else fitconverged = 0
      if fitconverged then print,'Sky fit converged. (Special case for FORS2)'
   endif
   if (instr eq 'lrisb' or instr eq 'lrisr') then begin
      if (n_notconverged lt 800) then fitconverged = 1 else fitconverged = 0
      if fitconverged then print,'Sky fit converged. (Special case for LRIS)'
   endif

   ; Force a check of the trace if we are running out of iterations
   if (traceobj and not(object_traced) and $
       niter eq (floor(niter_max/2.) > 1)) then begin
      print,'Forcing check of object trace.'
      fitconverged = 1
      check_forced = 1
   endif else begin
      check_forced = 0
   endelse

   ; If the fit has converged, check the object profile and trace
   ; (optional).  Tracing the object with respect to split position
   ; helps account for atmospheric dispersion).  'lr_traceobject'
   ; returns an object-centric slit grid, in which the object is at a
   ; constant position over all wavelengths.  These will be used only
   ; for masking and extracting the object, not for fitting the sky.
   ; Object profile and trace are not updated if an offset from
   ; a standard star trace has been provided.

   ;if (fitconverged and not(keyword_set(noobject)) and $
   ;    not(keyword_set(stdoffset)) and not(keyword_set(force))) then begin

   if (fitconverged and check_profile) then begin
      ; ****** May want to make spectrum extraction at this *******
      ; ****** point a keyword option.                      *******
      ; Get a rough spectrum extraction.  
      spec    = [0]   ; Default is no extracted spectrum
      speclam = [0]
      specvar = [0]
      extract_ls = where(abs(obj_slitorder-objposition) le 5*objsig and $
                         badorder gt 0,n_extract)
      if (n_extract eq 0) then begin
         print,' '
         print,'No pixels within object trace currently available.'
         print,' '
         temp_best_profile = 0
      endif else begin
         ; If a trace has been made,
         ; then use the x-intercept of the trace as an estimate of
         ; the object position.     
         if (total(trace_terms) ne -1) then $
            temp_objposition = trace_terms(0) else $
            temp_objposition = objposition
         temp_profile_params = [objheight,temp_objposition,objsig]
         porder = lr_probability_fn(obj_slitorder,temp_profile_params)
         temp_ls           = extract_ls
         temp_noskyorder   = noskyorder(temp_ls)
         temp_waveorder    = waveorder(temp_ls)
         temp_invvarorder  = invvarorder(temp_ls)
         temp_disporder    = disporder(temp_ls)
         temp_porder       = porder(temp_ls)
         print,'Getting rough extraction...'
         spec = opt_extract(temp_noskyorder,temp_waveorder,$
                            temp_invvarorder,temp_porder,speclam=speclam,$
                            bin=2*disp,maxcount=satlevel,$
                            modelcounts=temp_objfitorder,$
                            specvar=specvar,disp=temp_disporder,clip=5)
         ; Use extracted spectrum to get best Gaussian profile
         temp_best_profile = 1
         if not(no_show) then begin
            chan,10
            loadct,39,/silent
            n_specpix = n_elements(spec)
            buffer = round(n_specpix/20.)
            x_min = min(speclam)
            x_max = max(speclam)
            y_min = -0.5*median(spec) < 0
            y_max = max(spec(buffer:n_specpix-buffer)) > $
                    5*median(sqrt(specvar))
            plot,speclam,spec,yrange=[y_min,y_max],/xstyle
            oplot,speclam,sqrt(specvar),color=250
         endif
      endelse
      ; Refine object profile
      print,' '
      old_objposition = objposition
      old_objsig      = objsig
      param_est       = [objheight,objposition,objsig]
      porder  = lr_gaussprofile(obj_slitorder,noskyorder,waveorder,$
                                invvarorder,profile_params,$
                                waverange=[wavemin,wavemax],$
                                estimates=param_est,noshow=no_show,$
                                scaledcounts=profile_noskyorder,$
                                fitlist=profile_list,$
                                result=fit_result,$
                                bestprofile=temp_best_profile,$
                                spectrum=spec,specwave=speclam,$
                                specerr=sqrt(specvar),$
                                dispersion=disporder) 
      print,'New object profile terms: ',profile_params
      objheight   = profile_params(0)
      objposition = profile_params(1)
      objsig      = abs(profile_params(2))
      ; If the fit looks suspicious, substitue the user-specified
      ; FWHM and position.
      if (fit_result eq 0) then begin
         print,'Profile fit looks suspicious.'
         ; User-specified position (optional)
         if keyword_set(slitpos) then begin
            objposition = 1.*slitpos
            print,'Using specified object position = '+strtrim(slitpos,2)
         endif
         ; User-specified object width (optional)
         if keyword_set(fwhm) then begin
            objsig = fwhm / (2.*sqrt(2*alog(2)))
            print,'Using specified object FWHM = '+strtrim(fwhm,2)
         endif else begin
            objsig = old_objsig
            print,'Keeping previous profile width (sigma) = '+strtrim(objsig,2)
         endelse
         objheight = abs(objheight)
      endif
      if fitraw then objheight = abs(objheight)
      ; Check - did profile change much?
      if (abs(objposition-old_objposition) gt 0.1 or $
          abs(objsig-old_objsig) gt 0.1*old_objsig) then fitconverged=0 $
         else print,'Object profile looks stable.'
      ; If returning the slit position of the object, then store
      ; the current value
      if return_slitpos then slitpos = objposition
      ; Refine object trace
      if traceobj then begin
         old_obj_slitorder = obj_slitorder
         if keyword_set(bright) then obj_type='bright' else begin
            if keyword_set(faint) then obj_type='faint' else $
               obj_type='ordinary'
         endelse
         n_side = 20 < (max(slitorder) - min(slitorder))/3.
         obj_slitgrid = lr_traceobject(nosky,slitgrid,$
                               objposition,dispdir,range=chip_range,$
                               slitrange=[slitmin,slitmax],$
                               tracerange=trace_range,objtype=obj_type,$
                               spectrum=spec,speclam=speclam,specvar=specvar,$
                               wavegrid=wavegrid,pixtr=pixtr,slittr=slittr,$
                               slitfit=slitfit,$;waverange=[wavemin,wavemax],$
                               nside=n_side,noshow=no_show,$
                               bound=chip_bound,$
                               fitterms=trace_terms)
         trace_ok = 1
         if (total(obj_slitgrid) eq -1) then begin
            print,'WARNING: Object trace failed.  Using input slit positions.'
            obj_slitgrid = slitgrid
            trace_ok = 0
         endif
         obj_slitorder = obj_slitgrid(exposedorder)
         object_traced = 1
         ; Check - did trace change much?
         max_trace_change = max(abs(obj_slitorder - old_obj_slitorder))
         print,'Maximum change in slit position = ' + $
               strtrim(max_trace_change,2)
         slit_tolerance   = 0.1  ; Slit position tolerance, in units of
                                 ; the input slit grid (usually pixels).
         if (max_trace_change gt slit_tolerance) then fitconverged = 0 $
            else print,'Object trace looks stable.'
      endif 
      if (fitconverged eq 0) then $
         print,'Redoing fit with new object masking...'
   endif

   ; If a check of the object trace was forced, then resume
   ; iterating sky fit
   if check_forced then begin
      print,'Redoing fit with new object masking after forced check...'
      fitconverged = 0
   endif

   ; If the fit still seems to have converged, check for regions where
   ; the fit is poor and add additional break points in those regions.
   ; Skip this step if no additional break points can be added.
   if fitconverged then check_fit=1 else check_fit=0
   if (keyword_set(set_bkpts) or min_bk_space le bk_space) then check_fit=0

   if check_fit then begin
      ; Mask object and outlying pixels
      if (n_objpix ne 0) then invvarorder(objls) = 0.
      if (n_refpix ne 0) then invvarorder(refls) = 0.
      outlierls = where(abs(skyorder-skyfitorder) gt 5*sqrt(varorder) or $
                        skyorder ge satlevel or skyorder le -satlevel,$
                        n_outlier)
      if (n_outlier ne 0) then invvarorder(outlierls) = 0.
      print,' '
      print,'Checking fit...'
      old_bkpoints     = bkpoints
      n_bkpoints       = n_elements(bkpoints)
      red_chisq        = fltarr(n_bkpoints)
      boosted_varorder = varorder
      ; Speed up search by breaking arrays into regions that contain
      ; several break points
      ; *** Should rewrite this section using HISTOGRAM!! ***
      bkpoints_per_region = 1.*round(sqrt(n_bkpoints))
      num_regions = ceil(1.*n_bkpoints / bkpoints_per_region)
      for j=0,num_regions-1 do begin
         k_lo = j*bkpoints_per_region
         k_hi = ((j+1)*bkpoints_per_region - 1) < (n_bkpoints - 1)
         if (j eq 0) then $
            region_wave_lo = bkpoints(k_lo)-0.5*bk_space else $
            region_wave_lo = 0.5*(bkpoints(k_lo-1)+bkpoints(k_lo))
         if (j eq num_regions-1) then $
            region_wave_hi = bkpoints(k_hi)+0.5*bk_space else $
            region_wave_hi = 0.5*(bkpoints(k_hi)+bkpoints(k_hi+1))
         region = where(waveorder ge region_wave_lo and $
                        waveorder lt region_wave_hi and $
                        invvarorder ne 0,n_region)
         if (n_region gt 0) then begin
            region_wave   = waveorder(region)
            region_sky    = skyorder(region)
            region_skyfit = skyfitorder(region)
            region_var    = varorder(region)
            region_invvar = invvarorder(region)
            for k=k_lo,k_hi do begin
               if (k eq 0) then $
                  wave_lo = bkpoints(k)-0.5*bk_space else $
                  wave_lo = 0.5*(bkpoints(k-1)+bkpoints(k))
               if (k eq n_bkpoints-1) then $
                  wave_hi = bkpoints(k)+0.5*bk_space else $
                  wave_hi = 0.5*(bkpoints(k)+bkpoints(k+1))
               bk_bin = where(region_wave ge wave_lo and $
                              region_wave lt wave_hi,n_bk_bin)
               if (n_bk_bin ne 0) then begin
                  temp_resid   = region_sky(bk_bin) - region_skyfit(bk_bin)
                  temp_var     = region_var(bk_bin)
                  temp_chisq   = total((temp_resid^2)/temp_var)
                  red_chisq(k) = temp_chisq / n_bk_bin
                  ; If a large fraction of the pixels in the current bin
                  ; fall well outside the fit, treat this as a suspicious
                  ; case.  Boost the variance at those pixels, for use
                  ; only if additional breakpoints are added.
                  badfit_ls = where(abs(temp_resid) gt 5*sqrt(temp_var),$
                                    n_badfit)
                  if (n_badfit gt 0.75*n_bk_bin) then begin
                     print,'Supicious fit.  Boosting temporary variance...'
                     boosted_varorder(region(bk_bin)) = $
                                    ((1./3)^2)*median(temp_resid)^2
                     fitconverged = 0
                  endif
               endif
            endfor
         endif
      endfor
      ; Attempt to improve the fit by adding break points in regions
      ; where the standard deviation is significantly greater than
      ; expected.  Skip this step if extra break points were initially
      ; set around sky lines.
      if (not(keyword_set(set_bkpts)) and min_bk_space lt bk_space) then begin
         bkpoints = lr_refine_bkpoints(bkpoints,red_chisq,$
                                       min_bkspace=min_bk_space,$
                                       status=refine_status)
         ; If refinement was needed, set flag and continue to iterate fit,
         ; possibly using boosted variance.
         if (refine_status eq 1) then begin
            refined_bkpoints = 1
            fitconverged     = 0
            varorder = boosted_varorder
         endif
      endif else begin
         print,'Not adding additional break points.'
      endelse
   endif

   ; Quit after too many iterations
   if (niter ge niter_max) then begin
      print,' '
      print,'Too many iterations!'
      fitconverged = 1
   endif

   ; If the current fit has converged and break points do not need to
   ; be refined, then proceed to the ref-subtracted frame or quit.
   if fitconverged then begin
      if (fitraw and keyword_set(skyref)) then begin
         ; Proceed to ref-sky-subtracted
         fitraw           = 0
         fitconverged     = 0
         niter            = 0
         refined_bkpoints = 0
      endif else begin
         print,' '
         print,'Sky subtraction done.'
         print,' '
         done = 1  
      endelse
   endif

   ; Save the current sky fit to compare to the next fit.
   old_skyfitorder = skyfitorder

endwhile

; Re-calculate the inverse variance array.
invvarorder = 1./varorder
; Mask bad pixels
if (n_bad ne 0) then invvarorder(badls) = 0.

; Current profile paramters
profile_params = [objheight,objposition,objsig]

;;;
;;; Wavelength calibrate
;;;

; Perform wavelength fitting to the sky model without the slit
; illimination function.
sset_flat = sset_wave
if (illum_poly le 1) then begin
   wave_x2 = 0
endif else begin
   for i=2,illum_poly do sset_flat.coeff(i-1,*) = 0
   wave_x2 = x_2
endelse
if keyword_set(noskysub) then begin
   waveskyfitorder = 0
endif else begin
   waveskyfitorder = bspline_valu(waveorder,sset_flat,x2=x_2)
endelse

; Default is no skylines identified
xid    = -1
waveid = -1

; Default is to use the input wavelengths
wavecoeffs=[0.,1.]

if keyword_set(waveref) then begin

   ; Use a pre-made wavelength fit (optional)

   if (size(waveref,/tname) eq 'STRING') then begin
      checkexists = findfile(waveref,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Wave reference file not found.'
         return
      endif else begin
          wavereffile = 1
          waverefarr = readfits(waveref,hwaveref)
          sxaddpar,hobj,'WAVEREF',filename(waveref),'Wavelength reference.'
      endelse
   endif else begin
      wavereffile = 0
      waverefarr  = waveref
   endelse
   ; Possible that this is a subarray.  If smaller than the input
   ; arrays, try patching into a full array.
   if (n_elements(waverefarr) lt sz(0)) then begin
      waveref_sz = size(waverefarr)
      if (waveref_sz(1) ne subarray_x(1)-subarray_x(0)+1 or $
          waveref_sz(2) ne subarray_y(1)-subarray_y(0)+1) then begin
         print,'Error - Unexpected size for wavelength reference array.'
         return
      endif
      ; Fill remainder of wavelength array with dummy (high) values
      tmp_waverefarr = 1e9 + 0*rawobj
      tmp_waverefarr(subarray_x(0):subarray_x(1),subarray_y(0):subarray_y(1)) $
         = waverefarr
      waverefarr = tmp_waverefarr
   endif
   wavefitorder = waverefarr(exposedorder)
   ; If heliocentric velocity information is available, then undo the
   ; heliocentric velocity correction applied to the reference
   ; wavelength array, unless the /nohelio keyword is set.
   if (wavereffile) then begin
      vhelio_ref = sxpar(hwaveref,'HELIOVEL',count=heliovelcount)
      if (heliovelcount eq 1) then begin
         if keyword_set(nohelio) then begin
            ; Do not undo heliocentric velocity correction.  Record in
            ; header the heliocentric velocity correction applied to
            ; reference wavelength array.  (Standards may have
            ; coordinates and be taken at times similar to the
            ; science objects serving as wavelength references, therefore
            ; heliocentric velocities may be very similar.)
            sxaddpar,hobj,'HELIOCNT','Heliocentric correction from waveref.'
            sxaddpar,hobj,'HELIOVEL',vhelio_ref,$
                     'lam_hcnt = lam_obs*(1+heliovel/c)'
            print,'Heliocentric velocity borrowed from reference wavelength'
            print,'   array.'
            vhelio = 0.
         endif else begin
            ; Undo heliocentric velocity correction applied to the 
            ; reference wavelength array.
            wavefitorder = wavefitorder / (1. + vhelio_ref/c)
            print,'Heliocentric velocity correction from reference wavelength'
            print,'   array undone.' 
         endelse
      endif else begin
         print,'No heliocentric velocity information available for reference'
         print,'   wavelength array.'
         vhelio_ref = 0.
      endelse
   endif else begin
      print,'No heliocentric velocity information available for reference'
      print,'   wavelength array.'
      vhelio_ref = 0.
   endelse

endif else begin

   ; Perform skyline fit.  The model wavelength fit is different from
   ; the reference wavelength array above.  It is a 3-column ASCII
   ; table or an nx3 array that gives the intensity and wavelength
   ; over a range in wavelength-proxies (the values so far used as
   ; "wavelengths").  This fit is used as an estimate to compute a
   ; wavelength fit for the current array.

   if keyword_set(linelist) then line_file = linelist 

   ; Set initial order of polynomial to be used to fit wavelengths.
   if keyword_set(waveshift) then wc_order = 1 else wc_order = 3

   ;; Calculate the dispersion in units of input wavelength/pixel
   ;; Avoid pre/overscan regions.  Already done above!
   ;; NOTE: Pixels that are at the edge of the illuminated regions or
   ;; between orders may have the wrong dispersion calculated.
   ;print,'Calculating input dispersion...'
   ;disparr = dblarr(xsz,ysz)
   ;for thischip=0,nchip-1 do begin
   ;   xlo = chip_xmin(thischip) > 0
   ;   xhi = chip_xmax(thischip) < (xsz-1)
   ;   ylo = chip_ymin(thischip) > 0
   ;   yhi = chip_ymax(thischip) < (ysz-1)
   ;   if (dispdir eq 'x') then begin
   ;      for x=xlo,xhi-1 do begin
   ;         disparr(x,ylo:yhi) = abs(wavegrid(x+1,ylo:yhi) - $
   ;                                  wavegrid(x,ylo:yhi))
   ;      endfor
   ;      disparr(xhi,ylo:yhi) = disparr(xhi-1,ylo:yhi)
   ;   endif else begin
   ;      for y=ylo,yhi-1 do begin
   ;         disparr(xlo:xhi,y) = abs(wavegrid(xlo:xhi,y+1) - $
   ;                                  wavegrid(xlo:xhi,y))
   ;      endfor
   ;      disparr(xlo:xhi,yhi) = disparr(xlo:xhi,yhi-1)
   ;   endelse
   ;endfor
   ;disporder = disparr(exposedorder)
   ;; Median dispersion.  Note: If multiple orders are present on the
   ;; same ccd, then the dispersion as calculated above will be
   ;; spurious at the boundary between orders.  Therefore, take the
   ;; median instead of the mean.
   ;; NOTE: This is not the dispersion as wave_calib defines it.
   ;disp = median(disporder)

   ; If the dispersion at this point is significantly different from
   ; 1, then assume that the input 'wavearr' is already in physical
   ; wavelength units.  In that case, use fake dispersion and
   ; resolution inputs to get wave_calib to look for lines over the
   ; right range.
   if (disp lt 0.7 or disp gt 1.5) then begin
      wc_disp = 1
      wc_res  = 3*disp
   endif else begin
      wc_disp = 0
      wc_res  = 0
   endelse

   ; Use a fake resolution so that wave_calib looks for and fits
   ; lines over the right number of pixels.
   fake_resoln = 3*disp

   ; Automatic line identification?
   if (keyword_set(wavemodel) or keyword_set(autoid)) then begin
      auto_id   = 1
      wc_noshow = 1
   endif else begin
      auto_id   = 0
      wc_noshow = 0
   endelse

   ; Reject lines that appear to be noise (should be done only for
   ; spectra with well separated lines).
   if (instr eq 'hires' or instr eq 'mike' or $
       instr eq 'n_echelle' or instr eq 'mage') then $
       reject_noise = 1 else reject_noise = 0

   ; Use a model sky fit?
   if keyword_set(wavemodel) then wave_model = wavemodel $
      else wave_model = 0

   print,'Identifying skylines...'

   wave_calib,waveorder,waveskyfitorder,wavefitorder,reference=wave_model,$
           linefile=line_file,coeffs=temp_wavecoeffs,/nopersist,$
           xid=xid,waveid=waveid,order=wc_order,$
           resolution=wc_res,auto=auto_id,dispersion=wc_disp,$
           xoffset=xoff,rejectnoise=reject_noise,noshow=wc_noshow

   if (total(temp_wavecoeffs) ne 0) then wavecoeffs = temp_wavecoeffs
   print,'Wavelength fit coefficients:',wavecoeffs

   ; If the input wavearr was already calibrated, then reject large
   ; shifts.
   ;if (keyword_set(autoid) and n_elements(wavecoeffs) gt 1) then begin 
   ;   if (wavecoeffs(1) lt 0.5 or wavecoeffs(1) gt 1.5) then begin
   ;      print,'Suspicious wavelength fit - rejecting.'
   ;      xid    = -1
   ;      waveid = -1
   ;   endif
   ;endif
   if (keyword_set(autoid) and n_elements(uniq(xid)) gt 1) then begin 
      if (median(abs(wavefitorder-waveorder)) gt 12.0) then begin
         print,'Suspicious wavelength fit - rejecting.'
         xid    = -1
         waveid = -1
      endif
   endif

   ; If one and only one line is found, apply a simple bulk offset (in
   ; pixels) to the input wavelength array.  If no skylines are found,
   ; apply the user-supplied bulk shift (in wavelength units), if any.
   ; Otherwise, use the input wavelengths.
   if (n_elements(uniq(xid)) eq 1) then begin
      if (total(xid) eq -1) then begin
         print,'No skylines found.'
         if keyword_set(meanshift) then begin
            bulkoffset = meanshift
            print,'Applying user-supplied bulk offset of '+strtrim(bulkoffset,2)
            wavefitorder = waveorder + bulkoffset
         endif else begin
            print,'Using input wavelengths'
            wavefitorder = waveorder
            bulkoffset = 0.
         endelse
      endif else begin
         ; Print calculate offset in pixels, and apply to whole
         ; spectrum according to dispersion.
         print,'Only one skyline found.'
         wav_offset = waveid(0) - xid(0)
         local_disp = interpol(disporder,waveorder,xid(0))
         pix_offset = wav_offset / local_disp
         print,'Applying wavelength offset of '+strtrim(pix_offset,2)+' pixels'
         wavefitorder = waveorder + pix_offset*disporder
         ; Rough offset, in wavelength units
         bulkoffset = median(pix_offset*disporder)
      endelse
   endif else bulkoffset = 0

   ; Record the mean wavelength shift
   meanshift = mean(wavefitorder - waveorder)

endelse

; Apply heliocentric velocity correction
if not(keyword_set(nohelio)) then begin
   ; Get UT modified Julian date and convert to UT Julian date.
   ; Note: could also read UT date and time and convert to Julian
   ; date using JULDAY.
   mjday = sxpar(hobj,'MJD-OBS',count=mjd_count)
   if (mjd_count eq 0) then mjday = sxpar(hobj,'MJD',count=mjd_count) 
   mjday = double(mjday)
   jday  = mjday + 2400000.5
   if (mjd_count eq 0) then begin
      ; Read in UT date and time and convert to Julian date.  
      utdate = sxpar(hobj,'UT-DATE',count=utdate_count) ; Magellan
      if (utdate_count eq 0) then begin
         utdate = sxpar(hobj,'DATE',count=utdate_count) ; Palomar
      endif
      uttime = sxpar(hobj,'UT-TIME',count=uttime_count) ; Magellan
      if (uttime_count eq 0) then begin
         uttime = sxpar(hobj,'UT',count=uttime_count) ; Palomar
      endif
      if (utdate_count ne 0) then begin
         ; Assume format is YYYY-MM-DD'
         utdate = strtrim(utdate,2)
         year   = long(strmid(utdate,0,4))
         month  = long(strmid(utdate,5,2))
         day    = long(strmid(utdate,8,2))
         if (uttime_count ne 0) then begin
            ; Assume format is HH:MM:SS
            uttime = strtrim(uttime,2)
            hour   = float(strmid(uttime,0,2))
            minute = float(strmid(uttime,3,2))
            second = float(strmid(uttime,6,2))
         endif else begin
            print,'WARNING: UT time unknown.  Using 00:00:00.'
            hour = 0  & minute = 0  &  second = 0
         endelse
         frac_hour = hour + minute/60. + second/3600.
         jdcnv,year,month,day,frac_hour,jday
         print,'Calculated Julian day = '+strtrim(jday,2)
         mjd_count = 1
      endif
   endif
   ; Get telescope coordinates and epoch
   objra  = sxpar(hobj,'RA',count=ra_count)   ; Degrees, for NIRSPEC
   objdec = sxpar(hobj,'DEC',count=dec_count) ; Degrees, for NIRPSEC
   if (instr eq 'mage') then begin
      objra  = sxpar(hobj,'RA-D',count=ra_count)   ; Degrees, for MagE
      objdec = sxpar(hobj,'DEC-D',count=dec_count) ; Degrees, for MagE      
   endif
;   if (size(objra,/tname) eq 'STRING') then begin
;      decimal_coords = radecdecimal(objra,objdec)
;      objra  = decimal_coords(0)
;      objdec = decimal_coords(1)
;   endif
   if (instr eq 'luci') then begin
       objra  = sxpar(hobj,'OBJRA',count=ra_count)
       objdec = sxpar(hobj,'OBJDEC',count=dec_count)
       objra  = hms2dec(objra)*15.
       objdec = hms2dec(objdec)*1.
   endif

 
   epoch  = float(sxpar(hobj,'EQUINOX',count=eq_count))
   if (eq_count eq 0) then begin
      print,'Assuming epoch = J2000'
      epoch = 2000.
   endif
   if (mjd_count * ra_count * dec_count ne 0) then begin
      ; Get heliocentric velocity correction
      vcorr = heliocentric(objra,objdec,epoch,jd=jday,long=360.-longitude,$
                           lat=latitude,altitude=altitude)
      ; Change sign such that the velocity is that towards the object
      ; with respect to the sun
      vhelio = -1*vcorr
      ; Adjust wavelengths
      wavefitorder = wavefitorder * (1. + vhelio/c)
      ; Add keywords similar to those added by MAKEE.
      sxaddpar,hobj,'HELIOCNT','Heliocentric correction applied.'
      sxaddpar,hobj,'HELIOVEL',vhelio,'lam_hcnt = lam_obs*(1+heliovel/c)'
      print,'Heliocentric velocity correction applied = '$
            +strtrim(vhelio,2)+' km/s.'
   endif else begin
      print,'Could not get enough information from header to '+$
            'compute heliocentric velocity.'
      vhelio = 0.
   endelse
endif else vhelio = 0.

; Create 2D wavelength array
if keyword_set(waveref) then begin
   wavefit = waverefarr * (1. - vhelio_ref/c)
endif else begin
   ; Perform a bulk wavelength shift only?
   if (bulkoffset ne 0) then begin
      wavefit = wavegrid + bulkoffset
   endif else begin
      ; New wavelength solution based on skylines?
      if (total(xid) ne -1) then begin
         wavefit = 0.*wavegrid
         for i=0,n_elements(wavecoeffs)-1 do $
            wavefit = wavefit + wavecoeffs(i) * wavegrid^i
      endif else begin
         wavefit = wavegrid
      endelse
   endelse
endelse
; Heliocentric correction, if any
if (vhelio ne 0) then wavefit = wavefit * (1. + vhelio/c)

; Output wavelength calibration information
waveinfo = [wavecoeffs,vhelio]

; Shift to vacuum wavelengths?
if keyword_set(vacuum) then begin
   print,'Converting wavelengths to vacuum.'
   airtovac,wavefitorder
   airtovac,wavefit
   sxaddpar,hobj,'VACUUM','Wavelengths converted to vacuum.'
endif

; Calculate the dispersion in units of fitted wavelength/pixel
; Avoid pre/overscan regions.
; NOTE: Pixels that are at the edge of the illuminated regions or
; between orders may have the wrong dispersion calculated.
print,'Calculating dispersion...'
wavefitterms = poly_fit(dispindex,wavefitorder,disp_porder)
dispfitorder = dblarr(n_exposed)
for i=1,disp_porder do $
   dispfitorder = dispfitorder + i*wavefitterms(i)*(dispindex^(i-1))
dispfitorder = abs(dispfitorder)

; NOTE: This method gets the dispersion wrong at the 
;       edges, but accounted for multiple CCDs.  May have to
;       modify the above approach to work with LRIS-B
;dispfitarr = dblarr(xsz,ysz)
;for thischip=0,nchip-1 do begin
;   xlo = chip_xmin(thischip) > 0
;   xhi = chip_xmax(thischip) < (xsz-1)
;   ylo = chip_ymin(thischip) > 0
;   yhi = chip_ymax(thischip) < (ysz-1)
;   if (dispdir eq 'x') then begin
;      for x=xlo,xhi-1 do begin
;         dispfitarr(x,ylo:yhi) = abs(wavefit(x+1,ylo:yhi) - $
;                                     wavefit(x,ylo:yhi))
;      endfor
;      dispfitarr(xhi,ylo:yhi) = dispfitarr(xhi-1,ylo:yhi)
;   endif else begin
;      for y=ylo,yhi-1 do begin
;         dispfitarr(xlo:xhi,y) = abs(wavefit(xlo:xhi,y+1) - $
;                                     wavefit(xlo:xhi,y))
;      endfor
;      dispfitarr(xlo:xhi,yhi) = dispfitarr(xlo:xhi,yhi-1)
;   endelse
;endfor
;dispfitorder = dispfitarr(exposedorder)

;;;
;;; Extract object spectrum
;;;

; Make sure the object trace is covered
if not(keyword_set(noobject)) then begin
   extract_ls = where(abs(obj_slitorder-objposition) le 5*objsig and $
                      badorder gt 0,n_extract)
   if (n_extract eq 0) then begin
      print,' '
      print,'No pixels within object trace available.'
      print,' '
   endif
endif else n_extract = 0

; Adjust the wavelength interval over which to locate the object to
; account for wavelength calibration.
wavemin = interpol(wavefitorder,waveorder,wavemin)
wavemax = interpol(wavefitorder,waveorder,wavemax)

; User-specified binning (optional)
if keyword_set(bin)    then wavebin = bin
if keyword_set(loglin) then log_lin = 1 else log_lin = 0

if (n_extract ne 0 and not(keyword_set(noobject))) then begin
    
   ; Get a rough idea of the spectrum.  If using a spline profile
   ; ultimately, use boxcar extraction (the object should be bright
   ; enough to not be overly affected by cosmic rays).  If extracting
   ; using a Gaussian profile, use the current profile.
   if (keyword_set(sprofile) or keyword_set(bright) or $
       keyword_set(boxcar)) then begin
      boxcar_ls        = extract_ls
      n_boxcar         = n_elements(boxcar_ls)
      box_noskyorder   = noskyorder(boxcar_ls)
      box_wavefitorder = wavefitorder(boxcar_ls)
      box_invvarorder  = invvarorder(boxcar_ls)
      box_porder       = 1 + fltarr(n_boxcar)
      box_dispfitorder = dispfitorder(boxcar_ls)
      print,'Running boxcar extraction...'
      spec = opt_extract(box_noskyorder,box_wavefitorder,box_invvarorder,$
                         box_porder,speclam=speclam,bin=wavebin,loglin=log_lin,$
                         maxcount=satlevel,modelcounts=box_objfitorder,$
                         specvar=specvar,boxcar=1,disp=box_dispfitorder,clip=5)
   endif else begin
      porder = lr_probability_fn(obj_slitorder,profile_params)
      temp_ls           = extract_ls
      temp_noskyorder   = noskyorder(temp_ls)
      temp_wavefitorder = wavefitorder(temp_ls)
      temp_invvarorder  = invvarorder(temp_ls)
      temp_dispfitorder = dispfitorder(temp_ls)
      temp_porder       = porder(temp_ls)
      print,'Running initial extraction...'
      spec = opt_extract(temp_noskyorder,temp_wavefitorder,temp_invvarorder,$
                         temp_porder,speclam=speclam,bin=wavebin,$
                         loglin=log_lin,maxcount=satlevel,$
                         modelcounts=temp_objfitorder,specvar=specvar,$
                         disp=temp_dispfitorder,clip=5)
   endelse
   ;; Adjust temporary spectrum error estimates assuming Poisson stats
   ;specvar = specvar + abs(spec)/median(gainorder)

   ; If doing optimal extraction, refine the object profile and get
   ; probability distribution function.
   if not(keyword_set(boxcar)) then begin
      ; Compute profile over the entire spectrum?  Default is to model
      ; profile over a limited range in wavelength.
      best_profile = 0
      p_wavemin = wavemin
      p_wavemax = wavemax
      if (bestprofile eq 1) then begin
         best_profile = 1 
         ; Compute the min and max slit position covered in each
         ; wavelength bin.
         spec_slitmin = 0*spec
         spec_slitmax = 0*spec
         ; Convert wavelengths to integers corresponding to the
         ; index in the spectrum in which the wavelength would fall.
         if log_lin then begin
            w0 = alog10(speclam(0))
            dw = alog10(wavebin/c + 1)
            waveindex = (alog10(wavefitorder) - w0) / dw
         endif else begin
            w0 = speclam(0)
            dw = wavebin
            waveindex = (wavefitorder - w0) / dw
         endelse
         wave_hist = histogram(waveindex,min=-0.5,nbins=n_elements(spec),$
                               rev=wave_ri)
         for j=0,n_elements(wave_hist)-1 do begin
            if (wave_ri(j+1) gt wave_ri(j)) then begin
               temp_slit = obj_slitorder(wave_ri(wave_ri(j):wave_ri(j+1)-1))
               spec_slitmin(j) = min(temp_slit)
               spec_slitmax(j) = max(temp_slit)
            endif
         endfor
         ; Fit profile over wavelength range where the full 
         ; profile is sampled.
         fullprofile_ls = where(spec_slitmin le objposition-3*objsig and $
                                spec_slitmax ge objposition+3*objsig,$
                                n_fullprofile)
         if (n_fullprofile ne 0) then begin
           best_profile = 1
           p_wavemin = min(speclam(fullprofile_ls))
           p_wavemax = max(speclam(fullprofile_ls))
         endif
      endif
      if (keyword_set(bright) or keyword_set(sprofile)) then begin
         ; Spline profile
         print,'Using spline object profile..'
         ; For a bright object, allow profile to vary in the
         ; wavelength direction
         ;if keyword_set(bright) then variable=1 else variable=0
         ; Temporary!
         variable = 0
         porder = lr_splineprofile(obj_slitorder,noskyorder,wavefitorder,$
                                   invvarorder,$
                                   waverange=[p_wavemin,p_wavemax],$
                                   bestprofile=best_profile,spectrum=spec,$
                                   specwave=speclam,specerr=sqrt(specvar),$
                                   dispersion=dispfitorder,$
                                   scaledcounts=profile_noskyorder,$
                                   fitlist=profile_list,$
                                   bad=badorder,variable=variable,$
                                   noshow=no_show)
      endif else begin
         ; Attempt to fit a Gaussian profile to the object.
         param_est = [objheight,objposition,objsig]
         temp_porder =  lr_gaussprofile(obj_slitorder,noskyorder,$
                                     wavefitorder,$
                                     invvarorder,temp_params,$
                                     waverange=[p_wavemin,p_wavemax],$
                                     bestprofile=best_profile,spectrum=spec,$
                                     specwave=speclam,$
                                     specerr=sqrt(specvar),$
                                     dispersion=dispfitorder,$
                                     estimates=param_est,$
                                     scaledcounts=profile_noskyorder,$
                                     fitlist=profile_list,$
                                     result=fit_result,noshow=no_show)
         ; If the fit looks suspicious, use the user-provided profile
         ; parameters (if available).
         if (fit_result eq 0) then begin 
            print,'Profile looks suspicious - rejecting.'
            profile_params = param_est
         endif else begin
            objheight = abs(temp_params(0))
            objposition = temp_params(1)
            objsig = abs(temp_params(2))
            profile_params = [objheight,objposition,objsig]
         endelse
         ; Force the use of input profile parameters?
         if keyword_set(force) then begin
            print,'Using input profile parameters'
            if keyword_set(fwhm) then objsig = fwhm / (2*sqrt(2*alog(2)))
            if keyword_set(slitpos) then objposition = slitpos
            ; Override slitpos kewyord with stdoffset keyword
            if keyword_set(stdoffset) then $
               objposition = stdposition + stdoffset
            profile_params = [objheight,objposition,objsig]  
         endif
         ; Return slit position of object?
         if return_slitpos then slitpos = objposition
         print,'Final Gaussian profile fit terms: ',profile_params
         porder = lr_probability_fn(obj_slitorder,profile_params)
      endelse
      ; First pass optimal extraction.  Avoid potentially
      ; non-Gaussian wings.
      nsig_ex=3 
      optimal_ls = where(abs(obj_slitorder-objposition) le nsig_ex*objsig and $
                         porder ne 0 and badorder gt 0,n_optimal)
      if (n_optimal eq 0) then begin
         print,'LONGSLIT_REDUCE: Error - No points within range for extraction.'
         return
      endif
      opt_noskyorder   = noskyorder(optimal_ls)
      opt_wavefitorder = wavefitorder(optimal_ls)
      opt_invvarorder  = invvarorder(optimal_ls)
      opt_porder       = porder(optimal_ls)
      opt_dispfitorder = dispfitorder(optimal_ls)
      print,'Running optimal extraction...'
      spec = opt_extract(opt_noskyorder,opt_wavefitorder,opt_invvarorder,$
                         opt_porder,speclam=speclam,bin=wavebin,loglin=log_lin,$
                         maxcount=satlevel,modelcounts=opt_objfitorder,$
                         specvar=specvar,disp=opt_dispfitorder,clip=5)
      objfitorder = 0.*wavefitorder
      objfitorder(optimal_ls) = opt_objfitorder
      ; Recalculate expected variance 
      print,'Recalculating variance based on expected object counts.'
      varorder = lr_makevariance(rawobjorder,object=objfitorder,$
                                 flatfield=flatorder,$
                                 sky=skyfactor*rawskyfitorder,$
                                 flatvar=flatvarorder,$
                                 darkframe=darkfactor*darkmodel,$
                                 readnoise=rdnoiseorder,$
                                 objcoadds=skyfactor*ncoadds,$
                                 darkcoadds=darkfactor*dkcoadds,$
                                 gain=gainorder,$
                                 biasnoise=biasfactor*biasnoise_adu)
      invvarorder = 1./varorder
      if (n_bad ne 0) then invvarorder(badls) = 0.      
      opt_invvarorder = invvarorder(optimal_ls)
      ; Second pass optimal extraction
      print,'Re-running optimal extraction...'
      spec = opt_extract(opt_noskyorder,opt_wavefitorder,opt_invvarorder,$
                         opt_porder,speclam=speclam,bin=wavebin,loglin=log_lin,$
                         maxcount=satlevel,modelcounts=opt_objfitorder,$
                         specvar=specvar,disp=opt_dispfitorder,$
                         clip=5)
      objfitorder = 0.*wavefitorder
      objfitorder(optimal_ls) = opt_objfitorder
      ; Calculate final variance
      print,'Calculating final variance.'
      varorder = lr_makevariance(rawobjorder,object=objfitorder,$
                                 flatfield=flatorder,$
                                 sky=skyfactor*rawskyfitorder,$
                                 flatvar=flatvarorder,$
                                 darkframe=darkfactor*darkmodel,$
                                 readnoise=rdnoiseorder,$
                                 objcoadds=skyfactor*ncoadds,$
                                 darkcoadds=darkfactor*dkcoadds,$
                                 gain=gainorder,$
                                 biasnoise=biasfactor*biasnoise_adu)
      invvarorder = 1./varorder
      if (n_bad ne 0) then invvarorder(badls) = 0.  
   endif   

   ; Display extracted spectrum (optional)
   if not(no_show) then begin
      chan,10
      loadct,39,/silent
      n_specpix = n_elements(spec)
      buffer = round(n_specpix/20.)
      x_min = min(speclam)
      x_max = max(speclam)
      y_min = -0.5*median(spec) < 0
      y_max = max(spec(buffer:n_specpix-buffer)) > 5*median(sqrt(specvar))
      plot,speclam,spec,yrange=[y_min,y_max],/xstyle
      oplot,speclam,sqrt(specvar),color=220
      if fitraw then whichone = 'No reference sky subtraction.' $
         else whichone = 'Reference sky subtracted.'
      xyouts,0.9*x_max,0.9*y_max,whichone,/data,align=1.0
      xyouts,0.9*x_max,0.7*y_max,'illum_poly = '+strtrim(illum_poly,2),$
         /data,align=1.0
      xyouts,0.9*x_max,0.8*y_max,'niter      = '+strtrim(niter,2),$
         /data,align=1.0
   endif

   objfitmade = 1

endif

;;;
;;; Write outputs.
;;;

; NOTE: 1D spectrum and error are converted to single-precision floating
;       point.

; Plot in color
loadct,39,/silent

; Optional - Extract subarrays for output.
if getsubarray then begin
   ex_xmin = subarray_x(0)
   ex_xmax = subarray_x(1)
   ex_ymin = subarray_y(0)
   ex_ymax = subarray_y(1)
endif

if keyword_set(outfile) then rootname = outfile else rootname = 'longslit'

; 2D Subtracted frame
print,'Writing subtracetd frame...'
if getsubarray then begin
   hextract,nosky,hobj,smnosky,smhobj,ex_xmin,ex_xmax,ex_ymin,ex_ymax
   writefits,rootname+'_sub.fits',smnosky,smhobj
endif else writefits,rootname+'_sub.fits',nosky,hobj
; 2D Sky model
print,'Writing sky model frame...'
if getsubarray then begin
   hextract,skymodel,hobj,smskymodel,smhobj,ex_xmin,ex_xmax,ex_ymin,ex_ymax
   writefits,rootname+'_sky.fits',smskymodel,smhobj
endif else writefits,rootname+'_sky.fits',skymodel,hobj
; 2D Variance array
print,'Writing variance frame...'
noskyvar = 0.*rawobj
noskyvar(exposedorder) = varorder
if getsubarray then begin
   hextract,noskyvar,hobj,smnoskyvar,smhobj,ex_xmin,ex_xmax,ex_ymin,ex_ymax
   writefits,rootname+'_var.fits',smnoskyvar,smhobj
endif else writefits,rootname+'_var.fits',noskyvar,hobj
; 2D Wavelength array
print,'Writing wavelength frame...'
if getsubarray then begin
   hextract,wavefit,hobj,smwavefit,smhobj,ex_xmin,ex_xmax,ex_ymin,ex_ymax
   writefits,rootname+'_wav.fits',smwavefit,smhobj
endif else writefits,rootname+'_wav.fits',wavefit,hobj
; ASCII file containing wavelength fit (optional)
if keyword_set(writewave) then begin
   ; Create a list of wavelength-proxy values, sampling at least
   ; 2 times per pixel.
   ; Kludge for HIRES
   if (instr eq 'hires') then wavefitgridstep = 0.02 else $
      wavefitgridstep    = (1/2.) * min(dispfitorder)
   nwavesteps         = ceil((max(wavefitorder) - min(wavefitorder)) / $
                             wavefitgridstep)
   wavegridstep       = (max(wavefitorder) - min(wavefitorder)) / nwavesteps
   tmp_waveorder      = min(waveorder) + wavegridstep*findgen(nwavesteps)
   tmp_wavefitorder   = interpol(wavefitorder,waveorder,tmp_waveorder)
   tmp_waveskyfitorder = interpol(waveskyfitorder,waveorder,tmp_waveorder)
   print,'Writing wavelength fit to ASCII file...'
   openw,wav_unit,rootname+'_wav.dat',/get_lun
   for i=0,nwavesteps-1 do $
      printf,wav_unit,tmp_waveorder(i),tmp_waveskyfitorder(i),$
                      tmp_wavefitorder(i)
   free_lun,wav_unit
endif
; PS plot of wavelength fit
if (n_elements(uniq(xid)) gt 1) then begin   
   ; Reconstruct wavelength fit on coarsely spaced grid
   tmp_nbin    = 1000.
   tmp_binsz   = (max(waveorder)-min(waveorder)) / tmp_nbin
   tmp_inwave  = min(waveorder) + tmp_binsz*dindgen(tmp_nbin)
   tmp_outwave = dblarr(tmp_nbin)
   tmp_nterms  = n_elements(wavecoeffs)
   for k=0,tmp_nterms-1 do $
      tmp_outwave = tmp_outwave + wavecoeffs(k)*(tmp_inwave^k)
   ps_openfile,rootname+'_lineid',/color
   !p.multi = [0,1,2]
   ; Plot fit
   tmp_xr = [min(tmp_inwave),max(tmp_inwave)]
   tmp_yr = [min(waveid) < min(tmp_outwave),max(waveid) > max(tmp_outwave)]
   plot,xid,waveid,xr=tmp_xr,yr=tmp_yr,xs=3,ys=3,psym=1,$
        title='Wavelength Fit (poly order = '+strtrim(tmp_nterms-1,2)+')',$
        xtitle='Input Wavelength',ytitle='Output Wavlength'
   oplot,tmp_inwave,tmp_outwave
   ; Plot residuals
   tmp_fitwave = 0.*waveid
   for k=0,tmp_nterms-1 do $
      tmp_fitwave = tmp_fitwave + wavecoeffs(k)*(xid^k)
   tmp_resid = waveid - tmp_fitwave
   plot,xid,tmp_resid,psym=1,xr=tmp_xr,xs=3,ys=3,title='Residuals',$
        xtitle='Input Wavelength',ytitle='Fit residuals'
   oplot,[-1e6,1e6],[0,0],linestyle=1
   !p.multi=[0,1,1]
   ps_closefile
endif
; Outputs for extracted spectrum, if any
if (total(spec) ne 0 and not(keyword_set(boxcar) and keyword_set(force))) $
    then begin
   ; 1D Extracted spectrum - NAXIS keywords updated automatically.
   hspec = hobj
   if keyword_set(loglin) then begin
      sxaddpar,hspec,'CRVAL1',alog10(min(speclam))
      sxaddpar,hspec,'CDELT1',alog10(wavebin/c + 1)
      sxaddpar,hspec,'CD1_1',alog10(wavebin/c + 1)
      sxaddpar,hspec,'CRPIX1',1
      sxaddpar,hspec,'CTYPE1','LINEAR'    ; MAKEE convention
      sxaddpar,hspec,'CTYPE2','PIXEL'
   endif else begin
      sxaddpar,hspec,'CRVAL1',min(speclam)
      sxaddpar,hspec,'CDELT1',wavebin
      sxaddpar,hspec,'CD1_1',wavebin
      sxaddpar,hspec,'CRPIX1',1
      sxaddpar,hspec,'CTYPE1','LAMBDA'
      sxaddpar,hspec,'CTYPE2','PIXEL'
   endelse
   spec = float(spec)
   writefits,rootname+'_spec.fits',spec,hspec
   ; 1D Error array
   specerr = sqrt(specvar)
   specerr = float(specerr)
   writefits,rootname+'_err.fits',specerr,hspec
   ; Trace and profile information
   openw,tr_unit,rootname+'_trace.dat',/get_lun
   printf,tr_unit,'# Trace fit coeffients'
   for i=0,n_elements(trace_terms)-1 do printf,tr_unit,trace_terms(i)
   printf,tr_unit,'# Object centroid'
   printf,tr_unit,objposition
   free_lun,tr_unit
   ; PS Plot - 1D Spectrum + Error
   ps_openfile,rootname+'_spec',/color
   n_specpix = n_elements(spec)
   buffer = round(n_specpix/20.)
   y_min = -0.5*median(spec) < 0
   y_max = max(spec(buffer:n_specpix-buffer-1)) > 5*median(sqrt(specvar))
   plot,speclam,spec,yrange=[y_min,y_max],xstyle=1,$
        title=rootname,xtitle='wavelength (angstroms)',ytitle='flux'
   oplot,speclam,sqrt(specvar),color=220
   ps_closefile
   ; PS Plot - Object Profile
   ; If not finding the bestprofile and not using a spline profile,
   ; then reconstruct the non-normalized Gaussian fit.
   if (not(bestprofile eq 1 or keyword_set(sprofile)) or $
       keyword_set(boxcar)) then begin
      temp_z      = (obj_slitorder-profile_params(1))/profile_params(2)
      porder_plot = profile_params(0) * exp(-(temp_z^2)/2.)
   endif else porder_plot = porder
   x_min   = min(obj_slitorder)
   x_max   = max(obj_slitorder)
   y_min   = -1.5*max(porder_plot(extract_ls))
   y_max   = 3*max(porder_plot(extract_ls))
   y_range = [y_min,y_max]
   ; Only plot up to a maximum of 500 points per unit slit length
   delta_slit   = max(obj_slitorder) - min(obj_slitorder)
   max_npoints  = n_elements(profile_list)
   if (max_npoints eq 0) then begin
      max_npoints  = n_elements(obj_slitorder)
      profile_list = lindgen(max_npoints)
   endif
   npoints = (500*delta_slit) < max_npoints
   if (npoints ne 0) then begin
      ps_openfile,rootname+'_profile',/color
      plotls       = ((max_npoints-1.)/npoints)*findgen(npoints)
      slit_plot    = obj_slitorder(profile_list(plotls))
      nosky_plot   = profile_noskyorder(profile_list(plotls))
      profile_plot = porder_plot(profile_list(plotls))
      tmp_order    = sort(slit_plot)
      plot,slit_plot(tmp_order),nosky_plot(tmp_order),psym=3,xstyle=1,$
         yrange=y_range,ystyle=2,title='Object Profile - '+rootname,$
         xtitle = 'slit position'
      oplot,slit_plot(tmp_order),profile_plot(tmp_order),color=250,thick=5
      if not(keyword_set(sprofile)) then begin
         objfwhm = 2.*sqrt(2.*alog(2))*objsig
         xyouts,x_min+0.1*(x_max-x_min),0.95*y_max,$
            "Centroid = "+strtrim(string(objposition),2),/data
         xyouts,x_min+0.1*(x_max-x_min),0.85*y_max,$
            "FWHM = "+strtrim(string(objfwhm),2),/data
      endif
      ps_closefile
   endif else begin
      print,'Unable to plot profile.'
   endelse
   ; PS Plot of object trace (optional)
   if traceobj and trace_ok then begin
      ps_openfile,rootname+'_trace',/color
      plot,pixtr,slittr,psym=3,/ystyle,title='Object Trace - '+rootname,$
         xtitle='pixel',ytitle='slit position',/xs
      oplot,pixtr,slitfit,psym=3
      
      ps_closefile
   endif
endif

; Set object profile keyword return values.
fwhm = 2.*sqrt(2*alog(2))*objsig
if keyword_set(stdtrace) then begin
   stdoffset = objposition - stdposition
endif else stdoffset = 0.

; Header keyword output
hdr = hobj

; Restore system variables
!p.multi=currentmulti

return
end

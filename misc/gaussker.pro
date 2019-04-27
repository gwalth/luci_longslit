Function gaussker,sig,fwhm=fwhm,clip=clip,size=size,$
            center=center,peak=peak
;+----------------------------------------------------------------
;
; GAUSSKER	   3/2002 GDB
;
; Returns a 1-D unit area Gaussian kernel.  The value in each bin is
; the integral of a Gaussian over the pixel, with the entire array
; normalized to unity (unless the peak keyword is set).
; 
; INPUTS
;	sig	- Sigma of distribution (pixels)
;
; KEYWORDS
;	fwhm	- FWHM of distribution (pixels) (overides sig 
;		  parameter).
;	clip	- Min value of Gaussian kernel.  Default is 0.001.
;		  Note: Sets size of kernel if size keyword is not
;		  set.
;	size	- Total size, in pixels, of kernel.  Overides
;		  clip keyword.  Must be odd!
;	center  - Pixel index of peak.  Defualt is (size-1)/2.
;		  Setting = 9999 centers peak at zero. 
;	peak	- Returns a non-unit area Gaussian with this value
;		  as its maximum
;
; HISTORY
;	Written 3/20/2002 GDB
;	Keyword PEAK added 3/10/2003 GDB
;       Modified to use ERF (faster) 4/20/2006 GDB
;-------------------------------------------------------------------

if (n_params() eq 0 and not(keyword_set(fwhm))) then begin
   print,'CALLING SEQUENCE: result=gaussker(sig,fwhm=fwhm,clip=clip,'
   print,'                           size=size,center=center,peak=peak)'
   return,-1
endif

if (n_params(0) gt 0) then sigma=sig

if keyword_set(fwhm) then begin
   ratio = double(2*sqrt(2*alog(2)))
   sigma = fwhm / ratio
endif

limit=0.001

if keyword_set(clip) then limit=clip

sz = 2*ceil(sigma*sqrt(-2.0*alog(limit))) + 1

if keyword_set(size) then begin
   if (long(size/2.) eq size/2.) then begin
      print,'Size must be odd!'
      return,-1
   endif
   sz=size
endif

;;; Create kernel.

w=findgen(sz)
cent=(sz-1)/2.
if keyword_set(center) then begin
   if (center eq 9999) then cent=0 else cent=center
endif

; Integrate over a simple Gaussian with limits for each pixel
; scaled to reflect sigma and the line center.
intlimit_lo = (w-cent-0.5)/(sqrt(2)*sigma)
intlimit_hi = (w-cent+0.5)/(sqrt(2)*sigma)
kernel      = 0.5*(erf(intlimit_hi)-erf(intlimit_lo))

; Normalize
kernel=kernel/total(kernel)

if keyword_set(peak) then begin
   oldpeak = kernel(cent)
   kernel = kernel*peak/oldpeak
endif

return,kernel

end

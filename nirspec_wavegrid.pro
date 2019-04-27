Pro nirspec_wavegrid,exposure,slitgrid,wavegrid,slitrange=slitrange,$
                     list=list,outlist=outlist

;+-------------------------------------------------------------------------
;
; NIRSPEC_WAVEGRID         4/2005
;
;
; Trace skylines or arclines in a 2D spectrum and create a 2D array of
; uncalibrated wavelengths.  The uncalibrated wavelength at a pixel is
; the x-intercept of the curve of constant wavelength that would pass
; through the center of that pixel.
;
; The user may either select lines by hand or, optionally,
; identify lines in an input file (see below).
;
; INPUTS
;      exposure    - 2D array or FITS filename, raw object frame with
;                    pronounced skylines, or arc lamp exposure.
;      slitgrid    - 2D array or FITS filename, frame with slit
;                    position at each pixel.  Typically, the slit
;                    position is the x-intercept of a line of constant
;                    slit position that passes through the center of
;                    that pixel.
;      
; KEYWORDS
;      slitrange   - 2-element array, minimum and maximum slit
;                    positions over which to trace sky or arc lines.
;                    Default is [520,800], which covers a typical
;                    low-resolution NIRSPEC exposure.
;      list        - String, name of 2-column ascii file containing x
;                    and y positions of lines to trace.  Default is
;                    to allow the user to select lines from a cut
;                    across the input exposure.
;      outlist     - String, name of 2-column ascii file to which to
;                    write the x and y positions of the identified
;                    lines.  If the keyword is set without any value,
;                    line positions will be written to 'list.dat'.
;
; OUTPUTS
;      wavegrid    - 2D array or FITS filename, uncalibrated wavelengths.
;                    If given as a string, the output will be written
;                    to a file with that name.
;
; HISTORY
;      Written 4/27/2005 GDB
;--------------------------------------------------------------------------

if (n_params() lt 3) then begin
 print,'CALLING SEQUENCE: '
 print,'   nirspec_wavegrid,exposure,slitgrid,wavegrid,slitrange=slitrange,'
 print,'                    list=list,outlist=outlist'
 return
endif

;;; Read in arrays, if necessary.

if (size(exposure,/tname) eq 'STRING') then begin
   checkexists = findfile(exposure,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Exposure file not found.'
      return
   endif else begin
      obj = readfits(exposure)
   endelse
endif else obj = exposure

if (size(slitgrid,/tname) eq 'STRING') then begin
   checkexists = findfile(slitgrid,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Slitgrid file not found.'
      return
   endif else begin
      slit = readfits(slitgrid)
   endelse
endif else slit = slitgrid

;;; Take a cut along the sky and identify lines, or read in line
;;; positions from a file.
if keyword_set(list) then begin
   templ =  {   $
          version: 1.0, $
          datastart: 0L, $
          delimiter: ' ', $
          missingvalue: -99.0, $
          commentsymbol: '#', $
          fieldcount: 2L, $
          fieldtypes: [4L,4L], $
          fieldnames: ['x','y'], $
          fieldlocations: [6L,23L], $
          fieldgroups: [0L,1L] $
        }
   position = read_ascii(list,template=templ)
   linex    = position.x
   liney    = position.y
   nlines   = n_elements(linex)
endif else begin
   cut = median(obj(698:702,*),dim=1)
   y   = findgen(1024)
   device,get_screen_size=ssz
   window,10,xs=ssz(0),ys=(ssz(1) < 600)
   !p.multi = [0,1,1]
   plot,y,cut,/xstyle,title='Select lines to trace',xtitle='y'
   print,'Left-click on lines to trace.  Lines should be fairly well spaced.'
   print,'Right-click to exit.'
   getpoints,liney,linecounts,/plotx
   nlines = n_elements(liney)
   linex  = 700 + fltarr(nlines)
endelse

;;; Optional - write line positions to a file.
if keyword_set(outlist) then begin
   if (size(outlist,/tname) eq 'STRING') then out_list = outlist $
      else out_list = 'list.dat'
   openw,list_lun,out_list,/get_lun
   for i=0,nlines-1 do printf,list_lun,linex(i),liney(i)
   close,list_lun
endif

;;; Set range of slit positions over which to trace lines.
if keyword_set(slitrange) then range=slitrange else range=[520,800]

;;; Trace lines
chan,11
traceline,obj,linex,liney,xtr,ytr,slittr,/xdir,nside=10,/bygauss,$
   median=3,coord=slit,cbound=range,merit=merit

;;; Fit polynomials in x.  Reject outliers and iterate until 
;;; fit converges.
trace_order = 2
b           = fltarr(trace_order+1,nlines)
yfit        = 0*ytr
nacross     = floor(sqrt(nlines))
ndown       = ceil((1.*nlines)/nacross)
chan,11
!p.multi    = [0,nacross,ndown]
for i=0,nlines-1 do begin
   temp_ls = where(merit(*,i) ne 0,n_temp)
   temp_x = xtr(temp_ls,i)
   temp_y = ytr(temp_ls,i)
   fit_converged = 0
   old_yfit      = 0
   n_iter        = 0
   good_ls       = indgen(n_temp)
   while not(fit_converged) do begin
      temp_b = poly_fit(temp_x(good_ls),temp_y(good_ls),trace_order,$
                        yerror=yerr)
      temp_yfit   = 0*temp_y
      for j=0,trace_order do temp_yfit = temp_yfit + temp_b(j)*temp_x^j
      fit_diff = temp_yfit - old_yfit
      if (max(abs(fit_diff)) eq 0) then begin
         fit_converged = 1
      endif else begin
         good_ls = where(abs(temp_y-temp_yfit) le 3*yerr,n_good)
         if (n_good lt trace_order+1) then begin
            ; Flag bad fit
            print,'Bad fit - rejecting.'
            temp_b(*)     = -9999
            fit_converged = 1
         endif
      endelse
      n_iter = n_iter + 1
      if (n_iter ge 10) then fit_converged = 1
   endwhile
   b(*,i) = temp_b
   ; Plot residuals
   temp_yresid = temp_y - temp_yfit
   plot,temp_x,temp_yresid,/xstyle,/ystyle,psym=3,$
      title='x = '+strtrim(linex(i),2)+', y = '+strtrim(liney(i),2),$
      xtitle='x',ytitle='y residual',charsize=2
endfor

;;; Fit coefficients as a function of y-intercept
; Reject bad fits
fitls = where(b(0,*) ne -9999 and b(0,*) ge -500 and b(0,*) le 1500)
b0 = b(0,fitls)
q  = fltarr(trace_order,4)
window,10,xs=(700 < ssz(0)),ys=(700 < ssz(1))
!p.multi = [0,1,1]
for i=0,trace_order-1 do begin
   print,' '
   print,'Fitting term '+strtrim(i+1,2)
   bi = b(i+1,fitls)
   qi = inter_poly_fit(b0,bi,1,include=inc)
   q(i,0:n_elements(qi)-1) = qi
endfor

;;; Construct the wavelength grid.
x = indgen(1024)
y = indgen(1024)
wavearr = xy2b(x,y,q)

;;; Optional - write to file.
if (size(wavegrid,/tname) eq 'STRING') then $
   writefits,wavegrid,wavearr else wavegrid = wavearr

return
end

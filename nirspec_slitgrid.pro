Pro nirspec_slitgrid,exposurelist,slitgrid

;+-------------------------------------------------------------------------
;
; NIRSPEC_SLITGRID         4/2005
;
; Trace objects in 2D spectra to create a 2D array of slit positions.
; The slit position at a pixel is the y-intercept of the curve of
; constant wavelength that would pass through the center of that
; pixel.
;
; Standard stars or reasonably bright objects are best for tracing.
; In order to get a good solution, the objects should be spaced out
; over the length of the slit.  A single object may used, but the slit
; positions will only be good around that one slit position.
;
; INPUTS
;      exposurelist - String, name of single-column ASCII file containing 
;                     a list of objects to trace..
; OUTPUTS
;      slitgrid     - 2D array or FITS filename, slit positions for
;                     each pixel.  If given as a string, the output
;                     will be written to a file with that name.
;
; HISTORY
;      Written 4/28/2005 GDB
;--------------------------------------------------------------------------

if (n_params() lt 2) then begin
 print,'CALLING SEQUENCE: '
 print,'   nirspec_slitgrid,exposurelist,slitgrid'
 return
endif

;;; Read in exposure list.
templ =  {   $
          version: 1.0, $
          datastart: 0L, $
          delimiter: ' ', $
          missingvalue: -99.0, $
          commentsymbol: '#', $
          fieldcount: 1L, $
          fieldtypes: [7L], $
          fieldnames: ['filename'], $
          fieldlocations: [6L], $
          fieldgroups: [0L] $
        }
list    = read_ascii(exposurelist,template=templ)
objname = list.filename
nobj    = n_elements(objname)

;;; For each file, display the array and allow the user to
;;; identify the object, trace the object and fit the traced
;;; points.
trace_order = 2
a           = fltarr(trace_order+1,nobj)
for i=0,nobj-1 do begin
   ; Read in and display the frame
   temp_objname = objname(i)
   print,' '
   print,'Now reading in '+temp_objname
   obj = readfits(temp_objname)
   help, obj
   chan,10
;   doit,obj(450:980,250:750),0,5000
   doit, obj(450:980,700:1300),80,2000
   ; Identify object
   print,'Click anywhere along the object trace'
   cursor,x0,y0,/device,/down
   x0 = x0 + 450
   y0 = y0 + 700
;   y0 = y0 + 250
 
   ; Trace object
   chan,11
   traceline,obj,x0,y0,xtr,ytr,/xdir,nside=10,/bygauss,median=3,merit=merit
   ; Fit trace points.  Reject outliers and refit until fit converges.
   fit_converged = 0
   old_xfit      = 0
   n_iter        = 0
   good_ls       = indgen(n_elements(xtr))
   while not(fit_converged) do begin
      temp_a   = poly_fit(ytr(good_ls),xtr(good_ls),trace_order,$
                        yerror=xerr)
      xfit     = 0*xtr
      for j=0,trace_order do xfit = xfit + temp_a(j)*ytr^j
      fit_diff = xfit - old_xfit
      if (max(abs(fit_diff)) eq 0) then begin
         fit_converged = 1
      endif else begin
         good_ls = where(abs(xtr-xfit) le 3*xerr,n_good)
         if (n_good lt trace_order+1) then begin
            ; Flag bad fit
            print,'Bad fit - rejecting.'
            temp_a(*)     = -9999
            fit_converged = 1
         endif
      endelse
      n_iter = n_iter + 1
      if (n_iter ge 10) then fit_converged = 1
   endwhile
   a(*,i) = temp_a
   ; Plot residuals
   xresid = xtr - xfit
   xresid_sdev = robust_sigma(xresid)
   chan,12
   plot,xresid,ytr,/xstyle,/ystyle,psym=3,title=temp_objname,$
      xtitle='x residual',ytitle='y',xrange=[-5*xresid_sdev,5*xresid_sdev]
endfor

;;; If more than one object used, then fit coefficients as a function
;;; of x-intercept.
p  = fltarr(trace_order,4)
if (nobj gt 1) then begin
   ; Reject bad fits
   fitls = where(a(0,*) ne -9999 and a(0,*) ge -500 and a(0,*) le 1500)
   a_0 = a(0,fitls)
   chan,10
   !p.multi = [0,1,1]
   for i=0,trace_order-1 do begin
      print,' '
      print,'Fitting term '+strtrim(i+1,2)
      a_i = a(i+1,fitls)
      p_i = inter_poly_fit(a_0,a_i,1,include=inc)
      p(i,0:n_elements(p_i)-1) = p_i
   endfor
endif else begin
   p(*,0) = a(1:trace_order)
endelse

;;; Construct the slit position grid.
;x = indgen(1024)
;y = indgen(1024)
help, p
x = indgen(2048)
y = indgen(2048)
slitarr = xy2a(x,y,p)

;;; Optional - write to file.
if (size(slitgrid,/tname) eq 'STRING') then $
   writefits,slitgrid,slitarr else slitgrid = slitarr

return
end

Function lris_subbias,raw,blue=blue,red=red,npoly=npoly,flat=flat,$
                      plot=plot

;+---------------------------------------------------------------------
;
; LRIS_SUBBIAS      7/2004
;
; Subtract the bias from a raw LRIS frame by fitting a low-order
; polynomial in the dispersion direction to the overscan region for
; each amplifier.
;
; Must specify eiter blue or red using keywords.
;
; INPUTS
;          raw     - Raw 2D array, can be a named variable or the name
;                    of a FITS file.
;
; KEYWORDS
;          blue    - LRIS-Blue.  Currently assumes four amplifiers,
;                    but only subtracts the bias for the second and
;                    third amplifiers (since this is the region
;                    illuminated by long-slit spectra).
;          red     - LRIS-Red.  Currently assumes two amplifiers.
;          npoly   - Order of polynomial used to fit to bias in the
;                    dispersion direction.  Default is 1 (linear).
;          flat    - Raw object is a flat-field, therefore, to avoid
;                    bleed problems from high counts,  use only
;                    the non-exposed rows to determine the bias level
;                    (npoly = 0 or 1 recommended).  Only effects
;                    red-side exposures.
;          plot    - Show plots of fits and residuals.  User will be 
;                    prompted to continue after each amplifier.
;
; OUTPUTS
;          <result> - 2D bias-subtracted array.
;
; HISTORY
;          Written 7/22/2004 GDB
;---------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE: result = lris_subbias(raw,/blue,/red,npoly=npoly,'
   print,'                                          /flat,/plot)'
   return,-1
endif

; Read in FITS file, if necessary
if (size(raw,/tname) eq 'STRING') then begin
   checkexists = findfile(raw,count=nfiles)
   if (nfiles eq 0) then begin
      print,'LRIS_SUBBIAS: Raw file not found.'
      return,-1
   endif else begin
       rawframe = readfits(raw,hrawframe)
       fromfitsfile = 1 
   endelse
endif else rawframe = raw


; Make sure that either the blue or the red keyword is set
if ( not(keyword_set(blue) or keyword_set(red)) or $
     (keyword_set(blue) and keyword_set(red)) ) then begin
   print,'LRIS_SUBBIAS: Must specify either blue or red side.'
   return,-1
endif

; Create x-grid and y-grid
sz  = size(rawframe)
xsz = sz(1)
YSZ = SZ(2)
xgrid = intarr(xsz,ysz)
ygrid = intarr(xsz,ysz)
for i=0,xsz-1 do xgrid(i,*) = i
for i=0,ysz-1 do ygrid(*,i) = i

; Specify blue or red side
if keyword_set(blue) then begin
   blue = 1
   red  = 0
endif else begin
   blue = 0
   red  = 1
endelse

; Set number of amps
if blue then namps=4 else namps=2

; Set starting x-values and sizes of data and overscan regions 
if blue then begin
   case xsz of
      4620: begin ; 1x1 binning
               datastart  = [204,1228,2252,3276]
               oscanstart = [4300,4380,4460,4540]
               dataxsize  = 1024   
               oscanxsize = 80
            end
      2308: begin ; 2x2 binnin
               datastart  = [100,612,1124,1636]
               oscanstart = [2148,2188,2228,2268]
               dataxsize  = 512
               oscanxsize = 40
            end
      else: begin
               print,'LRIS_SUBBIAS: Error - Array size not recognized.'
               return,-1
            end
   endcase
endif else begin
   case xsz of
      2248: begin
               datastart  = [40,1064]
               oscanstart = [2088,2168]
               dataxsize  = 1024   
               oscanxsize = 80
            end
      2048: begin
               print,'LRIS_SUBBIAS: Overscan region has been removed.'
               return,-1
            end
      else: begin
               print,'LRIS_SUBBIAS: Error - Array size not recognized.'
               return,-1
            end
   endcase
endelse

; If red-side flat field, omit exposed rows from fit.
if (red and keyword_set(flat)) then begin
   badymin = 85
   badymax = 944
endif

; Polynomial order
if keyword_set(npoly) then polyorder= npoly else polyorder=1

; Bias-subtracted frame
subframe = float(rawframe)

; For each amp, fit a polynomial to the overscan region and subtract
; the fit from the corresponding data region

for i=0,namps-1 do begin

   ; Select overscan region
   oscan_temp = rawframe(oscanstart(i):oscanstart(i)+oscanxsize-1,*)

   ; Compute median value for each row
   oscan_1d = median(oscan_temp,dim=1,/even)
   y = findgen(ysz)

   ; Omit exposed rows?
   if (red and keyword_set(flat)) then begin
      goodyls = where(y lt badymin or y gt badymax)
      oscan_good = oscan_1d(goodyls)
      y_good = y(goodyls)
   endif else begin
      y_good = y
      oscan_good = oscan_1d
   endelse

   ; Fit polynomial and evaluate for all y
   pfit = poly_fit(y_good,oscan_good,polyorder,yfit=oscan_good_fit)
   oscan_fit = fltarr(ysz)
   for n=0,polyorder do oscan_fit = oscan_fit + pfit(n) * y^n

   ; Plot fit and residuals
   if keyword_set(plot) then begin
      currentmulti = !p.multi
      !p.multi = [0,1,2]
      plot,y_good,oscan_good,psym=3,title = 'Fit',$
         yrange=[min(oscan_good),max(oscan_good)]
      oplot,y,oscan_fit
      plot,y_good,oscan_good - oscan_good_fit,psym=3,title = 'Residuals'
      !p.multi = currentmulti
      print,'Press any key to continue...'
      key = get_kbrd(1)
   endif

   ; Subtract from overscan and data regions
   for j=0,ysz-1 do begin
      subframe(oscanstart(i):oscanstart(i)+oscanxsize-1,j) = $
         subframe(oscanstart(i):oscanstart(i)+oscanxsize-1,j) - oscan_fit(j)
      subframe(datastart(i):datastart(i)+dataxsize-1,j) = $
         subframe(datastart(i):datastart(i)+dataxsize-1,j) - oscan_fit(j)
   endfor

endfor

return,subframe

end

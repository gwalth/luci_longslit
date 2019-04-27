Function xy2a,x,y,p,lrisb=lrisb,xbin=xbin,ybin=ybin

;+--------------------------------------------------------------------
;
; XY2A    February 2005
;
; Convert from x,y coordinates to "a_0" coordinates, where loci of
; constant a_0 are polynomials in y.  The coefficients of those
; polynomials have themselves been fit as a function of a_0 by
; polynomials with coeffiecients p.  The a_0 value at a given pixel is
; the x-intercept of the polynomial in y that passes through that
; pixel.
;
; This program is meant for converting from x,y pixel corrdinates in a
; 2D spectrum to proxies for slit position or wavelenth.  
;
; Prior to running this program, lines (e.g., arc/sky lines or flat
; field edges) that run along the y-direction have been fit as
;           x = a_0 + a_1*y + a_2*y^2 + a_3*y^3 + ... + a_m*y^m
; where a_i is parametrized as a polynomial in a_0:
;           a_i = p_i(0) + p_i(1)*a_0 + p_i(2)*a_0^2 + p_i(3)*a_0^3
;
; Input x,y,p ---> Get out a_0 which is a  proxy for (e.g.) wavelength
; or slit position.
;
; The coefficients that describe a_i as a function of a_0 are input as
; an m x 4 array, such that p_i(j) = p(i,j).  Some of the terms may be
; 0 to allow lower-order fits.
;
; The values of a_0 at a given x,y is the solution to the equation
;
;           0 = s0 + s1*a_0 + s2*a_0^2 + s3*a_0^3
; where
;           s0 = -x + Sum[p(i,0)*x^(i+1)]
;           s1 = 1 + Sum[p(i,1)*x^(i+1)]
;           s2 = Sum[p(i,2)*x^(i+1)]
;           s3 = Sum[p(i,3)*x^(i+1)]
;
; INPUTS
;          x,y   - Pixel coordinates.  If x is an nx-element vector
;                  and y is an ny element vector, the result will
;                  be an nx x ny element array.
;          p     - m x 4 coeffient array describing line traces that run
;                  in the y-direction, where m can be any integer.
;
; KEYWORDS
;          lrisb - Take into account the gap between the chips
;                  in LRIS-B.  If set, the offset and tilt
;                  of the right-hand chip are accounted for, assuming
;                  that the "true" x and y axes are those defined by
;                  the left-hand chip.  NOTE: The program assumes
;                  that the prescan regions have NOT been removed,
;                  and so the correction is applied for all input
;                  x-values >= 2252 (IDL zero-centric coordinates).
;                  Also assumes 1x1 binning.
;          xbin  - Binning in x-direction (LRIS-B only)
;          ybin  - Binning in y-direction (LRIS-B only)
;
; OUTPUTS
;          <result> = nx by ny array, where nx is the
;                     number of x coordinates and ny is the number
;                     of y coordinates, containing the value of
;                     a_0 at each point.
;
; HISTORY
;          XY2AB written 4/2/2004 GDB
;          XY2AB modified to allow an aribitrary number of terms
;            in p and q  8/20/2004 GDB
;          Adapted from xy2ab.pro 2/4/2005 GDB
;          XBIN and YBIN keywords added 9/14
;-----------------------------------------------------------------

if (n_params() lt 3) then begin
   print,' CALLING SEQUCE result = xy2a(x,y,p,/lrisb,xbin=xbin,ybin=ybin)'
   return,-1
endif

; Convert inputs to double precision

xsave = x
ysave = y
psave = p

x = double([x])
y = double([y])
p = double(p)

; Make sure p is a (some number) x 4 array.

psz = size(p)
if not(psz(0) eq 2 and psz(2) eq 4) then begin
   print,'p must be a (some number) x 4 array.'
   return,-1
endif
nyterms = psz(1)

; Determine the degree of the fits in a_0:

lsp3 = where(p(*,3) ne 0) 
lsp2 = where(p(*,2) ne 0)

if (total(lsp3) eq -1) then begin
   if (total(lsp2) eq -1) then afit = 'linear' else afit = 'quadratic'
endif else afit = 'cubic'

; Solve for a_0 at each x,y

nx = n_elements(x)
ny = n_elements(y)

xarr = dblarr(nx,ny)
yarr = dblarr(nx,ny)
for i=0,nx-1 do xarr(i,*) = x(i)
for j=0,ny-1 do yarr(*,j) = y(j)

; If this is LRIS-B data, perform shift for pixels in the
; right-hand chip.  Transformations provided by Alice Shapley.
if keyword_set(lrisb) then begin
   ; Default it 1x1 binning
   if keyword_set(xbin) then x_bin = xbin else x_bin = 1.
   if keyword_set(ybin) then y_bin = ybin else y_bin = 1.
   case x_bin of
     1: rh_ls = where(xarr ge 2252,n_rh)
     2: rh_ls = where(xarr ge 1124,n_rh)      
     else : begin
                print,'XY2A: Error - Not sure how to handle binning.'
                return,-1
            end
   endcase
   if (n_rh ne 0) then begin
      xoffset = 103.192 / x_bin
      yoffset = -2.982 / y_bin
      xscale  = 0.00162 * (y_bin/x_bin)
      x_rh = xarr(rh_ls)
      y_rh = yarr(rh_ls)
      xarr(rh_ls) = x_rh + xoffset + xscale*(y_rh + 1) ; IDL coords
      yarr(rh_ls) = y_rh + yoffset
   endif
endif

; Analytic solutions for a

a_array = dblarr(nx,ny)

; Avoid numerical problems when y=0
yeq0_ls = where(yarr eq 0,n_yeq0,complement=yne0_ls,ncomplement=n_yne0)

if (n_yne0 ne 0) then begin

   s0 = -xarr
   for k=0,nyterms-1 do s0 = s0 + p(k,0)*yarr^(k+1)

   s1 = 1.
   for k=0,nyterms-1 do s1 = s1 + p(k,1)*yarr^(k+1)

   case afit of
      'linear'    : begin
                       a_array = -s0 / s1
                    end
      'quadratic' : begin
                       s2 = 0.
                       for k=0,nyterms-1 do s2 = s2 + p(k,2)*yarr^(k+1)
                       a_array = (-s1 + sqrt(s1^2 - 4*s2*s0)) / (2*s2)
                    end
      'cubic'     : begin
                       s2 = 0.
                       for k=0,nyterms-1 do s2 = s2 + p(k,2)*yarr^(k+1)
                       s3 = 0.
                       for k=0,nyterms-1 do s3 = s3 + p(k,3)*yarr^(k+1)
                       ; Slow method.
                       for i=0,nx-1 do begin
                          for j=0,ny-1 do begin
                             aroots = cubic_roots([s0(i,j),s1(i,j),$
                                                   s2(i,j),s3(i,j)])
                             if (n_elements(aroots) eq 3) then $
                                a_array(i,j) = aroots(2) else $
                                a_array(i,j) = aroots(0)
                          endfor
                       endfor
                       ;Rewrite cubic_roots to accept/return arrays?
                       ;aroots = cubic_roots([s0,s1,s2,s3])
                    end
      else        : begin
                       print,'Error - unknown fit to a coeffs.'
                       return,-1
                    end
   endcase

endif

if (n_yeq0 ne 0) then a_array(yeq0_ls) = xarr(yeq0_ls)

x = xsave
y = ysave
p = psave

return,a_array

end

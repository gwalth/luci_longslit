Function xy2b,x,y,q,lrisb=lrisb,xbin=xbin,ybin=ybin

;+--------------------------------------------------------------------
;
; XY2B    February 2005
;
; Convert from x,y coordinates to "b_0" coordinates, where loci of
; constant b_0 are polynomials in x.  The coefficients of those
; polynomials have themselves been fit as a function of b_0 by
; polynomials with coeffiecients q.  The b_0 value at a given pixel is
; the y-intercept of the polynomial in x that passes through that
; pixel.
;
; This program is meant for converting from x,y pixel corrdinates in a
; 2D spectrum to proxies for slit position or wavelenth.  
;
; Prior to running this program, lines (e.g., arc/sky lines or flat
; field edges) that run along the x-direction have been fit as
;           y = b_0 + b_1*x + b_2*x^2 + b_3*x^3 + ... + b_n*x^n
; where b_i is parametrized as a polynomial in b_0:
;           b_i = q_i(0) + q_i(1)*b_0 + q_i(2)*b_0^2 + q_i(3)*b_0^3
; 
; Input x,y ---> Get out b_0 which is a  proxy for (e.g.) wavelength
; or slit position.
;
; The coefficients that describe b_i as a function of b_0 are input as
; an n x 4 array, such that q_i(j) = q(i,j).  Some of the terms may be
; 0 to allow lower-order fits.
;
; The value of b_0 at a given x,y is the solution to the equation
;
;           0 = t0 + t1*b_0 + t2*b_0^2 + t3*b_0^3
; where
;           t0 = -y + Sum[q(i,0)*y^(i+1)]
;           t1 = 1 + Sum[q(i,1)*y^(i+1)]
;           t2 = Sum[q(i,2)*y^(i+1)]
;           t3 = Sum[q(i,3)*y^(i+1)]
;
; INPUTS
;          x,y   - Pixel coordinates.  If x is an nx-element vector
;                  and y is an ny element vector, the result will
;                  be an nx x ny element array.
;          q     - m x 4 coeffient array describing line traces that run
;                  in the x-direction, where m can be any integer.
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
;          xbin  - Binning in x-direction (LRIS-B only)
;          ybin  - Binning in y-direction (LRIS-B only)
;
; OUTPUTS
;          <result> = nx by ny array, where nx is the
;                     number of x coordinates and ny is the number
;                     of y coordinates, containing the value of
;                     b_0 at each point.
;
; HISTORY
;          XY2AB written 4/2/2004 GDB
;          XY2AB modified to allow an aribitrary number of terms
;          in p and q  8/20/2004 GDB
;          Adapted from xy2ab.pro 2/4/2005 GDB
;          XBIN and YBIN keywords added 9/13/2006 GDB
;-----------------------------------------------------------------

if (n_params() lt 3) then begin
   print,' CALLING SEQUCE result = xy2b(x,y,q,/lrisb,xbin=xbin,ybin=ybin)'
   return,-1
endif

; Convert inputs to double precision

xsave = x
ysave = y
qsave = q

x = double([x])
y = double([y])
q = double(q)

; Make sure q is a (some number) x 4 array.

qsz = size(q)
if not(qsz(0) eq 2 and qsz(2) eq 4) then begin
   print,'q must be a (some number) x 4 array.'
   return,-1
endif
nxterms = qsz(1)

; Determine the degree of the fits in b_0:

lsq3 = where(q(*,3) ne 0) 
lsq2 = where(q(*,2) ne 0)

if (total(lsq3) eq -1) then begin
   if (total(lsq2) eq -1) then bfit = 'linear' else bfit = 'quadratic'
endif else bfit = 'cubic'

; Solve for a_0,b_0 at each x,y

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
                print,'XY2B: Error - Not sure how to handle binning.'
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

; Analytic solutions for b

b_array = dblarr(nx,ny)

; Identify pixles where x=0
xeq0_ls = where(xarr eq 0,n_xeq0,complement=xne0_ls,ncomplement=n_xne0)

if (n_xne0 ne 0) then begin

   t0 = -yarr
   for k=0,nxterms-1 do t0 = t0 + q(k,0)*xarr^(k+1)

   t1 = 1.
   for k=0,nxterms-1 do t1 = t1 + q(k,1)*xarr^(k+1)
 
   case bfit of
      'linear'    : begin
                             b_array = -t0 / t1
                    end
      'quadratic' : begin
                       t2 = 0.
                       for k=0,nxterms-1 do t2 = t2 + q(k,2)*xarr^(k+1)
                       b_array = (-t1 + sqrt(t1^2 - 4*t2*t0)) / (2*t2)
                    end
      'cubic'     : begin
                       t2 = 0.
                       for k=0,nxterms-1 do t2 = t2 + q(k,2)*xarr^(k+1)
                       t3 = 0.
                       for k=0,nxterms-1 do t3 = t3 + q(k,3)*xarr^(k+1)
                       ; Slow method
                       for i=0,nx-1 do begin
                          for j=0,ny-1 do begin
                             broots = cubic_roots([t0(i,j),t1(i,j),$
                                                   t2(i,j),t3(i,j)])
                             if (n_elements(broots) eq 3) then $
                                b_array(i,j) = broots(2) else $
                                b_array(i,j) = broots(0)
                          endfor
                       endfor
                       ;Rewrite cubic_roots to accept/return arrays?
                       ;broots = cubic_roots([t0,t1,t2,t3])
                    end
      else        : begin
                       print,'Error - unknown fit to b coeffs.'
                       return,-1
                    end
   endcase

endif

; Avoid numerical problems when x = 0
if (n_xeq0 ne 0) then b_array(xeq0_ls) = yarr(xeq0_ls)

x = xsave
y = ysave
q = qsave

return,b_array

end

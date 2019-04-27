;+
; NAME:
;   astrom_tweak
;
; PURPOSE:
;   Tweak astrometric solution, given a good initial guess
;
; CALLING SEQUENCE:
;   gsa_out = astrom_tweak(cat, im, maxsep, gsa_in, $
;    [ nminmatch=, /radial, maxrad=, errflag=, nmatch= ]
;
; INPUTS:
;   cat        - Structure (with fields .ra, .dec) of catalogue positions
;   im         - Structure (with fields .x, .y) of image star positions
;   maxsep     - Maximum allowed separation for good match (pixels)
;   gsa-in     - Initial guess for astrometric solution (struct)
;
; OPTIONAL KEYWORDS:
;   nminmatch  - Minimum number of stars for match; default to 5
;   radial     - include higher-order radial terms in fit
;   maxrad     - passed to ptstarmatch (???)
;  
; OUTPUTS:
;   gsa_out    - returned guess for astrometric solution (struct);
;                0 if solution failed
;
; OUTPUT OUTPUTS:
;   errflag    - set to 1 if fatal error occurs
;   nmatch     - Number of matched stars
;
; COMMENTS:
;
;   Uses preliminary solution given in astr structure to match image
;    and catalogue stars within maxsep pixels of each other.  These
;    are then used by astrom_warp to determine a new solution, returned
;    in astr.
; 
;   cat (.ra, .dec) will contain values from USNO SA2.0 catalogue or
;    other catalogue, possibly derived from PT frames with good astrometry
;
; BUGS:
;
; PROCEDURES CALLED:
;   poly_iter
;   gsssadxy
;   gsssxyad
;
; INTERNAL SUPPORT PROCEDURES:
;   astrom_starmatch
;   astrom_warp
;
; REVISION HISTORY:
;   26-Aug-2002  Written by D. Schlegel, Princeton.
;                Modified from D. Finkbeiner's PT_TWEAK_ASTR.
;-
;-----------------------------------------------------------------------------
; Returns the indices of the catalog stars that match with the image stars

pro astrom_starmatch, cat, im, maxsep, indcat, indobs, maxrad=maxrad

   ncat = n_elements(cat)
   nobs = n_elements(im)

   indcat = lonarr(ncat)
   indobs = lonarr(ncat)

   if (keyword_set(maxrad) EQ 0) then maxrad = 1E9
   mx = mean(im.x)
   my = mean(im.y)
   rdif = sqrt((im.x-mx)^2+(im.y-my)^2)
   in_rad = rdif LT maxrad

   k = 0L

   for i=0L, ncat-1 do begin

      dist = sqrt((cat[i].x-im.x)^2 + (cat[i].y-im.y)^2)

      whmin = where(dist LT maxsep, ct)
      if (ct EQ 1) then begin ; Only do if exactly 1 star is in range
         if (in_rad[whmin[0]]) then begin
            indcat[k] = i
            indobs[k] = whmin
            k = k+1
         endif
      endif
   endfor

   ; Trim extra zeros

   splog, k, ' matches between catalogue and image stars.'
   if (k EQ 0) then begin
      indcat = -1L
      indobs =  -1L
   endif else begin
      indcat = indcat[0:k-1]
      indobs = indobs[0:k-1]
   endelse

   return
end
;-----------------------------------------------------------------------------
pro astrom_warp, catmatch, immatch, deltax, deltay, rot, shift=shift

; determines the x and y offsets, rotation, and plate scale (first degree)

   degree = 1

; We subtract 1024.5 here so rotations are about
;  the center of image, not (0,0)

   if (NOT keyword_set(shift)) then shift = [1024.5, 1024.5]

   x_0 = immatch.x-shift[0]
   y_0 = immatch.y-shift[1]
   x_i = catmatch.x-shift[0]
   y_i = catmatch.y-shift[1]

; coefficients   

; (rot will always be 2x2)
   rot = dblarr(2, 2)

   polywarp, x_i, y_i, x_0, y_0, degree, k_x, k_y

; translational offset   

   deltax = k_x[0,0]
   deltay = k_y[0,0]

; rotation correction (radians)

   rot[0,0] = k_x[0,1]
   rot[0,1] = k_y[0,1]
   rot[1,0] = k_x[1,0]
   rot[1,1] = k_y[1,0]

; plate scale

   splog, 'Inital scale guess off by factor of', determ(rot)

; cross terms

   crossx = k_x[1,1]
   crossy = k_y[1,1]
   
   if (abs(crossx) > abs(crossy)) GT 1E-4 then begin
       splog, 'warning: Cross terms too big'
   endif 

   return
end 
;-----------------------------------------------------------------------------
function astrom_tweak, cat, im, maxsep, gsa_in, radial=radial, maxrad=maxrad, $
 nminmatch=nminmatch, errflag=errflag, catind=catind, obsind=obsind, $
 nmatch=nmatch

   ; Set default return values
   catind = -1L
   obsind = -1L
   nmatch = 0L

   if (keyword_set(errflag)) then return, 0
   if (NOT keyword_set(nminmatch)) then nminmatch = 5
  
   ; fill cat.x and .y fields using current astr structure

   gsa_out = gsa_in
   gsssadxy, gsa_out, cat.ra, cat.dec, catx, caty

   cat.x = catx
   cat.y = caty

   ; find matches between USNO catalogue and image stars; return index arrays
   astrom_starmatch, cat, im, maxsep, catind, obsind, maxrad=maxrad
   nmatch = n_elements(catind) * (catind[0] NE -1)
   if (nmatch LT nminmatch) then begin
      splog, 'Only', nmatch, ' stars found - skipping'
      nmatch = 0L
      errflag = 1
      return, 0
   endif 

   catmatch = cat[catind]        ; cat stars that match
   immatch = im[obsind]          ; image stars that match

   xcen = gsa_out.ppo3/gsa_out.xsz-.5d   ; 1023.5
   ycen = gsa_out.ppo6/gsa_out.ysz-.5d   ; 1023.5

   nord = 3
   if (keyword_set(radial)) then begin ; fit radial part
     
      cx = catmatch.x-xcen
      cy = catmatch.y-ycen
      ix = immatch.x-xcen
      iy = immatch.y-ycen
;      xsi = cdelt[0]*(cd[0,0]*ix + cd[0,1]*iy) ;no matrix notation, in
;      eta = cdelt[1]*(cd[1,0]*ix + cd[1,1]*iy) ;case X and Y are vectors
 
      catr = sqrt(cx^2+cy^2)
      imr  = sqrt(ix^2+iy^2)
      delr = imr-catr

      poly_iter, [imr, -imr], [delr, -delr], nord, 3.0, yfit, coeff=coeff
      cubecoeff  = -coeff[3]

      ; this is not strictly correct - astrom_tweak must be iterated a few times
      qd = (-cubecoeff*imr^2)
      catmatch.x = catmatch.x+ix*qd ; bend catalogue to image
      catmatch.y = catmatch.y+iy*qd
   endif
  
   ; astro_warp determines the translational, rotational, and plate scale coeffs
  
   astrom_warp, catmatch, immatch, deltax, deltay, rot, $
    shift=[xcen, ycen]

   ; This routine coverts X,Y to RA, dec using information in astr
   ; structure which contains the initial guess...

   gsssxyad, gsa_out, xcen+deltax, ycen+deltay, tru_ra, tru_dec

   err_ra  = tru_ra-gsa_out.crval[0]
   err_dec = tru_dec-gsa_out.crval[1]
   err = [err_ra, err_dec]*3600.   ; arcsec
   splog, 'Pointing error: ', err[0], err[1], ' arcseconds'

   ; update astr structure with results of ptwarp
   gsa_out.crval = [tru_ra, tru_dec]
   cd = [[gsa_out.amdx[0], gsa_out.amdy[1]], [gsa_out.amdx[1], gsa_out.amdy[0]]]/3600.
   sc = sqrt(abs(determ(cd)))*3600.
   cubelast = (gsa_out.amdx[10]+gsa_out.amdx[11])/(gsa_out.amdx[1]+gsa_out.amdx[0])*sc

   cd = transpose(rot)##cd
   gsa_out.amdx[0]  = cd[0, 0] *3600. ; x  coeff
   gsa_out.amdx[1]  = cd[0, 1] *3600. ; y  
   gsa_out.amdy[0]  = cd[1, 1] *3600. ; y  coeff (not a typo!)
   gsa_out.amdy[1]  = cd[1, 0] *3600. ; x
   sc = sqrt(abs(determ(cd)))*3600.

   if (n_elements(cubecoeff) EQ 0) then cubecoeff = 0
   cubecoeff = cubecoeff + cubelast

;   astr.projp1 = astr.projp1+cubecoeff
   ; Set amd coeffs corresponding to r^3 term 
   gsa_out.amdx[8]  = cubecoeff/sc * gsa_out.amdx[1] ; y x^2
   gsa_out.amdx[10] = cubecoeff/sc * gsa_out.amdx[1] ; y y^2
   gsa_out.amdx[11] = cubecoeff/sc * gsa_out.amdx[0] ; x r^2
   gsa_out.amdy[8]  = cubecoeff/sc * gsa_out.amdy[1] ; x y^2
   gsa_out.amdy[10] = cubecoeff/sc * gsa_out.amdy[1] ; x x^2
   gsa_out.amdy[11] = cubecoeff/sc * gsa_out.amdy[0] ; y r^2

   gsssadxy, gsa_out, cat.ra, cat.dec, catx, caty
   cat.x = catx & cat.y = caty

   return, gsa_out
end
;------------------------------------------------------------------------------

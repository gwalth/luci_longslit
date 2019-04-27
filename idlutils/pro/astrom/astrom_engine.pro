;+
; NAME:
;   astrom_engine
;
; PURPOSE:
;   Compute astrometric solution for a list of stars & catalogue stars
;
; CALLING SEQUENCE:
;   gsa_out = astrom_engine( xpos, ypos, catlon, catlat, gsa_in, $
;    [ search_rad=, search_scale=, search_angle=, $
;    poserr=, nmatch=, catind=, obsind=, /radial ] )
;
; INPUTS:
;   xpos       - X positions in CCD coordinates
;   ypos       - Y positions in CCD coordinates
;   catlon     - Catalog star longitudes in the same coordinate system as GSA_IN
;   catlat     - Catalog star latitutes in the same coordinate system as GSA_IN
;   gsa_in     - Input GSSS structure with initial guess for astrometric
;                solution
;   radial     - If set, then fit for radial distortion terms; default to 0
;
; OPTIONAL INPUTS:
;   search_rad   - Unused ???
;   search_scale - Unused ???
;   search_angle - ???
;   poserr       - ???
;
; OUTPUTS:
;   gsa_out    - Output GSSS structure with astrometric solution;
;                return 0 if astrometry failed
;
; OPTIONAL OUTPUTS:
;   nmatch     - Number of matched objects with the input catalog.
;   catind     - Indices of CATLON,CATLAT for matched objects.
;   obsind     - Indices of XPOS,YPOS for matched objects.
;
; COMMENTS:
;   We assume that we know the scale and rotation well enough, then solve
;   for the X,Y offsets by correlating with catalog stars.
;
; BUGS:
;
; PROCEDURES CALLED:
;   angle_from_pairs()
;   astrom_tweak
;   gsssadxy
;   gsssxyad
;   offset_from_pairs
;
; REVISION HISTORY:
;   10-Jun-2002  Written by D. Schlegel & D. Finkbeiner, Princeton.
;-
;------------------------------------------------------------------------------
function astrom_engine, xpos, ypos, catlon, catlat, gsa_in, $
 search_rad=search_rad, $
 search_scale=search_scale, search_angle=search_angle, $
 poserr=poserr, nmatch=nmatch, catind=catind, $
 obsind=obsind, radial=radial

   if (n_params() LT 5) then begin
      doc_library, 'astrom_engine'
      return, 0
   endif
   if (size(gsa_in,/tname) NE 'STRUCT') then $
    message, 'Must pass gsa structure with initial guess'

   if (NOT arg_present(radial)) then radial = 0B

   ;----------
   ; Set default return values

   catind = -1L
   obsind = -1L
   nmatch = 0L

   ;----------
   ; Set criterion for match on first pass through tweak_astrom [pix]
   ; The factor 0.6 deg is pulled out of thin air, but should be based
   ; upon the step size used to search for the position angle error.

   xsep = (max(xpos) - min(xpos)) > poserr
   ysep = (max(ypos) - min(ypos)) > poserr
   maxsep = poserr + 0.6 * sqrt(xsep^2 + ysep^2) / !radeg
   binsz = 0.5 * maxsep
   dmax = xsep > ysep

; These are some typical values ???
;maxsep = 25
;binsz = 5L
;dmax = 2000L

   gsa1 = gsa_in
   gsssadxy, gsa1, catlon, catlat, catx, caty

   ;----------
   ; Search for the angle between the catalogue stars and image stars

   if (search_angle GT 1.) then begin 
      ang = angle_from_pairs(catx, caty, xpos, ypos, $
       dmax=dmax, binsz=binsz, bestsig=bestsig, angrange=[-1., 1]*search_angle)
    
      if (bestsig gt 12) then begin 
         splog, 'Best Angle: ', ang, '  sigma: ', bestsig
         cd = fltarr(2, 2)
         cd[0, 0] = gsa1.amdx[0]
         cd[0, 1] = gsa1.amdx[1]
         cd[1, 1] = gsa1.amdy[0]
         cd[1, 0] = gsa1.amdy[1]
         angrad = ang * !pi / 180.
         mm = [[cos(angrad), sin(angrad)], [-sin(angrad), cos(angrad)]]
         cd = cd # mm
         print, cd
         print
         print, mm
         gsa1.amdx[0] = cd[0, 0]
         gsa1.amdx[1] = cd[0, 1]  
         gsa1.amdy[0] = cd[1, 1]
         gsa1.amdy[1] = cd[1, 0]
         gsssadxy, gsa1, catlon, catlat, catx, caty
      endif else begin
         splog, 'Warning: I think I am lost, but I will try anyway...'
      endelse
   endif 
  
   ;----------
   ; Search for the X,Y offset between the catalogue stars and image stars

   xyshift = offset_from_pairs(catx, caty, xpos, ypos, $
    dmax=dmax, binsz=binsz, errflag=errflag, bestsig=bestsig)

   if (errflag NE 0) then begin
      splog, 'XY shift FAILED in astrom_engine'
      return, 0
   endif
  
   splog, 'XYSHIFT: ', xyshift * binsz, ' pix'
  
   xcen = gsa1.ppo3 / gsa1.xsz - 0.5d   ; 1023.5
   ycen = gsa1.ppo6 / gsa1.ysz - 0.5d   ; 1023.5

   ; NOTE: FITS crpix is 1-indexed but argument of xy2ad is 0-indexed
   refpix = [xcen, ycen] - xyshift
   gsssxyad, gsa1, refpix[0], refpix[1], racen, deccen

   ; Update astrometry structure with new CRVALs
   gsa1.crval = [racen, deccen]

   ; Update catalogue .x and .y fields
   gsssadxy, gsa1, catlon, catlat, catx, caty

   im_template = {im_specs, $
                  x:   0.0, $
                  y:   0.0, $
                  ra:  0.d, $
                  dec: 0.d  }

   im = replicate(im_template, n_elements(xpos))
   im.x = xpos
   im.y = ypos

   ;----------
   ; Tweak astrometry structure with cat (ra,dec) and im (x,y) comparison

   ; First pass
   gsssadxy, gsa1, catlon, catlat, catx, caty

   cat = replicate(im_template, n_elements(catlat))
   cat.ra = catlon
   cat.dec = catlat
   cat.x = catx
   cat.y = caty

   gsa1 = astrom_tweak(cat, im, maxsep, gsa1, errflag=errflag, nmatch=nmatch)

   ; Hardwired ???
   maxrad = [1000, 1200, 2000]

   if (NOT keyword_set(errflag)) then $
    gsa1 = astrom_tweak(cat, im, maxsep/2, gsa1, radial=radial, $
     maxrad=maxrad[0], errflag=errflag, nmatch=nmatch)

   if (NOT keyword_set(errflag)) then $
    gsa1 = astrom_tweak(cat, im, maxsep/4., gsa1, radial=radial, $
     maxrad=maxrad[1], errflag=errflag, nmatch=nmatch)

   if (NOT keyword_set(errflag)) then $
    gsa1 = astrom_tweak(cat, im, maxsep/8., gsa1, radial=radial, $
     errflag=errflag, nmatch=nmatch)

   if (NOT keyword_set(errflag)) then $
    gsa1 = astrom_tweak(cat, im, maxsep/8., gsa1, radial=radial, $
     errflag=errflag, catind=catind, obsind=obsind, nmatch=nmatch)

   if (errflag NE 0) then begin
      splog, 'XY shift FAILED in astrom_tweak'
      return, 0
   endif

   ; Compute rotation angle w.r.t. initial guess
   xvec0 = gsa_in.amdx
   xvec1 = gsa1.amdx
   angerr = acos((transpose(xvec0[0:1]) # xvec1[0:1]) $
    / sqrt(total(xvec0[0:1]^2)*total(xvec1[0:1]^2))) / !dtor
   splog, 'Initial guess rotated by: ', angerr, ' deg'

   ; Compute the offset of the corner pixel (0,0)
   gsssxyad, gsa_in, 0.0, 0.0, ra_in, dec_in
   gsssxyad, gsa1, 0.0, 0.0, ra_out, dec_out
   gsssadxy, gsa1, ra_in, dec_in, xoff, yoff
   splog, 'RA offset = ', (ra_out - ra_in) * 3600. / cos(dec_out/!radeg), ' arcsec'
   splog, 'DEC offset = ', (dec_out - dec_in) * 3600., ' arcsec'
   splog, 'X offset = ', xoff, ' pix'
   splog, 'Y offset = ', yoff, ' pix'

   return, gsa1
end
;------------------------------------------------------------------------------

Pro nirspec_divstd,specfile,spectype,outfile

;+-------------------------------------------------------------------
;
; NIRSPEC_DIVSTD     6/2004 GDB
;
; Divide 1D stellar spectrum by normalized standard flux values.  Uses
; Meyers (1998) H-band standards.  Note that Meyers standards cover
; only 15154.8 - 17860.0 AA.  Any additional blue wavelengths covered
; by input spectrum are assumed to be unity in the standard.  The is a
; very good approximation for A0 stars but a relatively poor one for F
; & G stars.
;
; Program shifts the standard in wavelength to match the input
; spectrum by performing a crosscorelation.  (The input spectrum may
; have been adjusted for a heliocentric velocity shift.  However, the
; relative difference, in angstroms, between the blue and red ends of
; the spectrum for a heliocentric velocity correction is ~< 0.3 AA, or
; ~1/10 of a pixel.  Therefore, only a translation, and not a
; stretching, is used to match up the standard.)  Meyer standards have
; been found to be up to ~10 AA off.
;
; Program also degrades resolution of standards from R=3000 to
; (currently) R=1600 prior to dividing.
;
; Produces fits file of divided spectrum and PS plots of fit and
; result.
;
; INPUTS
;       specfile   - Filename of (preferably normalized) 1D stellar
;                    spectrum.  Header must contain wavelength
;                    information.
;       spectype   - Spectral type.  Currently accepted values are
;                    'A0V', 'F3V', 'F4V', 'F5V', 'F6V', and 'G1V'.
;      
; OUTPUTS
;       outfile    - Root filename for divided spectrum.  '.ps' and
;                    '.fits' will be appended.
;
; HISTORY
;       Written 6/24/2004 GDB
;----------------------------------------------------------------------

if (n_params() lt 3) then begin
   print,'CALLING SEQUENCE: nirspec_divstd,specfile,spectype,outfile'
   print," Accepted values for spectype: 'A0V', 'F3V', 'F4V', 'F5V', 'F6V', and 'G1V'"
   return
endif

;;; Read in user stellar spectrum
spec = rflam(specfile,hspec,wspec)

;;; Read in normalized std fluxes
stddir = '/scr2/gdb/keck/nirspec/stds/hband/'  ; On tallis
case spectype of
    'A0V': stdfile = 'HR7001.101.av.azc.h.nm.xc.fit'
    'F3V': stdfile = 'HR1279.55.av.azc.h.nm.xc.fit'
    'F4V': stdfile = 'HR1279_HR2943_mean.fits'
    'F5V': stdfile = 'HR2943.67.av.azc.h.nm.xc.fit'
    'F6V': stdfile = 'HR1538.11.av.azc.h.nm.xc.fit'
    'G1V': stdfile = 'HR483.85.av.azc.h.nm.xc.fit'
    else : begin
              print,'Unknown spectral type (use capital letters).'
              return
           end
endcase
std = rflam(stddir+stdfile,hstd,wstd,/invcm)

;;; Degrade resolution of standard from R=3000 to R=1600.  Meyer
;;; pixels are ~2 pixels/res element (2.66 AA/pix at 16000 AA).
kerfwhm = 3.2  ; standard pixels
ker   = gaussker(fwhm=kerfwhm)
smstd = convol(std,ker,/edge_truncate)

;;; Interpolate input spec and standard onto a fine wavelength grid
;;; and compute crosscorrelation.  Use a region of the spectrum where
;;; stellar features should domnidate over atmospheric absorption.
wlo = 1.60e4 > max([min(wspec),min(wstd)])
whi = 1.75e4 < min([max(wspec),max(wstd)])
gridstep = 0.1  ; angstroms
ngridsteps = (whi-wlo)/gridstep + 1
wgrid = wlo + gridstep*findgen(ngridsteps)
specgrid  = interpol(spec,wspec,wgrid)
smstdgrid = interpol(smstd,wstd,wgrid)
lag = findgen(3001) - 1500   ; Assume offset is less than 30 angstroms
cc  = c_correlate(smstdgrid,specgrid,lag)
bestlag = lag(where(cc eq max(cc)))
if (n_elements(bestlag) gt 1) then bestlag = mean(bestlag) else $
   bestlag = bestlag(0)
woffset = gridstep*bestlag

;;; Interpolate offset smoothed standard onto input spec wavelengths.
;;; Set values in standard at lambda ~< 15200 to unity.
stdmatch = interpol(smstd,wstd+woffset,wspec)
ls = where(wspec lt 15215+woffset)
if (total(ls) ne -1) then stdmatch(ls) = 1

;;; Divide stellar spec by standard
divspec = spec/stdmatch

;;; Output PS plots
currentmulti = !p.multi
!p.multi = [0,1,2]
ps_openfile,outfile
plot,wspec,spec,xstyle=1,yrange=[0,1.2],title=outfile+' - '+spectype
oplot,wspec,stdmatch,color=100
oplot,wspec,1+0*wspec,linestyle=1
plot,wspec,divspec,xstyle=1,yrange=[0,1.2],xtitle='wavelength (angstroms)'
oplot,wspec,1+0*wspec,linestyle=1
ps_closefile
!p.multi = currentmulti

;;; Output fits file
hout = hspec
sxaddpar,hout,'STDTYPE',spectype,'Standard star spectral type'
sxaddpar,hout,'STDSTAR',stdfile,'Filename of standard star'
writefits,outfile+'.fits',divspec,hout

return
end

Pro nirspec_flatfixer,flatlist,wavearr,slitarr,combflat,combflatvar,$
                      fixgain=fixgain,gain=gain,darklist=darklist,illum=illum

;+-------------------------------------------------------------------
;
; NIRSPEC_FLATFIXER    5/2004
;
; Fix and combine NIRSPEC longslit flats.
;
; In each flat, streaks caused by scratches or dust on the dewar
; window are identified and removed.  User identifies streaks with the
; cursor, which are then fit with a Gaussian profile and divided out.
;
; The program also optionally dark subtracts, adjusts the relative
; gains of the blue and red amplifiers for the illuminated areas of
; the detector, and normalizes each flat by the slit illumination
; function and the spectral shape.
;
; Processed flats are median combined, and the formal variance in the
; combined flat is computed assuming that the dominant source of
; error is Poisson noise.
;
; INPUTS
;       flatlist  - String, name of single-column ASCII file containing
;                   the FITS filenames of flats to be fixed and combined.
;                   Flatfield filenames should end in '.fits'.
;       wavearr   - 2D array or fits filename, same size as flats, where
;                   the pixel values are uncalibrated wavelengths.
;       slitarr   - 2D array or fits filename, same size as flats, where
;                   the pixel values are position along the slit.
;
; KEYWORDS
;       fixgain   - If set, divide the bluer half of the flat by the
;                   ratio of the blue and red gains.  Note that this
;                   accounts for the relative gains of only the two
;                   illuminated detector quadrants for longslit data.
;       gain      - Two-element array, [bluegain, redgain], giving the
;                   gains in the bluer and redder illuminated
;                   quadrants.  Default is [5*1.0343,5], which was
;                   measured from November 2003 1.5 sec x 7 coadd
;                   flats.  Most important is to get the ratio between
;                   the two gains correct.  Non-illuminated quadrants
;                   are ignored.
;       darklist  - String, name of single-column ASCII file containing
;                   the FITS filenames of dark exposures with the same
;                   exposure time as the flats.  If the same number of
;                   darks as flats are listed, then the darks are subtracted
;                   in the order listed.  If only one dark is listed, then
;                   it will be subtracted from all flat frames.
;       illum     - Fit and correct for the illumination function, both
;                   in the wavelength direction and slit direction.
;
;  OUTPUTS
;       For each file <flatfile>.fits listed in the input flatlist, the
;       following files will be written:
;
;          <flatfile>_model.fits - 2D model of the flat field. If the
;                                  /iilum keyword is not set, the
;                                  model will be all 1's, except for
;                                  regions where defects have been
;                                  removed.
;          <flatfile>_fix.fits   - 2D flat field with streaks removed,
;                                  and optionally dark-subtracted and
;                                  normalized such that
;                                  fixedflat = (flat-dark) / model .
;
;       combflat    - Median-combined flat field.  May either be a
;                     variable or a fits filename to which array is to
;                     be written.
;       combflatvar - Formal variance in the combined flat.  May
;                     either be a variable or a fits filename to which
;                     array is to be written.
;
;
; HISTORY
;       Written 5/12/2004 GDB
;       Converted from FLATFIXER 4/30/2005 GDB
;---------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE:  '
   print,' nirspec_flatfixer,flatlist,wavearr,slitarr,combflat,combflatvar,'
   print,'                   /fixgain,gain=gain,darklist=darklist,/illum'
   return
endif

currentpmulti = !p.multi

;;; Read in flat filenames
templ = {   $
       version: 1.0, $
       datastart: 0L, $
       delimiter: ' ', $
       missingvalue: -99.0, $
       commentsymbol: '#', $
       fieldcount: 1L, $
       fieldtypes: [7L], $
       fieldnames: ['name'], $
       fieldlocations: [6L], $
       fieldgroups: [0L] $
     }
list     = read_ascii(flatlist,template=templ)
flatfile = list.name
nflat    = n_elements(flatfile)

;;; Optional - Read in dark filenames
if keyword_set(darklist) then begin
   list     = read_ascii(darklist,template=templ)
   darkfile = list.name
   ndark    = n_elements(darkfile)
   if not(ndark eq nflat or ndark eq 1) then begin
      print,"'darklist' must either contain either a single filename or"
      print,"   the same number of filenames as 'flatlist'."
      return
   endif
   darksub = 1
endif else darksub = 0

;;; Read in wave and slit grids.

if (size(wavearr,/tname) eq 'STRING') then begin
   checkexists = findfile(wavearr,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Wavelength file not found.'
      return
   endif else begin
      wavegrid = readfits(wavearr)
   endelse
endif else wavegrid = wavearr

if (size(slitarr,/tname) eq 'STRING') then begin
   checkexists = findfile(slitarr,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Slit position file not found.'
      return
   endif else begin
      slitgrid = readfits(slitarr)
   endelse
endif else slitgrid = slitarr

;;;
;;; Process each flat.
;;;

; Arrays to store processed flats.
bigflat    = fltarr(1024,1024,nflat)
bigflatvar = fltarr(1024,1024,nflat)

; User-specified edges
xedge = fltarr(4)

for k=0,nflat-1 do begin

   ;;; Read in flat
   tmp_flatfile = flatfile(k)
   checkexists  = findfile(tmp_flatfile,count=nfiles)
   if (nfiles eq 0) then begin
      print,'Flat file '+tmp_flatfile+' not found.'
      return
   endif else begin
      flatfix = readfits(tmp_flatfile,hflatfix)
      sxaddpar,hflatfix,'ORIGFLAT',tmp_flatfile,$
         'Filename of original flat exposure'
   endelse

   ; Promote to floating point
   flatfix   = 1. * flatfix       ; Fixed-up flat field
   flatmodel = 1 + 0. * flatfix   ; Model of flat-field

   ; Add header information
   sxaddpar,hflatfix,'WAVEGRID',wavearr,'Filename of wavelength grid'
   sxaddpar,hflatfix,'SLITGRID',slitarr,'Filename of slit-position grid'

   ;;; Dark-subtract (optional)
   if darksub then begin
      if (ndark gt 1) then tmp_darkfile = darkfile(k) else $
         tmp_darkfile = darkfile(0)
      checkexists = findfile(tmp_darkfile,count=nfiles)
      if (nfiles eq 0) then begin
         print,'Dark file '+tmp_darkfile+' not found.'
         return
      endif else begin
         flatfix = flatfix - readfits(tmp_darkfile)
         sxaddpar,hflatfix,'DARK',tmp_darkfile,$
                  'Dark frame used for dark subtraction'
      endelse
   endif
   
   ; Roughly define illuminated regions of the array
   exposed    = where(slitgrid ge 500 and slitgrid le 800)
   notexposed = where(slitgrid lt 500 or slitgrid gt 800)
   
   ;;; Correct for difference in blue vs. red gain (optional)

   if keyword_set(gain) then begin
      bluegain = gain(0)
      redgain  = gain(1)
   endif else begin
      bluegain = 5. * 1.0343  ; Measured Nov 2003
      redgain  = 5.
   endelse
   rel_gain = bluegain / redgain

   if keyword_set(fixgain) then begin
      flatfix(*,0:511) = flatfix(*,0:511) / rel_gain
      if keyword_set(illum) then $
         flatmodel(*,0:511) = flatmodel(*,0:511) * rel_gain
      sxaddpar,hflatfix,'BLUEGAIN',bluegain,'Blue quadrant gain.'
      sxaddpar,hflatfix,'REDGAIN',redgain,'Red quadrant gain.'
   endif
   
   ;;; Correct illumination function (optional)
   
   if keyword_set(illum) then begin
   
      print,'Correcting slit illumination...'
   
      ; Select a region in the middle in the slit and produce
      ; a linear fit.  Divide out over entire wavelength range.
      centerls = where(wavegrid ge 360 and wavegrid le 380 and $
                       slitgrid ge 550 and slitgrid le 775)
      centerslit = slitgrid(centerls)
      centerflat = flatfix(centerls)
      !p.multi=[0,1,1]
      chan,13
      plot,centerslit,centerflat,psym=3,title='Slit illumination function',$
         xtitle='slit position',ytitle='counts'
      pcoeffs = poly_fit(centerslit,centerflat,1,$
                         yfit=pfit,yerror=pfiterr)
      tmp_ls=where(abs(centerflat-pfit) le 3*pfiterr)
      pcoeffs = poly_fit(centerslit(tmp_ls),centerflat(tmp_ls),1,$
                         yfit=pfit,yerror=pfiterr)
      oplot,centerslit(tmp_ls),pfit,color=100,thick=2
      tmp_arr = pcoeffs(0) + slitgrid*pcoeffs(1)
      tmp_arr = float(tmp_arr)
      flatfix = flatfix / tmp_arr
      flatmodel = flatmodel * tmp_arr
   
      ; Fit a 1D slowly-varying bspline in the dispersion direction
      ; At present, sampling is fine enough to remove medium-scale
      ; features.
      specls = where(slitgrid ge 650 and slitgrid le 670)
      specwave = wavegrid(specls)
      specflat = flatfix(specls)   
      order = sort(specwave)
      specwave = specwave(order)
      specflat = specflat(order)
      chan,13
      plot,specwave,specflat,psym=3,title='Spectral shape',$
           xtitle='uncalibrated wavelength'
      tmp_ls = where(specwave ge 300 and specwave le 400)
      tmp_stddev = stddev(specflat(tmp_ls))
      tmp_invvar = 0.*specflat
      tmp_invvar(*) = 1/(tmp_stddev^2)
      sset = bspline_iterfit(specwave,specflat,invvar=tmp_invvar,bkspace=25,$
                             yfit=specfit,upper=5,lower=5)
      oplot,specwave,specfit,psym=3,color=100,thick=2
      ; Evaluate bspline fit over entire array and divide out
      tmp_arr = 1. + 0.*tmp_arr
      tmp_arr(exposed) = bspline_valu(wavegrid(exposed),sset)
      flatfix = flatfix / tmp_arr   
      flatmodel = flatmodel * tmp_arr
   
      ; Fit a slowly-varying bspline across the left and right
      ; edges of slit in order to boost under-illuminated areas.
      ; Allow user to pick ranges to be fit from the first flat..  
      for i=0,1 do begin
         if (i eq 0) then left=1 else left=0   ; Flag which side
         if left then begin
            slitlo = 470
            slithi = 570
         endif else begin
            slitlo = 750
            slithi = 850
         endelse
         showls = where(wavegrid ge 100 and wavegrid le 300 and $
                        slitgrid ge slitlo and slitgrid le slithi)
         chan,13
         if (k eq 0) then begin
            ttl = 'Slit edge - Click twice to select range'
            print,'Select edge region to be fitted (click twice).'
         endif else begin
            ttl = 'Slit edge'
         endelse
         plot,slitgrid(showls),flatfix(showls),psym=3,title=ttl,$
              xtitle='slit position'
         oplot,[-1e6,1e6],[1,1],linestyle=2,color=100
         tmp_xedge = fltarr(2)
         for j=0,1 do begin
            cursor,tmp_x,tmp_y,/data,/down
            tmp_xedge(j)=tmp_x
            oplot,[tmp_xedge(j),tmp_xedge(j)],[-1e6,1e6],linestyle=2,$
                  color=100
         endfor
         xorder = sort(tmp_xedge) ; Sort entries.
         if left then xedge(0:1) = tmp_xedge(xorder) else $
            xedge(2:3) = tmp_xedge(xorder)
         if left then begin
            xedge_lo = xedge(0)
            xedge_hi = xedge(1)
         endif else begin
            xedge_lo = xedge(2)
            xedge_hi = xedge(3)
         endelse
         edgels = where(wavegrid ge 100 and wavegrid le 300 and $
                        slitgrid ge xedge_lo and slitgrid le xedge_hi)
         edgeslit = slitgrid(edgels)
         edgeflat = flatfix(edgels)
         order    = sort(edgeslit)
         edgeslit = edgeslit(order)
         edgeflat = edgeflat(order)
         chan,13
         plot,edgeslit,edgeflat,psym=3,title='Edge illumination function fit',$
              xtitle='slit position'
         tmp_invvar = 0.*edgeflat
         tmp_invvar(*) = 1/(tmp_stddev^2)
         sset = bspline_iterfit(edgeslit,edgeflat,invvar=tmp_invvar,bkspace=3,$
                             yfit=edgefit,upper=5,lower=5)
         oplot,edgeslit,edgefit,color=100,thick=2
         oplot,[-1e6,1e6],[1,1],linestyle=2,color=100
         ; Only correct the edge up until the fit equals 1
         if left then begin
            pixlimit = min(where(edgefit ge 1))
            if (pixlimit eq -1) then pixlimit = n_elements(edgeslit)-1
            edgelimit = edgeslit(pixlimit)
            entireedge = where(slitgrid ge xedge(0) and slitgrid le edgelimit)
         endif else begin
            pixlimit = max(where(edgefit ge 1))
            if (pixlimit eq -1) then pixlimit = 0
            edgelimit = edgeslit(pixlimit)
            entireedge = where(slitgrid ge edgelimit and slitgrid le xedge(3))
         endelse
         tmp_arr = 1. + 0.*tmp_arr
         tmp_arr(entireedge) = bspline_valu(slitgrid(entireedge),sset)
         flatfix = flatfix / tmp_arr
         flatmodel = flatmodel * tmp_arr
      endfor
      ; Set region beyond edge in flatmodel to 1
      covered_ls = where(slitgrid ge xedge(0) and slitgrid le xedge(3),$
                         complement=beyondedge)
      flatmodel(beyondedge)=1.
   
   endif
   
   ;;; Enter into a loop to indentify and remove defects
   
   ; Display range
   dispmin=0.8*median(flatfix(650:750,450:550))
   dispmax=1.2*median(flatfix(650:750,450:550))

   ; Displayed region
   device,get_screen_size=ssz
   screen_xsz = ssz(0)
   screen_ysz = ssz(1)
   margin     = 50
   disp_xlo   = round(512 - screen_xsz/2. + margin) > 0
   disp_xhi   = round(512 + screen_xsz/2. - margin) < 1023
   disp_ylo   = round(512 - screen_ysz/2. + margin) > 0
   disp_yhi   = round(512 + screen_ysz/2. - margin) < 1023
   
   ; Select first defect
   
   print,'Left click on a defect to remove.  Right click to exit'
   chan,10
   doit,flatfix(disp_xlo:disp_xhi,disp_ylo:disp_yhi),dispmin,dispmax
   cursor,x,y,/device,/down
   x = x + disp_xlo
   y = y + disp_ylo
   
   while(!mouse.button ne 4) do begin
   
      ; Extract a long section in wavelength and slit position
      ; around the selected point
      waveselect = wavegrid(x,y)
      slitselect = slitgrid(x,y)
      region  = where(wavegrid ge waveselect-50 and wavegrid le waveselect+50 $
                   and slitgrid ge slitselect-7 and slitgrid le slitselect+7)
      regionwave = wavegrid(region)
      regionslit = slitgrid(region)
      regionflat = flatfix(region)
   
      ;chan,11
      ;!p.multi=[0,1,2]
      ;plot,regionwave,regionflat,psym=3
   
      ; Fit a low-order polynomial in wavelength and divide out
      pcoeffs = poly_fit(regionwave,regionflat,2,yfit=pfit,yerror=pfiterr)
      tmp_ls  = where(abs(regionflat-pfit) le 3*pfiterr) ; Reject outliers
      pcoeffs = poly_fit(regionwave(tmp_ls),regionflat(tmp_ls),2,$
                         yfit=pfit,yerror=pfiterr)
      ;oplot,regionwave(tmp_ls),pfit,psym=3
      for i=long(0),long(n_elements(region)-1) do $
         regionflat(i) = regionflat(i) / ( pcoeffs(0) + $
                                           pcoeffs(1)*regionwave(i) + $
                                           pcoeffs(2)*regionwave(i)^2 )
      ;plot,regionwave,regionflat,psym=3
   
      ;chan,12
      ;plot,regionslit,regionflat,psym=3
   
      ; Fit a Gaussian in slit position
      gfit = gaussfit(regionslit,regionflat,gterms,nterms=4,$
                      estimates=[-0.02,slitselect,1,1],yerror=gfiterr)
      tmp_ls = where(abs(regionflat-gfit) le 3*gfiterr)
      gfit = gaussfit(regionslit(tmp_ls),regionflat(tmp_ls),gterms,nterms=4,$
                      estimates=gterms,yerror=gfiterr)
      ;oplot,regionslit(tmp_ls),gfit,psym=3
      ;plot,regionslit(tmp_ls),regionflat(tmp_ls)/gfit,psym=3
      
      ; Evaluate fit
      goodfit = 1   ; Default is an acceptable fit
      gterms(2) = abs(gterms(2))
      if (gterms(2) lt 0. or $
          gterms(2) gt 10. or $
          (keyword_set(illum) and abs(gterms(0)) gt 0.1)) then goodfit = 0
   
      ; Divide through by Gaussian fit over entire flat
      if goodfit then begin
         fixls = where( slitgrid ge gterms(1)-5*gterms(2) and $ 
                          slitgrid le gterms(1)+5*gterms(2) )
         for i=long(0),long(n_elements(fixls)-1) do begin
            slitval = slitgrid(fixls(i))
            gval = 1 + gterms(0)*exp(-0.5*((slitval-gterms(1))/gterms(2))^2)
            flatfix(fixls(i)) = flatfix(fixls(i)) / gval
            flatmodel(fixls(i)) = flatmodel(fixls(i)) * gval
         endfor
         ; Print out parameters of defect removed
         print,gterms(0),gterms(1),gterms(2)
      endif else begin
         print,'Bad fit - no defect removed.'
      endelse
   
      chan,10
      doit,flatfix(disp_xlo:disp_xhi,disp_ylo:disp_yhi),dispmin,dispmax
   
      ; Select next defect or exit
      print,'Left click on a defect to remove.  Right click to exit'
      chan,10
      cursor,x,y,/device,/down
      x = x + disp_xlo
      y = y + disp_ylo
   
   endwhile

   ;;; Estimate variance in fixed flat
   if keyword_set(illum) then begin
      ; Variance in the normalized counts
      flatfixvar = 1. / flatmodel
   endif else begin
      ; Variance in the raw counts
      flatfixvar = flatfix 
   endelse
   flatfixvar(*,0:511)    = flatfixvar(*,0:511) / bluegain
   flatfixvar(*,512:1023) = flatfixvar(*,512:1023) / redgain

   ;;; Write files for the individual flat
   lastdot = strpos(tmp_flatfile,'.',/reverse_search)
   if (lastdot eq -1) then root = tmp_flatfile else $
      root = strmid(tmp_flatfile,0,lastdot)
   writefits,root+'_fix.fits',flatfix,hflatfix
   writefits,root+'_model.fits',flatmodel,hflatfix

   ;;; Flag non-illuminated pixels
   flatfix(beyondedge) = -99999

   ;;; Store fixed flat 
   bigflat(*,*,k)    = flatfix
   bigflatvar(*,*,k) = flatfixvar

   ;;; Write files for the individual flat
   lastdot = strpos(tmp_flatfile,'.',/reverse_search)
   if (lastdot eq -1) then root = tmp_flatfile else $
      root = strmid(tmp_flatfile,0,lastdot)
   writefits,root+'_fix.fits',flatfix,hflatfix
   writefits,root+'_model.fits',flatmodel,hflatfix

endfor

;;; Combine flats

print,'Median combining flats...'
if (nflat eq 1) then begin
   medflat    = flatfix
   medflatvar = flatfixvar
endif else begin
   medflat    = fltarr(1024,1024) + 0.001
   medflatvar = fltarr(1024,1024) + 10000.
   for i=470,930 do begin
      for j=0,1023 do begin
         ls = where(bigflat(i,j,*) ne -99999,n_ls)
         if (n_ls ne 0) then begin
            medflat(i,j) = median(bigflat(i,j,ls),/even)
            ; Variance in the median
            medflatvar(i,j) = (!pi/2) * total(bigflatvar(i,j,ls)) / (n_ls^2)
         endif
      endfor
   endfor
endelse

;;; Set output arrays or write fits files

hcombflat = hflatfix
sxaddhist,'NIRSPEC_FLATFIXER: Median combination of the following:',hcombflat
for i=0,nflat-1 do sxaddhist,'                   '+flatfile(i),hcombflat

if (size(combflat,/tname) eq 'STRING') then begin
   print,'Writing '+combflat
   writefits,combflat,medflat,hcombflat
endif else combflat = medflat

if (size(combflatvar,/tname) eq 'STRING') then begin
   print,'Writing '+combflatvar
   writefits,combflatvar,medflatvar,hcombflatvar
endif else combflatvar = medflatvar

print,'Done.'

!p.multi=currentpmulti

return
end

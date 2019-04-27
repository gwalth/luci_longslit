LUCIFER data reduction package is based on
G Becker's NIRSPEC spectral data reduction.

Fuyan Bian@Steward Obervatory 2010

pipeline
1. Creat slit position array
IDL> luci_slitgrid, 'trace_object.list', 'slitgrid_J.fits'
trace_object.list is file contain the name of the standard stars.
slitgrid_J.fits is the output slit position arrray


2. Create uncalibrated wavelength array
IDL> luci_wavegrid,  'lucispec.fits', 'slitgrid_J.fits', 'uncalib_wavegrid_J.fits'
slitgrid_J.fits is the slit position array generated in step 1
and lucispec.fits is either science spectra or lamps spectra
output 

3. Normalize and combine flat fields.


IDL> luci_flatfixer,'flats_J.list','uncalib_wavegrid_J.fits',$
IDL>    'slitgrid_J.fits','normflat_J.fits','normflatvar_J.fits',$
IDL>    darklist='flatdarks_J.list',/illum
slitgrid_J.fits and uncalib_wavegrid_H.fits are generated from step 1 and 2.




4.Combine dark frames with median
As we use the a-b images to subtract the sky backrgound, 
therefore the dark is not crucial.

 IDL> median_combine, 'darks.list', 'median_dark_300s.fits'
darks.list is the list of the dark frames.

median_dark_300s is the output of the median dark image.


5. run the longslit_reduce for the first time to generate
the wavelength calibration model file for further input


longslit_reduce,obj1+'.fits',wave,slit,nosky,noskyvar,wavefit,$
skymodel,/luci, skyref=obj2+'.fits',$
darkframe=dark,texp=300, $
flatfield=flat,flatvar=flatvar,objcoadds=1,$
outfile=obj1,/writewave, /nobias, /nodarksub, /noobject


6. run the longslit_reduce to get the sky subtracted 2d images

readcol, 'objlist', obj, ref,format='a,a'
;readcol, 'qsolist', obj, ref, format='a,a'
;for k = 17,n_elements(obj)-1 do begin
for k = 0, n_elements(obj)-1 do begin
wave    = 'uncalib_wavegrid_J.fits'
slit    = 'slitgrid_J.fits'
flat    = 'normflat_J.fits'
flatvar = 'normflatvar_J.fits'
dark    = 'median_dark_300s.fits'
obj1    =  obj[k]                ;'luci01'     ; Frame to be reduced, here w/o the '.fits'
obj2    =  ref[k]                ;'luci02'     ; Nod exposure
wave_model = 'specJ_wav.dat'
longslit_reduce,obj1+'.fits',wave,slit,nosky,noskyvar,wavefit,$
skymodel,/luci, skyref=obj2+'.fits',$
darkframe=dark,$
flatfield=flat,flatvar=flatvar,objcoadds=2,$
outfile=obj1,/noobject,/maskall,wavemodel=wave_model,/nobias, $
/nodarksub,/noshow,/set_bkpts,texp=300, rdnoise=12, gn=4.1 

endfor


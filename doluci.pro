pro doluci
luci_slitgrid, 'trace_objects_N.list', 'slitgrid_N.fits'
luci_wavegrid, 'spec1.fits', 'slitgrid_N.fits', 'uncalib_wavegrid_N.fits'
luci_flatfixer,'flats_N.list','uncalib_wavegrid_N.fits',$
   'slitgrid_N.fits','normflat_N.fits','normflatvar_N.fits',$
    darklist='flatdarks_N.list',/illum

obj1 = 'luci.20101103.0041'
obj2 = 'luci.20101103.0042'

wave    = 'uncalib_wavegrid_N.fits'
slit    = 'slitgrid_N.fits'
flat    = 'normflat_N.fits'
flatvar = 'normflatvar_N.fits'
dark    = 'median_dark_300s.fits'

longslit_reduce,obj1+'.fits',wave,slit,nosky,noskyvar,wavefit,$
skymodel,/luci, skyref=obj2+'.fits',$
darkframe=dark,texp=300, $
flatfield=flat,flatvar=flatvar,objcoadds=1,$
outfile=obj1,/writewave, /nobias, /nodarksub, /noobject

;return



readcol, 'objlist', obj, ref,format='a,a'
;readcol, 'qsolist', obj, ref, format='a,a'
;for k = 17,n_elements(obj)-1 do begin
for k = 0, n_elements(obj)-1 do begin
wave    = 'uncalib_wavegrid_N.fits'
slit    = 'slitgrid_N.fits'
flat    = 'normflat_N.fits'
flatvar = 'normflatvar_N.fits'
dark    = 'median_dark_300s.fits'
obj1    =  obj[k]                ;'luci01'     ; Frame to be reduced, here w/o the '.fits'
obj2    =  ref[k]                ;'luci02'     ; Nod exposure
wave_model = 'spec_wav.dat'
longslit_reduce,obj1+'.fits',wave,slit,nosky,noskyvar,wavefit,$
skymodel,/luci, skyref=obj2+'.fits',$
darkframe=dark,$
flatfield=flat,flatvar=flatvar,objcoadds=2,$
outfile=obj1,/noobject,/maskall,wavemodel=wave_model,/nobias, $
/nodarksub,/noshow,/set_bkpts,texp=300, rdnoise=12, gn=4.1 

endfor
end


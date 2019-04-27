Pro nirspec_info,filenames,outfile=outfile

;+----------------------------------------------------------------------
;
; NIRSPEC_INFO      5/2004
;
; Extract useful information from a NIRSPEC fits file header,
; including target name, integration time per coadd, number of
; coadds, airmass, etc..
;
; INPUTS
;        file     - String array of fits filenames.
;
; KEYWORDS
;        outfile  - File to which to direct output (one line per input
;                   fits file)
;
; PROGRAMS CALLED
;        strsetlen
;
; HISTORY
;        Written 5/17/2004 GDB
;-----------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE:  nirspec_info,filenames,outfile=outfile'
   return
endif

; Output to a file
tofile = 0    ; Default is standard output
if keyword_set(outfile) then begin
   openw,1,outfile
   tofile = 1
endif

;;; Loop through input files and extract header information

filenames = [filenames]
nfiles = n_elements(filenames)

for i=0,nfiles-1 do begin
   
   ; Read header
   file = filenames(i)
   sxhread,file,hdr

   ; Determine whether this is an object or calibration exposure.
   calmpos = sxpar(hdr,'CALMPOS')   ; Cal. mirror position
   calcpos = sxpar(hdr,'CALCPOS')   ; Cal. cover position
   calppos = sxpar(hdr,'CALPPOS')   ; Cal. pinhole position
   calibration = 0 ; Default is not a calibration exposure
   if (calmpos eq 1 or calcpos eq 1 or calppos eq 1) then $
      calibration = 1

   ; If a calibration exposure, determine which type, otherwise
   ; get target name
   objname = ''
   if (calibration) then begin
       if (calppos eq 1) then objname = objname+'Pinhole'
       if (sxpar(hdr,'NEON') eq 1) then objname = objname+'Ne'
       if (sxpar(hdr,'ARGON') eq 1) then objname = objname+'Ar'
       if (sxpar(hdr,'KRYPTON') eq 1) then objname = objname+'Kr'
       if (sxpar(hdr,'XENON') eq 1) then objname = objname+'Xe'
       if (sxpar(hdr,'ETALON') eq 1) then objname = objname+'Etalon'
       if (sxpar(hdr,'Flat') eq 1) then objname = objname+'Flat'
       if (objname eq '') then objname = 'Dark'
   endif else begin
       objname = sxpar(hdr,'TARGNAME')
   endelse
   ; Make object name 18 cahracters long
   objname = strsetlen(objname,18)

   ; Get additional keywords
   filname  = strsetlen(sxpar(hdr,'FILNAME'),12)
   slitname = strsetlen(sxpar(hdr,'SLITNAME'),11)
   itime    = strsetlen(sxpar(hdr,'ITIME'),10)
   coadds   = strsetlen(sxpar(hdr,'COADDS'),4)
   sampmode = sxpar(hdr,'SAMPMODE')
   case sampmode of
      1: sampmode = strsetlen('Single',7)
      2: sampmode = strsetlen('CDS   ',7)
      3: sampmode = strsetlen('MCDS  ',7)
   endcase
   multispe = strsetlen(sxpar(hdr,'MULTISPE'),5) ; Number of multiple reads
   echlpos  = strsetlen(sxpar(hdr,'ECHLPOS'),10)
   disppos  = strsetlen(sxpar(hdr,'DISPPOS'),10)
   utc      = strsetlen(sxpar(hdr,'UTC'),14)
   airmass  = strsetlen(sxpar(hdr,'AIRMASS'),10)

   ; Print out results
   fileprint = strsetlen(file,18)
   entry = fileprint+objname+filname+slitname+itime+coadds+sampmode+$
               multispe+echlpos+disppos+utc+airmass
   if (tofile) then printf,1,entry else print,entry

endfor

if (tofile) then close,1

return
end



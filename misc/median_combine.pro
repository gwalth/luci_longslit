Pro median_combine,list,comb

;+-------------------------------------------------------------------------
;
; MEDIAN_COMBINE       May 2005 GDB
;
; Median combine 2D FITS files from a list.
;
; INPUTS
;        list      - String, name of single-column ASCII file containig
;                    the names of 2D FITS files to be median combined.
; 
; OUTPUTS
;        comb      - String or variable name, 2D median combined array.
;                    If a string, the array will be written to a FITS
;                    file with that name.
;
; HISTORY
;        Written 5/1/2005 GDB
;--------------------------------------------------------------------------

if (n_params() lt 2) then begin
   print,'CALLING SEQUENCE: median_combine,list,comb'
   return
endif

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
file = read_ascii(list,template=templ)

nfiles = n_elements(file.name)

for i=0,nfiles-1 do begin
   tmp_filename = file.name(i)
   tmp_arr = readfits(tmp_filename,hdr)
   if (i eq 0) then begin
      sz = size(tmp_arr)
      case sz(3) of
          2:    bigarr = intarr(sz(1),sz(2),nfiles)
          3:    bigarr = lonarr(sz(1),sz(2),nfiles)
          4:    bigarr = fltarr(sz(1),sz(2),nfiles)
          5:    bigarr = dblarr(sz(1),sz(2),nfiles)
          12:   bigarr = uintarr(sz(1),sz(2),nfiles)
          else: bigarr = lonarr(sz(1),sz(2),nfiles)
      endcase
      hcomb = hdr
      sxaddhist,'Median combination of the following files:',hcomb
   endif
   bigarr(*,*,i) = tmp_arr
   sxaddhist,'   '+tmp_filename,hcomb
endfor

combarr = median(bigarr,dim=3,/even)

if (size(comb,/tname) eq 'STRING') then begin
   print,'Writing '+comb
   writefits,comb,combarr,hcomb
endif else comb = combarr

return
end

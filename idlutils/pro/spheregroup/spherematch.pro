;+
; NAME:
;   spherematch
;
; PURPOSE:
;   Take two sets of ra/dec coords and efficiently match them. It 
;   returns all matches between the sets of objects, and the distances
;   between the objects. The matches are returned sorted by increasing
;   distance. A parameter "maxmatch" can be set to the number of matches 
;   returned for each object in either list. Thus, maxmatch=1 (the default)
;   returns the closest possible set of matches. maxmatch=0 means to
;   return all matches
;
; CALLING SEQUENCE:
;   spherematch( ra1, dec1, ra2, dec2, matchlength, match1, match2, $
;     		       distance12, [maxmatch=] )
;
; INPUTS:
;   ra1         - ra coordinates in degrees (N-dimensional array)
;   dec1        - dec coordinates in degrees (N-dimensional array)
;   ra2         - ra coordinates in degrees (N-dimensional array)
;   dec2        - dec coordinates in degrees (N-dimensional array)
;   matchlength - distance which defines a match (degrees)
;
; OPTIONAL INPUTS:
;   maxmatch    - Return only maxmatch matches for each object, at
;                 most. Defaults to maxmatch=1 (only the closest
;                 match for each object). maxmatch=0 returns all
;                 matches.
;
; OUTPUTS:
;   match1     - List of indices of matches in list 1 
;   match2     - List of indices of matches in list 2 
;   distance12 - Distance of matches
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The code breaks the survey region into chunks of size
;   4*matchlength. Matches are then performed by considering only 
;   objects in a given chunk and neighboring chunks. This makes the
;   code fast.
;
;   The matches are returned sorted by distance.
;
; EXAMPLES:
;
; BUGS:
;   Behavior at poles not well tested.
;
; PROCEDURES CALLED:
;   Dynamic link to spherematch.c
;
; REVISION HISTORY:
;   20-Jul-2001  Written by Mike Blanton, Fermiland
;-
;------------------------------------------------------------------------------
pro spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $
                 distance12, maxmatch=maxmatch, chunksize=chunksize

   ; Need at least 3 parameters
   if (N_params() LT 8) then begin
      print, 'Syntax - spherematch, ra1, dec1, ra2, dec2, matchlength, match1, match2, $'
      print, ' distance12, [maxmatch=]'
      return
  endif

   if (n_elements(maxmatch) eq 0) then begin
       maxmatch=1l
   end else begin
       if (maxmatch lt 0l) then begin
           print,'illegal maxmatch value: '+maxmatch 
           return
       endif
   endelse 

   if(NOT keyword_set(chunksize)) then chunksize=max([4.*matchlength,0.1])

   npoints1 = N_elements(ra1)
   if (npoints1 le 0l) then begin
       print, 'Need array with > 0 elements'
       return
   endif
   if (npoints1 ne N_elements(dec1)) then begin
       print, 'ra1 and dec1 must have same length'
       return
   endif
   if (N_elements(ra2) le 0l) then begin
       print, 'Need array with > 0 elements'
       return
   endif
   npoints2=N_elements(ra2)
   if (npoints2 ne N_elements(dec2)) then begin
       print, 'ra2 and dec2 must have same length'
       return
   endif

   if (matchlength le 0l) then begin
       print, 'Need matchlength > 0'
       return
   endif

   ; Call matching software, to get lengths
   soname = filepath('libspheregroup.so', $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
   onmatch=0l
   omatch1=lonarr(1)
   omatch2=lonarr(1)
   odistance12=dblarr(1)
   retval = call_external(soname, 'spherematch', long(npoints1), double(ra1), $
                          double(dec1), long(npoints2), double(ra2), $
                          double(dec2), double(matchlength), $
                          double(chunksize), long(omatch1),long(omatch2), $
                          double(odistance12), long(onmatch))

   ; Call matching software for real
   if (onmatch gt 0l) then begin
       omatch1=lonarr(onmatch)
       omatch2=lonarr(onmatch)
       odistance12=dblarr(onmatch)
       retval = call_external(soname, 'spherematch', long(npoints1), $
                              double(ra1), double(dec1), long(npoints2), $
                              double(ra2), double(dec2), $
                              double(matchlength), double(chunksize), $
                              long(omatch1),long(omatch2), $
                              double(odistance12), long(onmatch))
   endif else begin
       print, 'no matches'
       return
   endelse

   ; Retain only desired matches
   sorted=sort(odistance12)
   if (maxmatch gt 0l) then begin
       gotten1=lonarr(npoints1)
       gotten2=lonarr(npoints2)
       nmatch=0l
       for i = 0l, onmatch-1l do begin
           if ((gotten1[omatch1[sorted[i]]] lt maxmatch) and $
               (gotten2[omatch2[sorted[i]]] lt maxmatch)) then begin
               gotten1[omatch1[sorted[i]]]=gotten1[omatch1[sorted[i]]]+1l
               gotten2[omatch2[sorted[i]]]=gotten2[omatch2[sorted[i]]]+1l
               nmatch=nmatch+1l
           end
       end
       gotten1=lonarr(npoints1)
       gotten2=lonarr(npoints2)
       match1=lonarr(nmatch)
       match2=lonarr(nmatch)
       distance12=dblarr(nmatch)
       nmatch=0l
       for i = 0l, onmatch-1l do begin
           if ((gotten1[omatch1[sorted[i]]] lt maxmatch) and $
               (gotten2[omatch2[sorted[i]]] lt maxmatch)) then begin
               gotten1[omatch1[sorted[i]]]=gotten1[omatch1[sorted[i]]]+1l
               gotten2[omatch2[sorted[i]]]=gotten2[omatch2[sorted[i]]]+1l
               match1[nmatch]=omatch1[sorted[i]]
               match2[nmatch]=omatch2[sorted[i]]
               distance12[nmatch]=odistance12[sorted[i]]
               nmatch=nmatch+1l
           end
       end
   endif else begin
       nmatch=onmatch[sorted]
       match1=omatch1[sorted]
       match2=omatch2[sorted]
       distance12=odistance12[sorted]
   endelse

   return
end
;------------------------------------------------------------------------------

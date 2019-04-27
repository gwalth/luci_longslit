;-----------------------------------------------------------------------
;+
; NAME:
;   mwrfits_chunks
;
; PURPOSE:
;   Write a FITS binary table in chunks.
;
; CALLING SEQUENCE:
;   mwrfits_chunks, input, filename, [ header, chunksize=, /silent, $
;    _EXTRA=KeywordsForMwrfits ]
;
; INPUTS:
;   input     - Structure to be written.
;
; OPTIONAL INPUTS:
;   chunksize - Number of rows to write in each sub-write.
;   silent    - If set, then suppress informational messages.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This procedure can be used instead of MWRFITS to write a large
;   FITS binary table without doubling the memory usage.
;
; BUGS:
;   This routine is written to always create a new file.  It should
;   be easy to modify to allow the creation of a new HDU to an existing file.
;
; PROCEDURES CALLED:
;   fxpar()
;   fxposit()
;   host_to_ieee
;   ieee_to_host
;   mwrfits
;
; INTERNAL SUPPORT PROCEDURES:
;   mwrbin_append
;
; REVISION HISTORY:
;   05-Jun-2002  Written by David Schlegel, Princeton, based upon
;                e-mail from Tom McGlynn.
;-
;-----------------------------------------------------------------------
; The following procedure was written by Tom McGlynn on 05-Jun-2002.
pro mwrbin_append, file, ext, struct

    ; Get the beginning of the extension.
    lun   = fxposit(file, ext)
    point_lun, -lun, header_start

    ; Read the header and find where the data begins
    mrd_hread, lun, header
    point_lun, -lun, data_start

    ; Update the NAXIS2  -- this has a fixed location relative to the start
    n = fxpar(header, "NAXIS2")
    point_lun, lun, header_start+330

    z = string(n+n_elements(struct), format="(i20)")
    writeu, lun, z

    ; Write the new data
    data_size =  n*fxpar(header, "NAXIS1")

    point_lun, lun, data_start+data_size
    
    host_to_ieee, struct
    writeu,lun, struct
    ieee_to_host, struct

    ; Write any needed padding.
    data_size = (n+n_elements(struct))*fxpar(header, "NAXIS1")
    blanks = data_size mod 2880
    if blanks ne 0 then begin
        buf = bytarr(2880-blanks)
        writeu, lun, buf
    endif

    free_lun, lun
end 

;-----------------------------------------------------------------------
pro mwrfits_chunks, input, filename, header, chunksize=chunksize, $
 silent=silent, _EXTRA=KeywordsForMwrfits

   ;----------
   ; Need at least 2 parameters

   if (n_params() LT 2) then begin
      doc_library, 'mwrfits_chunks'
      return
   endif

   if (size(input, /tname) NE 'STRUCT') then begin
      print, 'INPUT must be a structure'
      return 
   endif

   if (size(filename, /tname) NE 'STRING') then begin
      print, 'FILENAME must be a string'
      return 
   endif

   if (keyword_set(header) AND size(header, /tname) NE 'STRING') then begin
      print, 'HEADER must be a string array'
      return 
   endif

   ;----------
   ; Pre-condition to FITS structure to have same-length strings
   ; (for any given tag name) by concatenating spaces.

   ntag = n_tags(input)
   for itag=0, ntag-1 do begin
      if (size(input[0].(itag), /tname) EQ 'STRING') then begin
         taglen = strlen(input.(itag))
         maxlen = max(taglen)
         padspace = string('', format='(a'+string(maxlen)+')')
         input.(itag) = strmid(input.(itag) + padspace, 0, maxlen)
      endif
   endfor

   ;----------
   ; Write the structure in chunks.

   nrow = n_elements(input)
   if (NOT keyword_set(chunksize)) then chunksize = nrow
   nchunk = nrow / chunksize + ((nrow MOD chunksize) NE 0)
   for ichunk=0L, nchunk-1 do begin
      if (NOT keyword_set(silent)) then $
       print,'Writing chunk number ', ichunk

      i1 = ichunk * chunksize
      i2 = ((ichunk+1) * chunksize - 1) < (nrow-1)
      indx = i1 + lindgen(i2-i1+1)

      if (ichunk EQ 0) then begin
         ; Create a new file, and write this structure as HDU #1.
         mwrfits, input[indx], filename, header, /create, $
          silent=silent, _EXTRA=KeywordsForMwrfits
         extnum = 1
         lun = fxposit(filename, extnum)
         point_lun, -lun, header_start
      endif else begin
         mwrbin_append, filename, extnum, input[indx]
      endelse
   endfor

   return
end 
;-----------------------------------------------------------------------

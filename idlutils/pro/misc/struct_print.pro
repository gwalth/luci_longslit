;+
; NAME:
;   struct_print
;
; PURPOSE:
;   Formatted print of a structure to standard out, a file, or an array.
;
; CALLING SEQUENCE:
;   struct_print, struct, [ lun=, filename=, tarray=, /no_head, /html, $
;    fdigit=, ddigit=, alias= ]
;
; INPUTS:
;   struct     - Structure
;
; OPTIONAL INPUTS:
;   filename   - Output file name; open and close this file
;   lun        - LUN number for an output file if one is already open
;   no_head    - Do not print the header lines that label the columns
;   html       - If set, then output as an HTML table
;   fdigit     - Number of digits for type FLOAT numbers; default to 5.
;   ddigit     - Number of digits for type DOUBLE numbers; default to 7.
;   alias      - Set up aliases to convert from the IDL structure
;                to the FITS column name.  The value should be
;                a STRARR(2,*) value where the first element of
;                each pair of values corresponds to a column
;                in the structure and the second is the name
;                to be used in the FITS file.
;                The order of the alias keyword is compatible with
;                use in MWRFITS,MRDFITS.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   tarray     - String array for output
;
; COMMENTS:
;   If neither FILENAME or LUN is set and TARRAY is not returned,
;   then write to the standard output.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;   struct_checktype()
;
; REVISION HISTORY:
;   01-Nov-2000  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
; This code is copied from MWR_CHECKTYPE from within "mwrfits.pro".
function struct_checktype, tag, alias=alias
   if not keyword_set(alias) then return, tag

   sz = size(alias)
   ; 1 or 2 D string array with first dimension of 2
   if (sz[0] eq 1 or sz[1] eq 2) and sz[1] eq 2 and sz[sz[0]+1] eq 7 then begin
      w = where(tag eq alias[0,*])
      if (w[0] eq -1) then begin
         return, tag
      endif else begin
         return, alias[1,w[0]]
      endelse
   endif else begin
      print, 'Warning: Alias values not strarr(2) or strarr(2,*)'
   endelse

   return, tag
end
;------------------------------------------------------------------------------
pro struct_print, struct, filename=filename, lun=lun, tarray=tarray, $
 no_head=no_head, html=html, fdigit=fdigit, ddigit=ddigit, alias=alias

   if (keyword_set(filename)) then $
    openw, lun, filename, /get_lun

   if (NOT keyword_set(lun) AND NOT arg_present(tarray)) then lun = -1

   if (NOT keyword_set(fdigit)) then fdigit = 5
   if (NOT keyword_set(ddigit)) then ddigit = 7

   if (size(struct,/tname) NE 'STRUCT') then return
   nrow = n_elements(struct)
   if (nrow EQ 0) then return

   tags = tag_names(struct)
   ntag = n_elements(tags)

   if (keyword_set(html)) then begin
      htmhdr = '<TABLE BORDER=1 CELLPADDING=1>'
      hdr1 = ''
      hdr2 = '<TR>'
      rowsep = '"<TR>",'
      colsep = '"<TD ALIGN=RIGHT>",'
      hdrsep = '<TH>'
      lastline = '</TABLE>'
   endif else begin
      hdr1 = ''
      hdr2 = ''
      hdrsep = ''
      rowsep = ''
      colsep = ''
   endelse

   ;----------
   ; Construct the header lines and format string

   for itag=0L, ntag-1 do begin
      narr = n_elements(struct[0].(itag))
      for iarr=0L, narr-1 do begin

         if (itag EQ 0 AND iarr EQ 0) then begin
            format = '(' + rowsep
         endif else begin
            hdr1 = hdr1 + ' '
            hdr2 = hdr2 + ' '
            format = format + ',1x,'
         endelse

         thisname = struct_checktype(tags[itag], alias=alias)
         if (narr GT 1) then thisname = thisname + strtrim(string(iarr),2)

         tname = size(struct[0].(itag),/tname)
         if (tname EQ 'BYTE' OR tname EQ 'INT' OR tname EQ 'LONG' $
          OR tname EQ 'LONG64' OR tname EQ 'UINT' OR tname EQ 'ULONG' $
          OR tname EQ 'ULONG64') then begin
            minval = min( struct.(itag)[iarr] )
            maxval = max( struct.(itag)[iarr] )

            nchar = strlen(strtrim(string(minval),2)) $
             > strlen(strtrim(string(maxval),2))
            nchar = nchar > strlen(thisname)
            thiscode = 'I' + strtrim(string(nchar),2)
         endif else if (tname EQ 'FLOAT') then begin
            minval = min( struct.(itag)[iarr] )
            if (minval LT 0) then nchar = fdigit + 7 $
             else nchar = fdigit + 6
            nchar = nchar > strlen(thisname)
            thiscode = 'G' + strtrim(string(nchar),2) + '.' $
             + strtrim(string(fdigit),2)
         endif else if (tname EQ 'DOUBLE') then begin
            minval = min( struct.(itag)[iarr] )
            if (minval LT 0) then nchar = ddigit + 7 $
             else nchar = ddigit + 6
            nchar = nchar > strlen(thisname)
            thiscode = 'G' + strtrim(string(nchar),2) + '.' $
             + strtrim(string(ddigit),2)
         endif else if (tname EQ 'STRING') then begin
            nchar = max(strlen( struct.(itag)[iarr] )) > strlen(thisname)
            thiscode = 'A' + strtrim(string(nchar),2)
         endif else begin
            message, 'Unsupported type code: ' + tname
         endelse

         schar = strtrim(string(nchar),2)
         hdr1 = hdr1 + hdrsep + string(thisname, format='(a' + schar + ')')
         if (NOT keyword_set(html)) then $
          hdr2 = hdr2 + string(replicate('-',nchar), $
           format='(' + schar + 'a)')
         format = format + colsep + thiscode

      endfor
   endfor

   format = format + ')'

   ;----------
   ; Now print one row at a time

   if (keyword_set(lun)) then begin
      if (keyword_set(htmhdr)) then printf, lun, htmhdr
      if (NOT keyword_set(no_head)) then begin
         printf, lun, hdr1
         printf, lun, hdr2
      endif
      for irow=0L, nrow-1 do begin
         printf, lun, struct[irow], format=format
      endfor
      if (keyword_set(lastline)) then printf, lun, lastline
      if (keyword_set(filename)) then close, lun
   endif
   if (arg_present(tarray)) then begin
      tarray = strarr(nrow)
      for irow=0L, nrow-1 do begin
         tarray[irow] = string(struct[irow], format=format)
      endfor
      if (NOT keyword_set(no_head)) then $
       tarray = [hdr1, hdr2, tarray]
      if (keyword_set(htmhdr)) then tarray = [htmhdr, tarray]
      if (keyword_set(lastline)) then tarray = [tarray, lastline]
   endif

   return
end
;------------------------------------------------------------------------------

;+
; NAME:
;   yanny_write
;
; PURPOSE:
;   Write a Yanny parameter file from an IDL structure.
;
; CALLING SEQUENCE:
;   yanny_write, filename, [ pdata, hdr=, enums=, structs=, stnames= ]
;
; INPUTS:
;   filename   - Output file name for Yanny parameter file
;
; OPTIONAL INPUTS:
;   pdata      - Array of pointers to all structures write.  The i-th data
;                structure is then referenced with "*pdata[i]"
;   hdr        - Header lines in Yanny file, which are usually keyword pairs.
;   enums      - All "typedef enum" structures.
;   structs    - All "typedef struct" structures, which define the form
;                for all the PDATA structures.
;   stnames    - Structure names, overriding the IDL structure names.
;                Typically, you will want to use this if a structure was
;                read using the /ANONYMOUS keyword.
;
; OUTPUT:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Read and write variables that are denoted INT in the Yanny file
;   as IDL-type LONG, and LONG as IDL-type LONG64.  This is because
;   Yanny files assume type INT is a 4-byte integer, whereas in IDL
;   that type is only 2-byte.
;
; EXAMPLES:
;   Read a Yanny parameter file, then re-write as a different file with:
;     yanny_read, 'testin.par', pdata, comments=comments
;     yanny_write, 'testout.par', pdata, comments=comments
;
; BUGS:
;   There is no testing that STRUCTS is consistent with PDATA.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Sep-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
pro yanny_write, filename, pdata, hdr=hdr, enums=enums, structs=structs, $
 stnames=stnames

   if (N_params() LT 1) then begin
      print, 'Syntax - yanny_write, filename, [ pdata, hdr=, enums=, structs=, stnames= ]'
      return
   endif

   ; Read and write variables that are denoted INT in the Yanny file
   ; as IDL-type LONG, and LONG as IDL-type LONG64.  This is because
   ; Yanny files assume type INT is a 4-byte integer, whereas in IDL
   ; that type is only 2-byte.
;   tname = ['char', 'short', 'int', 'long', 'long', 'float', 'double']
   tname = ['char', 'short', 'int', 'int', 'long', 'float', 'double']
   idlname = ['STRING', 'BYTE', 'INT', 'LONG', 'LONG64', 'FLOAT', 'DOUBLE']

   get_lun, olun
   openw, olun, filename

   ; Write the header of the Yanny file
   if (keyword_set(hdr)) then begin
      for i=0, N_elements(hdr)-1 do begin
         if (hdr[i] NE '') then printf, olun, hdr[i]
      endfor
      printf, olun, ''
   endif

   ; Write the "typedef enum" lines of the Yanny file
   if (keyword_set(enums)) then begin
      for i=0, N_elements(enums)-1 do begin
         if (enums[i] NE '') then printf, olun, enums[i]
      endfor
      printf, olun, ''
   endif

   ; Write the "typedef struct" lines of the Yanny file
   if (keyword_set(structs)) then begin
      for i=0, N_elements(structs)-1 do begin
         if (structs[i] NE '') then printf, olun, structs[i]
      endfor
      printf, olun, ''
   endif else begin
      ; The "typedef struct" lines were not passed to this routine,
      ; so generate those lines consistent with the data structures passed.
      if (keyword_set(pdata)) then begin
         for idat=0, N_elements(pdata)-1 do begin
            ntag = N_tags( *pdata[idat] )
            tags = strlowcase(tag_names( *pdata[idat] ))
            if (keyword_set(stnames)) then stname1 = strupcase(stnames[idat]) $
             else stname1 = tag_names( *pdata[idat], /structure_name)
            if (stname1 EQ '') then $
             stname1 = 'STRUCT' + strtrim(string(idat+1),2)

            printf, olun, 'typedef struct {'

            for itag=0, ntag-1 do begin          ; Loop through each variable
               tt = size( (*pdata[idat])[0].(itag), /tname )
               dims = size( (*pdata[idat])[0].(itag), /dimens )
               ndim = size( (*pdata[idat])[0].(itag), /n_dimen )
               tagname = tname[(where(idlname EQ tt))[0]]

               sline = ' ' + tagname + ' ' + tags[itag]
               for j=0, ndim-1 do begin
                  sline = sline + '[' + strtrim(string(dims[j]),2) + ']'
               endfor
               if (tagname EQ 'char') then $
                sline = sline + '[' $
                + strtrim(string(max(strlen((*pdata[idat]).(itag)))+1),2) + ']'
               sline = sline + ';'

               printf, olun, sline
            endfor

            printf, olun, '} ' + stname1 + ';'
            printf, olun, ''

         endfor
      endif
   endelse

   ; Write the data in the Yanny file
   if (keyword_set(pdata)) then begin
      for idat=0, N_elements(pdata)-1 do begin
         printf, olun, ''

         ntag = N_tags( *pdata[idat] )
         if (keyword_set(stnames)) then stname1 = stnames[idat] $
          else stname1 = tag_names( *pdata[idat], /structure_name)
         if (stname1 EQ '') then stname1 = 'STRUCT' + strtrim(string(idat+1),2)

         for iel=0, N_elements( *pdata[idat] )-1 do begin ; Loop thru each row

            sline = stname1

            for itag=0, ntag-1 do begin          ; Loop through each variable
               words = (*pdata[idat])[iel].(itag)
               nword = N_elements(words)

               ; If WORDS is type STRING, then check for white-space
               ; in any of its elements.  If there is white space, then
               ; put double-quotes around that element.  Or, use quotes
               ; if the string is completely empty.
               if (size(words,/tname) EQ 'STRING') then begin
                  for iw=0, nword-1 do $
                   if (strpos(words[iw],' ') NE -1 $
                    OR strtrim(words[iw],2) EQ '') then $
                    words[iw] = '"' + words[iw] + '"'
               endif else begin
                  words = string(words)
               endelse

               if (nword EQ 1) then begin
                  sline = sline + ' ' + words
               endif else begin
                  sline = sline + ' {'
                  for i=0, N_elements( (*pdata[idat])[iel].(itag) )-1 do $
                   sline = sline + ' ' + words[i]
                  sline = sline + ' }'
               endelse
            endfor
            sline = strcompress(sline)
            printf, olun, sline
         endfor

      endfor
   endif

   close, olun
   free_lun, olun

   return
end
;------------------------------------------------------------------------------

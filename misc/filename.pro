Function filename,fullpathname,dirname=dirname,noextension=noextension

;+----------------------------------------------------------------------
;
; FILENAME         4/2005
;
; Extract the bare filename from the full path name.
;
; E.g., name = filename('/home/machine/user/myfile.txt')
;       help,name
;       <Expression>    STRING    = 'myfile.txt'
;
; INPUTS
;         fullpathname  - String, filename preceeded by any number
;                         of directory names separated by '/'.
;
; KEYWORDS
;         dirname       - String, directory name, including final '/'.
;                         If no directory name is present, then an
;                         empty string will be returned
;         noextension   - If set, filename will be returned without
;                         an extension (the last dot and any following
;                         characters).
;
; HISTORY
;         Written 4/22/2005 GDB
;         DIRNAME and NOEXTENSION keyword added 11/5/2005 GDB
;-----------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE: name = filename(fullpathname,dirname=dirname,'
   print,'                                  /noextension)'
   return,-1
endif

; Remove leading and trailing spaces
trimpathname = strtrim(fullpathname,2)

; Get total string length
len = strlen(trimpathname)

; Locate last '/'
lastslash = strpos(trimpathname,'/',/reverse_search)  

; Loacte last '.'
lastdot = strpos(trimpathname,'.',/reverse_search)

; Extract filename
if (lastslash eq -1) then filename = trimpathname else $
   filename = strmid(trimpathname,lastslash+1)

; Remove extension (optional)
if keyword_set(noextension) then begin
   lastdot = strpos(filename,'.',/reverse_search)
   if (lastdot eq -1) then filename = filename else $
      filename = strmid(filename,0,lastdot)
endif

; Extract directory name
if (lastslash eq -1) then dirname = '' else $
   dirname = strmid(trimpathname,0,lastslash+1)

return,filename
end

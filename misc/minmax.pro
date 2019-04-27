Function minmax,arr,finite=finite

;+-------------------------------------------------------------------
; 
; MINMAX	6/2001
;
; Returns a 2-element vector giving the the min and max of an array
;
; KEYWORDS
;       finite   - Return the min/max of only the finite elements
;
; HISTORY
;	Written 6/25/2001 GDB
;       FINITE keyword added 8/1/2007 GDB
;---------------------------------------------------------------------

; Check whether array is defined
if (n_elements(arr) ge 1) then begin

   if keyword_set(finite) then begin

      ls = where(finite(arr) eq 1)
      minimum=min(arr(ls))
      maximum=max(arr(ls))

   endif else begin

      minimum=min(arr)
      maximum=max(arr)

   endelse

   return,[minimum,maximum]

endif else begin

   print,'MINMAX: Undefined array'
   return,-1

endelse

end

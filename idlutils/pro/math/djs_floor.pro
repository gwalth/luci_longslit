;------------------------------------------------------------------------------
;+
; NAME:
;   djs_floor
;
; PURPOSE:
;   Return largest integer not greater than xvalue.
;   This is identical to the C library function "floor()".
;
; CALLING SEQUENCE:
;   result = djs_floor(xvalue)
;
; INPUTS:
;   xvalue
;
; OUTPUTS:
;   result
;
; PROCEDURES CALLED:
;   fix()
;
; REVISION HISTORY:
;   Written D. Schlegel, 27 November 1996, Durham
;-
;------------------------------------------------------------------------------
function djs_floor, xvalue
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - result = djs_floor( xvalue )'
      return, -1
   endif

   result = fix(xvalue - fix(xvalue) + 1) + fix(xvalue) - 1

   return, result
end
;------------------------------------------------------------------------------

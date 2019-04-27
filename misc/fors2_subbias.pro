Function fors2_subbias,array,ccd=ccd

;+------------------------------------------------------------------
;
; FORS2_SUBBIAS      10/2006
;
; Compute bias level from pre/overscan region and subtract.
;
; Currently assumes 2x2 binning.
;
; INPUTS
;        array     - Raw 2D array
;
; KEYWORDS
;        ccd       - Scalar, CCD number (1 or 2).  Determines rows
;                    from which bias is calculated.  If not set,
;                    bias is calculated from row 0.
;
; OUTPUTS
;        <result>  - Raw 2D array with bias level subtracted.
;
; HISTORY
;        Written 10/31/2006 GDB
;-------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE: result = fors2_subbias(array,ccd=ccd)'
   return,-1
endif

; Determine whether pre/overscan regions are present
sz  = size(array)
ysz = sz(2)

if (ysz eq 1034) then begin

   ; Determine rows from which to calculate bias level
   if keyword_set(ccd) then begin
      case ccd of
         1:    begin
                  bias_ymin = 0
                  bias_ymax = 3
               end
         2 :   begin
                  bias_ymin = 1030
                  bias_ymax = 1033
               end
         else: begin
                  print,'Unknown CCD number.  Caculating bias from row 0.'
                  bias_ymin = 0
                  bias_ymax = 0
               end
      endcase
   endif else begin
      print,'CCD number not specified.  Caculating bias from row 0.'
      bias_ymin = 0
      bias_ymax = 0
   endelse
   
   ; Compute bias level and subtract
   bias_level = mean(array(*,bias_ymin:bias_ymax))
   print,'Subtracting bias level = '+strtrim(bias_level,2)
   sub_array  = array - bias_level

endif else begin

   print,'Pre/overscan regions not present.  No bias-subtraction done.'
   sub_array = array

endelse

return,sub_array
end

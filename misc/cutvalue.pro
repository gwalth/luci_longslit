Function cutvalue,array,p,incl=incl,above=above

;+---------------------------------------------------
;
; CUTVALUE	1/2003 GDB
;
; Returns the value below which a specified percent
; of array values reside.  For example, if the 
; desired percentile is 50, cutvalue returns the
; median (in the same sense as the IDL intrinsic 
; function MEDIAN, see below).
;
; In an ordered array with N entries, CUTVALUE will 
; return the value of the nth element where n has 
; the smallest value possible such that at least
; (p/100)*N elements have values smaller than the 
; nth element, where p is the percentile.
;
; Note that for array=[1,2] and p=50, the result
; is 2.  This is identical to the intrinsic IDL
; routine MEDIAN, but not the usual definition
; of the median of a set with an even number
; of members.
;
; CALLING SEQUENCE
;	result = cutvalue(array,p [,/incl,/above])
;
; INPUTS
;	array	- Integer, floating, or double array,
;		  need not be ordered
;	p	- Desired percentile of elements that 
;		  lie below the returned value
;
; OUPUTS
;	<result> - The value in the array with at least
;		  p percent of the elements having
;		  lower values.
;
; KEYWORDS
;	incl	- "Inclusive".  If set then the returned
;		  value is that such that at least 
;                 (p/100)*N elements have EQUAL OR lesser
;                 value.  If 'above' keyword is also set
;		  then the returned value is that such that
;		  at least (p/100)*N elements have EQUAL OR
;		  greater value.
;	above	- If set then the returned value is that
;		  value in the array with at least p
;		  percent of the elements having
;		  higher values
;
; HISTORY
;	Written 1/29/2003 GDB
;--------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE: result = cutvalue(array,p [,/incl,/above])'
   return,'Error'
endif

;;; Check that 'p' is a percentile

if (p lt 0 or p gt 100) then begin
   print,'p is a percentile and must be between 0 and 100.'
   return,'Error'
endif

n = n_elements(array)

;;; Create a sorted version of the array.  If 'above' keyword
;;; set, reverse order.

order = sort(array)

sorted = array(order)

if keyword_set(above) then begin
   temp = sorted
   for i=long(0),long(n-1) do sorted(i)=temp(n-1-i)
endif

;;; Return the value below (or equal to) which at least the 
;;; desired percentile of elements fall.

frac = float(p)/100.

index = ceil(frac*n) - 1

if keyword_set(incl) then begin
   if (index lt 0)   then index=0
   if (index gt n-1) then index=n-1
endif else begin
   if (ceil(frac*n) eq frac*n) then index=min([index+1,n-1])
   if (index lt 1 and p ne 0) then index=1
   if (index gt n-1) then index=n-1
endelse


return,sorted(index)

end

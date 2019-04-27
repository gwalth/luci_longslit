;+
; NAME:
;   hogg_usersym
; PURPOSE:
;   make an n-sided plotting point
; USAGE:
;   hogg_usersym, N [,scale=scale,/diagonal]
;   plot, x,y,psym=8
; INPUTS:
;   N           - number of sides on the polygon
; OPTIONAL INPUTS:
;   scale       - linear size
;   _extra      - keywords for usersym (see usersym help page)
; KEYWORDS
;   diagonal    - rotate symbol through 1/2 of 1/N turns
; REVISION HISTORY:
;   2002-04-09  written - Hogg
;-
pro hogg_usersym,N,diagonal=diagonal,scale=scale,_EXTRA=KeywordsForUsersym
  if keyword_set(diagonal) then delta= 0D0 else delta= 5D-1
  if NOT keyword_set(scale) then scale= 1D0
  theta= 2D0*!PI*(dindgen(N+2)+delta)/double(N)
  usersym, scale*sin(theta),-scale*cos(theta),_EXTRA=KeywordsForUsersym
end

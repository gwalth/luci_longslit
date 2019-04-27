Pro  doit,array,llow,hhigh,border=border,width=width

;+-----------------------------------------------------------------------------
;
; DOIT  6/18/2001 GDB
;
; Works like DOITPS in terms of the windowing.  For use in displaying to 
; screen.  Uses CTV instead of TV, therefore willsize window to image size.
;
; ADDITIONAL KEYWORDS
;
; 	border	- draws a border around the image specified width.  
;		  Warning: draws over image.
;	width	- width of border, in pixels.  Default is 2. 
;
; HISTORY
;   Written 6/01 GDB
;   Retro-adapted from DOITPS, itself an adaptation of DOIT written by RWO.
;------------------------------------------------------------------------------

if (n_params(0) eq 0) then begin
  print,'CALLING SEQUENCE: doitps,array,llow,hhigh,/border,width=width'
  return
 endif


topp=!d.n_colors - 2   ; default

low=llow
high=hhigh
   
barray=bytscl(array,low,high,top=topp)

if keyword_set(border) then begin	;;;Draw border

  if keyword_set(width) then w = width  else  w = 2

  sz = size(array)
  
  barray(0:w-1,*)=255
  barray(sz(1)-w:sz(1)-1,*)=255
  barray(*,0:w-1)=255
  barray(*,sz(2)-w:sz(2)-1)=255
  
 endif 

ctv,barray

return
end
Pro getpoints,xpoints,ypoints,device=device,plotx=plotx,ploty=ploty,$
              silent=silent,_extra=_extra

;+------------------------------------------------------------------
;
; GETPOINTS    March 2004
;
; Interactively read points off of a previously drawn x,y plot or 2D
; image and return vectors of selected x,y values.
;
; Use left mouse button to select points, right mouse button to exit.
;
; INPUTS
;       none
;
; KEYWORDS
;       device   - Read device coordinates (pixel coordinates).  
;                  Default is data coords.
;       plotx    - Plot a vertical dashed line through each selected
;                  point.  Not effective if the device keyword is set.
;       ploty    - PLot a horizontal dahsed line through each selected
;                  point.  Not effective if the device keyword is set.
;       silent   - Suppress messages.
;       _extra   - oplot keywords not otherwise used.
;
; OPTIONAL OUTPUTS 
;       xpoints,ypoints - N-element array of selected plot values or 
;                         pixel coordinates, where N is the number of 
;                         points clicked on.
;
; HISTORY
;       Written 3/31/2004 GDB
;       PLOTX,PLOTY keywords added 3/8/2005 GDB
;-------------------------------------------------------------------

if (n_params() eq 0) then begin
   print,'CALLING SEQUENCE: getpoints,xpoints,ypoints,/device,/plotx,/ploty,'
   print,'                            /silent,_extra=oplot_keywords'
   print,' '
   return
endif

if not(keyword_set(silent)) then $
   print,'Use left mouse button to select points, right button to exit.'

!mouse.button = 0
firstpoint = 1        ; Flag to indicate that this is the first point
                      ; to be selected

if keyword_set(device) then begin
   dev=1
   dat=0
endif else begin
   dev=0
   dat=1
endelse

; Default is no lines selected
xpoints = -1
ypoints = -1

while (!mouse.button ne 4) do begin
      cursor,xselect,yselect,/down,data=dat,device=dev
      if (!mouse.button ne 4) then begin
         print,xselect,yselect
         if (firstpoint) then begin
            xpoints = [xselect]
            ypoints = [yselect]
            firstpoint=0
         endif else begin
            xpoints = [xpoints,xselect]
            ypoints = [ypoints,yselect]
         endelse
         ; Plot lines though the current point?
         if (keyword_set(plotx) and not(keyword_set(device))) then begin
            oplot,[xselect,xselect],[-1e29,1e29],linestyle=1,_extra=_extra
            oplot,[xselect,xselect],[1e-29,1e29],linestyle=1,_extra=_extra
         endif
         if (keyword_set(ploty) and not(keyword_set(device))) then begin
            oplot,[-1e29,1e29],[yselect,yselect],linestyle=1,_extra=_extra  
            oplot,[1e-29,1e29],[yselect,yselect],linestyle=1,_extra=_extra  
         endif
      endif
endwhile

return
end

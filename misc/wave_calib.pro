;+---------------------------------------------------------------------
;
; WAVE_CALIB     4/2004
;
; General wavelength calibration utility.  Allows user to identity
; line features from a list and compute fit.
;
; In general, the widget buttons and options are self-explanatory.  To
; identify a line, first select it from the list at right and then
; left-click on the spectrum at a point near the line center.  The
; program will attempt to fit the upper part of the line profile.  To
; remove a line, center-click near it (will have to click 'Compute
; Fit' button to recompute fit).  Right-clicking displays either the
; x-value or fitted wavelength at the cursor position.  The program
; takes only the x-coordinate as mouse input, therefore one can click
; on either window.
;
; NOTE: The documentation assumes that the dispersion direction is
;       the x-direction.
;
; CALLING SEQUENCE
;
;      wave_calib,xin,yin,wave,linefile=linefile,reference=reference,
;              nopersist=nopersist,psout=psout,coeffs=coeffs,xid=xid,
;              waveid=waveid,order=order,dispersion=dispersion,
;              resolution=resolution,xoffset=xoffset,/auto,
;              /rejectnoise,/noshow,psym=psym
;
; INPUTS
;      xin        - 1D array of wavelength-ordered pixels
;      yin        - 1D array of wavelength-ordered intensities, may
;                   be arc lamp or sky exposure
;
; KEYWORDS
;      linefile   - Name of line list file, where the line wavelengths
;                   are in the first column.  See below for default.
;      guess      - Name of ASCII file or array containing estimated
;                   line x-positions.  For an ASCII file, the format
;                   is row = {x wavelength}.  For an array, the array
;                   dimensions should be n x 2, where guess(*,0)=x,
;                   guess(*,1)=wavelength.  If set, the program will
;                   attempt to automatically find these lines at the
;                   specified x-positions.  Note: The lines in guess
;                   must be included in the linefile.
;      reference  - Name of ASCII file or array containing a reference
;                   wavelength fit.  For an ASCII file, the format is
;                   row = {xref yref refwavelengthfit}.  For an array,
;                   the array dimensions should be n x 3, where
;                   reference(*,0)=xref, reference(*,1)=yref,
;                   refrence(*,3)=waveref.  If set, program will
;                   attempt to automatically find all lines within the
;                   present x-range based on their location in the
;                   reference fit.
;      nopersist  - Kill widget after auto-line identification.  Has an
;                   effect only if the reference keyword is set.
;      psout      - File name for output postscript plot of line
;                   ID's and residuals as a function of fitted
;                   wavelength.  Will automatically append
;                   *NOTE: PS output not currently functional!*
;      coeffs     - Array of fitted polynomial coefficients.
;      xid        - Array of x pixel values of identified lines
;      waveid     - Array of wavelength values of identified lines
;      order      - Initial order of polynomial fit.  Default is 3.
;                   If the user changes the polynomial order
;                   interactively, this keyword will contain
;                   the final value on output.
;      dispersion - Scalar approximate dispersion in (physical
;                   units)/(unit of xin), where the physical units are
;                   those in the line list (generally Angstroms) and
;                   the units of xin are generally pixels.  However,
;                   if xin is already calibrated (e.g., for doing a wavelength
;                   shift to match the skylines in a spectrum), then
;                   the dispersion may be = 1.  Affects how the
;                   program identifies and fits lines.  Default is 1.
;                   *NOTE: If setting either the dispersion or
;                   resolution keywords, it is recommended to set both.
;      resolution - Scalar approximate resolution in physical units
;                   used in the line list.  Affects how the program
;                   identifies and fits lines.  Default is
;                   3*dispersion.
;      xoffset    - Scalar bulk offset found between the reference
;                   wavelength solution and the input wavelength solution,
;                   found by crosscorrelating the input instensity with
;                   the reference intensities (positive if the reference
;                   x values are greater than the best-matching input
;                   x values).
;      auto       - Automatically identify lines.  If set without the
;                   reference keyword, xin values should themselves
;                   be in calibrated physical wavelength units, and
;                   the lines in the line list should appear near 
;                   the corresponding value in xin.  This is for small
;                   corrections to the wavelength solution only, such
;                   as making a small shift to match the skylines in 
;                   an individual exposure.
;      rejectnoise- Reject identifications that appear to be noise.  This
;                   keyword should be set only if the lines to be
;                   identified are well separated (e.g., high-resolution
;                   optical data, but not low-resolution IR data).
;      limit      - Scalar, reject lines with peak values above this value.
;                   For rejecting saturated lines.  Default is not to
;                   reject strong lines.
;      noshow     - Suppresses display.  Currently effective only
;                   if reference keyword is set.
;      psym       - Plot symbol to use for spectrum plot.  Default is 3
;                   (points).
;
; OUPUTS
;      wave       - 1D array of fitted wavelengths
;
; HISTORY
;      Written 4/2004 GDB
;      Dispersion, resolution keywords added 2/1/2005 GDB
;      LIMIT keyword added 5/14/2007 GDB
;----------------------------------------------------------------------

;----------------------------------------------------------------------
; Initialize common blocks, line lists, and any other settings
;----------------------------------------------------------------------

Pro wave_calib_initialize,linefile=linefile

; Read in line list

; Template for a single column of double-precision wavelengths
templ = {   $
          version: 1.0, $    
          datastart: 0L, $   
          delimiter: ' ', $        ; Fields separated by white space
          missingvalue: -99.0, $
          commentsymbol: '#', $
          fieldcount: 1L, $
          fieldtypes: 5L, $       ; Double precision
          fieldnames: 'wavelength', $
          fieldlocations: 0L, $
          fieldgroups: [0L] $
        }
if keyword_set(linefile) then line_file = linefile else $
   line_file = '~gdb/calib/nirspec/oh_lists/lowd_ir_ohlines_hband.lst'
linedata  = read_ascii(line_file,template=templ)
linewaves = linedata.wavelength
linenames = strtrim(string(linewaves),2)

nlinestot = n_elements(linenames) ; Total number of lines

lineswithid = fltarr(nlinestot) ; Flags to indicate which lines have been
                                ; marked. Default for all lines is no (0).

; Initialize common blocks

common wc_state, state
common wc_data, xinputarr, yinputarr, wavelengthfit, xplotarr

state = {              $
          show: 1L, $                 ; Flag to display plots
          whichwindow: 0L, $          ; Current display window
          specwindow: 0L, $           ; Spectrum plot window
          fitwindow: 0L, $            ; Wavelength fit window
          spec_xsz: 1000L, $          ; Screen x-size of spec window
          spec_ysz: 400L, $           ; Screen y-size of spec window
          fit_xsz: 1000L, $           ; Screen x-size of fit window
          fit_ysz: 400L, $            ; Screen y-size of fit window
          xmin: min(xinputarr), $     ; Min data x-value (min wavelength)
          xmax: max(xinputarr), $     ; Max data x-value (max wavelength)
          ymin: min(yinputarr), $     ; Min data y-value (min intensity)
          ymax: max(yinputarr), $     ; Max data y-value (max intensity)
          xdispmin: min(xplotarr), $  ; Min display x-value 
          xdispmax: max(xplotarr), $  ; Max display x-value
          xwave: 0L, $                ; Flag to plot vs. fitted wavelengths, default is no
          specydispmin: 1.0*min(yinputarr), $ ; Min display y-value 
          specydispmax: 1.2*min(yinputarr), $ ; Max display y-value
          fitydispmin: 0., $          ; Min display y-value 
          fitydispmax: 1., $          ; Max display y-value
          psym: 3, $                  ; Plot symbol for spectrum plot
          nlinestot: nlinestot, $     ; Total number of lines in list
          linenames: linenames, $     ; Line names
          linewaves: linewaves, $     ; Line wavelengths
          linefitwaves: 0*linewaves, $; Line wavelengths computed from fit
          lineresids: 0*linewaves, $  ; Line wavelength residuals
          lineswithid: 0*linewaves, $ ; Flag that lines are identified
          linesused: 0*linewaves, $   ; Flag that lines are included in fit
          linexvalues: 0*linewaves, $ ; X values of pixels at line centers
          lineyvalues: 0*linewaves, $ ; Y values of pixels at line centers
          fitorder: 3L, $             ; Order of wavelength polynomial fit
          fitcoeffs: dblarr(20),$     ; Wavelength polynomial fit coefficients
          fitcomputed: 0L, $          ; Flag that a wavelength fit has been computed
          fiterror: 0d0, $            ; Standard error between fitted and tabulated wavelengths
          coords: fltarr(3), $        ; Data coordinates of mouse click
          xcoord: 0., $               ; Data x-coordinate of mouse click
          ycoord: 0., $               ; Data y-coordinate of mouse click
          mbutton: '', $              ; Mouse button clicked
          currentwave: 0d0, $         ; Wavelength of current line
          dispersion: 1., $           ; Dispersion in physical units/units of xin
          resolution: 3., $           ; Approximate resolution, in physical units
          searchpixrad: 3., $         ; Distance in pixels around x-coord to search for line
          searchwaverad: 3., $        ; Distance in Angstroms around x-coord to serach for line
          gfitpixrad: 3., $           ; Distance in pixels +/- around line peak to fit Gaussian 
          gfitwaverad: 3., $          ; Distance in Angstroms +/- around line peak to fit Gaussian
          rejectnoise: 0, $           ; Flag to reject lines that appear to be noise
          rejectstrong: 0, $          ; Flag to reject strong (saturated) lines
          limit: 0., $                ; Limit for rejecting strong lines
          ticklength: 0.03, $         ; Tick length as a fraction of y display range
          plotresiduals: 0L, $        ; Flag to plot wavelength residuals 
          createps: 0L, $             ; Flag to produce ps plot on exit, default is no
          psfile: '', $               ; Filename for postscript output
          alldone: 0L, $              ; Flag that fitting is complete 
          base_id: 0L, $              ; Widget ids...
          leftbase_id: 0L, $          ;
          rightbase_id: 0L, $         ;
          specwindow_id: 0L, $        ;
          fitwindow_id: 0L, $         ;
          optionsbase_id: 0L, $       ;
          ymin_text_id: 0L, $         ;
          ymax_text_id: 0L, $         ;
          xmin_text_id: 0L, $         ;
          xmax_text_id: 0L, $         ;
          fullrange_button_id: 0L, $  ;
          xpixels_button_id: 0L, $    ;
          xwave_button_id: 0L, $      ;
          buttonbase_id: 0L, $        ;
          showwavelengths_id: 0L, $   ;
          showresiduals_id: 0L, $     ;
          polyorder_text_id: 0L, $    ;
          computefit_button_id: 0L, $ ;
          findmorelines_button_id: 0L, $ ;
          rejectoutliers_button_id: 0L, $;
          rejectmultipleids_button_id: 0L, $;
          done_button_id: 0L, $       ;
          linelist_id: 0L $           ;
        }

end

;-----------------------------------------------------------------------

Pro wave_calib_window_event,event

common wc_state
common wc_data

widget_control,event.id,get_uvalue=eventval

case eventval of
    'specwindow'  : begin
                       ; Convert from the device coordinates of the cursor
                       ; to the data coordinates.  Since we are using the 
                       ; same x-range for both plots and are drawing them 
                       ; on top of one another, it doesn't matter
                       ; for the x-coordinate which plot was drawn most
                       ; recently.
                       if (event.press gt 0) then begin
                          case event.press of
                             1: state.mbutton = 'left'
                             2: state.mbutton = 'center'
                             4: state.mbutton = 'right'
                             else: print,'Error - unknown mouse button.'
                          endcase
                          wset,state.specwindow
                          state.coords = convert_coord(event.x,event.y,/device,/to_data)
                          state.xcoord = state.coords(0)
                          state.ycoord = state.coords(1)
                          print,state.mbutton+' clicking in spec window at x = ',state.xcoord
                          if (state.mbutton eq 'left') then findline
                          if (state.mbutton eq 'center') then removeline
                       endif
                   end
    'fitwindow'   : begin
                       if (event.press gt 0) then begin
                          case event.press of
                             1: state.mbutton = 'left'
                             2: state.mbutton = 'center'
                             4: state.mbutton = 'right'
                             else: print,'Error - unknown mouse button.'
                          endcase
                          wset,state.fitwindow
                          state.coords = convert_coord(event.x,event.y,/device,/to_data)
                          state.xcoord = state.coords(0)
                          state.ycoord = state.coords(1)
                          print,state.mbutton+' clicking in fit window at x = ',state.xcoord
                          if (state.mbutton eq 'center') then removeline
                       endif
                    end
    'line'        : begin
                       state.currentwave = event.value
                       print,'Line Selected with wavelength ',state.currentwave
                    end
    'ymin_text'   : begin
                       print,'Changing y display min to ',event.value
                       state.specydispmin = event.value
                       ; If new ymin is greater than ymax, wait for
                       ; ymax to be changed before replotting
                       if (state.specydispmin lt state.specydispmax) then begin
                          draw_fitplot
                          draw_specplot,/norescale
                       endif
                    end
    'ymax_text'   : begin
                       print,'Changing y display max to ',event.value
                       state.specydispmax = event.value
                       ; If new ymax is less than ymin, wait for ymin
                       ; to be changed before replotting
                       if (state.specydispmin lt state.specydispmax) then begin
                          draw_fitplot
                          draw_specplot,/norescale
                       endif
                    end
    'xmin_text'   : begin
                       print,'Changing x display min to ',event.value
                       state.xdispmin = event.value
                       ;; If new xmin is greater than xmax, wait for
                       ;; xmax to be changed before replotting
                       ;if (state.xdispmin lt state.xdispmax) then begin
                          draw_fitplot
                          draw_specplot
                       ;endif
                    end
    'xmax_text'   : begin
                       print,'Changing x display max to ',event.value
                       state.xdispmax = event.value
                       ;; If new xmax is less than xmin, wait for xmin
                       ;; to be changed before replotting
                       ;if (state.xdispmin lt state.xdispmax) then begin
                          draw_fitplot
                          draw_specplot
                       ;endif
                    end
    'full_range'  : begin
                       state.xdispmin = min(xplotarr)
                       ;widget_control,state.xmin_text_id,set_value=state.xdispmin
                       state.xdispmax = max(xplotarr)                       
                       ;widget_control,state.xmax_text_id,set_value=state.xdispmax
                       draw_fitplot
                       draw_specplot
                    end
    'x_pixels'    : begin ; Plot vs. pixel (or other input) value
                       if (state.xwave eq 1) then begin ; Not already plotting vs pixels?
                          state.xwave = 0
                          swapxplot
                       endif
                    end
    'x_wave'      : begin ; Plot vs. fitted wavelength
                       if (state.xwave eq 0) then begin ; Not already plotting vs wavelength?
                          tmp_fitrange = max(wavelengthfit)-min(wavelengthfit)
                          if (tmp_fitrange gt 0) then begin
                             state.xwave = 1
                             swapxplot
                          endif else begin
                             print,'No wavelength fit computed'
                          endelse
                       endif
                    end
    'show_wavelengths'  : begin ; Show all selected lines in fit window
                       state.plotresiduals = 0
                       draw_fitplot
                       draw_specplot
                    end
    'show_residuals' : begin ; Show residuals for all selected lines in fit window
                          state.plotresiduals = 1
                          computeresiduals
                          draw_fitplot
                          draw_specplot
                       end
    'polyorder_text' : begin ; Change polynomial order
                          state.fitorder = event.value
                          computefit
                       end
    'compute_fit' : begin
                       print,'Computing fit...'
                       computefit
                    end
    'find_more_lines'    : find_more_lines
    'reject_outliers'    : reject_outliers
    'reject_multiple_ids': reject_multiple_ids
    'done'        : begin
                       state.alldone = 1
                       widget_control,event.top,/destroy
                    end
    else          : print,'Something else happened...'
endcase

; Reset mouse button value
state.mbutton = ''

end

;------------------------------------------------------------------------
; Initialize widgets
;------------------------------------------------------------------------

Pro create_wave_calib_window,group=group

common wc_state

state.base_id = widget_base(/row,title='Wavelength Calibration')

state.leftbase_id  = widget_base(state.base_id,/column,/base_align_right)

; Adjust plot windows to fit screen
device,get_screen_size=screen
state.spec_xsz = state.spec_xsz < (screen(0)-200)
state.spec_ysz = state.spec_ysz < (screen(1)-200)/2.
state.fit_xsz  = state.fit_xsz  < (screen(0)-200)
state.fit_ysz  = state.fit_ysz  < (screen(1)-200)/2.

state.specwindow_id = widget_draw(state.leftbase_id,$
                                  scr_xsize=state.spec_xsz,$
                                  scr_ysize=state.spec_ysz,$
                                  uvalue='specwindow',$
                                  /button_events)

state.fitwindow_id = widget_draw(state.leftbase_id,$
                                 scr_xsize=state.fit_xsz,$
                                 scr_ysize=state.fit_ysz,$
                                 uvalue='fitwindow',$
                                 /button_events)

state.optionsbase_id = widget_base(state.leftbase_id,/row)

state.ymin_text_id = cw_field(state.optionsbase_id, $
                              uvalue = 'ymin_text',/floating,$
                              title = 'Y Min =', $
                              value = state.specydispmin,  $
                              /return_events, $
                              xsize = 12)

state.ymax_text_id = cw_field(state.optionsbase_id, $
                              uvalue = 'ymax_text',/floating, $
                              title = 'Y Max =', $
                              value = state.specydispmax,  $
                              /return_events, $
                              xsize = 12)

state.xmin_text_id = cw_field(state.optionsbase_id, $
                              uvalue = 'xmin_text',/floating,$
                              title = 'X Min =', $
                              value = state.xdispmin,  $
                              /return_events, $
                              xsize = 12)

state.xmax_text_id = cw_field(state.optionsbase_id, $
                              uvalue = 'xmax_text',/floating, $
                              title = 'X Max =', $
                              value = state.xdispmax,  $
                              /return_events, $
                              xsize = 12)

state.fullrange_button_id = widget_button(state.optionsbase_id,$
                                          value='Full X Range',$
                                          uvalue='full_range')

state.xpixels_button_id = widget_button(state.optionsbase_id,$
                                        value='X Pixels',$
                                        uvalue='x_pixels')

state.xwave_button_id = widget_button(state.optionsbase_id,$
                                      value='X Wavelength',$
                                      uvalue='x_wave')

state.buttonbase_id = widget_base(state.leftbase_id,/row)

state.showwavelengths_id = widget_button(state.buttonbase_id,$
                                   value='Show Wavelengths',$
                                   uvalue='show_wavelengths')

state.showresiduals_id = widget_button(state.buttonbase_id,$
                                       value='Show Residuals',$
                                       uvalue='show_residuals')

state.polyorder_text_id = cw_field(state.buttonbase_id, $
                                   uvalue = 'polyorder_text',/long, $
                                   title = 'Polynomial Order =', $
                                   value = state.fitorder,  $
                                   /return_events, $
                                   xsize = 4)

state.findmorelines_button_id = widget_button(state.buttonbase_id,$
                                           value='Compute Fit',$
                                           uvalue='compute_fit')

state.computefit_button_id = widget_button(state.buttonbase_id,$
                                           value='Find More Lines',$
                                           uvalue='find_more_lines')

state.rejectoutliers_button_id = widget_button(state.buttonbase_id,$
                                           value='Reject Outliers',$
                                           uvalue='reject_outliers')

state.rejectmultipleids_button_id = widget_button(state.buttonbase_id,$
                                           value='Reject Multiple IDs',$
                                           uvalue='reject_multiple_ids')

state.done_button_id = widget_button(state.buttonbase_id,value='Done',$
                                     uvalue='done')

state.rightbase_id = widget_base(state.base_id,/column)

; Get geometry of window
base_geometry = widget_info(state.base_id,/geometry)
yscrollsz = base_geometry.scr_ysize - 2*base_geometry.margin - 40

state.linelist_id = cw_bgroup(state.rightbase_id,state.linenames,$
                              /column,/exclusive,/scroll,$
                              y_scroll_size=yscrollsz,uvalue='line',$
                              button_uvalue=state.linewaves,/no_release)

; Realize widgets
widget_control,state.base_id,/realize

; Determine window numbers
widget_control,state.specwindow_id, get_value = tmp_value
state.specwindow = tmp_value
widget_control,state.fitwindow_id, get_value = tmp_value
state.fitwindow = tmp_value

;xmanager,'wave_calib_window',state.base_id,group_leader=group,$
;   cleanup='wave_calib_cleanup',/no_block

end

;----------------------------------------------------------------------
; Convert from pixel x-value to fitted wavelength coordinate
;----------------------------------------------------------------------
Function xpix2wave,tmp_xpix

common wc_state

if (state.fitcomputed eq 0) then begin
   print,'xpix2wave: No fit computed'
   return,-1
endif

tmp_fitorders = where(state.fitcoeffs ne 0) ; Determine polynomial order
tmp_norders   = n_elements(tmp_fitorders)

tmp_fitwave = 0d0

for i=0,tmp_norders do begin
   tmp_coeff = state.fitcoeffs(i)
   tmp_fitwave = tmp_fitwave + tmp_coeff*(tmp_xpix^i)
endfor

return,tmp_fitwave

end

;---------------------------------------------------------------------
; Draw spectrum plot
;---------------------------------------------------------------------
Pro draw_specplot,norescale=norescale

common wc_state
common wc_data

thisspecregion = where(xplotarr ge min([state.xdispmin,state.xdispmax]) and $
                       xplotarr le max([state.xdispmin,state.xdispmax]),$
                       n_region)

; Ignore top and bottom 0.5% of points when scaling plot
yinput_region = yinputarr(thisspecregion)
region_order  = sort(yinput_region)
yinput_sorted = yinput_region(region_order)
dispmin_cut   = 0.005
dispmax_cut   = 0.995

if not(keyword_set(norescale)) then begin
   ;region_min = min(yinputarr(thisspecregion))
   ;region_max = max(yinputarr(thisspecregion))
   ;state.specydispmin = region_min
   ;state.specydispmax = region_max + 0.2*(region_max - region_min)
   region_min = yinput_sorted(ceil(dispmin_cut*n_region))
   region_max = yinput_sorted(floor(dispmax_cut*n_region))
   state.specydispmin = region_min - 0.1*(region_max - region_min)
   state.specydispmax = region_max + 0.5*(region_max - region_min)
endif

wset,state.specwindow
plot,xplotarr,yinputarr,xrange=[state.xdispmin,state.xdispmax],xstyle=1,$
   yrange=[state.specydispmin,state.specydispmax],ystyle=1,$
   psym=state.psym

marklines

; Refresh displays of plot range values
widget_control,state.ymin_text_id,set_value=state.specydispmin
widget_control,state.ymax_text_id,set_value=state.specydispmax
widget_control,state.xmin_text_id,set_value=state.xdispmin
widget_control,state.xmax_text_id,set_value=state.xdispmax

end

;---------------------------------------------------------------------
; Draw wavelength fit plot
;---------------------------------------------------------------------
Pro draw_fitplot

common wc_state
common wc_data

; Plot wavelength of any lines selected.  

selectedls = where(state.lineswithid eq 1)

wset,state.fitwindow

if (total(selectedls) ne -1) then begin

   if (state.xwave eq 1) then begin
      selectedlinexvalues = state.linefitwaves(selectedls)
   endif else begin
      selectedlinexvalues = state.linexvalues(selectedls)
   endelse
   selectedlinewaves   = state.linewaves(selectedls)

   linestoplot = $
       where(selectedlinexvalues ge min([state.xdispmin,state.xdispmax]) and $
             selectedlinexvalues le max([state.xdispmin,state.xdispmax]))

   if (state.plotresiduals eq 1) then begin    ; Plot residuals?

      ; If a fit has been computed, plot wavelength residuals vs.
      ; pixel values
      if (state.fitcomputed eq 1) then begin
         tmp_maxresid = max(abs(state.lineresids))
         plot,selectedlinexvalues,state.lineresids(selectedls),$
            xrange=[state.xdispmin,state.xdispmax],xstyle=1,$
            yrange=[-1*tmp_maxresid,tmp_maxresid],psym=1
         oplot,[state.xdispmin,state.xdispmax],[0,0],linestyle=1
         ; Overplot 1 and 3-"sigma" ranges
         oplot,[state.xdispmin,state.xdispmax],[state.fiterror,state.fiterror],linestyle=1
         oplot,[state.xdispmin,state.xdispmax],[-1*state.fiterror,-1*state.fiterror],linestyle=1
         oplot,[state.xdispmin,state.xdispmax],[3*state.fiterror,3*state.fiterror],linestyle=1
         oplot,[state.xdispmin,state.xdispmax],[-3*state.fiterror,-3*state.fiterror],linestyle=1
      endif else begin
         print,'Must perform fit before displaying residuals.'
         return
       endelse

   endif else begin ; Otherwise plot wavelength vs. pixel value

      ; Reset plotting y-range so long as some lines will be displayed
      if (total(linestoplot) ne -1) then begin
         if (n_elements(linestoplot) eq 1) then begin
            state.fitydispmin = 0.99999*min(selectedlinewaves(linestoplot))
            state.fitydispmax = 1.00001*max(selectedlinewaves(linestoplot))
         endif else begin
            ydispdiff = abs(min(selectedlinewaves(linestoplot)) $
                            - max(selectedlinewaves(linestoplot)))
            state.fitydispmin = min(selectedlinewaves(linestoplot)) - 0.1*ydispdiff
            state.fitydispmax = max(selectedlinewaves(linestoplot)) + 0.1*ydispdiff
         endelse
      endif

      ; If a wavelength fit has been computed, set display range to show
      ; fit for the selected region
      if (state.fitcomputed eq 1) then begin
         thisregion = where(xplotarr ge min([state.xdispmin,state.xdispmax]) $
                        and xplotarr le max([state.xdispmin,state.xdispmax]))
         alt_fitydispmin = min(wavelengthfit(thisregion)) > 0  ; Don't show neg wavelengths
         alt_fitydispmax = max(wavelengthfit(thisregion))
         ; Adopt these values if outside the current display range
         state.fitydispmin = state.fitydispmin < alt_fitydispmin
         state.fitydispmax = state.fitydispmax > alt_fitydispmax
      endif

      plot,selectedlinexvalues,selectedlinewaves,$
         xrange=[state.xdispmin,state.xdispmax],xstyle=1,$
         yrange=[state.fitydispmin,state.fitydispmax],ystyle=1,psym=1

      ; If a wavelength fit has been computed, overplot
      if (state.fitcomputed eq 1) then oplot,xplotarr,wavelengthfit,linestyle=1

   endelse

endif else begin  ;If no lines selected, produce a blank plot.

   plot,[0,0],[0,0],xrange=[state.xdispmin,state.xdispmax],xstyle=1,$
      yrange=[0,1]

endelse

thisspecregion = where(xplotarr ge min([state.xdispmin,state.xdispmax]) and $
                       xplotarr le max([state.xdispmin,state.xdispmax]))

end

;--------------------------------------------------------------------
; Mark lines selected for wavelength fitting
;--------------------------------------------------------------------
Pro marklines

common wc_state

; If a line has been selected, mark it with a tick
for i=0,state.nlinestot-1 do begin
   if (state.lineswithid(i) eq 1) then begin
      drawtick,i
   endif
endfor

end   

;---------------------------------------------------------------------
; Locate line
;---------------------------------------------------------------------
Pro findline

common wc_state
common wc_data

; Make sure a wavelength has been selected from the list
if (state.currentwave eq 0) then begin
   print,'Must choose a line from the list first!'
   return
endif

; Find the x-value of the maximum point near the current position.

; +/- distance in data units (pixels) around selected x-value to
; search for maximum.
if (state.xwave eq 1) then search_dist = abs(state.searchwaverad) $
   else search_dist = abs(state.searchpixrad)

xcent   = state.xcoord

chunkls = where(xplotarr ge xcent-search_dist and $
                xplotarr le xcent+search_dist,nchunkpix)

xchunk    = xplotarr(chunkls)  ; x values currently plotted
xpixchunk = xinputarr(chunkls) ; x pixel values (these may be the same)
ychunk    = yinputarr(chunkls)

; Ignore the top 2% of pixels when estimating the peak value
temp_order = sort(ychunk)
yorder     = ychunk(temp_order)
ychunk_max = yorder(floor(0.98*(nchunkpix-1)))
;ychunk_max = max(ychunk)
maxpix = where(ychunk eq ychunk_max)

; More than one pixel at max?
if (n_elements(maxpix) gt 1) then begin 
   print,'More than one pixel found at max, taking the lowest.'
   maxpix = min(maxpix)
endif

; Make sure that the max pixel is not all the way at one end of
; the selected region (no local maximum).
if (maxpix eq 0 or maxpix eq nchunkpix-1) then begin
   print,'No local maximum found.'
   return
endif

xpixchunk_max = (xpixchunk(maxpix))(0)  ; Force conversion to scalar

; Fit a Gaussian to the neighboring pixels (in pixel space).  Start by
; ignoring the top 2% of pixels.  Reject outliers and refit.
tmp_lox = xpixchunk_max - state.gfitpixrad
tmp_hix = xpixchunk_max + state.gfitpixrad
gausschunkls = where(xinputarr ge tmp_lox and xinputarr le tmp_hix,n_gauss)
xgausschunk = xinputarr(gausschunkls)
ygausschunk = yinputarr(gausschunkls)
gauss_order = sort(ygausschunk)
cut = (ygausschunk(gauss_order))(floor(0.98*(n_gauss-1)))
fit_ls = where(ygausschunk lt cut)
for i=0,2 do begin
   heightguess = max(ygausschunk(fit_ls)) - min(ygausschunk(fit_ls))
   centguess   = xpixchunk_max
   sigguess    = state.gfitpixrad / 3.
   constguess  = 0 ;min(ychunk)
   est = [heightguess,centguess,sigguess,constguess]
   tmp_gfit = gaussfit(xgausschunk(fit_ls),ygausschunk(fit_ls),$
                       tmp_terms,nterms=4,estimates=est,sigma=tmp_sigma,$
                       yerror=temp_yerr)
   z = (xgausschunk-tmp_terms(1))/tmp_terms(2)
   gfit = tmp_terms(0)*exp(-(z^2)/2.) + tmp_terms(3)
   fit_ls = where(abs(ygausschunk-gfit) lt 5*temp_yerr)
endfor
tmp_xpixcenter = tmp_terms(1) ; Line center
tmp_ypixmax    = tmp_terms(0) + tmp_terms(3) ; Line amplitude

; If the fitted center falls outside the current seach range, reject 
; it.
if (tmp_xpixcenter lt min(xgausschunk) or $
    tmp_xpixcenter gt max(xgausschunk)) then begin
   print,'Fitted line center falls outside of range.'
   return
endif

; Optional - If the fitted line appears to be a noise feature, reject it.
; Do this only when locating lines automatically (i.e., do not reject
; if the line was selected with the mouse).
if (state.rejectnoise and state.mbutton ne 'left') then begin
   compare_ls  = where((xplotarr ge xcent-10*search_dist and $
                        xplotarr le xcent-search_dist) or $
                       (xplotarr ge xcent+search_dist and $
                        xplotarr le xcent+10*search_dist),n_compare)
   local_ymedian = median(yinputarr(compare_ls))
   local_ysdev = robust_sigma(yinputarr(compare_ls))
   if (tmp_ypixmax lt local_ymedian + 3*local_ysdev) then begin
      print,'Line appears to be noise.'
      return
   endif
endif

; Optional - If the intensity of the line exceeds a given limit,
; reject it.
if state.rejectstrong then begin
   ; Use a narrow window around the fitted center
   check_ls  = where(xinputarr ge tmp_xpixcenter-tmp_terms(2) and $
                     xinputarr le tmp_xpixcenter+tmp_terms(2),n_check)
   if (n_check gt 0) then begin
      max_value = max(yinputarr(check_ls))
      if (max_value ge state.limit) then begin
         print,'Line intensity exceeds limit.'
         return
      endif
   endif else print,'No pixels near line center!'
endif

; Line list index of the current line
thiswavelength = where(state.linewaves eq state.currentwave)
; If there are duplicate entries in the list, take the first one.
if (n_elements(thiswavelength) gt 1) then thiswavelength = thiswavelength(0)

; If the line has already been selected, 'erase' the old tick mark
; before drawing a new one.
if (state.lineswithid(thiswavelength) eq 1) then begin
   drawtick,thiswavelength,/erase
endif

; Save this pixel x-value as belonging to the currently selected wavelength
state.linexvalues(thiswavelength) = tmp_xpixcenter
state.lineyvalues(thiswavelength) = tmp_ypixmax

; Compute the wavelength at this x-value under the current fit, if one
; exists 
if (state.fitcomputed eq 1) then $
   state.linefitwaves(thiswavelength) = xpix2wave(tmp_xpixcenter)

; Print line center
if (state.xwave eq 1) then begin  ; Displaying wavelengths
   tmp_print_center = xpix2wave(tmp_xpixcenter)
endif else begin
   tmp_print_center = tmp_xpixcenter
endelse

print,'Line with wavelength ',strtrim(state.currentwave,2),$
      ' found at x = ',strtrim(tmp_print_center,2)

; Flag this line as identified
state.lineswithid(thiswavelength) = 1. 

; Mark all lines found so far
marklines

end

;----------------------------------------------------------------------
; Remove line from wavelength fit
;----------------------------------------------------------------------
Pro removeline

common wc_state
common wc_data

; Make sure that there is at least one line currently selected for 
; fitting
selectedls = where(state.lineswithid eq 1)
if (total(selectedls) eq -1) then begin
   print,'No lines currently selected.'
   return
endif

; Compute min x difference between current x coord and any selected
; line and find the already-selected line closest to the current
; cursor position
if (state.xwave eq 1) then begin ; Fitted wavelength display mode
   minabsxdiff = min(abs(state.linefitwaves(selectedls)-state.xcoord(selectedls)))
   linetocut = where( abs(state.linefitwaves-state.xcoord) eq minabsxdiff and $
                      state.lineswithid eq 1 )
endif else begin ; Pixel display mode
   minabsxdiff = min(abs(state.linexvalues(selectedls)-state.xcoord(selectedls)))
   linetocut = where( abs(state.linexvalues-state.xcoord) eq minabsxdiff and $
                      state.lineswithid eq 1 )
endelse

; If more than one line has been assigned to the same pixel x-value, remove
; the first one in the list
if (n_elements(linetocut) gt 1) then linetocut=linetocut(0)

; Flag the line as not selected
state.lineswithid(linetocut) = 0.

; 'Erase' the current tick mark
drawtick,linetocut,/erase

print,'Line with wavelength ',state.linewaves(linetocut),' unselected from x = ',state.linexvalues(linetocut)

end

;---------------------------------------------------------------------
; Identify additional lines based on the current fit
;---------------------------------------------------------------------
Pro find_more_lines

common wc_state
common wc_data

; Make sure a fit has been computed
if (state.fitcomputed eq 0) then begin
   print,'Must compute a fit first!'
   return
endif

; For now, make sure plotting vs. pixels
if (state.xwave eq 1) then begin
   print,'Must be plotting vs. pixels (fix this).'
   return
endif

; Select lines within the current wavelength range that
; have not yet been identified.
minwave = min(wavelengthfit)
maxwave = max(wavelengthfit)
linestoidentify_ls = where(state.linewaves gt minwave and $
                           state.linewaves lt maxwave and $
                           state.lineswithid eq 0,n_linestoidentify)
if (n_linestoidentify eq 0) then begin
   print,'No more lines within wavelength range.'
   return
endif

; Attempt to identify each new line
for i=0L,long(n_linestoidentify-1) do begin
   state.currentwave = state.linewaves(linestoidentify_ls(i))
   line_xguess = interpol(xinputarr,wavelengthfit,state.currentwave)
   if (line_xguess gt state.xmin and $
       line_xguess lt state.xmax) then begin
      state.xcoord = line_xguess     
      findline
   endif
endfor

end

;---------------------------------------------------------------------
; Reject multiple identifications for the same lines
;---------------------------------------------------------------------
Pro reject_multiple_ids

common wc_state
common wc_data

; Get list of lines that have been identified from the list
id_ls = where(state.lineswithid eq 1,n_id)
if (n_id lt 2) then begin
   print,'Must identify at least two lines first!'
   return
endif

; Step through the list of identified lines.  Flag those that have
; very similar x-values.
x_tol = 0.1 * state.searchpixrad
for i=0,n_id-1 do begin
   tmp_x = state.linexvalues(id_ls(i))
   x_ls = where(state.linexvalues gt tmp_x - x_tol and $
                state.linexvalues lt tmp_x + x_tol,n_x)
   if (n_x gt 1) then begin
      state.lineswithid(x_ls) = 0
      ;tmp_resid = state.lineresids(x_ls)
      ;best = where(abs(tmp_resid) eq min(abs(tmp_resid)),$
      ;             complement=redundant_ls,ncomplement=n_redundant)
      ;if (n_redundant ne 0) then $
      ;   state.lineswithid(x_ls(redundant_ls)) = 0
   endif
endfor

; Reject redundant identifications
reject_ls = where(state.lineswithid eq 0,n_reject)
if (n_reject ne 0) then begin
   state.linefitwaves(reject_ls) = 0
   state.lineresids(reject_ls)   = 0
   state.linesused(reject_ls)    = 0
   state.linexvalues(reject_ls)  = 0
   state.lineyvalues(reject_ls)  = 0
endif

; Redraw plots
if state.show then begin
   draw_fitplot
   draw_specplot
endif

end

;---------------------------------------------------------------------
; Reject outliers from the computed fit
;---------------------------------------------------------------------
Pro reject_outliers

common wc_state

; Make sure a fit has been computed
if (state.fitcomputed eq 0) then begin
   print,'Must compute a fit first!'
   return
endif

id_ls = where(state.lineswithid eq 1,n_id)
if (n_id lt 2) then begin
   print,'Must identify at least two lines first!'
   return
endif

nsig_rej = 3

rej_ls = where(abs(state.lineresids) gt nsig_rej*state.fiterror,n_rej)
if (n_rej ne 0) then begin
   state.linefitwaves(rej_ls) = 0
   state.lineresids(rej_ls)   = 0
   state.lineswithid(rej_ls)  = 0
   state.linesused(rej_ls)    = 0
   state.linexvalues(rej_ls)  = 0
   state.lineyvalues(rej_ls)  = 0
endif

; Redraw plots
if state.show then begin
   draw_fitplot
   draw_specplot
endif

end

;---------------------------------------------------------------------
; Draw tick with wavelength
;---------------------------------------------------------------------
Pro drawtick,lineindex,erase=erase

common wc_state

if state.show then begin

   ; Draw tick with wavelenth in spec window

   if keyword_set(erase) then tcolor=0 else tcolor=250

   tl = state.ticklength * (state.specydispmax - state.specydispmin)

   xtickval = state.linexvalues(lineindex)
   ; If plotting vs fitted wavelength, convert x coordinate
   if (state.xwave eq 1) then xtickval = xpix2wave(xtickval)
   ytickval = state.lineyvalues(lineindex)
   if (xtickval ge min([state.xdispmin,state.xdispmax]) and $
       xtickval le max([state.xdispmin,state.xdispmax])) then begin
      oplot,[xtickval,xtickval],[ytickval+0.5*tl,ytickval+1.5*tl],color=tcolor
      xyouts,xtickval,ytickval+1.75*tl,state.linenames(lineindex),$
             orientation=90.,color=tcolor,/data
   endif

endif

end

;---------------------------------------------------------------------
; Compute wavelength fit
;---------------------------------------------------------------------
Pro computefit

common wc_state
common wc_data

; For now, include all lines that have been identified
state.linesused = state.lineswithid

; Make sure at least two lines have been selected for fitting
fitls = where(state.linesused eq 1,n_fit)
if (n_fit lt 1) then begin
   print,'Must select at least one line before computing fit.'
   state.fitcomputed = 0
   wavelengthfit(*) = 0.
   return
endif

;; Case where only one line has been identified
;if(n_fit eq 1) then begin
;   wavelengthfit(*)   = state.linewaves(fitls)
;   state.fitcoeffs(*) = 0.
;   state.fitcoeffs(0) = state.linewaves(fitls)
;   state.fitcomputed  = 1
;   return
;endif

; If number of lines with unique x-values is equal to or smaller than
; the specified polynomial fit order, use a lower order polynomial
unique_ls    = uniq(state.linexvalues(fitls))
n_unique     = n_elements(unique_ls)
tmp_fitorder = state.fitorder < (n_unique - 1)

; Case where only one unique line has been identified
if(n_fit eq 1) then begin
   wavelengthfit(*)   = mean(state.linewaves(fitls))
   state.fitcoeffs(*) = 0.
   state.fitcoeffs(0) = state.linewaves(fitls)
   state.fitcomputed  = 1
   return
endif

;if (n_unique le tmp_fitorder + 2) then begin
;   tmp_coeffs = poly_fit(state.linexvalues(fitls),state.linewaves(fitls),$
;                         tmp_fitorder,yerror=tmp_fitsdev,yfit=tmp_fit)
;endif else begin
;   tmp_coeffs = robust_poly_fit(state.linexvalues(fitls),$
;                                state.linewaves(fitls),tmp_fitorder,tmp_fit,$
;                                tmp_fitsdev)
;endelse

; TESTING: Try robust_poly_fit first (better at fitting points where range in
; x is much smaller than x values, as well as dealing with outliers).  
; If robust_poly_fit returns a poor result, use standard poly_fit
tmp_coeffs = robust_poly_fit(state.linexvalues(fitls),$
                             state.linewaves(fitls),tmp_fitorder,tmp_fit,$
                             tmp_fitsdev)
if (total(tmp_coeffs) eq 0) then begin
   tmp_coeffs = poly_fit(state.linexvalues(fitls),state.linewaves(fitls),$
                         tmp_fitorder,yerror=tmp_fitsdev,yfit=tmp_fit)
endif


; Set flag that a fit has been computed
state.fitcomputed = 1

; Save coefficients and standard error
state.fitcoeffs(*)=0.
state.fitcoeffs(0:tmp_fitorder) = tmp_coeffs
; Note: The standard error returned by the 'yerror' keyword in
; poly_fit is the standard deviation corrected for the degrees of
; freedom used up in computing the fit.
state.fiterror = tmp_fitsdev

; Compute fitted wavelength for all input pixels
wavelengthfit(*) = 0.
for i=0,tmp_fitorder do wavelengthfit = wavelengthfit $
                                         + tmp_coeffs(i)*(xinputarr^i)

; Update plotting array if plotting against x-values
if (state.xwave eq 1) then xplotarr = wavelengthfit

; Compute new fitted wavelengths and residuals for all identified lines
computeresiduals

; Compute the mean dispersion
dispersionfit = 0*wavelengthfit
for i=1,tmp_fitorder do dispersionfit = dispersionfit $
                                         + i*tmp_coeffs(i)*(xinputarr^(i-1))
state.dispersion = mean(dispersionfit)

; Adjust wavelength range over which to search for and fit lines
state.searchwaverad = abs(state.searchpixrad * state.dispersion)
state.gfitwaverad   = abs(state.gfitpixrad * state.dispersion)

; Redraw plots with fit
if state.show then begin
   draw_fitplot
   draw_specplot
endif

end

;----------------------------------------------------------------------
; Compute wavlength fit values and residuals
;----------------------------------------------------------------------
Pro computeresiduals

common wc_state

selectedls = where(state.lineswithid eq 1)

; Return if no lines selected
if (total(selectedls) eq -1) then begin
   print,'No lines selected'
   return
endif

n_selected = n_elements(selectedls)

; Compute fitted wavelengths

for i=0,n_selected-1 do begin
   state.linefitwaves(selectedls(i)) = xpix2wave(state.linexvalues(selectedls(i)))
endfor

; Compute residuals for all identified lines
state.lineresids(*)=0.
state.lineresids(selectedls) = state.linewaves(selectedls) $
                                 - state.linefitwaves(selectedls)

end

;----------------------------------------------------------------------
; Convert from plotting vs pixels to vs wavelength, or visa versa
;----------------------------------------------------------------------
Pro swapxplot

common wc_state
common wc_data

; Swap from plotting vs pixels to vs wavelength
if (state.xwave eq 1) then begin
   xplotarr = wavelengthfit
   pix_dispmin = state.xdispmin
   pix_dispmax = state.xdispmax
   wave_dispmin = xpix2wave(pix_dispmin)
   wave_dispmax = xpix2wave(pix_dispmax)
   ; If, by some strage fit, wave_dispmin > wave_dispmax, switch
   state.xdispmin = wave_dispmin < wave_dispmax
   state.xdispmax = wave_dispmin > wave_dispmax
endif

; Swap from plotting vs wavelength to vs pixels
if (state.xwave eq 0) then begin
   xplotarr = xinputarr
   ; Rather than invert polynomial, reset to full x range
   state.xdispmin = min(xplotarr)
   state.xdispmax = max(xplotarr)
endif
   
; Reset values in display
widget_control,state.xmin_text_id,set_value=state.xdispmin
widget_control,state.xmax_text_id,set_value=state.xdispmax

draw_fitplot
draw_specplot

end

;----------------------------------------------------------------------
; Cleanup program when base widget dies
;----------------------------------------------------------------------
Pro wave_calib_cleanup,widget_id

common wc_state

end
   
;----------------------------------------------------------------------
; Main program - See help at beginning of file
;----------------------------------------------------------------------
Pro wave_calib,xin,yin,wave,linefile=linefile,guess=guess,reference=reference,$
               nopersist=nopersist,psout=psout,coeffs=coeffs,xid=xid,$
               waveid=waveid,order=order,dispersion=dispersion,$
               resolution=resolution,xoffset=xoffset,auto=auto,$
               rejectnoise=rejectnoise,limit=limit,noshow=noshow,$
               psym=psym

if (n_params() lt 2) then begin
   print,'CALLING SEQUENCE:'
   print,'   wave_calib,xin,yin,wave,linefile=linefile,guess=guess,'
   print,'         reference=reference,nopersist=nopersist,psout=psout,'
   print,'         coeffs=coeffs,xid=xid,waveid=waveid,order=order,'
   print,'         dispersion=dispersion,resolution=resolution,xoffset=xoffset,'
   print,'         /auto,/rejectnoise,limit=limit,/noshow,psym=psym'
   return
endif

common wc_state
common wc_data

; Plot in color
if not(keyword_set(noshow)) then begin
   device,decomp=0
   loadct,39
endif

xinputarr = xin
yinputarr = yin

wavelengthfit = double(0*xin)

; Set input x array as inital x array for plotting
xplotarr = xinputarr

; Create window if one does not already exist
if not(xregistered('wave_calib_window')) then begin
   ; Initialize common blocks and settings
   if keyword_set(linefile) then begin
      wave_calib_initialize,linefile=linefile
   endif else begin
      wave_calib_initialize
   endelse
   ; Plot symbol?
   if keyword_set(psym) then state.psym = psym
   ; Display window?
   if keyword_set(noshow) then state.show = 0 else state.show = 1
   state.xdispmin = min(xplotarr)
   state.xdispmax = max(xplotarr)   
   if (state.show) then begin
      create_wave_calib_window
      widget_control,state.xmin_text_id,set_value=state.xdispmin  
      widget_control,state.xmax_text_id,set_value=state.xdispmax 
      draw_fitplot
      draw_specplot   
   endif
endif else begin
   ; Clear line identification and wavelength fit information
   state.linefitwaves(*) = 0
   state.lineresids(*)   = 0
   state.lineswithid(*)  = 0
   state.linesused(*)    = 0
   state.linexvalues(*)  = 0
   state.lineyvalues(*)  = 0
   state.fitcoeffs(*)    = 0
   state.fiterror        = 0
   state.fitcomputed     = 0
   state.fitorder        = 3
endelse

; Change default fit order?
if keyword_set(order) then state.fitorder = order
if state.show then $
   widget_control,state.polyorder_text_id,set_value=state.fitorder

; Establish min/max data values
state.xmin = min(xinputarr)
state.xmax = max(xinputarr)
state.ymin = min(yinputarr)
state.ymax = max(yinputarr)

;Let initial plot range span entire input range
state.xdispmin = min(xplotarr)
state.xdispmax = max(xplotarr)

; Produce PS plot when finished?
if keyword_set(psout) then begin
   state.createps = 1
   state.psfile = psout
endif

; Adjust parameters for finding/fitting lines (optional)?
if keyword_set(dispersion) then state.dispersion = 1.*dispersion
state.resolution = 3.*state.dispersion ; Default
if keyword_set(resolution) then state.resolution = 1.*resolution
state.searchpixrad  = abs(state.resolution / state.dispersion)
state.searchwaverad = abs(state.searchpixrad * state.dispersion)
state.gfitpixrad    = abs(state.resolution / state.dispersion)
state.gfitwaverad   = abs(state.gfitpixrad * state.dispersion)

; Default is no offset from a reference array
xoffset = 0

; Default is to identify lines
id_lines = 1

; Reject lines that appear to be noise?
if keyword_set(rejectnoise) then state.rejectnoise = 1

; Reject strong lines above a given limit?
if keyword_set(limit) then begin
   state.rejectstrong = 1
   state.limit        = limit
endif

; Automatic line identification (optional)
if (keyword_set(reference) or keyword_set(auto)) then begin
   ; Create widget without blocking
   if state.show then begin
      xmanager,'wave_calib_window',state.base_id,group_leader=group,$
         cleanup='wave_calib_cleanup',/no_block
      draw_fitplot
      draw_specplot
   endif
   ; Use a reference fit to get expected offsets (optional)
   if keyword_set(reference) then begin
      ; Read in reference array from a file or from an array.
      if (size(reference,/tname) eq 'STRING') then begin
         ; Template for three columns of double-precision entries (x,y,wave,
         ; which are pixel, intensity, wavelength)
         templ = {   $
                  version: 1.0, $    
                  datastart: 0L, $   
                  delimiter: ' ', $        ; Fields separated by white space
                  missingvalue: -99.0, $
                  commentsymbol: '#', $
                  fieldcount: 3L, $
                  fieldtypes: [5L,5L,5L], $       ; Double precision
                  fieldnames: ['xref','yref','waveref'], $
                  fieldlocations: [6L,23L,39L], $
                  fieldgroups: [0L,1L,2L] $
                 }
         ; Read in reference data
         print,'Reading in reference wavelength fit...'
         refdata = read_ascii(reference,template=templ)
         xref    = refdata.xref
         yref    = refdata.yref
         waveref = refdata.waveref
      endif else begin
         refarrsz = size(reference)
         if (refarrsz(0) ne 2 or refarrsz(2) ne 3) then begin
            print,'WAVE_CALIB: Reference wave solution must be an ascii file'
            print,'            or an n x 3 array.'
         endif
         xref    = reference(*,0)
         yref    = reference(*,1)
         waveref = reference(*,2)
      endelse
      ; Compute bulk x-offset between xin and xref using crosscorrelation
      ; Ignore the first and last 5% of the input array to avoid edge
      ; effects.
      print,'Computing offset to reference spectrum...'
      xchop = 0.05 * (max(xinputarr)-min(xinputarr))
      print,'WAVE_CALIB: Discarding ends of input spectrum to compute'
      print,'            cross-corr with reference spectrum.'
      ;print,'WAVE_CALIB: Warning - Using full wavelength range to compute '
      ;print,'                      cross-corr.  This may cause problems.'
      ; Make sure there is at least some overlap
      overlap_ls = where(xref ge min(xinputarr)+xchop and $
                         xref le max(xinputarr)-xchop,n_overlap)
      if (n_overlap eq 0) then begin
         print,'No overlap with reference spectrum.  Assuming zero offset.'
         xoffset = 0
         minwave = state.xmin
         maxwave = state.xmax
      endif else begin
         ; Sort reference spectrum
         ref_order = sort(xref)
         xref      = xref(ref_order)
         yref      = yref(ref_order)
         waveref   = waveref(ref_order)
         ; Interpolate input and reference arrays onto a common 
         ; wavelength grid.
         xlo   = (min(xinputarr)+xchop) > min(xref)
         xhi   = (max(xinputarr)-xchop) < max(xref)
         xstep = (max(xref)-min(xref)) / n_elements(xref)
         nbins = floor((xhi-xlo)/xstep)
         xgrid = fltarr(nbins)
         for i=0L,long(nbins-1) do xgrid(i) = xlo + i*xstep
         yingrid  = interpol(yinputarr,xinputarr,xgrid)
         yrefgrid = interpol(yref,xref,xgrid)
         ; Assume offset is within +/- 5% of x-range, but evaluate cc
         ; over a slightly larger range.  If the max in the cc occurs
         ; at the ends, then it is probably spurious.
         nlag = 0.12 * nbins
         lag = findgen(nlag) - round(nlag/2.)
         cc  = c_correlate(yingrid,yrefgrid,lag)
         best_lag = lag(where(cc eq max(cc)))
         if (abs(best_lag) ge (5./6)*max(lag) or $
             max(cc) lt median(cc)+5*robust_sigma(cc)) then begin
            print,'No good cross-correlation with reference wavelength fit found.'
            id_lines = 0
         endif else begin
            xoffset = xstep * lag(where(cc eq max(cc)))
            xoffset = xoffset(0)
            minwave = min(waveref)
            maxwave = max(waveref)
         endelse
      endelse
   endif else begin
      ; If auto-identifying lines without first matching to a 
      ; reference spectrum, assume that the input x-values are 
      ; themselves calibrated physical wavelengths.
      minwave = state.xmin
      maxwave = state.xmax 
   endelse
   ; Search for each line in the input list based on expected position
   ; Only look for lines within the input or reference wavelength 
   ; range.
   if id_lines then begin
      print,'Identifying lines...'
      linesinrange_ls = where(state.linewaves ge minwave and $
                              state.linewaves le maxwave,n_linesinrange)
      if (n_linesinrange eq 0) then begin
         print,'WAVE_CALIB: Line list does not contain any lines within the'
         print,'            wavelength range of the reference spectrum.'
      endif else begin
         for i=0L,long(n_linesinrange-1) do begin
            state.currentwave = state.linewaves(linesinrange_ls(i))
            if keyword_set(reference) then begin
               line_xref    = interpol(xref,waveref,state.currentwave)
               line_xguess  = line_xref - xoffset
            endif else line_xguess = state.currentwave
            if (line_xguess gt state.xmin and $
                line_xguess lt state.xmax) then begin
               state.xcoord = line_xguess     
               findline
            endif
         endfor
         ; Reject lines ID'd at the same x-values.
         reject_multiple_ids
         ; If lines have been identified, then compute fit
         id_ls = where(state.lineswithid ne 0,n_id)
         if (n_id ne 0) then begin
            if state.show then state.plotresiduals = 1
            computefit
            ; Reject multiple identifications and recompute fit
            reject_multiple_ids
            computefit
            ; Reject outliers & recompute fit until fit converges
            fit_converged = 0
            old_wavelengthfit = 0
            n_fit_iter = 0
            while not(fit_converged) do begin
               reject_outliers
               computefit
               fit_diff = wavelengthfit - old_wavelengthfit
               if (max(abs(fit_diff)) eq 0) then fit_converged = 1
               n_fit_iter = n_fit_iter + 1
               if (n_fit_iter ge 10) then fit_converged = 1
               old_wavelengthfit = wavelengthfit
            endwhile
            print,'Wavelength fit completed.'
         endif
      endelse
   endif

   ; Kill widget (optional)
   if state.show then begin
      if keyword_set(nopersist) then widget_control,state.base_id,/destroy
   endif

endif else begin

   ; Optional - Attempt to identify lines based on input guesses
   if keyword_set(guess) then begin
      ; Read in line position guesses from a file or from an array.
      if (size(guess,/tname) eq 'STRING') then begin
         ; Template for 2 columns of double-precision entries (x,wave)
         guess_templ = {   $
                  version: 1.0, $    
                  datastart: 0L, $   
                  delimiter: ' ', $        ; Fields separated by white space
                  missingvalue: -99.0, $
                  commentsymbol: '#', $
                  fieldcount: 2L, $
                  fieldtypes: [5L,5L], $       ; Double precision
                  fieldnames: ['xguess','waveguess'], $
                  fieldlocations: [6L,23L], $
                  fieldgroups: [0L,1L] $
                 }
         ; Read in reference data
         print,'Reading in line position guesses'
         guessdata  = read_ascii(guess,template=guess_templ)
         x_guess    = guessdata.xguess
         wave_guess = guessdata.waveguess
      endif else begin
         guessarrsz = size(guess)
         if (guessarrsz(0) ne 2) then begin
            print,'WAVE_CALIB: Guess identifications must be an ascii file'
            print,'            or an n x 2 array.'
            return
         endif
         x_guess    = guess(*,0)
         wave_guess = guess(*,1)
      endelse
      n_guess = n_elements(x_guess)
      ; Attempt to identify each guess
      for i=0,n_guess-1 do begin
         state.xcoord = x_guess(i)
         line_match   = where(abs(state.linewaves-wave_guess(i)) eq $
                              min(abs(state.linewaves-wave_guess(i))),n_match)
         if (n_match gt 1) then line_match = line_match(0)
         state.currentwave = state.linewaves(line_match)
         findline
      endfor
      ; Compute fit using these lines and show residuals
      state.plotresiduals = 1
      computefit
   endif

   ; If identifying lines by hand, create widget with blocking
   xmanager,'wave_calib_window',state.base_id,group_leader=group,$
      cleanup='wave_calib_cleanup'

endelse

; Output polynomial coefficients
tmp_coeffs = state.fitcoeffs
ls = where(tmp_coeffs ne 0)
if (total(ls) ne -1) then begin
   highestorder = max(ls)
   coeffs = tmp_coeffs(0:highestorder)
endif else coeffs=0.

; Output fitted wavelength array
wave = wavelengthfit

; Output the x pixel values and wavelengths of identified lines
; included in fit.
used_ls = where(state.linesused eq 1,n_used)
if (n_used ne 0) then begin
   xid    = state.linexvalues(used_ls)
   waveid = state.linewaves(used_ls)
endif else begin
   xid    = -1
   waveid = -1
endelse

; Output the order of the polynomial fit
order = state.fitorder

; Set color table to greyscale
if not(keyword_set(noshow)) then loadct,0

return
end

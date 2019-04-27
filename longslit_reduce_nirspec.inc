; Instrument-specific parameters for NIRSPEC
;
; Amplifiers are ordered from left to right, then from bottom to top.
;
; New relative gains, darkcurrent, and read noise measured Feb 2006.
; Assumed gain for lower-left amplifier = 5.0

nchip         = 1        ; Number of CCDs included in each array
chip_xmin     = [0]      ; Min x for each chip (not including overscan)
chip_xmax     = [2047]   ; Max x for each chip (not including overscan)
chip_ymin     = [0]      ; Min y for each chip (not including overscan)
chip_ymax     = [2047]   ; Max y for each chip (not including overscan)
namp          = 1        ; Number of amplifiers
gain          = [4.10] ; Gain for each amplifier (e-/ADU)
amp_xmin      = [0]        ; Min x for each amplifier (not including overscan)
amp_xmax      = [2047]  ; Max x for each amplifier (not including overscan)
amp_ymin      = [0]        ; Min y for each amplifier (not including overscan)
amp_ymax      = [2047]  ; Max y for each amplifier (not including overscan)
erdnoise      = [7.3] ; Read noise per coadd (e-), scalar or for each amplifier
satlevel      = 30000.   ; Saturation level (ADU) 
dispdir       = 'x'      ; Dispersion direction, 'x' or 'y'
prescan_x     = -1    ; [x_min,x_max] of prescan region (-1 = none)
prescan_y     = -1    ; [y_min,y_max] of prescan region (-1 = none)
oscan_x       = -1    ; [x_min,x_max] of overscan region (-1 = none)
oscan_y       = -1    ; [y_min,y_max] of overscan region (-1 = none)
trim          = 0     ; Trim pre/overscan regions before processing?
rotate        = 0     ; Rotate frame (90 deg clockwise) before processing?
darkcurrent   = [0]  ; Dark current (e-/sec/pixel), scalar or for each amplifier
badpixelfile  = -1    ; Filename of bad pixel map (-1 = none)
bk_space      = 0.7   ; Initial spacing of b-spline break points, in input wave units
min_bk_space  = 0.7   ; Minimum spacing of b-spline break points, in input wave units
slitedge1     = 547   ; First edge of slit, in slit position coords
slitedge2     = 1473   ; Second edge of slit, in slit position coords
wavemin       = 500   ; Min input wavelength over which to locate object
wavemax       = 600   ; Max input wavelength over which to locate object
slitmin       = 547   ; Min input slit position over which to locate object
slitmax       = 1473   ; Max input slit position over which to locate object
sub_width     = 100.  ; Range in slit position over which to perform sky subtraction
wavebin       = 3.    ; Binning for extracted spectrum, in angstroms 
illum_poly    = 2.    ; 1 + polynomial order for slit illumination fn (2 = linear)
x_display     = [0,2047] ; [xmin,xmax] for display of subtracted array
y_display     = [550,1500]  ; [ymin,ymax] for display of subtracted array
traceobj      = 0         ; Trace object with respect to slit position? (1 = yes)
bestprofile   = 1         ; Estimate object profile using entire spectrum? (1 = yes)
line_file     = '/Users/gwalth/idl/luci_fuyan/lowd_ir_ohlines.lst' ; Default sky line list
wave_model    = -1        ; Default tabulated sky model
getsubarray   = 0         ; Extract sub-arrays for output files? (1 = yes)
subarray_x    = [0,2047]  ; [xmin,xmax] of extracted subarray
subarray_y    = [0,2047]  ; [ymin,ymax] of extracted subarray
latitude      = 32.70   ; Observatory latitude (degrees; Keck = 19.8283; LBT = 32.7000)
longitude     = 109.89167  ; Observatory longitude (degrees; Keck = 155.478; LBT = 109.89167)
altitude      = 3200.     ; Observatory altitude (meters; Keck = 4160.; LBT = 3200)


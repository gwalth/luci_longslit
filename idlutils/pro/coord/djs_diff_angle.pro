;-----------------------------------------------------------------------
;+
; NAME:
;   djs_diff_angle
;
; PURPOSE:
;   Compute the angular distance between two points on a sphere.
;
;   Note that either (ra1,dec1) or (rap,decp) must be scalars.
;
; CALLING SEQUENCE:
;   adist = djs_diff_angle( ra, dec, ra0, dec0, [ units=units ] )
;
; INPUTS:
;   ra1:        RA of first point(s) in radians/degrees/hours
;   dec1:       DEC of first point(s) in radians/degrees
;   rap:        RA of second point(s) in radians/degrees/hours
;   decp:       DEC of second point(s) in radians/degrees
;
; OPTIONAL INPUTS:
;   units:      Set to
;                  degrees - All angles in degrees
;                  hrdeg - RA angles in hours, DEC angles and output in degrees
;                  radians - All angles in radians
;               Default to "degrees".
;
; OUTPUTS:
;   adist:      Angular distance(s) in radians/degrees
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   14-May-1997  Written by D. Schlegel, Durham
;-
;-----------------------------------------------------------------------
function djs_diff_angle, ra1, dec1, rap, decp, units=units

   DPIBY2 = 0.5d0 * !dpi

   ; Need 4 parameters
   if N_params() LT 4 then begin
      print, 'Syntax - adist = djs_diff_angle(ra1, dec1, rap, decp, [units=units] )'
      return, -1
   endif

   if (NOT keyword_set(units)) then units='degrees'

   case units of
      "hrdeg" : begin
         convRA = !dpi / 12.d0
         convDEC = !dpi / 180.d0
      end
      "radians" : begin
         convRA = 1.d0
         convDEC = 1.d0
      end
      else : begin
         convRA = !dpi / 180.d0
         convDEC = !dpi / 180.d0
      end
   endcase

   theta1 = dec1*convDEC + DPIBY2
   theta2 = decp*convDEC + DPIBY2
   cosgamma= sin(theta1) * sin(theta2) $
    * cos((ra1-rap)*convRA)  + cos(theta1) * cos(theta2)

   adist = 0.d0 * cosgamma
   ivalid = where( cosgamma LT 1 )
   if (ivalid[0] NE -1) then adist[ivalid] = acos(cosgamma[ivalid]) / convDEC

   return, adist
end
;-----------------------------------------------------------------------

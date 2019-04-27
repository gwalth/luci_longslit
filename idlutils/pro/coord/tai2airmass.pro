;+
; NAME:
;   tai2airmass
;
; PURPOSE:
;   Compute airmass.
;
; CALLING SEQUENCE:
;   airmass = tai2airmass( ra, dec, [ equinox, jd=, tai=, mjd=, $
;    longitude=, latitude=, altitude= ] )
;
; INPUTS:
;   ra             - Right ascension [degrees]
;   dec            - Declination [degrees]
;   equinox        - Equinox of observation for RA, DEC; default to 2000.
;
; OPTIONAL KEYWORDS:
;   jd             - Decimal Julian date.  Note this should probably be
;                    type DOUBLE.
;   tai            - Number of seconds since Nov 17 1858
;                    Note this should probably either be type DOUBLE or LONG64.
;   mjd            - Modified Julian date.
;   longitude      - Longitude of observatory;
;                    default to (360-105.820417) deg for APO
;   latitute       - Latitude of observatory; default to 32.780361 deg for APO
;   altitude       - Altitude of observatory; default to 2788 m for APO
;
; OUTPUTS:
;   airmass        - Airmass; 1.0 for zenith
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   TAI, JD, or MJD must be specified.
;
;   This routine only returns sec(z) for the airmass.
;   Formula from Smart, Spherical Astronomy.
;
; EXAMPLES:
;
; BUGS:
;   Outputs SLIGHTLY different airmasses from those computed by the PT
;     system.  We think that they may be going to second order.
;   Equinox does NOTHING!!!
;
; PROCEDURES CALLED:
;   ct2lst
;
; REVISION HISTORY:
;   10-May-2000  Written by D. Schlegel, Princeton, & D. Hogg, IAS
;   02-Jun-2000  Fixed minor bugs, Schlegel
;   05-Nov-2000  Added HA keyword
;-
;------------------------------------------------------------------------------
function tai2airmass, ra, dec, equinox, jd=jd, tai=tai, mjd=mjd, $
 longitude=longitude, latitude=latitude, altitude=altitude, ha=ha

   ; Default to location of Apache Point Observatory
   if (NOT keyword_set(longitude)) then longitude = 360. - 105.820417
   if (NOT keyword_set(latitude)) then latitude = 32.780361
   if (NOT keyword_set(altitude)) then altitude = 2788.

   if (NOT keyword_set(jd)) then begin
      if (keyword_set(tai)) then begin
         jd = 2400000.5D + tai / (24.D*3600.D)
      endif else if (keyword_set(mjd)) then begin
         jd = 2400000.5D + mjd
      endif
   endif

   if (NOT keyword_set(jd)) then begin
      message, 'Must specify TAI, JD or MJD', /cont
      return, 0
   endif

   DRADEG = 180.d0 / !DPI

   ;----------
   ; Compute the hour angle, HA, in degrees

   ct2lst, LST, longitude, junk, jd
   LST = 15. * LST ; convert from hours to degrees
   HA = LST - ra

   decrad = dec / DRADEG
   harad = HA / DRADEG
   latrad = latitude / DRADEG

   ; Compute airmass with spherical trig
   cosz = sin(decrad)*sin(latrad)+ cos(decrad)*cos(harad)*cos(latrad)
   airmass = 1.0/cosz

   return, airmass
end

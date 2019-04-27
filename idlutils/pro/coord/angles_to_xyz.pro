;-----------------------------------------------------------------------
;  Convert spherical coordinates (r,phi,theta) into Cartesion coordinates
;  (x,y,z).  The angles must be in the following ranges:
;    0 <= phi < 360
;    0 <= theta <= 180
;  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;  If you want to convert from RA and DEC, pass the following
;  arguments (in degrees):  RA, 90-DEC
 
pro angles_to_xyz,r,phi,theta,x,y,z
   DRADEG = 180.d0/!dpi
 
   stheta = sin(theta / DRADEG)
   x = r * cos(phi / DRADEG) * stheta
   y = r * sin(phi / DRADEG) * stheta
   z = r * cos(theta / DRADEG)
 
   return
end
;-----------------------------------------------------------------------

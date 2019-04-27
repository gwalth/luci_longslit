;-----------------------------------------------------------------------
;  Convert Cartesian coordinates (x,y,z) to spherical coordinates
;  (r,phi,theta).  The returned angles are in the following ranges:
;    0 <= phi < 360
;    0 <= theta <= 180
;  where theta=0 corresponds to the N pole, and theta=180 is the S pole.
;  Note that RA=phi and DEC=theta-90.

pro xyz_to_angles,x,y,z,r,phi,theta
   DRADEG = 180.d0/!dpi

   x2y2 = x*x + y*y
   r = sqrt(x2y2 + z*z)

   indx1 = where(x2y2 GT 0.d0)
   indx2 = where(r GT 0.d0)
   indx3 = where(r GT 0.d0 AND y LT 0.d0)

   phi = dblarr(N_elements(x))
   theta = dblarr(N_elements(x))

   if (indx1(0) NE -1) then $
    theta(indx1) = DRADEG * acos(z(indx1) / r(indx1))
   if (indx2(0) NE -1) then $
    phi(indx2) = DRADEG * acos(x(indx2) / sqrt(x2y2(indx2)))
   if (indx3(0) NE -1) then $
    phi(indx3) = phi(indx3) + 180.d0

   return
end
;-----------------------------------------------------------------------

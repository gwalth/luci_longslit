function angle_from_pairs, x1, y1, x2, y2, dmax=dmax, binsz=binsz, $
              bestsig=bestsig, angrange=angrange

   sigcut = 5.0 ; A detection is sigcut sigma above the mean
   if (NOT keyword_set(binsz)) then binsz = 1
   xyshift = [0,0]

   num1 = n_elements(x1)
   num2 = n_elements(x2)

   IF (num1 > num2) GT 32000 THEN print, '  WARNING:  Awful lot of stars...', num1, num2

   ;----------
   ; Construct an image and populate with the vector offsets for all
   ; distances between object positions X1,Y1 and X2,Y2.

;   nhalf = fix(2*dmax/binsz) + 1
   nhalf = fix(dmax/binsz)  ; I think this is what is intended - DPF
   nn = 2*nhalf + 1

   xy1 = [[x1], [y1]]
   
   if NOT keyword_set(angrange) then  angrange = [-5, 5]
   dangle = 1.
   nang = (max(angrange)-min(angrange))/dangle

   theta = (0.5+findgen(nang))/nang * (angrange[1]-angrange[0])+angrange[0]
   pk = fltarr(nang)
   for k=0, nang-1 do begin 
      img = fltarr(nn,nn)
      angrad = theta[k] *!dtor
      mm = [[cos(angrad), sin(angrad)], [-sin(angrad), cos(angrad)]]
      xyr = xy1 # mm

; Speed up change - 29 March 2001  - DPF
      xoff = fltarr(num1, num2, /nozero)
      yoff = fltarr(num1, num2, /nozero)
      for i=0L, num1-1 do begin
         xoff[i, *] = (x2 - xyr[i, 0])/float(binsz)
         yoff[i, *] = (y2 - xyr[i, 1])/float(binsz)
      endfor 
      xoff = reform(xoff, num1*num2) + (nhalf + 1)
      yoff = reform(yoff, num1*num2) + (nhalf + 1)

; IMPORTANT:  only call populate_image once - calling overhead is
;                                             large. 
;      t1 = systime(1)
;      print, 'Start: ', t1
      populate_image, img, xoff, yoff, assign='cic'
;      print, 'End: ', systime(1)-t1

      pk[k] = max(img)
      print, k, theta[k], pk[k]
   endfor

   bestsig = (max(pk, ind)-mean(img))/stddev(img)
   ang = theta[ind]

   return, ang
end


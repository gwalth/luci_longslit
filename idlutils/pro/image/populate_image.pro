;+
; NAME:
;   populate_image
;
; PURPOSE:
;   Populate a vector or image with weights at the specified positions.
;
; CALLING SEQUENCE:
;   populate_image, image, x, [y, weights=, assign=]
;
; INPUTS:
;   image      - Image vector or array
;   x          - X coordinate(s) of locations to populate, 0-indexed
;
; OPTIONAL INPUTS:
;   y          - Y coordinate(s) of locations to populate, 0-indexed
;   weights    - Weight(s) to add at each X or X,Y position
;   assign     - Assignment scheme:
;                'ngp': nearest grid point assignment; default
;                'cic': cloud-in-cell assignment
;
; OUTPUTS:
;   image      - (Modified)
;
; COMMENTS:
;   IMAGE, X, and Y are treated as floating-point values for the assignment.
;
; BUGS:
;
; PROCEDURES CALLED:
;   Dynamic link to pop_image.c
;
; REVISION HISTORY:
;   17-May-2000  Written by D. Schlegel, Princeton
;-
;------------------------------------------------------------------------------
pro populate_image, image, x, y, weights=weights, assign=assign

   npts = n_elements(x)
   if (keyword_set(y)) then $
    if (npts NE n_elements(y)) then $
     message, 'Dimensions of X and Y do not agree'
   if (keyword_set(weights)) then $
    if (npts NE n_elements(weights)) then $
     message, 'Dimensions of X and WEIGHTS do not agree'
   ndim = size(image,/n_dimen)
   if (ndim NE 1 AND ndim NE 2) then $
    message, 'Number of dimensions for IMAGE not supported'
   if (NOT keyword_set(assign)) then assign = 'ngp'
   iassign = (where(assign EQ ['ngp', 'cic']))[0]
   if (iassign EQ -1) then $
    message, 'Unknown value for ASSIGN'

   dims = size(image, /dimens)
   nx = dims[0]
   if (ndim EQ 1) then ny = 1L $
    else ny = dims[1]

   if (NOT keyword_set(y)) then y = fltarr(npts) + 0
   if (NOT keyword_set(weights)) then weights = fltarr(npts) + 1.0

   soname = filepath('libimage.so', $
    root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')

   if size(image, /tname) EQ 'FLOAT' then begin 
      retval = call_external(soname, 'pop_image', $
        npts, float(x), float(y), float(weights), nx, ny, image, iassign)
   endif else begin  
      fimage = float(image)
      retval = call_external(soname, 'pop_image', $
        npts, float(x), float(y), float(weights), nx, ny, fimage, iassign)
      image[*] = fimage[*]
   endelse 

end
;------------------------------------------------------------------------------

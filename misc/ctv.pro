PRO CTV,IMAGE,P1,P2,P3,P4, PIXMAP=PX,TITLE= title, FIELD = FLD, WINSIZE=WINSIZE
;+
; NAME:
;	CTV
; PURPOSE:
;	Load an array to the image display.  It differs from the native IDL 
;	procedure TV in the following ways:
;
;	(1)  The size and position of the image is stored in the common block
;		IMAGES.
;	(2)  The user is warned if the image array is larger than the screen 
;		size
;	(3)  The window is resized to the size of the image, unless the position
;		parameter is given.
;	(4)  Negative numbers may be used for the X and Y parameters to extract 
;		a subarray of a large image.
;	(5)  The optional TITLE and PIXMAP keywords are available
;
; CALLING SEQUENCE:
;	CTV, IMAGE
;	CTV, IMAGE, Header [,X,Y,Field]
;	CTV, IMAGE, Header [,Position]
;	CTV, IMAGE, X, Y [,Field]
;	CTV, IMage, Position
;
; INPUTS:
;	IM - Image array to be loaded to the TV screen
;
; OPTIONAL INPUT PARAMETERS:
;	USER can enter up to 4 input parameters: at most, 3 can be integer 
;	numbers and one can be a FITS header (String array).  The numerical 
;	parameters will be defined thus:
;
;	1 numerical parameter  ==> Position
;	2    ''     parameters ==> X and Y (in same order as given)
;	3    ''         ''     ==> X, Y and Field (in same order as given)
;
;	The placement of the FITS header, if there is one, can be arbitrary.
;	The FIELD parameter can also be given as a keyword parameter.
;
;	X - X position of window, starting from the lower left hand corner of
;		the screen. (If X but not Y is specified, X will be the 
;		POSITION)
;	Y - Y position of window, starting from the lower left hand corner of
;		the screen.
;		If X and Y are not specified, then the current window position
;		is used.
;		If X and Y are negative, then a subimage is extracted from IM 
;		beginning at position X,Y.    The size of the subimage is 
;		either taken from the WINSIZE keyword, or the maximum size 
;		that can fit on the screen. 
;
;	FIELD - Image window index to be used.  If omitted, the current image
;		channel (from the common block TV) is used.  FIELD can also be
;		 given as a Keyword Input (See below.)
;	HEADER   - FITS header.  Information from the header will be used to
;		create the window's title, unless a specific TITLE (keyword) is 
;		specified.
;	POSITION - Scalar giving position on screen where image will be
;		displayed.  Position 0 is in the upper left hand corner and
;		positions run from left to right, and from top to bottom.
;		See documentation for TV for more information.
;
; OPTIONAL KEYWORD INPUTS:
;	PIXMAP - A keyword parameter indicating whether the image should be
;		written to a pixmap instead of a visible window.  If set and 
;		non-zero,the image is written to a pixmap.  The PIXMAP keyword
;		may be used with any of the four calling sequences.
;	FIELD - Same as parameter FIELD above.
;	TITLE - A scalar string to be displayed as the title to the image 
;		window.
;	WINSIZE - Specified the size of the window to contain the image.
;		Useful only if (X,Y) are negative
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	The current image plane is read from the common block TV.  The lower
;	left hand corner and the window size are stored in the common block
;	IMAGES.
;
; SIDE-EFFECTS:
;	The current image plane is over-written.  
;
; RESTRICTIONS: 
;	The image array must be two-dimensional.
;
; MODIFICATION HISTORY:
;	Written, W. Landsman September, 1986
;	Adapted to workstations.  M. Greason, May 1990.
;	Subimage extraction re-installed.  M. Greason, June 1990.
;	Use GET_SCREEN_SIZE  to get max window size  W. Landsman March 1991
;	Added TITLE keyword         W. Landsman   Sep 1991
;	Added HEADER parameter.     J. D. Offenberg  Oct 1991
;	Added WINSIZE keyword       W. Landsman Dec 1995
;-
on_error,2                             ;Return to caller 
COMMON TV, chan, zoom, xroam, yroam
COMMON IMAGES, x00, y00, xsize, ysize
;
;			Check the image array for validity.
;
npar = n_params(0)
;
; Check to see which set of parameters has been called.
;
IF npar EQ 0 THEN BEGIN		;Image array supplied?
	print,'Syntax - CTV,image,[ x, y, field, /PIXMAP, TITLE=, FIELD = ]'
	RETURN
ENDIF

HDR_YN = 0		;All variables ending in _YN are booleans.
X_YN = 0		;All are initialized to zero.
Y_YN = 0
FIELD_YN = 0
POS_YN = 0
TRBL_YN = 0
NP = Intarr(4)
Index = 0

IF npar GE 2 then CTVParams,np,P1,header,hdr_yn,index,trbl_yn
IF npar GE 3 then CTVParams,np,P2,header,hdr_yn,index,trbl_yn
IF npar GE 4 then CTVParams,np,P3,header,hdr_yn,index,trbl_yn
IF npar GE 5 then CTVParams,np,P4,header,hdr_yn,index,trbl_yn
IF trbl_yn then begin	;Trouble from ParamType (Debugging routine)
	print,'Oops.  Problems in PARAMTYPE routine.'
	return
    endIF

Case INDEX of
0:	Begin		;NO PLOTTING PARAMETERS
	POS_YN = 0 & X_YN = 0 & Y_YN = 0 & Field_YN = 0
	X = 0 & Y = 0
	End
1:	Begin		;POSITION parameter
	POS_YN = 1 & X_YN = 0 & Y_YN = 0 & Field_YN = 0
	X = NP(0) & Y = 0
	end
2:	Begin		;X and Y parameter
	POS_YN = 0 & X_YN = 1 & Y_YN = 1 & Field_YN = 0
	X = NP(0) & Y = NP(1)
	end
3:	Begin		;X, Y and FIELD parameters
	POS_YN = 0 & X_YN = 1 & Y_YN = 1 & Field_YN = 1
	X = NP(0) & Y = NP(1) & FIELD = NP(2)
	end
4:	Begin		;Oops
	message,"TOO MANY PLOTTING PARAMETERS: 3 AT MOST!",/inf
	;return
	end
else:	begin
	message,"Encountered a problem.  Try again"
	endElse
endCASE
;
;Check FIELD keyword parameter
;
IF Keyword_set(FLD) then BEGIN
	FIELD = FLD
	FIELD_YN = 1 
    endIF
sz = size(image)		;Warn user if image isn't two dimensional.
IF sz(0) NE 2 THEN $
    message,'Image array (first parameter) is not 2 dimensional'
;

IF not(HDR_YN) then HEADER = ''

;
;			Check the PIXMAP keyword.
;
p = keyword_set(px)
;
;			If a window index has been supplied, make the
;			corresponding window the active window.
;
IF not(FIELD_YN) THEN field = chan
chan, field

IF (HDR_YN) then begin
	HTitle = Headerstring(Header)
   endIF $
ELSE	 $
HTitle = ' '

wtitle = 'IDL '+ strtrim(field,2)
IF keyword_set(TITLE) THEN 			$
	wtitle = wtitle + '   ' + title		$
ELSE						$
	wtitle = wtitle + '   ' + HTitle	
;
;			Initialize window variables.
;
device,get_screen_size = screen
wi = !d.window
xwsz = !d.x_vsize
ywsz = !d.y_vsize
if not keyword_set(WINSIZE) then winsize = min([sz(1),sz(2)]) < min(screen)
;
;			Check for a small window / screen.
;
rszflg = 0
IF (sz(1) GT xwsz) OR (sz(2) GT ywsz) THEN BEGIN
    IF ((sz(1) GT screen(0) ) OR (sz(2) GT screen(1))) AND $
        ((x GE 0) AND (y GE 0)) THEN $
        message,'Image array too large for the display screen'
    xwsz = sz(1)
    ywsz = sz(2)
    rszflg = 1
ENDIF

;    If POSITION was given, resize the window (if necessary)
;    and display the image using the POSITION parameter.

 IF POS_YN THEN BEGIN
;
;			Open window, display image.
;
IF Rszflg EQ 1 THEN BEGIN
    window, field, xsize=xwsz, ysize=ywsz,pixmap=p, TITLE = wtitle 
endIF  
    tv, image, x
;
;			Fill common block variables.
;
    nx = xwsz / sz(1)
    x00(chan) = sz(1) * (x MOD nx)
    y00(chan) = ywsz - (sz(2) * (1 + x / nx))
    xsize(chan) = sz(1)
    ysize(chan) = sz(2)
;
;			Otherwise, resize the window and display the image
;			without using the image parameter.
;
ENDIF ELSE BEGIN
;
;			If a negative position was specified, extract subimage.
;
    IF (x LT 0) OR (y LT 0) THEN BEGIN
        xst = abs(x) & xfn = (xst + winsize-1) < (sz(1) - 1)
        yst = abs(y) & yfn = (yst + winsize-1) < (sz(2) - 1)
        xwsz = xfn - xst + 1
        ywsz = yfn - yst + 1
        npar = -1
    ENDIF ELSE BEGIN
        xwsz = xwsz < sz(1)
        ywsz = ywsz < sz(2)
        xst = 0 & xfn = sz(1) - 1
        yst = 0 & yfn = sz(2) - 1
    ENDELSE
;
;			Open window, display image.
;
device,window_state=opnd, get_window_pos = XX, get_screen_size = Screen

if not(opnd(chan)) OR (!D.X_VSIZE NE xwsz) or (!D.Y_VSIZE NE ywsz) or $
   Keyword_set(Title) or (HDR_YN) or (X_YN) then begin

	IF not((X_YN) and (Y_YN) and (Y GE 0) and (X GE 0)) THEN BEGIN
		XZ = XX(0)
		YZ = XX(1)
	endIF ELSE begin
		XZ = X
		YZ = Y
	endELSE
 	If (XZ + xwsz) GT Screen(0) THEN XZ = Screen(0) - xwsz
	IF (YZ + ywsz) GT Screen(1) THEN YZ = Screen(1) - ywsz
	IF XZ LT 0 THEN XZ = 0
	IF YZ LT 0 THEN YZ = 0

    	window, chan, xsize=xwsz, ysize=ywsz, xpos=xz, ypos=yz,$
		pixmap=p, TITLE = wtitle 

endif
    tv, image(xst:xfn,yst:yfn)
;
;			Fill common block variables.
;
;;; Commented out 7/3/2002 GDB
;;;    x00(chan) = x < 0
;;;    y00(chan) = y < 0
;;;    zoom(chan) = 1
;;;    xsize(chan) = sz(1)
;;;    ysize(chan) = sz(2)
ENDELSE
;
RETURN
END

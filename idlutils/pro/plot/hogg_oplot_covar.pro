;-----------------------------------------------------------------------
;+
; NAME:
;   hogg_oplot_covar
; PURPOSE:
;   Wrapper for oplot that plots covariance ellipses.
; CALLING SEQUENCE:
;   hogg_oplot_covar, x,y,covar,...
; INPUT:
;   x,y     - [N] arrays of points
;   covar   - [2,2,N] array of symmetric covariance matrices
; OPTIONAL INPUTS:
;   nsigma  - plot n-sigma ellipse (default to 1)
; COMMENTS:
;   Allows color to be a [N] array.
; REVISION HISTORY:
;   2002-04-11  written by Hogg
;-
;-----------------------------------------------------------------------
pro hogg_oplot_covar,x,y,covar,color=color,nsigma=nsigma,_EXTRA=KeywordsForPlot

; set defaults, etc
  if NOT keyword_set(color) then color = !p.color
  ndata= n_elements(x)
  x= reform([x],ndata)
  y= reform([y],ndata)
  covar= reform(covar,2,2,ndata)
  ncolor= N_elements(color)
  if ncolor EQ 1 then icolor= intarr(ndata)+djs_icolor(color) $
  else icolor= djs_icolor(color)
  if NOT keyword_set(nsigma) then nsigma= 1D0

; set up useful vectors
  npoint= 100
  theta= 2D0*!PI*dindgen(npoint+2)/double(npoint)
  xx= cos(theta)
  yy= sin(theta)

; loop over points
  for ii=0L,ndata-1L do begin

; symmetrize and eigensolve covariance matrix
    tcovar= double(covar[*,*,ii])
    tcovar= 5d-1*(tcovar+transpose(tcovar))
    eval= sqrt(eigenql(tcovar,eigenvectors=evec,/double))

; make ellipses
    xxx= x[ii]+xx*nsigma*eval[0]*evec[0,0]+yy*nsigma*eval[1]*evec[0,1]
    yyy= y[ii]+xx*nsigma*eval[0]*evec[1,0]+yy*nsigma*eval[1]*evec[1,1]
    oplot, xxx,yyy,color=icolor[ii],_EXTRA=KeywordsForPlot
  endfor
end

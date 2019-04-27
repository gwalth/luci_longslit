pro s_plotdev, x, y, dev, maxsize, cap=cap, _EXTRA=KeywordsForPlot

   plot, x, y, /nodata, _EXTRA=KeywordsForPlot


   if (NOT keyword_set(cap)) then cap = max(abs(dev))

   size = abs(dev)/cap * maxsize

   t = findgen(21)*!Pi*2/20.0
   usersym, cos(t), sin(t)

   neg = where(dev LT 0)
   pos = where(dev GT 0)

   if (pos[0] NE -1) then begin
     for i=0, n_elements(pos)-1 do $
       djs_oplot, [x[pos[i]]], [y[pos[i]]], symsize=size[pos[i]], $
           ps=8,color='green'
   endif

   if (neg[0] NE -1) then begin
     for i=0, n_elements(neg)-1 do $
       djs_oplot, [x[neg[i]]], [y[neg[i]]], symsize=size[neg[i]], $
           ps=8, color='red'
   endif

   return
end

   


PRO GETQE,energy,QE

QE = fltarr(n_elements(energy))

for i=0, n_elements(energy)-1 do begin &$
  if energy(i) gt 1.0 AND energy(i) lt 8.0 then begin 
    QE(i) = 1.0
  endif else begin
    QE(i) = 0.0
  endelse 
   
endfor  
  
;print, QE
  
end
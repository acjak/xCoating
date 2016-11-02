PRO GETAXIONSPECTRUM,energy,aspectrum

aspectrum = fltarr(n_elements(energy))

for i=0, n_elements(energy)-1 do begin &$
  aspectrum(i) = (energy(i)^(2.481))*exp(-energy(i)/1.205)
endfor

amax = max(aspectrum)

for i=0, n_elements(energy)-1 do begin &$
  aspectrum(i) = aspectrum(i)/amax
endfor


END
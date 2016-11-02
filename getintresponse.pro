PRO GETINTRESPONSE,energy,RA,aspectrum,QE,dE,int_response,int_array

int_array = fltarr(n_elements(energy))
for i=0, n_elements(energy)-1 do begin &$
  
  int_array(i) = RA(i)^2 * aspectrum(i) * QE(i)
endfor

int_response = total(int_array)*dE

;  float(sum((self.idl.RA[::-1])**2 * A_spectrum * QE)/len(self.idl.energy))

END
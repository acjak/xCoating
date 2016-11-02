PRO GETNK,subs,mat1,mat2,emin,emax,evalues,nk_array,lambda,energy,dE

energy= vector(emin,emax,evalues+1)
hc = 4.1356692e-15                      ; Planck's constant, in eV
lightS = 299792458.                     ; Speed of light, in m/s
LAMBDA = ((hc*lightS)*1e7)/(energy)
dE = (emax-emin)/(evalues+1)

;mat_array = ['W_llnl_cxro','SiO_llnl_cxro','B4C_llnl_cxro','a-C_llnl_cxro','Mo_llnl_cxro','Ni.93V.07_llnl_cxro','Pt_llnl_cxro','Si_llnl_cxro']
mats = [subs,mat1,mat2]
nk_array = dcomplexarr(n_elements(mats),n_elements(LAMBDA))

for i=0, n_elements(mats)-1 do begin &$
  for j=0, n_elements(LAMBDA)-1 do begin &$
    nk_array(i,j) = IMD_NK(mats(i),LAMBDA(j))
  endfor
endfor

print,'Optical constants array done'
return
END


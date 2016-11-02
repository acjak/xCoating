PRO GETNKARRAY,MAT_ARRAY,LAMBDA,NK_ARRAY,NC

NC = dcomplexarr(n_elements(mat_array)+2,n_elements(LAMBDA))

for i=0, n_elements(LAMBDA)-1 do begin &$
  NC(0,i) = complex(1.0,0.0)
endfor

for i=0, n_elements(mat_array)-1 do begin &$
  if mat_array(i) eq 1 then begin
    NC(i+1,*) = NK_ARRAY(1,*)
    endif
  if mat_array(i) eq 2 then begin
    NC(i+1,*) = NK_ARRAY(2,*)
    endif
    
NC(n_elements(mat_array)+1,*) = NK_ARRAY(0,*)
    
endfor


end
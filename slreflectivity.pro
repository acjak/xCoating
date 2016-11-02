function getaxionspectrum,energy
	aspectrum = fltarr(n_elements(energy))

	for i=0, n_elements(energy)-1 do begin &$
	  aspectrum(i) = (energy(i)^(2.481))*exp(-energy(i)/1.205)
	endfor

	amax = max(aspectrum)

	for i=0, n_elements(energy)-1 do begin &$
	  aspectrum(i) = aspectrum(i)/amax
	endfor
	return,aspectrum
end


function getintresponse,energy,RA,aspectrum,QE,dE,int_response,int_array
	int_array = fltarr(n_elements(energy))
	for i=0, n_elements(energy)-1 do begin &$
  
	  int_array(i) = RA(i)^2 * aspectrum(i) * QE(i)
	endfor

	int_response = total(int_array)*dE
	;return,int_response,int_array
end

function getnkarray,MAT_ARRAY,LAMBDA,NK_ARRAY
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
	return, NC
end
	
function getqe,energy
	;QE = fltarr(n_elements(energy))

	;for i=0, n_elements(energy)-1 do begin &$
	;  if energy(i) gt 1.0 AND energy(i) lt 8.0 then begin 
	;    QE(i) = 1.0
	;  endif else begin
	;    QE(i) = 0.0
	;  endelse 
   
	;endfor  
	;return, qe
	readcol,'Eff_mM_No_Strongback.txt',F='F,F',file_energy_values,file_qe_values
	;datafromfile = fltarr(2,99)
    ;openr,lun,'MM_QE.txt',/get_lun
    ;readf,lun,datafromfile
    ;close,/all
        
                            ;print,datafromfile                                                                                                    
    ;file_energy_values = reform(datafromfile(0,*))
    ;file_qe_values = reform(datafromfile(1,*))

    QE = interpol(file_qe_values,file_energy_values,energy)


    return,qe
	
end

PRO SLREFLECTIVITY,alpha,D,matR,subsR,nk_array,lambda,energy,dE,RA,aspectrum,qe,int_array,int_response

	QE = getqe(energy)
	aspectrum = getaxionspectrum(energy)
	my_colors
	
	Z_array = fltarr(1)
	mat_Array = intarr(1)
	sigmaarray = fltarr(2)
	
	mat_array(0) = 1
	Z_array(0) = D
	sigmaarray(1) = matR
	sigmaarray(0) = subsR

	NC = getnkarray(mat_array,LAMBDA,nk_array)

	;;;;;;;;; FRESNEL CALCULATION ;;;;;;;;;;
	fresnel,alpha,LAMBDA,NC,Z_array,sigmaarray,1,RA=RA

	GETINTRESPONSE,energy,RA,aspectrum,QE,dE,int_response,int_array
END


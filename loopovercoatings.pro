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

function getqe,energy
	;QE = fltarr(n_elements(energy))

	;for i=0, n_elements(energy)-1 do begin &$
	 ; if energy(i) gt 1.0 AND energy(i) lt 8.0 then begin 
	 ;   QE(i) = 1.0
	 ; endif else begin
	 ;   QE(i) = 0.0
	 ; endelse 
   
	;endfor  
	;return,qe
	
	readcol,'Eff_mM_No_Strongback.txt',F='F,F',file_energy_values,file_qe_values
	
	;datafromfile = fltarr(2,99)
    ;openr,lun,'MM_QE.txt',/get_lun
	;openr,lun,'Eff_mM_No_Strongback.txt',/get_lun
    ;readf,lun,datafromfile
    ;close,/all
        
                            ;print,datafromfile                                                                                                    
    ;file_energy_values = reform(datafromfile(0,*))
    ;file_qe_values = reform(datafromfile(1,*))

    QE = interpol(file_qe_values,file_energy_values,energy)


    return,qe
	
end

function getnk,subs,mat1,mat2,emin,emax,evalues,nk_array,lambda,energy,dE
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
	    	nk_array(i,j) = IMD_NK(mats(i),LAMBDA(j));
	  	endfor
	endfor

	print,'Optical constants array done'
	;return,nk_array,lambda,energy,dE
end

function getintresponse,energy,RA,aspectrum,QE,dE,int_response,int_array
	int_array = fltarr(n_elements(energy))
	;for i=0, n_elements(energy)-1 do begin &$
  	 int_array = RA^2 * aspectrum * QE
		;int_array(i) = RA(i)^2 * aspectrum(i) * QE(i)
	;endfor

	int_response = total(int_array*dE)
	;return,int_response,int_array
end

function getnkarray,MAT_ARRAY,LAMBDA,NK_ARRAY
	NC = dcomplexarr(n_elements(mat_array)+2,n_elements(LAMBDA))

	;for i=0, n_elements(LAMBDA)-1 do begin &$
	NC(0,*) = complex(1.0,0.0)
	;endfor

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
				
function MAKEMLARRAY,dmin,dmax,gamma,n,mat1R,mat2R,subsR,Z_array,mat_array,sigmaarray
	Z_array = fltarr(n*2)
	mat_Array = intarr(n*2)
	dArray = fltarr(n)
	sigmaarray = fltarr(n*2+1)

	if n lt 2 then begin &$
		dArray = dmax
	endif else begin
		for i=0, n-1 do begin &$
			dArray(i) = dmax-i*(dmax-dmin)/(n-1)
		endfor
	endelse

	
	if n lt 2 then begin &$
		mat_array(0) = 2
		mat_array(1) = 1
		sigmaarray(0) = mat2R
		sigmaarray(1) = mat1R
		Z_array(0) = dArray*gamma
		Z_array(1) = dArray*(1-gamma)
	endif else begin
		for i=0, n-1 do begin &$
			mat_array(i*2) = 2	
			mat_array(i*2+1) = 1
			sigmaarray(i*2) = mat2R
			sigmaarray(i*2+1) = mat1R
			Z_array(i*2) = dArray(i)*gamma
			Z_array(i*2+1) = dArray(i)*(1-gamma)
		
		endfor
	endelse
	sigmaarray(-1) = subsR
	
	return,0
end

PRO LOOPOVERCOATINGS,N,GAM,DMIN,DMAX,alpha,nk_array,lambda,energy,dE,totalcalc,mat1R,mat2R,subsR,bestml_array,best_int_array,bestRA,aspectrum,qe,recipe_refl_array,recipe_ra_array,recipe_int_response_array
	
	;GETNK,subs,mat1,mat2,emin,emax,evalues,nk_array,lambda,energy,dE

	aspectrum = getaxionspectrum(energy)

	QE = getqe(energy)


	my_colors
	m = 0.0
	BestML = [[0.,0.,0.,0.,0.]]
	;nc_array_recipe = fltarr(200,n_elements(alpha))
	;z_array_recipe = fltarr(200,n_elements(alpha))
	;sigma_array_recipe = fltarr(200,n_elements(alpha))
	loadct = 5
	for i=0,n_elements(N)-1 do begin &$
		for j=0,n_elements(GAM)-1 do begin &$
			for k=0,n_elements(DMIN)-1 do begin &$
				for l=0,n_elements(DMAX)-1 do begin &$

					if DMAX(l) gt DMIN(k) OR DMAX(l) eq DMIN(k) then begin &$
						m += 1.0
						nul = MAKEMLARRAY(DMIN(k),DMAX(l),GAM(j),N(i),mat1R,mat2R,subsR,Z_array,mat_array,sigmaarray)
						NC = getnkarray(mat_array,LAMBDA,nk_array)
						
						;;;;;;;;; FRESNEL CALCULATION ;;;;;;;;;;
						fresnel,alpha(0),LAMBDA,NC,Z_array,sigmaarray,1,RA=RA

						GETINTRESPONSE,energy,RA,aspectrum,QE,dE,int_response,int_array
						
						if int_response gt BestML(0,-1) then begin &$
							BestML = [[BestML],[int_response,N(i),GAM(j),DMIN(k),DMAX(l)]]
							best_int_array = int_array
							bestRA = RA
							nc_array_best = NC
							z_array_best = Z_array
							sigma_array_best = sigmaarray

							print,Z_array

							print,sigmaarray
							
							;if n_elements(BestML) gt 50 then begin &$
							;print,m
							;print,"Best coating for the past",m," calculations:"
							;print,"Int-response: ",BestML(0,-1), " N: ",BestML(1,-1)
							;print,"Gamma: ",BestML(2,-1)," d_min: ",BestML(3,-1)," d_max: ",BestML(4,-1)
								;plot,energy,RA,col=!col.red
								;oplot,energy,int_array,col=!col.blue
								;oplot,energy,aspectrum,col=!col.green
								;oplot,energy,QE,col=!col.white
							;endif
						endif
						
						if m eq 1 then begin &$
							begintime = SYSTIME(/SECONDS)
						endif
							
						if m eq totalcalc-2 then begin &$
							endtime = SYSTIME(/SECONDS)
							;Time for 1 calculation:
							timepercalc = (endtime - begintime)/(totalcalc-2)
							;Time to calculate totalcalc:
							print,totalcalc*timepercalc/60.,' minutes'
							;print,z_array_best
							
						endif
						
						if m mod 20000 eq 0. then begin &$
						    print,m,'	/',totalcalc
							print,"N: ",BestML(1,-1)
							print,"Gamma: ",BestML(2,-1)
							print,"d_min: ",BestML(3,-1)
							print,"d_max: ",BestML(4,-1)
						endif
						
						
					endif
				endfor
			endfor
		endfor	
	endfor
	bestml_array = bestml(*,-1)
	
	recipe_refl_array = fltarr(n_elements(alpha),n_elements(energy))
	recipe_ra_array = fltarr(n_elements(alpha),n_elements(energy))
	recipe_int_response_array = fltarr(n_elements(alpha))


	for h=0,n_elements(alpha)-1 do begin &$

		fresnel,alpha(h),LAMBDA,nc_array_best,z_array_best,sigma_array_best,1,RA=RA
		
		GETINTRESPONSE,energy,RA,aspectrum,QE,dE,int_response,int_array
		
		recipe_refl_array(h,*) = int_array
		recipe_ra_array(h,*) = RA
		recipe_int_response_array(h) = int_response
	endfor
		

end




					




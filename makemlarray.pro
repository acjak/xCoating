PRO MAKEMLARRAY,dmin,dmax,gamma,n,mat1R,mat2R,subsR,Z_array,mat_array,sigmaarray
	
	Z_array = fltarr(n*2)
	mat_Array = intarr(n*2)
	dArray = fltarr(n)
	sigmaarray = fltarr(n*2+1)
	
	for i=0, n-1 do begin &$
		dArray(i) = dmin+i*(dmax-dmin)/(n-1)
	endfor
	
	
	
	for i=0, n-1 do begin &$
		mat_array(i*2) = 2	
		mat_array(i*2+1) = 1
		sigmaarray(i*2) = mat2R
		sigmaarray(i*2+1) = mat1R
		Z_array(i*2) = dArray(i)*gamma
		Z_array(i*2+1) = dArray(i)*(1-gamma)
		
	endfor
	sigmaarray(-1) = subsR
	

	
end
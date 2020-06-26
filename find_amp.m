function a = find_amp(A)
	a = 0;
	for  i = 1: length(A)
			a = a + A(i);
	end
	a = a/ length(A);

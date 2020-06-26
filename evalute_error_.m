% функция вычисления ошибки оценки угловой координаты для каждого выбранного метода
function [err1 var1 err2 var2 err3 var3 err4 var4 err5 var5 err6 var6 err7 var7 err8 var8] = evalute_error_(Err_method, Sko_method)

	N = length(Err_method);
	err1=0; var1 = 0;  err2 = 0; var2 =0;  err3 = 0; var3 = 0; err4 = 0; var4 = 0; err5 = 0; var5 = 0; err6 = 0; var6 = 0;
	err7 = 0;  var7 = 0;  err8 = 0; var8 = 0;

	for  i = 1:N
		 eval('base', ['err',int2str(i), '= Err_method(', int2str(i),');']);
		 eval('base', ['var',int2str(i), '= Sko_method(', int2str(i),');']);
	end
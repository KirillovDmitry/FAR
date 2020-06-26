
function evaluate 
	for met = 1:8
		if ~strcmp('NULL',  eval(['method',int2str(met)]))
			ogr_ = min(eval(['NumberSignal',int2str(met)]), ogr);
			for temp = 1:length(eval(['eps_method',int2str(met)]))
				evalin('base', [ 'Angle_', (['method',int2str(met)]),'_', int2str(temp),'(mm,zz) = eps_',...
					(['method',int2str(met)]), '(',int2str(temp), ');']);
			end

			for temp = 1:length(eval(['P',int2str(met)]))
				evalin('base', [ 'Power_', (['method',int2str(met)]),'_', int2str(temp),'(mm,zz) = P',...
					int2str(met), '(',int2str(temp), ');']);
			end

			string_temp = '';
			evalin('base', [ '[Err_', (['method',int2str(met)]), ' Sko_', (['method',int2str(met)]),...
				'] = eror(eps, eps_', (['method',int2str(met)]), '(1:ogr_), ogr_, alpha, regl);']);
			for temp = 1:length(eval(['P',int2str(met)]))
				string_temp = [string_temp, ' Err_', (['method',int2str(met)]), '_', int2str(temp),...
					'(mm,zz)', ' Sko_', (['method',int2str(met)]),  '_', int2str(temp),'(mm,zz)'];
			end

			string_temp = ['[', string_temp, ' ]'];
			evalin('base', [ string_temp, ' = evalute_eror(Err_', (['method',int2str(met)]), ', Sko_',...
				(['method',int2str(met)]),  ');']);
			evalin('base', [ 'Iter', int2str(met), '(mm,zz) = Iteration', int2str(met), ';']);
			evalin('base', [ 'NSignal', int2str(met), '(mm,zz) = NumberSignal', int2str(met), ';']);
		end
	end

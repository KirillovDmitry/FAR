% вспомогательное окно для задания исследуемых алгоритмов
function R = GetMethod


    Fig2 = figure('MenuBar', 'None',...
        'Name', 'Tune',...
        'NumberTitle', 'Off',...
        'Position', [20 180 145 290]);

    m1 = uicontrol('Style', 'RadioButton',   'String', 'MpOneMethod', 'Tag', 'MpOneMethod','Position', [10, 260, 130, 20], 'Value', 0);
    m2 = uicontrol('Style', 'RadioButton',   'String', 'MpOneIterMethod', 'Tag', 'MpOneIterMethod','Position', [10, 240, 130, 20], 'Value', 0);
    m3 = uicontrol('Style', 'RadioButton',   'String', 'MpMethodErm', 'Tag', 'MpMethodErm','Position', [10, 220, 130, 20], 'Value', 0);
    m4 = uicontrol('Style', 'RadioButton',   'String', 'MpIterMethodErm', 'Tag', 'MpIterMethodErm','Position', [10, 200, 130, 20], 'Value', 0);
    m5 = uicontrol('Style', 'RadioButton',   'String', 'LinMethodOdn', 'Tag', 'LinMethodOdn','Position', [10, 180, 130, 20], 'Value', 0);
    m6 = uicontrol('Style', 'RadioButton',   'String', 'LinMethod', 'Tag', 'LinMethod','Position', [10, 160, 130, 20], 'Value', 1);
    m7 = uicontrol('Style', 'RadioButton',   'String', 'PsiMethod', 'Tag', 'PsiMethod','Position', [10, 140, 130, 20], 'Value', 0);
    m8 = uicontrol('Style', 'RadioButton',   'String', 'Psi4Method', 'Tag', 'Psi4Method','Position', [10, 120, 130, 20], 'Value', 0);
    m9 = uicontrol('Style', 'RadioButton',   'String', 'PSI4+', 'Tag', 'PSI4','Position', [10, 100, 130, 20], 'Value', 0);
    m10 = uicontrol('Style', 'RadioButton',  'String', 'PolpOne', 'Tag', 'PolpOne','Position', [10, 80, 130, 20], 'Value', 0);
    m11 = uicontrol('Style', 'RadioButton',  'String', 'PolpOneIter', 'Tag', 'PolpOneIter','Position', [10, 60, 130, 20], 'Value', 0);
    m12 = uicontrol('Style', 'RadioButton',  'String', 'KeponOne', 'Tag', 'KeponOne','Position', [10, 40, 130, 20], 'Value', 0);
    
   
    b1    = get(m1, 'Value');
    b2    = get(m2, 'Value');
    b3    = get(m3, 'Value');
    b4    = get(m4, 'Value');
    b5    = get(m5, 'Value');
    b6    = get(m6, 'Value');
    b7    = get(m7, 'Value');
    b8    = get(m8, 'Value');
    b9    = get(m9, 'Value');
    b10   = get(m10, 'Value');
    b11   = get(m11, 'Value');
    b12   = get(m12, 'Value');

    R = [b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 0];
    handles.data = R;
    uicontrol('Style', 'PushButton',  'String', 'all', 'Tag', 'all','Position', [10, 20, 130, 20], 'CallBack', @all)
    uicontrol('Style', 'PushButton',  'String', 'continue', 'Tag', 'continue','Position', [10, 0, 130, 20], 'CallBack', @asd)

    uiwait;
    close(Fig2);
    R = handles.data;

	function asd(src,evt)
		b1 = get(m1, 'Value');
		b2 = get(m2, 'Value');
		b3 = get(m3, 'Value');
		b4 = get(m4, 'Value');
		b5 = get(m5, 'Value');
		b6 = get(m6, 'Value');
		b7 = get(m7, 'Value');
		b8 = get(m8, 'Value');
		b9 = get(m9, 'Value');
		b10 = get(m10, 'Value');
		b11 = get(m11, 'Value');
		b12 = get(m12, 'Value');
		R = [b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 0];
		handles = guidata(src);
		handles.data = R;
		%guidata(handles, src);
		uiresume;
	end

	function all(src,evt)
		b1 = ~get(m1, 'Value');
		set(m1, 'Value',b1);
		set(m2, 'Value',b1);
		set(m3, 'Value',b1);
		set(m4, 'Value',b1);
		set(m5, 'Value',b1);
		set(m6, 'Value',b1);
		set(m7, 'Value',b1);
		set(m8, 'Value',b1);
		set(m9, 'Value',b1);
		set(m10, 'Value',b1);
		set(m11, 'Value',b1);
		set(m12, 'Value',b1);

		R = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 0];
		handles = guidata(src);
		handles.data = R;
	end
end


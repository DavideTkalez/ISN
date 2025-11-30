function matrix = apri_dati(txt_file)

data = table2array(readtable(txt_file, VariableNamingRule='preserve'));
matrix = zeros(size(data));
matrix(:, 6) = data(:, 1);
matrix(:, 1:5) = data(:, 2:6);
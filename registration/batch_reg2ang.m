function [] = batch_reg2ang(root, trf_var)
%% batch_reg2ang: calculates & saves angles from registered videos
%
%   INPUT:
%       root    : root directory
%       trf_var : affine 2D variable name
%
%   OUTPUT:
%       -
%

[FILES, PATH] = uigetfile({'*.mat'},'Select registered videos', root, 'MultiSelect','on');
FILES = string(FILES);
n_file = length(FILES);

angdir = fullfile(PATH,'reg_angles');
mkdir(angdir)
for n = 1:n_file
    disp(FILES(n))
    data = load(fullfile(PATH,FILES(n)), trf_var);
    
    % Calculate angles
    [angles] = reg2ang(data.(trf_var));
    %pause
    
    % Save
    save(fullfile(angdir,FILES{n}), '-v7.3', 'angles')
end
disp('ALL DONE')

end
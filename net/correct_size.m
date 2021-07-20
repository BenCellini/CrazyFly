function [] = correct_size(root, new_sz)
%% correct_label: correct image size
%
%   INPUT:
%       root : root directory
%
%   OUTPUT:
%

[FILES,PATH] = uigetfile({'*.jpg'},'Select incorrect labeled data', root, 'MultiSelect','on');
FILES = cellstr(FILES);
n_file = length(FILES);

for n = 1:n_file
    %disp(n)
    fpath = fullfile(PATH, FILES{n});
    I = imread(fpath);
    Irz = imresize(I, new_sz);
  	imwrite(Irz, fpath)    
end

end
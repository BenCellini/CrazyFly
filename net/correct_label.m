function [] = correct_label(root)
%% correct_label: correct image headings manually
%
%   INPUT:
%       root    	:   root directory with "Up & "Down" directories
%
%   OUTPUT:
%

% Directories with images and flipped images
updir = fullfile(root, 'Up');
downdir = fullfile(root, 'Down');

[FILES,~] = uigetfile({'*.jpg'},'Select incorrect labeled data', updir, 'MultiSelect','on');
FILES = cellstr(FILES);
n_file = length(FILES);

for n = 1:n_file
    disp(n)
    uppath = fullfile(updir, FILES{n});
    downpath = fullfile(downdir, FILES{n});
    
    upI = imread(uppath);
    downI = imread(downpath);
    %[~,basename,~] = fileparts(FILES{n});
        
    imwrite(upI, downpath)
    imwrite(downI, uppath)    
end

end
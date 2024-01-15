function [] = batch_headtrack_roll(root_yaw, root_vid)
%% batch_headtrack_roll: runs head roll tracker for user selected video files
%
%   INPUT:
%       root_yaw   	: root directory with yaw angles
%       root_vid   	: root directory with videos corresponding to yaw angles
%
%   OUTPUT:
%       -
%

[FILES, path_yaw] = uigetfile({'*.mat', 'MAT-files'},'Select yaw angles', root_yaw, 'MultiSelect','on');
FILES = cellstr(FILES);
nfile = length(FILES);

if nargin == 1
    sub_fold = regexp(path_yaw, '\', 'split');
    path_vid = fullfile(sub_fold{1:end-2});
elseif nargin == 2
    path_vid = root_vid;
end

eye_ratio = 0.3;
roll_cal = 36.33;

headdir = fullfile(path_yaw,'tracked_head_roll');
mkdir(headdir)
for f = 1:nfile
    disp(FILES(f))
    disp('---------------------------------------')
    load(fullfile(path_yaw,FILES{f}),'head_data','head_mask')
    load(fullfile(path_vid,FILES{f}),'regvid')

    [roll, roll_idx] = track_head_roll(regvid, head_data.angle, head_mask.move_points.rot, ...
        eye_ratio, roll_cal, false);
    
    figure (100) ; clf ; hold on
        plot(roll, 'k')
        plot(hampel(1:length(roll), roll, 20), 'r')
        
	x = input('Save?: ');
    switch x
        case 1
            disp('Saving...')
            save(fullfile(headdir,FILES{f}),'-v7.3','roll', 'roll_idx')
        otherwise
            disp('Not saving...')
    end
                                                                       
end
disp('ALL DONE')
end
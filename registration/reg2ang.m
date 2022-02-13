function [ang] = reg2ang(trf)
%% reg2ang: computes angle from registration affine 2D object
%
%   INPUT:
%       trf     : affine 2D object 
%
%   OUTPUT:
%       ang     : angle
%

T11 = cellfun(@(x) x.T(1,1), trf);
T12 = cellfun(@(x) x.T(1,2), trf);
ang = rad2deg(unwrap(atan2(T12, T11)));

end
function pth = cortical_obs_rootpath()
% ROOTPATH - returns the path to the uppermost containing directory of this
% project    
    pth=which('rootpath');    
    
    pth=fileparts(fileparts(pth));
end
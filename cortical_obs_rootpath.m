function pth = cortical_obs_rootpath()
% ROOTPATH - returns the path to the uppermost containing directory of this
% project    
    pth=which('cortical_obs_rootpath');    
    
    pth=fileparts(pth);
end
function [r, v] = pcp_get_rv(bodyName, et)
% 
% Author: Sourav Ghosh, Intelligent Space Systems Laboratory, The University of Tokyo.
% 
% Return heliocentric (centered at SUN) state in J2000 (km, km/s).
name = upper(strtrim(bodyName));
try
    st = cspice_spkezr(name, et, 'J2000', 'NONE', 'SUN');
catch ME
    % If the center body isn't in the SPK (common with DE ephemerides),
    % retry using the planetary barycenter name.
    if contains(ME.message,'SPKINSUFFDATA') && ~contains(name,'BARYCENTER')
        st = cspice_spkezr([name ' BARYCENTER'], et, 'J2000', 'NONE', 'SUN');
    else
        rethrow(ME);
    end
end
r = st(1:3);
v = st(4:6);
end

function kernels

% Leapsecond kernel;
cspice_furnsh('kernels/naif0012.tls');
% Gravity constant kernel;
cspice_furnsh('kernels/gm_de440.tpc');
% Planet constant kernel;
cspice_furnsh('kernels/pck00011.tpc');
cspice_furnsh('kernels/de440.bsp');

end
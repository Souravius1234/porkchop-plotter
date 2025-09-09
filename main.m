kernels;  % your loader (naif0012.tls, gm_de440.tpc, pck00011.tpc, de440.bsp)

% G = porkchop_plotter('EARTH','MARS',...
%       '2070-01-01','2088-01-01',...   % departures
%       '2080-01-01','2092-01-01',...   % arrivals
%       'Ndep',500,'Narr',500,'LW','both','TOFminDays',50,'TOFmaxDays',5000);

% kernels;  % your loader (tls, tpc, pck, bsp already handled here)

G = porkchop_plotter('EARTH','MARS', ...
      '2005-04-30','2005-10-17', ... % Departure Dates range
      '2005-11-16','2006-12-21', ... % Arrival Dates range
      'Ndep',200,'Narr',200,'LW','both','TOFminDays',50,'TOFmaxDays',500);

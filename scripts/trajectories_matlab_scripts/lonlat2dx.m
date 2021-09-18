% $Header: /tphs1/user2/mlosch/cvsroot/FEMSECT/src/lonlat2dx.m,v 1.2 2003/11/06 14:25:00 mlosch Exp $
% $Name:  $

function dx = lonlat2dx(lon,lat)
%function dx = lonlat2dx(lon,lat)
%
% returns distance (in meters) between postions pairs (lon,lat) on the earth
% (sphere). length(dx) = length(lon)-1
  
  %earth=6371000;
  earth=6367.5e3;  %feom parameter
  deg2rad=pi/180;
  nx=length(lon);
  lt = deg2rad*lat;
  ln = deg2rad*lon;
  alpha(1:nx-1) = ...
     acos( cos(lt(1:nx-1)).*cos(lt(2:nx)).*cos(ln(1:nx-1)-ln(2:nx)) ...
	 + sin(lt(1:nx-1)).*sin(lt(2:nx)) );
  dx = earth*abs(alpha');

  return
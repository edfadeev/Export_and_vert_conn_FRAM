function newmap = bluewhite(n)

%Define colormap
c1=[1 1 1]; %W
c2=[0 0 1]; %B
% c3=[1 0 0]; %R

cmap=[linspace(c1(1),c2(1),n);linspace(c1(2),c2(2),n);linspace(c1(3),c2(3),n)];

newmap=cmap';

% in this script, we compute:
% (1) the bathymetry at mooring positions in lon and lat directions
% (2) the ice concentration (mean Jul 2016) at mooring positions in lon and lat directions 
% needed for plot of vertical distribution of particles (grey line and grey bar)


close all
clear all

path(path,'/home/ollie/cwekerle/bin/m_map/')

datapath_traj='/home/ollie/cwekerle/code_post/trajectories_fram_back_eddie2/';

% get mooring locations
blon_EG=ncread([datapath_traj,'result/drifter_Arc22_start2009_EG_bottom_speed60md.nc'],'blon');
blat_EG=ncread([datapath_traj,'result/drifter_Arc22_start2009_EG_bottom_speed60md.nc'],'blat');
x_moo_EG=blon_EG(1,1);
y_moo_EG=blat_EG(1,1);

blon_HG=ncread([datapath_traj,'result/drifter_Arc22_start2009_HG_bottom_speed60md.nc'],'blon');
blat_HG=ncread([datapath_traj,'result/drifter_Arc22_start2009_HG_bottom_speed60md.nc'],'blat');
x_moo_HG=blon_HG(1,1);
y_moo_HG=blat_HG(1,1);

blon_N=ncread([datapath_traj,'result/drifter_Arc22_start2009_N_bottom_speed60md.nc'],'blon');
blat_N=ncread([datapath_traj,'result/drifter_Arc22_start2009_N_bottom_speed60md.nc'],'blat');
x_moo_N=blon_N(1,1);
y_moo_N=blat_N(1,1);


bathy = 0
ice = 1

if bathy

% bathymetry from Schaffer at al. 2016
datapath = '/work/ollie/clidyn/topography/';
topofile = [datapath,'RTopo-2.0.1_30sec_bedrock_topography_2016-12-13_fram_subset.nc']; 
topo = ncread(topofile,'bedrock_topography');
lon = ncread(topofile,'lon');
lat = ncread(topofile,'lat');
[lon_m,lat_m]=meshgrid(lon,lat); lon_m = lon_m'; lat_m = lat_m';



figure
pcolor(lon_m,lat_m,topo); shading flat;
hold on
plot(x_moo_EG,y_moo_EG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(x_moo_HG,y_moo_HG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(x_moo_N,y_moo_N,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')




[minval,ind_lon_EG] = min(abs(lon-x_moo_EG))
[minval,ind_lat_EG] = min(abs(lat-y_moo_EG))

dep_lon_EG = double([lon_m(:,ind_lat_EG) topo(:,ind_lat_EG)]); 
dep_lat_EG = [lat_m(ind_lon_EG,:)' topo(ind_lon_EG,:)']; 

figure
plot(dep_lon_EG(:,1),dep_lon_EG(:,2))
figure
plot(dep_lat_EG(:,1),dep_lat_EG(:,2))


% HG

[minval,ind_lon_HG] = min(abs(lon-x_moo_HG))
[minval,ind_lat_HG] = min(abs(lat-y_moo_HG))

dep_lon_HG = [lon_m(:,ind_lat_HG) topo(:,ind_lat_HG)]; 
dep_lat_HG = [lat_m(ind_lon_HG,:)' topo(ind_lon_HG,:)']; 

figure
plot(dep_lon_HG(:,1),dep_lon_HG(:,2))
figure
plot(dep_lat_HG(:,1),dep_lat_HG(:,2))

% N

[minval,ind_lon_N] = min(abs(lon-x_moo_N))
[minval,ind_lat_N] = min(abs(lat-y_moo_N))

dep_lon_N = [lon_m(:,ind_lat_N) topo(:,ind_lat_N)]; 
dep_lat_N = [lat_m(ind_lon_N,:)' topo(ind_lon_N,:)']; 

figure
plot(dep_lon_N(:,1),dep_lon_N(:,2))
figure
plot(dep_lat_N(:,1),dep_lat_N(:,2))

% save
% 
% save('data/dep_lon_EG.mat','dep_lon_EG','dep_lat_EG',...
%     'dep_lon_HG','dep_lat_HG',...
%     'dep_lon_N','dep_lat_N')

end

%%% sea ice

if ice

datapath_ice='/work/ollie/cwekerle/data/data_seaice_ifremer/';
icefile = [datapath_ice,'20160701-20160731.nc'];  % mean July 2016
ice_conc = ncread(icefile,'concentration');
lon_ice_aux = ncread([datapath_ice,'grid_north_12km.nc'],'longitude');
lat_ice = ncread([datapath_ice,'grid_north_12km.nc'],'latitude');
lon_ice=mod((lon_ice_aux+180),360)-180;


figure
m_proj('stereographic','lat',90,'long',30,'radius',25);
[X,Y]=m_ll2xy(lon_ice,lat_ice);
pcolor(X,Y,ice_conc); shading flat;
m_grid;
m_coast;
% hold on
% plot(x_moo_EG,y_moo_EG,'--gs',...
%     'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
% hold on
% plot(x_moo_HG,y_moo_HG,'--gs',...
%     'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
% hold on
% plot(x_moo_N,y_moo_N,'--gs',...
%     'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')

% EG

lon_EG_x = [-180:0.05:180];%[0:0.05:360];
lon_EG_y = zeros(size(lon_EG_x)); lon_EG_y(:)=y_moo_EG;

lat_EG_y = [75:0.05:85];
lat_EG_x = zeros(size(lat_EG_y)); lat_EG_x(:)=x_moo_EG;

ice_lon_EG = griddata(double(reshape(lon_ice,1,[])),double(reshape(lat_ice,1,[])),...
    double(reshape(ice_conc,1,[])),lon_EG_x,lon_EG_y);
ice_lat_EG = griddata(double(reshape(lon_ice,1,[])),double(reshape(lat_ice,1,[])),...
    double(reshape(ice_conc,1,[])),lat_EG_x,lat_EG_y);

ice_lon_EG_all = [lon_EG_x' lon_EG_y' ice_lon_EG' ];
ice_lat_EG_all = [lat_EG_x' lat_EG_y' ice_lat_EG' ];

% HG
lon_HG_x = [-180:0.05:180];%[0:0.05:360];
lon_HG_y = zeros(size(lon_HG_x)); lon_HG_y(:)=y_moo_HG;

lat_HG_y = [75:0.05:85];
lat_HG_x = zeros(size(lat_HG_y)); lat_HG_x(:)=x_moo_HG;

ice_lon_HG = griddata(double(reshape(lon_ice,1,[])),double(reshape(lat_ice,1,[])),...
    double(reshape(ice_conc,1,[])),lon_HG_x,lon_HG_y);
ice_lat_HG = griddata(double(reshape(lon_ice,1,[])),double(reshape(lat_ice,1,[])),...
    double(reshape(ice_conc,1,[])),lat_HG_x,lat_HG_y);

ice_lon_HG_all = [lon_HG_x' lon_HG_y' ice_lon_HG' ];
ice_lat_HG_all = [lat_HG_x' lat_HG_y' ice_lat_HG' ];

% N
lon_N_x = [-180:0.05:180];%[0:0.05:360];
lon_N_y = zeros(size(lon_N_x)); lon_N_y(:)=y_moo_N;

lat_N_y = [75:0.05:85];
lat_N_x = zeros(size(lat_N_y)); lat_N_x(:)=x_moo_N;

ice_lon_N = griddata(double(reshape(lon_ice,1,[])),double(reshape(lat_ice,1,[])),...
    double(reshape(ice_conc,1,[])),lon_N_x,lon_N_y);
ice_lat_N = griddata(double(reshape(lon_ice,1,[])),double(reshape(lat_ice,1,[])),...
    double(reshape(ice_conc,1,[])),lat_N_x,lat_N_y);

ice_lon_N_all = [lon_N_x' lon_N_y' ice_lon_N' ];
ice_lat_N_all = [lat_N_x' lat_N_y' ice_lat_N' ];

figure
m_proj('stereographic','lat',90,'long',30,'radius',25);
[X,Y]=m_ll2xy(lon_ice,lat_ice);
pcolor(X,Y,ice_conc); shading flat;
m_grid;
m_coast;
caxis([0 100])
colorbar
% saveas2('plots/ice_con_072016.png',300)


figure
a = 10;
[XX,YY]=m_ll2xy(lon_EG_x,lon_EG_y);
scatter(XX,YY,a,ice_lon_EG,'filled')
hold on
[XX,YY]=m_ll2xy(lat_EG_x,lat_EG_y);
scatter(XX,YY,a,ice_lat_EG,'filled')

hold on
[XX,YY]=m_ll2xy(lon_HG_x,lon_HG_y);
scatter(XX,YY,a,ice_lon_HG,'filled')
hold on
[XX,YY]=m_ll2xy(lat_HG_x,lat_HG_y);
scatter(XX,YY,a,ice_lat_HG,'filled')

hold on
[XX,YY]=m_ll2xy(lon_N_x,lon_N_y);
scatter(XX,YY,a,ice_lon_N,'filled')
hold on
[XX,YY]=m_ll2xy(lat_N_x,lat_N_y);
scatter(XX,YY,a,ice_lat_N,'filled')

caxis([0 100])
colorbar
m_grid;
m_coast;
% saveas2('plots/ice_con_072016_EG_HG_N.png',300)


save('data/ice_EG_HG_N.mat','ice_lon_EG_all','ice_lat_EG_all',...
    'ice_lon_HG_all','ice_lat_HG_all',...
    'ice_lon_N_all','ice_lat_N_all')

end





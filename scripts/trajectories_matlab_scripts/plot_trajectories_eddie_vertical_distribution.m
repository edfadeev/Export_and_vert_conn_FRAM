%%% this script computes the vertical distribution of particles
%%% both in longitude and latitude direction


close all
clear all

%path(path,'/home/ollie/cwekerle/bin/m_map/')

%datapath='/home/ollie/cwekerle/code_post/trajectories_fram_back_eddie2/';
datapath='';



%%% time period: Mar - Jul 2016 at the surface
st_mar1=60; % day of year march 1
en_jul31=212;  % day of year july 31
num_days=en_jul31-st_mar1+1;
%%% then for each location, we plot only the days when the
%%% particle is at the surface between days 60 and 212

%%% make a grid (distance x depth)
% two grids: lon x depth and lat x depth
xmin_lon=-20;
xmax_lon=20;
xmin_lat=72;
xmax_lat=86;
ymin=-4000;
ymax=0;
xres=0.05;
yres=25.;
x_lon=[xmin_lon:xres:xmax_lon];
x_lat=[xmin_lat:xres:xmax_lat];
y=[ymin:yres:ymax];

[X_lon,Y_lon]=meshgrid(x_lon,y);
[X_lat,Y_lat]=meshgrid(x_lat,y);


% initialize arrays
distr_EG_lon=zeros(size(X_lon));
distr_EG_lat=zeros(size(X_lat));

distr_HG_lon=zeros(size(X_lon));
distr_HG_lat=zeros(size(X_lat));

distr_N_lon=zeros(size(X_lon));
distr_N_lat=zeros(size(X_lat));


% start reading the data (EG, HG, N; bottom)
loc={'EG';'HG';'N'}
dep={'bottom'}
depid=char(dep(1,:))


%%% EG
locid=char(loc(1,:))
sp={'52'}
spid=char(sp(1,:))

% read only data that is at the surface between 1 mar and 31 jul
A=load([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md_ice_dist_len.txt']);
% find out the indices that correspond to 1 mar - 31 jul
ind_time_A=find( A(:,3)>=st_mar1 & A(:,3)<=en_jul31);
ind_time_A_st=min(ind_time_A)
ind_time_A_en=max(ind_time_A)
aux=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday');
len_traj=length(aux);
start=[1 ind_time_A_st];
count=[len_traj ind_time_A_en-ind_time_A_st+1];

blon_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blon',start,count);
blat_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blat',start,count);
btemp_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'btemp',start,count);
bsalt_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bsalt',start,count);
bday_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday',start,count);
bdepth_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bdepth',start,count);
byear_EG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'byear',start,count);
% coordinates and depth of mooring EG
x_moo_EG=blon_EG(1,1);
y_moo_EG=blat_EG(1,1);
dep_moo_EG=bdepth_EG(1,1);



%%% HG
locid=char(loc(2,:))
sp={'29'}
spid=char(sp(1,:))

% read only data that is at the surface between 1 mar and 31 jul
B=load([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md_ice_dist_len.txt']);
ind_time_B=find( B(:,3)>=st_mar1 & B(:,3)<=en_jul31);
ind_time_B_st=min(ind_time_B)
ind_time_B_en=max(ind_time_B)
aux=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday');
len_traj=length(aux);
start=[1 ind_time_B_st];
count=[len_traj ind_time_B_en-ind_time_B_st+1];

blon_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blon',start,count);
blat_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blat',start,count);
btemp_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'btemp',start,count);
bsalt_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bsalt',start,count);
bday_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday',start,count);
bdepth_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bdepth',start,count);
byear_HG=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'byear',start,count);
% coordinates and depth of mooring HG
x_moo_HG=blon_HG(1,1);
y_moo_HG=blat_HG(1,1);
dep_moo_HG=bdepth_HG(1,1);

%%% N
locid=char(loc(3,:))
sp={'52'}
spid=char(sp(1,:))

% read only data that is at the surface between 1 mar and 31 jul
C=load([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md_ice_dist_len.txt']);
ind_time_C=find( C(:,3)>=st_mar1 & C(:,3)<=en_jul31);
ind_time_C_st=min(ind_time_C)
ind_time_C_en=max(ind_time_C)
aux=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday');
len_traj=length(aux);
start=[1 ind_time_C_st];
count=[len_traj ind_time_C_en-ind_time_C_st+1];

blon_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blon',start,count);
blat_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blat',start,count);
btemp_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'btemp',start,count);
bsalt_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bsalt',start,count);
bday_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday',start,count);
bdepth_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bdepth',start,count);
byear_N=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'byear',start,count);
% coordinates and depth of mooring N
x_moo_N=blon_N(1,1);
y_moo_N=blat_N(1,1);
dep_moo_N=bdepth_N(1,1);

% finished reading trajectory data...


% compute the vertical distribution in lon and lat direction
st=1
en=num_days %153
for ii=st:en % loop over all trajectories that are at the surface between 1 mar and 31 jul
    
    aux_EG_lon=zeros(size(X_lon));
    aux_HG_lon=zeros(size(X_lon));
    aux_N_lon=zeros(size(X_lon));
    aux_EG_lat=zeros(size(X_lat));
    aux_HG_lat=zeros(size(X_lat));
    aux_N_lat=zeros(size(X_lat));
    
    
    ind_EG=max(find(byear_EG(:,ii)>2000))-1; % length of the trajectory
    for jj=1:ind_EG % loop over the length of trajectory
        
        % for every point of the trajectory (lon/depth or lat/depth) we check in which grid box we are
        % this is done by function getgrid.m
        point_x_lon=blon_EG(jj,ii);
        point_x_lat=blat_EG(jj,ii);
        point_y=-bdepth_EG(jj,ii);
        [u v]=getgrid(point_x_lon, point_y, xmin_lon, xmax_lon, xres, ymin, ymax, yres);
        aux_EG_lon(v,u)=aux_EG_lon(v,u)+1;
        [u v]=getgrid(point_x_lat, point_y, xmin_lat, xmax_lat, xres, ymin, ymax, yres);
        aux_EG_lat(v,u)=aux_EG_lat(v,u)+1;
    end
    
    ind_HG=max(find(byear_HG(:,ii)>2000))-1;
    for jj=1:ind_HG
        point_x_lon=blon_HG(jj,ii);
        point_x_lat=blat_HG(jj,ii);
        point_y=-bdepth_HG(jj,ii);
        [u v]=getgrid(point_x_lon, point_y, xmin_lon, xmax_lon, xres, ymin, ymax, yres);
        aux_HG_lon(v,u)=aux_HG_lon(v,u)+1;
        [u v]=getgrid(point_x_lat, point_y, xmin_lat, xmax_lat, xres, ymin, ymax, yres);
        aux_HG_lat(v,u)=aux_HG_lat(v,u)+1;
    end
    
    ind_S=max(find(byear_N(:,ii)>2000))-1;
    for jj=1:ind_S
        point_x_lon=blon_N(jj,ii);
        point_x_lat=blat_N(jj,ii);
        point_y=-bdepth_N(jj,ii);
        [u v]=getgrid(point_x_lon, point_y, xmin_lon, xmax_lon, xres, ymin, ymax, yres);
        aux_N_lon(v,u)=aux_N_lon(v,u)+1;
        [u v]=getgrid(point_x_lat, point_y, xmin_lat, xmax_lat, xres, ymin, ymax, yres);
        aux_N_lat(v,u)=aux_N_lat(v,u)+1;
    end
    
    % the value can only be 1 or 0 (1: particle passed through grid box; 0: particle didn't pass through grid box)
    % this is the distribution for one trajectory:
    aux_EG_lon(aux_EG_lon>1)=1;
    aux_EG_lat(aux_EG_lat>1)=1;
    % this is the the distribution for all 153 (num_days) trajectories
    % (when we plot this field, we divide every grid box by num_days to get the relative number of particles)
    distr_EG_lon=distr_EG_lon+aux_EG_lon;
    distr_EG_lat=distr_EG_lat+aux_EG_lat;
    
    aux_HG_lon(aux_HG_lon>1)=1;
    aux_HG_lat(aux_HG_lat>1)=1;
    distr_HG_lon=distr_HG_lon+aux_HG_lon;
    distr_HG_lat=distr_HG_lat+aux_HG_lat;
    
    aux_N_lon(aux_N_lon>1)=1;
    aux_N_lat(aux_N_lat>1)=1;
    distr_N_lon=distr_N_lon+aux_N_lon;
    distr_N_lat=distr_N_lat+aux_N_lat;
end


% plot mooring locations
figure
plot(x_moo_EG,y_moo_EG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(x_moo_HG,y_moo_HG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(x_moo_N,y_moo_N,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')


% load the water depth at the longitude of the mooring position (grey line in plot)
% computed with script prepare_bathy_sections.m
load data/dep_lon_EG.mat

% load ice concentration (mean July 2016 from Ifremer) (grey bar in plot)
load data/ice_EG_HG_N.mat


a = 40; % size of scatter plot

X_lat=X_lat-xres/2;
X_lon=X_lon-xres/2;



hf=figure;
set(hf, 'Position', [100 100 1400 1000])

%%% plot 1
%%% -----------------------------------------------------------------------
subplot(2,3,1)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[0.7 1 1 1]) % stretch its width and height
lat_min=77;
lat_max=80.5;

% on top of plot, add distance in km from mooring position
dist1=-lonlat2dx([0,0],[lat_min,y_moo_EG])/1000;
dist2=lonlat2dx([0,0],[lat_max,y_moo_EG])/1000;

% plot vertical distribution of particles
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
p1 = pcolor(ax1,X_lat,Y_lat,distr_EG_lat./num_days); shading flat;
caxis(ax1,[0 0.3])
colormap(ax1,bluewhite(50))

% plot sea ice bar (only if concentration is >15%)
hold on
ind=find(ice_lat_EG_all(:,3)>15);
plot(ice_lat_EG_all(ind,2),zeros(size(ice_lat_EG_all(ind,2))),'--gs',...
    'MarkerSize',6,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4])

% plot mooring position (yellow square)
hold on
plot(y_moo_EG,-dep_moo_EG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')

% plot water depth along the latitude
hold on
plot(dep_lat_EG(:,1),dep_lat_EG(:,2),'-','Linewidth',2,'color',[.6 .6 .6])

% set axis
set(ax1,'xlim',[lat_min lat_max])
set(ax1,'ylim',[-dep_moo_EG 0])
text(77.3,-200,'Station EG','fontsize',16)
text(77.3,-400,'52 m/d','fontsize',16)
set(gca,'fontsize',16)
% add distance in km on top of plot
ax3 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
set(ax3,'xlim',[dist1 dist2])
set(ax3,'ylim',[-dep_moo_EG 0])
set(ax3,'yticklabel','')
xlabel(ax1,'Latitude')
ylabel(ax1,'Depth (m)')
xlabel(ax3,'Distance (km)')
set(gca,'fontsize',16)


%%% plot 2
%%% -----------------------------------------------------------------------
subplot(2,3,2)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[0.87 1 1 1]) % stretch its width and height
lat_min=77;  
lat_max=80.5;    
dist1=-lonlat2dx([0,0],[lat_min,y_moo_HG])/1000;
dist2=lonlat2dx([0,0],[lat_max,y_moo_HG])/1000;

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
p1 = pcolor(ax1,X_lat,Y_lat,distr_HG_lat./num_days); shading flat;
caxis(ax1,[0 0.3])
colormap(bluewhite(50))

hold on
ind=find(ice_lat_HG_all(:,3)>15);
plot(ice_lat_HG_all(ind,2),zeros(size(ice_lat_HG_all(ind,2))),'--gs',...
    'MarkerSize',6,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4])
hold on
plot(y_moo_HG,-dep_moo_HG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(dep_lat_HG(:,1),dep_lat_HG(:,2),'-','Linewidth',2,'color',[.6 .6 .6])

set(ax1,'xlim',[lat_min lat_max])
set(ax1,'ylim',[-dep_moo_HG 0])
text(77.3,-200,'Station HG','fontsize',16)
text(77.3,-400,'29 m/d','fontsize',16)
set(gca,'fontsize',16)
ax3 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
set(ax3,'xlim',[dist1 dist2])
set(ax3,'ylim',[-dep_moo_HG 0])
set(ax3,'yticklabel','')
xlabel(ax1,'Latitude')
xlabel(ax3,'Distance (km)')
set(gca,'fontsize',16)


%%% plot 3
%%% -----------------------------------------------------------------------
subplot(2,3,3)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[0.9 1 1 1]) % stretch its width and height
lat_min=77;  
lat_max=80.5;    
dist1=-lonlat2dx([0,0],[lat_min,y_moo_N])/1000;
dist2=lonlat2dx([0,0],[lat_max,y_moo_N])/1000;

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
p1 = pcolor(ax1,X_lat,Y_lat,distr_N_lat./num_days); shading flat;
caxis(ax1,[0 0.3])
colormap(bluewhite(50))

% add colorbar
cb1=colorbar(ax1,'Position', [ax1_pos(1)+ax1_pos(3)+0.03  ax1_pos(2)-0.1  0.02 0.4]);
%'Position', [left bottom width height])
xlabel(cb1,'Rel. number of particles','fontsize',16)
cb1.Ticks = [0:0.05:1];
cb1.TickLabels = [0:0.05:1];
%cbarrow('up')

hold on
ind=find(ice_lat_N_all(:,3)>15);
plot(ice_lat_N_all(ind,2),zeros(size(ice_lat_N_all(ind,2))),'--gs',...
    'MarkerSize',6,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4])
hold on
plot(y_moo_N,-dep_moo_N,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(dep_lat_N(:,1),dep_lat_N(:,2),'-','Linewidth',2,'color',[.6 .6 .6])

set(ax1,'xlim',[lat_min lat_max])
set(ax1,'ylim',[-dep_moo_N 0])
text(77.3,-200,'Station N','fontsize',16)
text(77.3,-400,'52 m/d','fontsize',16)

set(gca,'fontsize',16)
ax3 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
set(ax3,'xlim',[dist1 dist2])
set(ax3,'ylim',[-dep_moo_N 0])
set(ax3,'yticklabel','')

xlabel(ax1,'Latitude')
xlabel(ax3,'Distance (km)')
set(gca,'fontsize',16)


%%% plot 4
%%% -----------------------------------------------------------------------
subplot(2,3,4)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[0.7 1 1 1]) % stretch its width and height
lon_min=-6;  
lon_max=6;    
dist1=-lonlat2dx([lon_min,x_moo_EG],[79,79])/1000;
dist2=lonlat2dx([lon_max,x_moo_EG],[79,79])/1000;

p1 = pcolor(X_lon,Y_lon,distr_EG_lon./num_days); shading flat;
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
set(gca,'fontsize',16)
colormap(ax1,bluewhite(50))
caxis(ax1,[0 0.3])

hold on
ind=find(ice_lon_EG_all(:,3)>15);
plot(ice_lon_EG_all(ind,1),zeros(size(ice_lon_EG_all(ind,2))),'--gs',...
    'MarkerSize',6,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4])
hold on
plot(x_moo_EG,-dep_moo_EG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(dep_lon_EG(:,1),dep_lon_EG(:,2),'-','Linewidth',2,'color',[.6 .6 .6])

set(ax1,'xlim',[lon_min lon_max])
set(ax1,'ylim',[-dep_moo_EG 0])
text(-5,-200,'Station EG','fontsize',16)
text(-5,-400,'52 m/d','fontsize',16)
set(gca,'fontsize',16)

ax3 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
set(ax3,'xlim',[dist1 dist2])
set(ax3,'ylim',[-dep_moo_EG 0])
set(ax3,'yticklabel','')
xlabel(ax1,'Longitude')
ylabel(ax1,'Depth (m)')
xlabel(ax3,'Distance (km)')
set(gca,'fontsize',16)

%%% plot 5
%%% -----------------------------------------------------------------------
subplot(2,3,5)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[0.87 1 1 1]) % stretch its width and height
lon_min=-2;  
lon_max=10;   
dist1=-lonlat2dx([lon_min,x_moo_HG],[79,79])/1000;
dist2=lonlat2dx([lon_max,x_moo_HG],[79,79])/1000;

p1 = pcolor(X_lon,Y_lon,distr_HG_lon./num_days); shading flat;
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
set(gca,'fontsize',16)
colormap(ax1,bluewhite(50))
caxis(ax1,[0 0.3])

hold on
ind=find(ice_lon_HG_all(:,3)>15);
plot(ice_lon_HG_all(ind,1),zeros(size(ice_lon_HG_all(ind,1))),'--gs',...
    'MarkerSize',6,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4])
hold on
plot(x_moo_HG,-dep_moo_HG,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(dep_lon_HG(:,1),dep_lon_HG(:,2),'-','Linewidth',2,'color',[.6 .6 .6])

set(ax1,'xlim',[lon_min lon_max])
set(ax1,'ylim',[-dep_moo_HG 0])
text(-1,-200,'Station HG','fontsize',16)
text(-1,-400,'29 m/d','fontsize',16)
set(gca,'fontsize',16)

ax3 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
set(ax3,'xlim',[dist1 dist2])
set(ax3,'ylim',[-dep_moo_HG 0])
set(ax3,'yticklabel','')
xlabel(ax1,'Longitude')
xlabel(ax3,'Distance (km)')
set(gca,'fontsize',16)

%%% plot 6
%%% -----------------------------------------------------------------------
subplot(2,3,6)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[0.9 1 1 1]) % stretch its width and height
lon_min=-2;  
lon_max=10;    
dist1=-lonlat2dx([lon_min,x_moo_N],[79,79])/1000;
dist2=lonlat2dx([lon_max,x_moo_N],[79,79])/1000;

p1 = pcolor(X_lon,Y_lon,distr_N_lon./num_days); shading flat;
ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
set(gca,'fontsize',16)
colormap(ax1,bluewhite(50))
caxis(ax1,[0 0.3])

hold on
ind=find(ice_lon_N_all(:,3)>15);
plot(ice_lon_N_all(ind,1),zeros(size(ice_lon_N_all(ind,2))),'--gs',...
    'MarkerSize',6,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4])
hold on
plot(x_moo_N,-dep_moo_N,'--gs',...
    'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
hold on
plot(dep_lon_N(:,1),dep_lon_N(:,2),'-','Linewidth',2,'color',[.6 .6 .6])

set(ax1,'xlim',[lon_min lon_max])
set(ax1,'ylim',[-dep_moo_N 0])
text(-1,-200,'Station N','fontsize',16)
text(-1,-400,'52 m/d','fontsize',16)
set(gca,'fontsize',16)

ax3 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
set(ax3,'xlim',[dist1 dist2])
set(ax3,'ylim',[-dep_moo_N 0])
set(ax3,'yticklabel','')
xlabel(ax1,'Longitude')
xlabel(ax3,'Distance (km)')
set(gca,'fontsize',16)

% saveas2('plots/vertical_distr_EG_N_HG_bottom_2016.png',400)
% saveas2('plots/vertical_distr_EG_N_HG_bottom_2016.pdf',400)

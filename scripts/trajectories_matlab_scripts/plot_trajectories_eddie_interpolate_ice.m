%%% this script computes the following for each of the 9 experiments:
% (1) lon at surface
% (2) lat at surface
% (3) day at surface
% (4) year at surface
% (5) day at surface
% (6) year at surface
% (7) ice concentration at surface
% (8) 1: ice>15%, 0: ice<15%
% (9) mean distance to mooring loc
% (10) trajectory path length

%%% 9 experiments:
% moorings EG, HG, N
% depth: bottom
% speed: 60 m/d, 20 m/d and 52 m/d (EG and N) and 29 m/d (HG)

%%% trajectories started on every day in 2016

%%% plot should only contain trajectories where particles are at the surface
% between 1 Mar and 31 Jul 2016
%%% datevec2doy(datevec('3/1/2016'))
%%% datevec2doy(datevec('7/31/2016'))
% st=60;
% en=212;

close all
clear all

%path(path,'/home/ollie/cwekerle/bin/m_map/')
datapath='/home/ollie/cwekerle/code_post/trajectories_fram_back_eddie2/';


prepare_data=0
plot_data=1


% load ice concentration / lon / lat  (daily data)
datapath_ice='/work/ollie/cwekerle/data/data_seaice_ifremer/';
lon_ice_aux=ncread([datapath_ice,'grid_north_12km.nc'],'longitude');
lon_ice=mod((lon_ice_aux+180),360)-180;
lat_ice=ncread([datapath_ice,'grid_north_12km.nc'],'latitude');
% make dataset smaller -> interpolation much faster...
lon_v=reshape(lon_ice,1,[]);
lat_v=reshape(lat_ice,1,[]);
iind=find(lat_v>72 & lat_v<86 & lon_v<25 & lon_v>-25);
lon_v_small=double(lon_v(iind));
lat_v_small=double(lat_v(iind));
llon=size(lon_ice,1)
llat=size(lat_ice,2)
len_ice=length(lon_v)


loc={'EG';'HG';'N'}
sp={'60';'20';'xx'}
dep={'bottom'}
% EG: speed xx=52 m/d
% HG: speed xx=29 m/d
% N:  speed xx=52 m/d

%%% mooring location
for ii=1:3
    locid=char(loc(ii,:))
    blon=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_bottom_speed60md.nc'],'blon');
    blat=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_bottom_speed60md.nc'],'blat');
    x_moo(ii)=blon(1,1);
    y_moo(ii)=blat(1,1);
end


if prepare_data
    
    
    st=1
    en=365
    num_days=en-st+1;
    
    for ii=1:3 %%% location
        for jj=1:3  %%% speed
            locid=char(loc(ii,:))
            spid=char(sp(jj,:))
            if (spid=='xx' & locid=='EG')
                spid='52'
            elseif (spid=='xx' & locid=='HG')
                spid='29'
            elseif (spid=='xx' & locid=='N')
                spid='52'
            end
            %             if (spid=='xx' & locid=='EG')
            %                 spid='47'
            %             elseif (spid=='xx' & locid=='HG')
            %                 spid='30'
            %             elseif (spid=='xx' & locid=='N')
            %                 spid='77'
            %             end
            
            depid=char(dep(1,:))
            
            blon=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blon');
            blat=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'blat');
            bday=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bday');
            bdepth=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'bdepth');
            byear=ncread([datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md.nc'],'byear');
            
            
            
            % surface position / day / year
            ind=max(find(byear(:,1)>2000))-1;
            x_surf=blon(ind,:);
            y_surf=blat(ind,:);
            day_surf=double(bday(ind,:));
            year_surf=byear(ind,:);
            day_bot=double(bday(1,:));
            year_bot=byear(1,:);
            
            traj_dist=zeros(365,1);
            traj_len=zeros(365,1);
            ice_conc=zeros(365,1);
            ice_ind=zeros(365,1);
            
            for kk=1:365
                kk
                % distance from mooring position to particle position at
                % surface
                traj_dist(kk)=deg2km(distance(y_moo(ii),x_moo(ii),y_surf(kk),x_surf(kk)));
                
                % trajectory path length
                traj_len(kk)=sum(deg2km(distance(blat(1:ind-1,kk)',blon(1:ind-1,kk)',blat(2:ind,kk)',blon(2:ind,kk)')));
                
                
                % read ice data
                ice=ncread([datapath_ice,'daily/',num2str(year_surf(kk)),datestr(day_surf(kk),'mmdd'),'.nc'],'concentration');
                ice_v=reshape(ice,1,[]);
                ice_v_small=double(ice_v(iind));
                % interpolate ice data
                ice_conc(kk) = griddata(lon_v_small,lat_v_small,ice_v_small,x_surf(kk),y_surf(kk));
                
                if ice_conc(kk)>0.15
                    ice_ind(kk)=1;
                end
                
                
            end
            
            
            % what to save?
            % (1) lon at surface
            % (2) lat at surface
            % (3) day at surface
            % (4) year at surface
            % (5) day at surface
            % (6) year at surface
            % (7) ice concentration at surface
            % (8) 1: ice>15%, 0: ice<15%
            % (9) mean distance to mooring loc
            % (10) trajectory path length
            filename=[datapath,'result/drifter_Arc40_start2016_',locid,'_',depid,'_speed',spid,'md_ice_dist_len.txt'];
            fid=fopen(filename,'wt');
            for n=1:365
                fprintf(fid,'%f %f %d %d %d %d %f %d %f %f \n', ...
                    x_surf(n),y_surf(n),day_surf(n),year_surf(n),day_bot(n),year_bot(n),...
                    ice_conc(n),ice_ind(n),traj_dist(n),traj_len(n));
            end
            fclose(fid);
            
            
        end
        
    end
    
end

if plot_data
    
    
    
    %%% time period: Mar - Jul 2016 at the surface
    st_mar1=60; % day of year march 1
    en_jul31=212;  % day of year july 31
    %%% then for each location, we plot only the days when the
    %%% particle is at the surface between days 60 and 212
    
    %     A=load([datapath,'result/drifter_Arc22_start2009_EG_bottom_speed60md_ice_dist_len.txt']);
    %     B=load([datapath,'result/drifter_Arc22_start2009_HG_bottom_speed60md_ice_dist_len.txt']);
    %     C=load([datapath,'result/drifter_Arc22_start2009_N_bottom_speed60md_ice_dist_len.txt']);
    
    % plot: sinking speeds from Morten
    A=load([datapath,'result/drifter_Arc40_start2016_EG_bottom_speed52md_ice_dist_len.txt']);
    B=load([datapath,'result/drifter_Arc40_start2016_HG_bottom_speed29md_ice_dist_len.txt']);
    C=load([datapath,'result/drifter_Arc40_start2016_N_bottom_speed52md_ice_dist_len.txt']);
    
    
    % load menthly mean ice concentration
    ice_mar = ncread([datapath_ice,'20160301-20160331.nc'],'concentration');
    ice_apr = ncread([datapath_ice,'20160401-20160430.nc'],'concentration');
    ice_may = ncread([datapath_ice,'20160501-20160531.nc'],'concentration');
    ice_jun = ncread([datapath_ice,'20160601-20160630.nc'],'concentration');
    ice_jul = ncread([datapath_ice,'20160701-20160731.nc'],'concentration');
    
    
    % ice concentration
%     m_proj('lambert','long',[-8 15],'lat',[75 81]);
%     hf=figure(1);
%     set(hf, 'Position', [96 455 800 700])
%     [X,Y]=m_ll2xy(lon_ice,lat_ice);
%     pcolor(X,Y,ice_mar);
%     colorbar
%     caxis([0 100])
%     shading flat
%     hold on
%     contour(X,Y,ice_mar,[15 999],'linewidth',2,'color',[.5 .5 .5]);
%     m_grid('tickdir','out','xtick',[-20:5:20],'ytick',[50:2:82],...
%         'linest','none','fontsize',12);
%     m_gshhs_i('patch',[.6 .6 .6],'edgecolor',[.6 .6 .6]);
%     m_etopo2('contour',[-4000:1000:-1000],'color',[.5 .5 .5],'linewidth',1);
%     m_coast;
    
    % make plot
    m_proj('lambert','long',[-8 15],'lat',[75 81]);
    hf=figure(2);
    set(hf, 'Position', [96 455 800 700])
    
    
    % plot ice edge
    [X,Y]=m_ll2xy(lon_ice,lat_ice);
    contour(X,Y,ice_mar,[15 999],'linewidth',2,'color',[255 255 204]./255);
    hold on
    contour(X,Y,ice_apr,[15 999],'linewidth',2,'color',[255 255 153]./255);
    hold on
    contour(X,Y,ice_may,[15 999],'linewidth',2,'color',[255 255 102]./255);
    hold on
    contour(X,Y,ice_jun,[15 999],'linewidth',2,'color',[255 255 0]./255);
    hold on
    contour(X,Y,ice_jul,[15 999],'linewidth',2,'color',[204 204 0]./255);
    c1=plot(nan,'Color',[255 255 204]./255,'linewidth',2);
    c2=plot(nan,'Color',[255 255 153]./255,'linewidth',2);
    c3=plot(nan,'Color',[255 255 102]./255,'linewidth',2);
    c4=plot(nan,'Color',[255 255 0]./255,'linewidth',2);
    c5=plot(nan,'Color',[204 204 0]./255,'linewidth',2);
    
    % plot surface pos of trajectories (different color for ice / no ice)
    sz=15;
    % EG 60 m/d
    
    ind_time_A=find( A(:,3)>=st_mar1 & A(:,3)<=en_jul31);
    [x,y]=m_ll2xy(A(ind_time_A,1),A(ind_time_A,2));
    ind= find(A(ind_time_A,8)==0); % no ice (dark color)
    s1=scatter(x(ind),y(ind),sz,'MarkerEdgeColor',[51 51 255]./255,...
        'MarkerFaceColor',[51 51 255]./255,...
        'LineWidth',1.5);
    hold on
    ind= find(A(ind_time_A,8)==1); % with ice (light color)
    s2=scatter(x(ind),y(ind),sz,'MarkerEdgeColor',[204 204 255]./255,...
        'MarkerFaceColor',[204 204 255]./255,...
        'LineWidth',1.5);
    hold on
    
    % HG 60 m/day
    ind_time_B=find( B(:,3)>=st_mar1 & B(:,3)<=en_jul31);
    [x,y]=m_ll2xy(B(ind_time_B,1),B(ind_time_B,2));
    ind= find(B(ind_time_B,8)==0); % no ice (dark color)
    s3=scatter(x(ind),y(ind),sz,'MarkerEdgeColor',[255 51 51]./255,...
        'MarkerFaceColor',[255 51 51]./255,...
        'LineWidth',1.5);
    hold on
    ind= find(B(ind_time_B,8)==1); % with ice (light color)
    s4=scatter(x(ind),y(ind),sz,'MarkerEdgeColor',[255 204 204]./255,...
        'MarkerFaceColor',[255 204 204]./255,...
        'LineWidth',1.5);
    hold on
    
    % N
    ind_time_C=find( C(:,3)>=st_mar1 & C(:,3)<=en_jul31);
    [x,y]=m_ll2xy(C(ind_time_C,1),C(ind_time_C,2));
    ind= find(C(ind_time_C,8)==0); % no ice (dark color)
    s5=scatter(x(ind),y(ind),sz,'MarkerEdgeColor',[.4 .4 .4],...
        'MarkerFaceColor',[.4 .4 .4],...
        'LineWidth',1.5);
    hold on
    ind= find(C(ind_time_C,8)==1); % with ice (light color)
    s6=scatter(x(ind),y(ind),sz,'MarkerEdgeColor',[.8 .8 .8],...
        'MarkerFaceColor',[.8 .8 .8],...
        'LineWidth',1.5);
    
    
    
    m_grid('tickdir','out','xtick',[-20:5:20],'ytick',[50:2:82],...
        'linest','none','fontsize',12);
    m_gshhs_i('patch','k','edgecolor','k');
    m_etopo2('contour',[-4000:1000:-1000],'color',[.5 .5 .5],'linewidth',1);
    
    
    
    % plot mooring locations
    sz=120;
    hold on
    [x,y]=m_ll2xy(x_moo(1),y_moo(1));
    ss1=scatter(x,y,sz,'MarkerEdgeColor','k',...
        'MarkerFaceColor','b',...
        'LineWidth',1.0)   ;
    hold on
    [x,y]=m_ll2xy(x_moo(2),y_moo(2));
    ss2=scatter(x,y,sz,'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'LineWidth',1.0)    ;
    hold on
    [x,y]=m_ll2xy(x_moo(3),y_moo(3));
    ss3=scatter(x,y,sz,'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'LineWidth',1.0)   ;
    
    
    leg=legend([c1,c2,c3,c4,c5,s1,s2,s3,s4,s5,s6],...
        'ice edge Mar','ice edge Apr','ice edge May','ice edge Jun','ice edge Jul',...
        'EG ice free','EG ice covered',...
        'HG ice free','HG ice covered',...
        'N ice free','N ice covered');
    set(leg,'AutoUpdate','off');
    
    
    % (first element = top, last element = bottom),
    % chi=get(gca, 'Children')
    % chi_new=[chi(4:end);chi(1:3)]
    % %set(gca, 'Children',chi_new)
    % set(gca,'Children',[chi(1:3);chi(4:end)])
    
    set(gca,'fontsize',16)
    
    % % saveas2('plots/surface_positions_EG_HG_N_2016.png',300)
    
    %set(gcf,'Renderer','OpenGL')
    %  saveas2('plots/test.png',300)
    
    
            % (8) 1: ice>15%, 0: ice<15%
            % (9) mean distance to mooring loc
            % (10) trajectory path length
    

    disp('EG percentage of particles in sea ice:' )
    ind1=find(A(ind_time_A,8)==1);
    length(ind1)/length(ind_time_A)*100
    
    disp('HG percentage of particles in sea ice:' )
    ind1=find(B(ind_time_B,8)==1);
    length(ind1)/length(ind_time_B)*100
    
    disp('N percentage of particles in sea ice:' )
    ind1=find(C(ind_time_C,8)==1);
    length(ind1)/length(ind_time_C)*100
    
    
    disp('EG mean / std mean distance to mooring:' )
    round(mean(A(ind_time_A,9)))
    round(std(A(ind_time_A,9)))
    
    disp('HG mean / std mean distance to mooring:' )
    round(mean(B(ind_time_B,9)))
    round(std(B(ind_time_B,9)))
    
    disp('N mean / std mean distance to mooring:' )
    round(mean(C(ind_time_C,9)))
    round(std(C(ind_time_C,9)))
    
    disp('EG mean / std mean travel distance:' )
    round(mean(A(ind_time_A,10)))
    round(std(A(ind_time_A,10)))
    
    disp('HG mean / std mean travel distance:' )
    round(mean(B(ind_time_B,10)))
    round(std(B(ind_time_B,10)))
    
    disp('N mean / std mean travel distance:' )
    round(mean(C(ind_time_C,10)))
    round(std(C(ind_time_C,10)))
        
    
    
end




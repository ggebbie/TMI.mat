%% tmi_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
% load the master file: choose either woce2deg or woce4deg 
% (2 deg == 2 x 2 horiz resolution, 33 vertical levels)
% (4 deg == 4 x 4 horiz resolution, 33 vertical levels)
% (woce == trained on the WOCE Global Hydrographic Climatology,
%          Gouretski and Koltermann 2004)

load A_woce4deg  %% choose between the 4 deg and 2 deg MAT-files

addpath ~/mcode/toolbox/seawater
addpath ~/mcode/toolbox/general
addpath ~/mcode/toolbox/tmi

load ~/matlab/box_method/c13/work_ref_11feb2012.mat A_ref x_ref cobs ...
    atbotside grd
load ~/matlab/box_method/woce_4deg/stats_4deg_33lev_woce_4jan.mat inmixlyr

Tobs = x_ref.Tmod;
Sobs = x_ref.Smod;
O18obs = x_ref.O18mod;
Pobs = x_ref.Pmod;
Nobs = x_ref.Nmod;
Oobs = x_ref.Omod;
POobs = 170.*Pobs + Oobs;
NOobs = (170/15.5).*Nobs +Oobs;
NOjobs = 10.625.*Nobs +Oobs; % johnson definition
A = A_ref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the data a la Gordon 2005
%  does my model fit the data? 

%% do 115-130 E.
%% do  2N -  9 S.
it = i; jt = j; kt = k;
itf = find( LON(it) > 115 & LON(it) < 130 & LAT(jt)<2 & LAT(jt)>-10 ...
           & ~atbotside);

figure % looks similar to Gordon 2005 (Fig. 4)
plot(Tobs(itf),-DEPTH(kt(itf)),'o')
hold on
plot(cobs.Tmod(itf),-DEPTH(kt(itf)),'r.')

figure % looks similar to Gordon 2005 (Fig. 4)
plot(Sobs(itf),-DEPTH(kt(itf)),'o')
hold on
plot(cobs.Smod(itf),-DEPTH(kt(itf)),'r.')

figure % looks similar to Gordon 2005 (Fig. 4)
plot(Sobs(itf),Tobs(itf),'o')
hold on
plot(cobs.Smod(itf),cobs.Tmod(itf),'r.')
axis([34 35.4  3  32])
 

%% get temperature and salinity on sigma0=25.5
sigma0orig = sw_dens(Sobs,Tobs,0)-1000;
sigma0 = sigma0orig;
sigma0(inmixlyr) = nan;
%sigma0(atbotside) = nan;

%% find the depth of the sigma0 isopycnal.
sigma0fld = vector_to_field(sigma0,it,jt,kt);
%Error using griddedInterpolant
%The point coordinates are not sequenced in strict monotonic order.
%Error in interp1>Interp1D (line 335)
%    F = griddedInterpolant(X,V(:,1),method);
%Error in interp1 (line 220)
%    Vq = Interp1D(X,V,Xq,method); 

sigdepth = nan(45,90);
sigk     = nan(45,90);
for xx = 1:90
  for yy = 1:45
     jjj = find(it==xx & jt==yy & ~inmixlyr);
     %tmp = sq(sigma0fld(:,yy,xx));
     tmp = sq(sigma0(jjj));
     tmp(diff(tmp)==0) = nan;
     dtmp = DEPTH(kt(jjj));;

     ktmp = kt(jjj);     
     dtmp(isnan(tmp)) = [];
     ktmp(isnan(tmp)) = [];
     tmp(isnan(tmp)) = [];
     if  length(tmp)>1
       sigdepth(yy,xx) = interp1(tmp,dtmp,25.5);
       %% get the K-level.
       sigk(yy,xx) = interp1(tmp,ktmp,25.5);
     end
   end
end
 
sigma0 = sw_dens(Sobs,Tobs,0)-1000;

%%using sigk, map 3d field onto 2d. 
Ssig = vector_to_field_2dinterp(Sobs,it,jt,kt,sigk);
Tsig = vector_to_field_2dinterp(Tobs,it,jt,kt,sigk);
npacsig = vector_to_field_2dinterp(sq(c(:,4)),it,jt,kt, ...
                                   sigk);
spacsig = vector_to_field_2dinterp(sq(c(:,6)),it,jt,kt,sigk);



figure
%pcolor(LON,LAT,sq(Ssig))
contourf(LON,LAT,sq(Ssig),34.4:0.1:36.5)
%image(LON-2,LAT-2,sq(Ssig))
hold on
plot(coast(:,1)+360,coast(:,2),'LineWidth',1,'Color','k')
plot(coast(:,1),coast(:,2),'LineWidth',1,'Color','k')
axis equal
axis([35  210  -30  30])
title('WOCE salinity: \sigma_0=25.5')
caxis([34.4 36.5])
colorbar('horiz')

figure
%pcolor(LON,LAT,sq(Tsig))
contourf(LON,LAT,sq(Tsig),15:0.5:23)
hold on
plot(coast(:,1)+360,coast(:,2),'LineWidth',1,'Color','k')
plot(coast(:,1),coast(:,2),'LineWidth',1,'Color','k')
axis equal
axis([35  210  -30  30])
title('WOCE \theta: \sigma_0=25.5')
caxis([15 23])
colorbar('horiz')

figure
contourf(LON,LAT,sq(sigdepth),0:25:500)
%pcolor(LON-2,LAT,sq(sigdepth))
hold on
plot(coast(:,1)+360,coast(:,2),'LineWidth',1,'Color','k')
plot(coast(:,1),coast(:,2),'LineWidth',1,'Color','k')
axis equal
axis([35  210  -30  30])
title('WOCE depth: \sigma_0=25.5')
%caxis([17 23])
colorbar('horiz')

figure
contourf(LON,LAT,100.*sq(npacsig),[0:2:20 30:10:100])
%pcolor(LON-2,LAT,sq(100.*npacsig))
hold on
plot(coast(:,1)+360,coast(:,2),'LineWidth',1,'Color','k')
plot(coast(:,1),coast(:,2),'LineWidth',1,'Color','k')
axis equal
axis([35  250  -30  60])
caxis([0 50])
title('TMI NPAC % by volume: \sigma_0=25.5')
%caxis([17 23])
colorbar('horiz')

figure
contourf(LON,LAT,100.*sq(spacsig),[0:2:20 30:10:100])
%pcolor(LON-2,LAT,sq(100.*spacsig))
hold on
plot(coast(:,1)+360,coast(:,2),'LineWidth',1,'Color','k')
plot(coast(:,1),coast(:,2),'LineWidth',1,'Color','k')
axis equal
axis([35  210  -30  30])
caxis([0 100])
title('TMI SPAC % by volume: \sigma_0=25.5')
%caxis([17 23])
colorbar('horiz')

% variables in A.mat %%%
%
% A: pathways matrix
% i,j,k: vectors that contain the spatial coordinates 
%        consistent with the A matrix.
% LAT,LON,DEPTH: the latitude, longitude and depth of the 3-d grid. 
% dP = amount of remineralized phosphate in each box.

%% Most examples are computationally more efficient with the output of a LU decomposition.
% Use an LU decomposition so that the inverse does not have
% to be stored explicitly.
 [L, U, P, Q, R] = lu (A) ;

 %% HERE: TRACK MAJOR WATER-MASSES TO ITF
 %% REGION. (AABW,AAIW,NPIW,NATL,N.TROP,S.TROP)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1: Track the pathways of a user-defined water mass.         %
% Steps: (a) define the water mass 1). by a rectangular surface patch %
%            dyed with passive tracer concentration of 1,             % 
%            or 2. load a pre-defined surface patch in d_all.mat.     %
%        (b) propagate the dye with the matrix A, with the result     % 
%            being the fraction of water originating from the         %
%            surface region.                                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % timor = find(LAT(jt)==-9 & LON(it)==123);

 %iii = find(kt > 24);

 % 1) GLOBAL, 2) ANT, 3 NATL, 4 SUBANT, 5 NPAC, 
 % 6) ARC, 7 MED, 8) TROP, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
 % 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
 % 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP

load d_all_4deg

jt = j;
it = i;
kt = k;
dntrop = d_all(:,8).*( LAT(jt)> 0 );
dstrop = d_all(:,8).*( LAT(jt)< 0 );
%dntrop = d_all(:,8).*(LAT(jt)==0)./2; % split the equator.
%dstrop = d_all(:,8).*(LAT(jt)==0)./2; % split the equator.

%atl = d_all(:,3) + d_all(:,6) + d_all(:,7) + d_all(:,17);

%% water-masses defined in this case. 
%% 1. tropindian. 2. NH-pactrop+npac+arc+med+natl, 3. SH-pactrop,
%% 4. pacific gate (under subant boundary). 5. indian gate (under subant boundary)

%% Fern has asked to combine the Pacific and Indian Gates.
%% first define Indian ocean. 

itf = (LON(it)>=95 & LON(it) < 125 & LAT(jt) >= -10 & LAT(jt) < ...
       30 & kt ==1);

itf = itf + (LON(it)>=95 & LON(it) <140 & LAT(jt) <= -4 & LAT(jt) >-14 ...
                 & itf ==0 & kt==1);

%satltrop = 0.5.*LAT(jt)==0 + LAT(jt)<0).* d_all(:,17);
indian = d_all(:,19).*(itf==0 & LAT(jt) > -28);% + satltrop ; 

%indian = ( LON(it)>20 & LON(it)< 95 & LAT(jt) < 30 & kt == 1 & itf==0);
%indian = indian + ( LAT(jt) < -5 & LON(it) >=95 & LON(it) <=148 & ...
%                    kt ==1 & itf==0); 

%nh = (LAT(jt)>0 & kt ==1 & indian == 0 & itf == 0);
%nh = nh + 0.5.*(LAT(jt)==0 & kt==1 & indian == 0 & itf==0);


pacgate = (LAT(jt) ==-28 & LON(it) > 130 & LON(it) <= 300);
indgate = (LAT(jt) ==-28 & pacgate==0);

%% npac and trop npac
%nh = (LAT(jt)>0 & (d_all(:,5)==1 | d_all(:,18)==1) & indian == 0 & ...
%      itf == 0);
nh = (LAT(jt)>0 & indian == 0 & itf == 0 & kt ==1);
nh = nh + 0.5.*(LAT(jt)==0 & indian == 0 & itf == 0 & kt ==1);

% trop pac southern hemisphere.
sh =  (LAT(jt)<0 &  indian==0 & itf==0 & pacgate==0 & ...
       indgate == 0 & kt ==1 ); 
sh = sh + 0.5.*(LAT(jt)==0 & kt==1 & indian == 0 & itf == 0 & pacgate ...
                == 0 & indgate ==0 & kt ==1 );
%                LON(it) >=95  & LON(it) < 280);

%% IT SEEMS LIKE THERE ARE WAYS THROUGH THIS WALL.
%% I WILL DYE THE ENTIRE SUBANT BOX INSTEAD.
%% pacgate = max latitude of any subantpac d_all(:,15).
% plus make it a wall.
% $$$ for xx = 1:90
% $$$   jpacgate(xx) = max( jt .* d_all(:,15).*(it==xx));
% $$$ end
% $$$ for xx = 1:90
% $$$   jindgate(xx) = max( jt .* d_all(:,16).*(it==xx));
% $$$ end
% $$$ ipacgate = 1:90;
% $$$ ipacgate(jpacgate==0) = [];
% $$$ jpacgate(jpacgate==0) = [];
% $$$ 
% $$$ iindgate = 1:90;
% $$$ iindgate(jindgate==0) = [];
% $$$ jindgate(jindgate==0) = [];
% $$$ 
% $$$ 
% $$$ % at the end of each gate go north until hit land.
% $$$ jtmp = jt(it ==ipacgate(1) & jt> jpacgate(1) & kt==1)
% $$$ tmp = find(diff(jtmp)>1);
% $$$ jtmp(tmp(1)+1:end) = []; % must be a continuous wall.
% $$$ ipacgate = [ ipacgate(1).*ones(length(jtmp),1)' ipacgate];
% $$$ jpacgate = [ jtmp' jpacgate];
% $$$ 
% $$$ jtmp = jt(it ==ipacgate(end) & jt> jpacgate(end) & kt==1)
% $$$ tmp = find(diff(jtmp)>1);
% $$$ jtmp(tmp(1)+1:end) = []; % must be a continuous wall.
% $$$ ipacgate = [ ipacgate(end).*ones(length(jtmp),1)' ipacgate];
% $$$ jpacgate = [ jtmp' jpacgate];
% $$$ 
% $$$ 
% $$$ % at the end of each gate go north until hit land.
% $$$ jtmp = jt(it ==iindgate(1) & jt> jindgate(1) & kt==1)
% $$$ tmp = find(diff(jtmp)>1);
% $$$ jtmp(tmp(1)+1:end) = []; % must be a continuous wall.
% $$$ iindgate = [ iindgate(1).*ones(length(jtmp),1)' iindgate];
% $$$ jindgate = [ jtmp' jindgate];
% $$$ 
% $$$ jtmp = jt(it ==iindgate(end) & jt> jindgate(end) & kt==1)
% $$$ tmp = find(diff(jtmp)>1);
% $$$ jtmp(tmp(1)+1:end) = []; % must be a continuous wall.
% $$$ iindgate = [ iindgate(end).*ones(length(jtmp),1)' iindgate];
% $$$ jindgate = [ jtmp' jindgate];
% $$$ 
% $$$ %pacgate = ( LAT(jt)==-36 & LON(it)>=135 & LON(it) < 290 & kt > 1);
% $$$ %indgate = (LAT(jt)==-36 & pacgate==0 & kt > 1);
% $$$ 
% $$$ % $$$ pacgate = zeros(Nfield,1);
% $$$ % $$$ for np = 1:length(ipacgate)
% $$$ % $$$   pacgate( it==ipacgate(np) & jt == jpacgate(np)) = 1;
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ indgate = zeros(Nfield,1);
% $$$ % $$$ for np = 1:length(iindgate)
% $$$ % $$$   indgate( it==iindgate(np) & jt == jindgate(np)) = 1;
% $$$ % $$$ end
% $$$ 
% $$$ indgate = zeros(Nfield,1);
% $$$ for nv = 1:Nfield
% $$$   nv
% $$$   indgate(nv) = sum(iindgate==it(nv) & jindgate == jt(nv) & ~inmixlyr(nv));
% $$$ end
% $$$ 
% $$$ pacgate = zeros(Nfield,1);
% $$$ for nv = 1:Nfield
% $$$   nv
% $$$   pacgate(nv) = sum(ipacgate==it(nv) & jpacgate == jt(nv) & ~inmixlyr(nv));
% $$$ end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANOTHER MEHTOD === DYE THE ENTIRE SUBANT BOX FROM TOP TO BOTTOM
% $$$ indgate = zeros(Nfield,1);
% $$$ isubantind = find(d_all(:,16)==1);
% $$$ 
% $$$ for nv = 1:Nfield
% $$$  nv
% $$$  indgate(nv) = sum( it(isubantind) == it(nv) & jt(isubantind)==jt(nv));
% $$$ end
% $$$ 
% $$$ pacgate = zeros(Nfield,1);
% $$$ isubantpac = find(d_all(:,15)==1);
% $$$ for nv = 1:Nfield
% $$$  nv
% $$$  pacgate(nv) = sum( it(isubantpac) == it(nv) & jt(isubantpac)==jt(nv));
% $$$ end

dITF(:,1) = itf;
dITF(:,2) = indian;
dITF(:,3) = nh;
dITF(:,4) = sh;
dITF(:,5) = pacgate;
dITF(:,6) = indgate;
%dITF(:,7) = atl; % could include with indian gate.

Agate = A;

ipacgate = find(pacgate);
iindgate= find(indgate);

for nn = 1:length(pacgate)
   nn
   Agate(ipacgate(nn),:) = 0;
   Agate(ipacgate(nn),ipacgate(nn)) = 1;
end
for nn = 1:length(iindgate)
   nn
   Agate(iindgate(nn),:) = 0;
   Agate(iindgate(nn),iindgate(nn)) = 1;
end

%dITF = [d_all(:,2:5) dntrop dstrop];

 % to get quantity of dyed water throughout ocean:
 [L, U, P, Q, R] = lu (Agate) ;
 c = Q * (U \ (L \ (P * (R \ dITF))));

 iii = find( it == 32 & jt==22 ); % ITF
 figure
 plot(100.*c(iii,:),-DEPTH(kt(iii)))
 hold on
 legend('ITF','IND','NPAC','SPAC','DEEPPAC','DEEPIND','ATL')
  % legend('ANT','NATL','SUBANT','NPAC','NTROP','STROP')
 ylabel('DEPTH [m]')
 xlabel('% volume')  
 grid
 title('126^{\circ}E, 4^{\circ}S')

 %% Fern: combine deeppac and deepind.
 c(:,5) = c(:,5)+c(:,6);
 c(:,6) = [];

 set(0,'DefaultLineLineWidth',2.0)
 set(0,'DefaultTextFontSize',15)
 set(0,'DefaultAxesLineWidth',1.6)
 set(0,'DefaultAxesFontSize',15)
 
 figure
 subplot('position',[.2 .1 .45 .8])
 plot(100.*c(iii,:),DEPTH(kt(iii)))
 set(gca,'ydir','reverse')
 hold on
 legend('ITF','IND','NPAC','SPAC','ANT','FontSize',12)
 ylabel('DEPTH [m]')
 xlabel('% volume')  
 grid
 title('126^{\circ}E, 4^{\circ}S')
 axis([0 100 0 2000])
 set(gca,'XTick',0:20:100)
 set(gca,'YTick',0:500:2000)
 
%%  Same thing as above but for ITF outflow.
%   outflow

%% 5 Apr 2013: Move outflow location
outflow = find(LON(it) == 122 & LAT(jt)==-12 ); 
%outflow = find(LON(it) == 114 & LAT(jt)==-12 ); 
 figure
  subplot('position',[.2 .1 .45 .8])
 plot(100.*c(outflow,:),DEPTH(kt(outflow)))
 hold on
 ylabel('DEPTH [m]')
 xlabel('% volume')  
 set(gca,'ydir','reverse')
 legend('ITF','IND','NPAC','SPAC','ANT','FontSize',12)
 grid
 title('122^{\circ}E, 12^{\circ}S')
 axis([0 100 0 2000])
 set(gca,'XTick',0:20:100)
 set(gca,'YTick',0:500:2000)

 
%   inflow
%inflow = find(LON(it) == 134 & LAT(jt)==4 ); 
inflow = find(LON(it) == 126 & LAT(jt)==4 ); 
figure
subplot('position',[.2 .1 .45 .8])
plot(100.*c(inflow,:),DEPTH(kt(inflow)))
hold on
ylabel('DEPTH [m]')
xlabel('% volume')  
set(gca,'ydir','reverse')
grid
title('126^{\circ}E, 4^{\circ}N')
 legend('ITF','IND','NPAC','SPAC','ANT','FontSize',12)
 axis([0 100 0 2000])
 set(gca,'XTick',0:20:100)
 set(gca,'YTick',0:500:2000)

 %% make map of surface origins for Fern 
LONplot = [LON(11:end-1); LON(1:10)+360];
wm = sq(c(isfc,1) + 2.*c(isfc,2) + 3.*c(isfc,3) + ...
     4.*c(isfc,4) + 5.*c(isfc,5)); 
wm( LAT(jt(isfc))<-28) = 5;
wmfld = vector_to_field_2d(wm,it(isfc),jt(isfc));
pvol_plot = sq(wmfld(:,[lonafr:end-1 1:lonafr-1],1));
figure
m_proj('miller','longitudes',[40 200],'latitudes',[-40 40])
m_grid('box','fancy','tickdir','in','ytick',[-60 -40 0 40 60])%,'YTickLabel',80:-20:-80);
hold on
m_pcolor(LONplot,LAT,pvol_plot)
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');

gtext('IND','Color',[1 1 1],'Fontweight','bold')
gtext('ITF','Color',[1 1 1],'Fontweight','bold')
gtext('NPAC','Color',[0 0 0],'Fontweight','bold')
gtext('SPAC','Color',[0 0 0],'Fontweight','bold')
gtext('ANT','Color',[1 1 1],'Fontweight','bold')


[xx,yy]=m_ll2xy(126,-4)
plot(xx,yy,'wo','MarkerSize',7,'MarkerFaceColor',[1 1 1])

[xx,yy]=m_ll2xy(114,-12)
plot(xx,yy,'wo','MarkerSize',7,'MarkerFaceColor',[1 1 1])

[xx,yy]=m_ll2xy(134,4)
plot(xx,yy,'ok','MarkerSize',8,'MarkerFaceColor',[1  1 1])

%m_contourf(LONplot,LAT,sq(wm))%,[.2:.2:1.6])
%m_contourf(LONplot,LAT,sq(pDPb(19,:,:)),[0:.05:1.6])
%clabel(c,h)

 
 figure
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultTextFontSize',6)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',6)

subplot(221)
pvol_plot = 100.*sq(Cmax(:,[lonafr:end-1 1:lonafr-1],1));
hold on
m_proj('miller','longitudes',[40 400],'latitudes',[-80 80])
[c,h]=m_contourf(LONplot,LAT,pvol_plot,0:5:100)
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');
hotreverse
caxis([0 100])
clabel(c,h,'manual')
title('ANT')
set(gca,'XTick',60:30:390)
set(gca,'XTickLabel',[60:30:180 -150:30:0])
set(gca,'YTick',-80:20:80)
set(gca,'YTickLabel',80:-20:-80)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% track the water mass core (maximum values).

for nn = 1:size(c,2)
  C(:,:,:,nn) = vector_to_field(sq(c(:,nn)),i,j,k);  % g is the 3-d field of water-mass concentration.
end

Cmax = sq(max(C,[],1));
contourf(LON,LAT,sq(Cmax(:,:,1)),0:0.05:1) 

%% make a figure of core values.
lonafr = find(LON<40);
lonafr = lonafr(end);
LONplot = [LON(lonafr:end-1); LON(1:lonafr-1)+360];

figure
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultTextFontSize',6)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',6)

subplot(221)
pvol_plot = 100.*sq(Cmax(:,[lonafr:end-1 1:lonafr-1],1));
hold on
m_proj('miller','longitudes',[40 400],'latitudes',[-80 80])
[c,h]=m_contourf(LONplot,LAT,pvol_plot,0:5:100)
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');
hotreverse
caxis([0 100])
clabel(c,h,'manual')
title('ANT')
set(gca,'XTick',60:30:390)
set(gca,'XTickLabel',[60:30:180 -150:30:0])
set(gca,'YTick',-80:20:80)
set(gca,'YTickLabel',80:-20:-80)

subplot(222)
pvol_plot = 100.*sq(Cmax(:,[lonafr:end-1 1:lonafr-1],2));
hold on
m_proj('miller','longitudes',[40 400],'latitudes',[-80 80])
[c,h]=m_contourf(LONplot,LAT,pvol_plot,0:5:100)
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');
hotreverse
caxis([0 100])
clabel(c,h,'manual')
title('NATL')
set(gca,'XTick',60:30:390)
set(gca,'XTickLabel',[60:30:180 -150:30:0])
set(gca,'YTick',-80:20:80)
set(gca,'YTickLabel',80:-20:-80)

subplot(223)
pvol_plot = 100.*sq(Cmax(:,[lonafr:end-1 1:lonafr-1],3));
hold on
m_proj('miller','longitudes',[40 400],'latitudes',[-80 80])
[c,h]=m_contourf(LONplot,LAT,pvol_plot,0:5:100)
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');
hotreverse
caxis([0 100])
clabel(c,h,'manual')
title('SUBANT')
set(gca,'XTick',60:30:390)
set(gca,'XTickLabel',[60:30:180 -150:30:0])
set(gca,'YTick',-80:20:80)
set(gca,'YTickLabel',80:-20:-80)

subplot(224)
pvol_plot = 100.*sq(Cmax(:,[lonafr:end-1 1:lonafr-1],4));
hold on
m_proj('miller','longitudes',[40 400],'latitudes',[-80 80])
[c,h]=m_contourf(LONplot,LAT,pvol_plot,0:5:100)
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');
hotreverse
caxis([0 100])
clabel(c,h,'manual')
title('NPAC')
set(gca,'XTick',60:30:390)
set(gca,'XTickLabel',[60:30:180 -150:30:0])
set(gca,'YTick',-80:20:80)
set(gca,'YTickLabel',80:-20:-80)

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2: Find the ocean volume that has originated from each    %
%            surface box.                                           %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of each box %
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.                                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% $$$  NZ = max(k);
% $$$  NY = max(j);
% $$$  NX = max(i);
% $$$  
% $$$  % get vector of volumes.
% $$$  dx = diff(LON(1:2));
% $$$  dy = diff(LAT(1:2));
% $$$  for ny = 1:length(LAT)
% $$$    distx(ny) = sw_dist([LAT(ny) LAT(ny)],[0 dx],'km');
% $$$  end
% $$$  disty = sw_dist([0 dy],[0 0],'km')
% $$$  distx = reshape(distx,length(distx),1); % get dimensions correct.
% $$$  area = 1e6.*disty(ones(NY,NX)).*distx(:,ones(NX,1)); 
% $$$  zface= (DEPTH(1:end-1)+DEPTH(2:end))./2;
% $$$  dz = ([zface(1) ; diff(zface); 500]);
% $$$  vol = permute(area(:,:,ones(NZ,1)),[3 1 2]) .*dz(:,ones(NY,1),ones(NX,1));
% $$$  v   = field_to_vector(vol,i,j,k);

 v = zeros(Nfield,1);
 v(iii(10:end)) = 1;
 v = v./sum(v);
 
 % effectively take inverse of transpose A matrix.
 vtot = R' \ (P' * (L' \ (U' \ (Q' * v))));  % inverse of transpose
                                             % = transpose of inverse
 
 Vtot = vector_to_field(vtot,i,j,k);
 Vtot = squeeze(Vtot(1,:,:)); % scale the volume by surface area.
 contourf(LON,LAT,log(Vtot))


lonafr = find(LON<40);
lonafr = lonafr(end);
LONplot = [LON(lonafr:end-1); LON(1:lonafr-1)+360];
pvol_plot = log10(squeeze(Vtot(:,[lonafr:end-1 1:lonafr-1])))+2; %percent
figure
hold on
 set(0,'DefaultLineLineWidth',1.6)
 set(0,'DefaultTextFontSize',12)
 set(0,'DefaultAxesLineWidth',1.6)
 set(0,'DefaultAxesFontSize',12)
m_proj('miller','longitudes',[40 400],'latitudes',[-80 85])
m_contourf(LONplot,LAT,pvol_plot,-4:0.5:1)
shading flat
m_coast('linewidth',.5,'color','k')
m_coast('patch',[.5 .5 .5])
m_grid('box','fancy','tickdir','in');
hotreverse
caxis([-4 1])
colorbar
chandle = colorbar;%('FontSize',15)
set(chandle,'FontSize',12)
set(chandle,'YTick',-4:1)
set(chandle,'YTickLabel',[0.0001 0.001 0.01 0.1 1 10])
set(gca,'XTick',60:30:390)
set(gca,'XTickLabel',[60:30:180 -150:30:0])
set(gca,'YTick',-80:20:80)
set(gca,'YTickLabel',80:-20:-80)
xlabel('longitude','FontSize',8)
ylabel('latitude','FontSize',8)

title('Surface origin: 126^{\circ}E, 4^{\circ}S, z> 1500m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% track adjoint pathways?

vtot2 = vtot;
vtot2(~inmixlyr) = -vtot(~inmixlyr);

vtot2fld = vector_to_field(vtot2,it,jt,kt);

vtot3fld = max(vtot2fld,[],1);


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3: Find the surface origin of water for some interior box %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of a
% given interior box,
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.      
% Very similar mathematically to example 2.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


 NZ = max(k);
 NY = max(j);
 NX = max(i);
 Nfield = size(A,1);
 
 % choose an interior location X (Xlon[lon], Xlat [lat], Xdepth [m depth]).
% -7.38, 115.26E
 for np = 1:3
     if np ==1 % inflow
         Xlon = 126;
         Xlat = 4;
         Xdepth = 200;
     elseif np ==2 %outflow
         Xlon = 122;
         Xlat = -12;
         Xdepth = 200;
     elseif np ==3 % itf.
         Xlon = 126;
         Xlat = -4;
         Xdepth = 200;
     end
 end
 
 % Find the coordinate on the grid by linear interpolation. 
 LONtmp = [LON(1)-4 ; LON; LON(end)+4];  %% watch out for the
                                         %wraparound of the grid.
 iX = interp1(LONtmp,0:length(LON)+1,Xlon,'linear')
 jX = interp1(LAT,1:length(LAT),Xlat,'linear')
 kX = interp1(DEPTH,1:length(DEPTH),Xdepth,'linear')
 kX(Xdepth>5500) = 33; % if deeper than the bottom level, put it on
                       % the bottom level.

 % watch out for wraparound again.
 i2 = i;
 i2( i-iX > NX/2) = i2( i-iX > NX/2) - NX;
 i2( iX-i > NX/2) = i2( iX-i > NX/2) + NX;
  
 % find the closest gridpoints to the location of interest.
 [dist,loc] = sort((i2-iX).^2 + (j-jX).^2 + (k-kX).^2);
 dist = dist+0.1; % to eliminate singularity if point matches up
                  % exactly with the grid. 
 
 mask = zeros(Nfield,1); 
 mask(loc(1:6)) =(1./dist(1:6))./sum(1./dist(1:6)); % mask adds up
                                                    % to 1.
 
 vtot = R' \ (P' * (L' \ (U' \ (Q' * mask)))); 
 Vtot = vector_to_field(vtot,i,j,k);
 
 %Vtot = squeeze(Vtot(1,:,:)); % just keep the surface as that is
                              % the ultimate source region.

 %contourf(LON,LAT,Vtot,[0 .001 .005 .01 .05 .1 .5]) % one way to
                                                    % plot it.
 contourf(LON,-DEPTH,squeeze(log10(Vtot(:,16,:))),-8:.5:0)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4: Find the distribution of a tracer given:              %
%       (a) the pathways described by A,                          %
%       (b) interior sources and sinks given by dC,               % 
%           that best fits observations, Cobs,                    %
%   and (c) inequality constraints on the tracer concentration.   %
%                                                                  %
% Mathematically, minimize J = (C-Cobs)^T W (C-Cobs) subject to    %
%                         AC = d + Gamma u                         %
%  where u is the estimated change in surface concentration.       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
  % use temperature as an example, with Tobs the observed temperature vector.
  % load temperature in vector form, along with associated error.
  load tracerobs_4deg_33lev_woce  % could also choose 2deg

  
  dT = Tobs;      % set rhs vector to the observations.
  dT(k>1) = 0;    % no internal sinks or sources.

  % a first guess: observed surface boundary conditions are perfect.  
  Tmod1 =  Q * (U \ (L \ (P * (R \ dT))));
  Tmod1_field = vector_to_field(Tmod1,i,j,k);
 
  % get cost function (J) for first guess temperature field.
  % here the data-model misfit is weighted by the expected error. 
  Nfield = length(Tobs);
  Terr(ibotside) = Terr(ibotside).*50;
  %Terr(isfc) = 0;
  invWT = sparse(1:Nfield,1:Nfield,1./Terr.^2,Nfield,Nfield); % diagonal matrix.
  J1 = (Tmod1-Tobs)'*invWT*(Tmod1-Tobs)./Nfield % optimized J1 = 1
 
  % Get better fit to observations by modifying the surface boundary 
  % condition within its uncertainty.
  % mathematical form: Find uT such that:
  % J = (Tmod-Tobs)'*inv(WT)*(Tmod-Tobs) is minimized 
  % subject to:  A*Tmod = dT + Gamma * uT.
  Nsfc = sum(k==1);
  Gamma = sparse(Nfield,Nsfc);
  nu = 0;
  for nv = 1:Nfield
    if k(nv)==1
      nu = nu+1;
      Gamma(nv,nu) = 1;
    end
  end
  u0 = zeros(Nsfc,1); % first guess of change to surface boundary
                      % conditions.
  lbT = -2.*ones(Nsfc,1); % temperature lower bound: can not freeze.
  ubT = 40.*ones(Nsfc,1); % ad-hoc temperature upper bound: 40 C.

  %% 2 methods: a) compute full Hessian, b) use functional form of
  %Hessian. (b) is preferred.
  method = 3;
  if method==1
    H1 = Q * (U \ (L \ (P * (R \ dT)))) -Tobs; 
    H2 = Q * (U \ (L \ (P * (R \ Gamma))));   % lengthy calculation step: get full
                                 % Hessian. Can be avoided by specifying
				 % functional form of Hessian-vector product.

    HT = H2'*invWT*H2;
    fT = 2.*H1'*invWT*H2;
    uT = quadprog(HT,fT,[],[],[],[],lbT,[],u0); % solve quadratic
                                                 % programming problem for
                                                 % shifts (uT) to
                                                 % sfc. temp.
    
  elseif method==2
    [fval, exitflag, output, uT] = hessianprob(A,Gamma,invWT,Tobs,u0,dT,lbT,ubT);
    Tmod2 =  Q * (U \ (L \ (P * (R \ (dT+Gamma*uT))))) ; % best estimate
    J2 = (Tmod2-Tobs)'*invWT*(Tmod2-Tobs)./Nfield % expect J2 ~ 1, J2<J1

  elseif method==3 % slower because doesn't use Hessian.
      %%
     options = optimset('Algorithm','interior-point','Display','iter', ...
                  'GradObj','on','LargeScale','on');
     %options = optimset('Algorithm','trust-region-reflective','Display','iter', ...
     %        'GradObj','on','LargeScale','on');
     noncons = 0;
     uTtest= fmincon(@(x)objfun(x,A,invWT,Tobs,isfc,iint,noncons),Tobs(isfc),[],[],[],[],lbT,ubT,[],options);
     J2 = objfun(uTtest,A,invWT,Tobs,isfc,iint,noncons) %
                      % expect J2 ~ 1, J2<J1
     d2 = zeros(Nfield,1); d2(isfc) = uTtest;
     Tmod2 =  Q * (U \ (L \ (P * (R \ d2)))) ; % best estimate
  
  end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 5: How to update the TMI pathways and tracer distribution    %
%            to improve the fit to observations.                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% modify directions in Schlitzer 2007.

load ../woce_4deg/A0_woce_4deg_4jan
load ../woce_4deg/stats_4deg_33lev_woce_4jan.mat W*sfc *err
invWTsfc = WTsfc; clear WTsfc;
invWSsfc = WSsfc; clear WSsfc;
invWO18sfc = WO18sfc; clear WO18sfc;
invWPsfc = WPsfc; clear WPsfc;
invWNsfc = WNsfc; clear WNsfc;
invWOsfc = WOsfc; clear WOsfc;

%% Get all 6 tracers.
%% add something for ibotside?
Terr_nosfc = Terr; Terr_nosfc(isfc) = Inf;
invWTsm = sparse(1:Nfield,1:Nfield,1./Terr_nosfc.^2,Nfield,Nfield); % diagonal matrix.
invWTsm = invWTsm + Gamma*invWTsfc*Gamma'; % weighting matrix with
                                           % surface smoothing

Serr_nosfc = Serr; Serr_nosfc(isfc) = Inf;
invWSsm = sparse(1:Nfield,1:Nfield,1./Serr_nosfc.^2,Nfield,Nfield);
invWSsm = invWSsm + Gamma*invWSsfc*Gamma'; 

O18err_nosfc = O18err; O18err_nosfc(isfc) = Inf;
invWO18sm = sparse(1:Nfield,1:Nfield,1./O18err_nosfc.^2,Nfield,Nfield);
invWO18sm = invWO18sm + Gamma*invWO18sfc*Gamma';

Perr_nosfc = Perr; Perr_nosfc(isfc) = Inf;
invWPsm = sparse(1:Nfield,1:Nfield,1./Perr_nosfc.^2,Nfield,Nfield);
invWPsm = invWPsm + Gamma*invWPsfc*Gamma';

Nerr_nosfc = Nerr; Nerr_nosfc(isfc) = Inf;
invWNsm = sparse(1:Nfield,1:Nfield,1./Nerr_nosfc.^2,Nfield,Nfield);
invWNsm = invWNsm + Gamma*invWNsfc*Gamma';

Oerr_nosfc = Oerr; Oerr_nosfc(isfc) = Inf;
invWOsm = sparse(1:Nfield,1:Nfield,1./Oerr_nosfc.^2,Nfield,Nfield);
invWOsm = invWOsm + Gamma*invWOsfc*Gamma';

Wall{1} = invWTsm;
Wall{2} = invWSsm;
Wall{3} = invWO18sm;
Wall{4} = invWPsm;
Wall{5} = invWNsm;
Wall{6} = invWOsm;

Wpno{1} = invWPsm;
Wpno{2} = invWNsm;
Wpno{3} = invWOsm;

%% use objfun without changing paths -- but with 6 tracers and
%% nonconservative effects.



%%
m = A(iintdst+(iintsrc-1).*Nfield)';
x0 = [Tmod2(k==1); m];
Atemplate = A;
Atemplate(~inmixlyr,:) = 0;

D = sparse(iintdst,1:length(iintsrc),ones(Nm,1),Nfield,Nm);
D(inmixlyr,:) = [];

%%
%[Jtest1,gtest1] = objfun(x0,Atemplate,Wall,Tobs,isfc,iint,0, ...
%                         iintdst,iintsrc); 


% test case: take A as fixed. Compare J with optimized sfc + dp
% vs. without.

x0 = [Tobs(isfc);Sobs(isfc);O18obs(isfc);Pobs(isfc);Nobs(isfc); ...
      Oobs(isfc);dP(iint)];
cobs = [Tobs Sobs O18obs Pobs Nobs Oobs];
nclist = [0 0 0 1 15.5 -170];
[J0,g0] = objfun(x0,A,Wall,cobs,isfc,iint,nclist); 

delta = 1e-3; deltaloc = 50000;
x_path = x0;
x_path2 = x_path; x_path2(deltaloc) = x_path2(deltaloc) + delta;
[J2A,g2A] = objfun(x_path,A,Wall,cobs,isfc,iint,nclist);
[J2B,g2B] = objfun(x_path2,A,Wall,cobs,isfc,iint,nclist);
(g2B(deltaloc)+g2A(deltaloc))./2
(J2B-J2A)./delta

%% P,N, O optimization.
lb = [zeros(3*Nsfc,1); zeros(Nint,1)];
ub = [];

%ub = [35 + 0.*Tmod2(k==1); 1 + 0.*m];
%ub = [35 + 0.*Tmod2(k==1)];

while J2 > 1
   options = optimset('Algorithm','interior-point','Display','iter', ...
                   'GradObj','on','LargeScale','on','MaxIter',20, ...
                   'MaxFunEvals',40,'Hessian',{'lbfgs',4});
   %options = optimset('Algorithm','trust-region-reflective','Display','iter', ...
   %                   'GradObj','on','LargeScale','on','MaxIter',5, ...
   %                   'MaxFunEvals',10); % too much memory.

 [x_path,fval,exitflag,output]= fmincon(@(x)objfun(x,A, ...
     Wpno,sq(cobs(:,4:6)),isfc,iint,nclist(4:6)),x0(3*Nsfc+1:end), ...
     [],[],[],[],lb(3*Nsfc+1:end),ub,[],options);




   noncons = 0;
   [x_path,fval,exitflag,output]= fmincon(@(x)objfun(x,Atemplate,invWTsm,Tobs, ...
     isfc,iint,noncons,iintdst,iintsrc),x_path,[],[],[],[],lb,ub,[],options);

   %x_path= fmincon(@(x)objfun(x,Atemplate,invWT,Tobs, ...
   %  isfc,iint,noncons,iintdst,iintsrc),x0,[],[],[],[],lb,ub,[],options);

   [J2,g2] = objfun(x_path,Atemplate,invWTsm,Tobs,isfc,iint,noncons,iintdst,iintsrc); %
                      % expect J2 ~ 1, J2<J1

   % renormalize the "m" part of the state vector.
   diagsum = D*x_path(end-Nm+1:end);

   M = sparse(iintdst,iintsrc,x_path(end-Nm+1:end),Nfield,Nfield);
   Diagsum = sparse(1:Nint,1:Nint,1./diagsum,Nint,Nint);
   M(iint,:) = Diagsum*M(iint,:);
   
   m = M(iintdst+(iintsrc-1).*Nfield);
   x_path = [x_path(1:2806); m'];
   save work_1dec2011
end


%% check the gradient.                      
delta = 1e-3
deltaloc = 2000;

x_path2 = x_path; x_path2(deltaloc) = x_path2(deltaloc) + delta;
[J2A,g2A] = objfun(x_path,Atemplate,invWTsm,Tobs,isfc,iint,noncons,iintdst,iintsrc); 
[J2B,g2B] = objfun(x_path2,Atemplate,invWTsm,Tobs,isfc,iint, ...
                   noncons,iintdst,iintsrc); 
(g2B(deltaloc)+g2A(deltaloc))./2
(J2B-J2A)./delta

%%
load ../woce_4deg/grid_4deg_33lev_woce.mat

% Go through all interior locations to clean up indices.
iint = find(~inmixlyr);
imix = find(inmixlyr);

% clean up my destination and source info to not include mixed-layer
% destinations.
iintsrc(isnan(iintsrc))= [];
iintdst(isnan(iintdst))= [];
iintdst = transpose(iintdst);

for nm = 1:length(imix);
   
   nv = imix(nm)
   iintsrc(iintdst==nv)=[];
   iintdst(iintdst==nv)=[];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6: How to modify or reconstruct the TMI pathways,            %
%            i.e., the A matrix.                                       %
%                                                                      %
% Each row of A is determined by the solution of a non-negative        %
% least squares problem via the method of Lawson and Hanson.           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Pick an example box number "nv". Eventually must go through every box.  
% nv = 30000;

% If in mixed layer, then stop. The local inversion is not applicable.
if ~inmixlyr(nv)
  
  %% Find the neighboring gridcells.
  %  Account for longitudinal wraparound.
  it2 = it;
  it2( it-it(nv) > NX/2) = it2( it-it(nv) > NX/2) - NX;
  it2( it(nv)-it > NX/2) = it2( it(nv)-it > NX/2) + NX;
  boxsize =  1;

  %% six sides to a box.
  isource1=find( abs(it2-it2(nv))==boxsize & ...
              jt + boxsize > jt(nv)   & ...
              jt - boxsize < jt(nv)   & ...
              kt + boxsize > kt(nv)   & ...
              kt - boxsize < kt(nv)  );
  isource2=find( abs(jt-jt(nv))==boxsize & ...
              it2 + boxsize > it2(nv)   & ...
              it2 - boxsize < it2(nv)   & ...
              kt + boxsize > kt(nv)   & ...
              kt - boxsize < kt(nv)  );
  isource3=find( abs(kt-kt(nv))==boxsize & ...
              it2 + boxsize > it2(nv)   & ...
              it2 - boxsize < it2(nv)   & ...
              jt + boxsize > jt(nv)   & ...
              jt - boxsize < jt(nv)  );
  isource = [isource1; isource2; isource3];   
  NS = length(isource);

  %% set parameters of the least-squares problem.
  alpha2 = 5;
  obs = [Tobs(nv);Sobs(nv);O18obs(nv);Pobs(nv);Nobs(nv);Oobs(nv)];
  NOBS = length(obs);
     
  %% add an equation for mass conservation.
  obs = [ obs ; 1];

  NS = length(srctmp); % how many sources?

  %% add tapering to the solution.
  % But do it this funny way to fit the matlab lsqnonneg function.
  xbase = (1./NS).*ones(NS,1);
  xbase(end+1) = 0;
  xbase = reshape(xbase,length(xbase),1);
    
  %% Set up the local properties matrix, E.
  Etmp = [T(isource)'; S(isource)'; O18(isource)';P(isource)'; ...
          N(isource)';O(isource)'];
  Etmp(end+1,:) = ones(1,NS);

  %% added next 2 lines for nonconservative tracers.
  Etmp(:,end+1) = zeros(NOBS+1,1);
  Etmp(:,end) = [0 0 0 1 15.5 -170 0];

  obserr =[Terr(nv); Serr(nv); O18err(nv); Perr(nv); ...
                Nerr(nv); Oerr(dsttmp)];
  obserr = [obserr; .001]; % .001 is for mass conservation.
       
  W = diag(1./obserr); % equal to W^(-1/2) in TOCIP, Wunsch 1996.
       
  obsnew = W*obs;
  Enew = W*Etmp;
  if alpha2>0 % then tapering is on.
    Enew = [Enew ; (alpha2./NSlocal).*eye(NSlocal+1)];
    obsnew(end+1:end+length(xbase)) = (alpha2./NSlocal).*xbase;
  end

  x0 = zeros(NS+1,1); % first guess of solution, x ("m" in G&H, 2010)
  [x1,Jtmp,restmp,exflagnew] = lsqnonneg(Enew,obsnew,x0,options); % check to see if Jtmp is reasonable.

  % modify/overwrite the global pathways "A" matrix.
  A(nv,:) = 0;
  A(nv,nv) = -1;
  A(nv, isource) = x1(1:end-1);
  dP = x1(end); % nonconservative source.
end
   
   

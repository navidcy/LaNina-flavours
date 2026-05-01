
% Figure 2: Equatorial Pacific SST anomaly structure during CP and EP La Niña
% events as well as their difference 
%
% This script generates Figure 2 from:
% "Characteristics and impacts of La Niña diversity on Pacific
% teleconnections" by Freund et al.
%
% The figure shows zonal SST anomaly (SSTa) profiles averaged over
% 10°S–10°N for:
%   - Central Pacific (CP) La Niña events
%   - Eastern Pacific (EP) La Niña events
%   - CP minus EP differences
%
% SSTa profiles are calculated for:
%   - SON
%   - DJF
%   - MAM(+1)
%
% Three observational SST datasets are used:
%   - HadISST
%   - ERSST
%   - COBE
%
% Shading indicates the event-to-event range across individual
% La Niña composites.


%% Calculate ENSO SST composites for different types of La Nina events 

% create input for Figure  
% read in three SST datasets: HadISST, ERSST, COBE


% regrid data on common grid 

% 1) ERSST
load ERSST_1854_2023.mat;


% first recenter ERSST 
[Lat,Lon] = meshgrid(lat,lon);
[lat_new,lon_new,sst_recentered] = recenter(Lat,Lon,sst);


% regrid to match HadISST coordinates 
lon1 = [-179.5:1:179.5]';
lat1 = flipud([-89.5:1:89.5]');
[Lat1,Lon1] = meshgrid(lat1,lon1);


% interpolate the first two dimensions (lon and lat) and leave third
% dimension the same (time)
n = length(sst_recentered);

% the default is linear interpolation
XSST  = zeros(360,180,n);

for i = 1:n
    
   XSST(:,:,i) = interp2(lat_new, lon_new,sst_recentered(:,:,i),Lat1,Lon1); 
    
end

ERSST = XSST;

ERSST_resized = reshape(ERSST(:,:,1:2040),[360,180 12 170]);


%% 2) COBE

load COBE_1891_2023.mat;

[Lat,Lon] = meshgrid(lat,lon);
[lat_new,lon_new,sst_recentered] = recenter(Lat,Lon,sst);


% regrid to match HadISST coordinates 
lon1 = [-179.5:1:179.5]';
lat1 = flipud([-89.5:1:89.5]');
[Lat1,Lon1] = meshgrid(lat1,lon1);

% interpolate the first two dimensions (lon and lat) and leave third
% dimension the same (time)
n = length(sst_recentered);

% the default is linear interpolation
XSST  = zeros(360,180,n);

for i = 1:n
    
   XSST(:,:,i) = interp2(lat_new, lon_new,sst_recentered(:,:,i),Lat1,Lon1); 
    
end

COBE = XSST;

COBE_resized = reshape(COBE(:,:,1:1584),[360 180 12 132]);
COBE_resized(:,:,:,133)=NaN;
COBE_resized(:,:,1:6,133) = COBE(:,:,1585:1590);


%% 3) HadISST

 load HadISST_1870_2023.mat;

 HadISST = sst;
 HadISST_resized = SST_resized;

 %% shorten all datasets to post 1900
 % set time period for composites here: 1900 to 2023

 HadISST_res=HadISST_resized(:,:,:,31:154); 

 ERSST_res=ERSST_resized(:,:,:,47:170);

 COBE_res=COBE_resized(:,:,:,10:133);

 clear HadISST_resized ERSST_resized COBE_resized
 


%%
% Get La Nina years 

t1 = datetime(1900,1,15,0,0,0);
t2 = datetime(2023,3,15,0,0,0);
yearly_T = (t1:calendarDuration(1,0,0,0,0,0):t2)';
    
% set base period: 1981-2010
base_period = datetime(1981,1,15,0,0,0):calendarDuration(1,0,0,0,0,0):datetime(2010,12,15,0,0,0);
% find index of base period in the specific timeseries automatically using
% its index
idx = ismember(yearly_T,base_period);
base_idx= find(idx==1);

% La Nina years from Mandy --> mixed is when 3 datasets are inconsistent regarding
% CP/EP, e.g. 2 out of 3 datasets say CP, then we classify as mixed 
yr_CP = [1904;1910;1925;1934;1939;1943;1951;1955;1956;1965;1971;1972;1974;1975;1976;1985;1989;1999;2000;2001;2008;2011;2012;2021;2023];
yr_EP= [1907;1963;1968;2018;2022];
yr_Mixed= [1909;1911;1917;1918;1923;1950;1996];



%% calculate SST anomalies as seasonal averages 

% loop through all three datasets

for i = 1:3 



    if i==1
        dataset='HadISST'
        SST =HadISST_res;
    elseif i==2
        dataset='ERSST'
        SST =ERSST_res;
    else
        i==3
        dataset='COBE'
        SST = COBE_res;
    end 


length_years = size(SST,4);

% Seasons: MAM to JJA the next year
SST_Year = SST(:,:,3:12,:);
SST_Year(:,:,:,2:length_years) = SST_Year(:,:,:,1:length_years-1);
SST_Year(:,:,:,1) = NaN;
SST_Year(:,:,11:18,:)= SST(:,:,1:8,:);

% take seasonal averages

    SST_Seasons = squeeze(nanmean(SST_Year(:,:,1:3,:),3)); % MAM
    SST_Seasons(:,:,:,2)=squeeze(nanmean(SST_Year(:,:,4:6,:),3)); % JJA
    SST_Seasons(:,:,:,3)=squeeze(nanmean(SST_Year(:,:,7:9,:),3)); % SON
    SST_Seasons(:,:,:,4)=squeeze(nanmean(SST_Year(:,:,10:12,:),3)); % DJF
    SST_Seasons(:,:,:,5)=squeeze(nanmean(SST_Year(:,:,13:15,:),3)); % MAM+
    SST_Seasons(:,:,:,6)=squeeze(nanmean(SST_Year(:,:,16:18,:),3)); % JJA+

 % calculate climatology 

SST_Seasons_clim = squeeze(nanmean(SST_Seasons(:,:,base_idx,:),3));

i = [];

for s=1:size(SST_Seasons,4)
i = 1:1:length(yearly_T);
SST_Seasons_anomaly(:,:,i,s) = SST_Seasons(:,:,i,s) - SST_Seasons_clim(:,:,s); 
end 

i = [];   

% detrend SST anomalies before generating composites
for s = 1:size(SST_Seasons,4)
input_SST(:,:,:,s) = detrend3(SST_Seasons_anomaly(:,:,:,s),'omitnan');
end


tstart = datenum(yearly_T(1));
tend = datenum(yearly_T(end));



T = yearly_T.Year;

% this finds each year with an ENSO event and saves the index 
Nina_CP_idx_yr = zeros(1,length(yr_CP));
for i = 1:length(yr_CP)
Nina_CP_idx_yr (1,i) = find(T==yr_CP(i));
end

Nina_EP_idx_yr = zeros(1,length(yr_EP));
for i = 1:length(yr_EP)
Nina_EP_idx_yr (1,i) = find(T==yr_EP(i));
end


Nina_Mixed_idx_yr = zeros(1,length(yr_Mixed));
for i = 1:length(yr_Mixed)
Nina_Mixed_idx_yr (1,i) = find(T==yr_Mixed(i));
end

% Create composites

i = 1:1:length(Nina_CP_idx_yr);
Nina_composite_CP = input_SST(:,:,Nina_CP_idx_yr(i),:);

i = 1:1:length(Nina_EP_idx_yr);
Nina_composite_EP = input_SST(:,:,Nina_EP_idx_yr(i),:);

i = 1:1:length(Nina_Mixed_idx_yr);
Nina_composite_Mixed = input_SST(:,:,Nina_Mixed_idx_yr(i),:);


%set input data 

if strcmp("HadISST",dataset)

    Nina_composite_CP_HadISST=Nina_composite_CP; 
    Nina_composite_EP_HadISST=Nina_composite_EP; 
    Nina_composite_Mixed_HadISST=Nina_composite_Mixed;
    % combine all along third dimension to get all La Ninas in one matrix
    Nina_composite_All_HadISST=cat(3,Nina_composite_CP_HadISST,Nina_composite_EP_HadISST,Nina_composite_Mixed_HadISST); 
elseif strcmp("ERSST",dataset)
    Nina_composite_CP_ERSST=Nina_composite_CP; 
    Nina_composite_EP_ERSST=Nina_composite_EP; 
    Nina_composite_Mixed_ERSST=Nina_composite_Mixed;
    % combine all along third dimension to get all La Ninas in one matrix
    Nina_composite_All_ERSST=cat(3,Nina_composite_CP_ERSST,Nina_composite_EP_ERSST,Nina_composite_Mixed_ERSST);
else
    strcmp("COBE",dataset)
    Nina_composite_CP_COBE=Nina_composite_CP; 
    Nina_composite_EP_COBE=Nina_composite_EP; 
    Nina_composite_Mixed_COBE=Nina_composite_Mixed;
    % combine all along third dimension to get all La Ninas in one matrix
    Nina_composite_All_COBE=cat(3,Nina_composite_CP_COBE,Nina_composite_EP_COBE,Nina_composite_Mixed_COBE);

end 


end 
  



%% Prepare figure 


% get zonal and meridional region 
slice_lon_1=1:90; %(Dateline to 90 degrees west (eastern Pacific)
slice_lon_2= 301:360; %120 E to the dateline
slice_lon=[slice_lon_2 slice_lon_1];
slice_lat=80:101; % 10 degrees North and South

slice = 'zonal'


% 1) generate input data for each season and save - SON DJF and MAM+1
% needed 

for season = 3:5 % these are the indices for SON, DJF and MAM %3=SON 4=DJF 5=MAM


SSTa_slice=zeros(length(slice_lon),length(slice_lat),9);
SSTa_slice_all=cell(1,9);


for i=1:9
    
  
    
        if i==1 
        phase = 'La Nina';
        type  = 'CP';
        ENSO_Index=1;
        Dataset = 1;
        input_plot_all = Nina_composite_CP_HadISST;  

        elseif i==2
        phase = 'La Nina';
        type  = 'CP';
        ENSO_Index=1;
        Dataset = 2;
        input_plot_all = Nina_composite_CP_ERSST;

        elseif i==3
        phase = 'La Nina';
        type  = 'CP';
        ENSO_Index=1;
        Dataset = 3;
        input_plot_all = Nina_composite_CP_COBE;
       
        
        elseif i==4 
        phase = 'La Nina';
        type  = 'EP';
        ENSO_Index=2 ;
        Dataset = 1;
        input_plot_all = Nina_composite_EP_HadISST;  

        
        elseif i==5 
        phase = 'La Nina';
        type  = 'EP';
        ENSO_Index=2 ;
        Dataset = 2;
        input_plot_all = Nina_composite_EP_ERSST;

        
        elseif i==6 
        phase = 'La Nina';
        type  = 'EP';
        ENSO_Index=2 ;
        Dataset = 3;
        input_plot_all = Nina_composite_EP_COBE;

        
        elseif i==7 
        phase = 'La Nina';
        type = 'Mixed';
        ENSO_Index=3;
        Dataset = 1;
        input_plot_1 = squeeze(nanmean(Nina_composite_CP_HadISST,3));
        input_plot_2 = squeeze(nanmean(Nina_composite_EP_HadISST,3));
        

        elseif i==8 
        phase = 'La Nina';
        type = 'Mixed';
        ENSO_Index=3;
        Dataset = 2; 
        input_plot_1 = squeeze(nanmean(Nina_composite_CP_ERSST,3));
        input_plot_2 = squeeze(nanmean(Nina_composite_EP_ERSST,3));

        else 
        i==9 
        phase = 'La Nina';
        type = 'Mixed';
        ENSO_Index=3;
        Dataset = 3; 
        input_plot_1 = squeeze(nanmean(Nina_composite_CP_COBE,3));
        input_plot_2 = squeeze(nanmean(Nina_composite_EP_COBE,3));


        
        end 
    
 
input_plot = []

if i<=6
% average for composite
input_plot = squeeze(nanmean(input_plot_all,3));
SSTa_slice(:,:,i) = input_plot(slice_lon,slice_lat,season);

    SSTa_slice_all{i} = squeeze(nanmean(input_plot_all(slice_lon,slice_lat,:,season),2));
else
    i>=7;
input_plot= input_plot_1 - input_plot_2; 
SSTa_slice(:,:,i) = input_plot(slice_lon,slice_lat,season);

end 



end 


   SSTa_slice = squeeze(nanmean(SSTa_slice,2));

   if season == 3
        SSTa_slice_SON = SSTa_slice;
        SSTa_slice_all_SON = SSTa_slice_all;

    elseif season == 4
        SSTa_slice_DJF = SSTa_slice;
        SSTa_slice_all_DJF = SSTa_slice_all;
    
    elseif season == 5
        SSTa_slice_MAM = SSTa_slice;
        SSTa_slice_all_MAM = SSTa_slice_all;
    end


    
end 

%% Plot   
    

figure('pos',[10 10 1200 600])

input_lat=glat(slice_lat);
input_lon=glon(slice_lon);

plot_coords = input_lon;


xlabelstr='Longitude (^{\circ})';
xticklabelstr = {'120^{\circ} E','150^{\circ} E','180^{\circ} E','150^{\circ} W','120^{\circ} W','90^{\circ} W'};


for plot_num = 1:9 % 9 groups in total

if plot_num == 1 || plot_num ==2 || plot_num == 3 
    phase = 'La Nina';
    type  = 'CP';
elseif plot_num == 4 || plot_num == 5 || plot_num == 6
    phase = 'La Nina';
    type  = 'EP';
else
    plot_num == 7 || plot_num == 8 || plot_num == 9
    phase = 'La Nina';
    type  = 'Mixed';
end
    
  

for Datasets = 1:3

    
if strcmp('La Nina',phase) && strcmp('CP',type)
    event_index = Datasets;
elseif strcmp('La Nina',phase) && strcmp('EP',type)
     event_index = Datasets+3;
else
    strcmp('La Nina',phase) && strcmp('Mixed',type)
    event_index = Datasets+6;
end 
    


hold on 

SSTa_slice_all = [];
SSTa_slice = [];


% get the right seasons:
if plot_num == 1 || plot_num == 4 
    SSTa_slice_all = SSTa_slice_all_SON; 
    SSTa_slice = SSTa_slice_SON;
    
elseif plot_num == 2 || plot_num == 5 
    SSTa_slice_all = SSTa_slice_all_DJF; 
    SSTa_slice = SSTa_slice_DJF;
elseif plot_num == 3 || plot_num == 6 
    SSTa_slice_all = SSTa_slice_all_MAM; 
    SSTa_slice = SSTa_slice_MAM;
elseif plot_num == 7
    SSTa_slice = SSTa_slice_SON;
elseif plot_num ==8
    SSTa_slice = SSTa_slice_DJF;
else
    plot_num=9;
    SSTa_slice = SSTa_slice_MAM;

end 


if plot_num<=6
All_ENSO = SSTa_slice_all{event_index};
else

end 


if plot_num ==1 
    a_1=subplot(3,3,plot_num);
elseif plot_num ==2
    a_2 = subplot(3,3,plot_num);
elseif plot_num ==3
    a_3 = subplot(3,3,plot_num) ;
    
elseif plot_num ==4
    c_1 = subplot(3,3,plot_num);
elseif plot_num ==5
    c_2 = subplot(3,3,plot_num) ;
elseif plot_num ==6
    c_3 = subplot(3,3,plot_num) ;
    
elseif plot_num ==7
    e_1 = subplot(3,3,plot_num);
elseif plot_num ==8
    e_2 = subplot(3,3,plot_num) ;
else
    plot_num=9
    e_3 = subplot(3,3,plot_num) ;

end 





if Datasets==1
    p1=plot(SSTa_slice(:,event_index),'color','b','LineWidth',1.5);
elseif Datasets==2
    SSTa_slice(61,event_index)=SSTa_slice(60,event_index);
    SSTa_slice(62,event_index)=SSTa_slice(63,event_index);
    p2=plot(SSTa_slice(:,event_index),'color','r','LineWidth',1.5);
else
    Datasets==3
    p3=plot(SSTa_slice(:,event_index),'color','k','LineWidth',1.5);
end 

xticks([0 30 60 90 120 150])
xticklabels([])




if plot_num<=6
% the below colours in the range of SSTa during individual ENSO events
% below and above the average (as a polygon), with light pink for the SSTa
% event range for positive IPO and light blue for negative IPO 
bounds_positive = max(All_ENSO,[],2);
bounds_negative = min(All_ENSO,[],2);

% get rid of NaNs/replace with nearest coordinate values 
bounds_positive(61,1) = bounds_positive(60,1);
bounds_positive(62,1) = bounds_positive(63,1);

bounds_negative(61,1) = bounds_negative(60,1);
bounds_negative(62,1) = bounds_negative(63,1);

x_1=1:150;
x_2=flip(x_1);
x2 = [x_1 x_2];


if Datasets==1
    hold on 
    pol1=plot_ts_patch(x_1,bounds_positive,bounds_negative,rgb('light blue'));
    pol1.FaceAlpha = 0.7;
elseif Datasets==2
    hold on 
    pol2=plot_ts_patch(x_1,bounds_positive,bounds_negative,rgb('light pink'));
    pol2.FaceAlpha = 0.6;
else
    Datasets=3
    hold on 
    pol3=plot_ts_patch(x_1,bounds_positive,bounds_negative,rgb('light yellow'));
    pol3.FaceAlpha = 0.4;
end 




else

end 



hold on 
hline(0,'--');


hold on 


end 





for Datasets = 1:3
    
    if strcmp('La Nina',phase) && strcmp('CP',type)
    event_index = Datasets;
elseif strcmp('La Nina',phase) && strcmp('EP',type)
     event_index = Datasets+3;
else
    strcmp('La Nina',phase) && strcmp('Mixed',type)
    event_index = Datasets+6;
    end 

if Datasets==1
    p1=plot(SSTa_slice(:,event_index),'color','b','LineWidth',2);

elseif Datasets==2
    SSTa_slice(61,event_index)=SSTa_slice(60,event_index);
    SSTa_slice(62,event_index)=SSTa_slice(63,event_index);
    p2=plot(SSTa_slice(:,event_index),'color','r','LineWidth',2);

else
    Datasets==3;
    p3=plot(SSTa_slice(:,event_index),'color','y','LineWidth',2);


end 
end



% set ylim 
if strcmp('La Nina',phase) && strcmp('Mixed',type) 
    ylim([-1 0.5])
         
elseif strcmp('La Nina',phase) && strcmp('EP',type)
    ylim([-2 1])
    
else 
    strcmp('La Nina',phase) && strcmp('CP',type)
     ylim([-2 1])

     
end 

if plot_num == 9
lgd=legend([p1 p2 p3],'HadISST','ERSST','COBE','Location','southeast')
lgd.FontSize = 12;
lgd.Position = [0.88 0.13 0.1 0.05]
end 



if plot_num == 1 || plot_num == 2 || plot_num == 3 
    ylabel('^{\circ}C','FontSize',12)
elseif plot_num == 7

else

end 




if plot_num == 1
    title('CP La Niña','FontSize',14,'FontWeight','bold')
elseif plot_num == 4
    title('EP La Niña','FontSize',14,'FontWeight','bold')
elseif plot_num == 7
    title('CP-EP','FontSize',14,'FontWeight','bold')
else

end 


% get rid of longitude labels as they are shown in the below plot
xticklabels([])


if plot_num == 6 
    xlabel(xlabelstr,'FontSize',14)
end 



if plot_num == 4 || plot_num == 5 || plot_num == 6
    yticklabels([])
elseif plot_num == 7 || plot_num == 8 || plot_num == 9
    yticklabels([])
else

end 

if plot_num == 1 || plot_num == 2 || plot_num == 3|| plot_num == 4 || plot_num == 5|| plot_num == 6 
    xticklabels([])


end 


if plot_num == 3 || plot_num == 6  
xticklabels({'120^{\circ} E','150^{\circ} E','180^{\circ} E','150^{\circ} W','120^{\circ} W','             '})
elseif plot_num == 9
   xticklabels({'120^{\circ} E','150^{\circ} E','180^{\circ} E','150^{\circ} W','120^{\circ} W','90^{\circ} W'}) 
end 

end 



%% Set labels, sizes of plots etc 



% SON row
set(a_1,'Position',[0.09 0.70 0.28 0.25]); % CP, SON
set(c_1,'Position',[0.39 0.70 0.28 0.25]); % EP, SON
set(e_1,'Position',[0.69 0.70 0.28 0.25]); % CP-EP, SON

% DJF row
set(a_2,'Position',[0.09 0.40 0.28 0.25]); % CP, DJF
set(c_2,'Position',[0.39 0.40 0.28 0.25]); % EP, DJF
set(e_2,'Position',[0.69 0.40 0.28 0.25]); % CP-EP, DJF

% MAM row
set(a_3,'Position',[0.09 0.10 0.28 0.25]); % CP, MAM
set(c_3,'Position',[0.39 0.10 0.28 0.25]); % EP, MAM
set(e_3,'Position',[0.69 0.10 0.28 0.25]); % CP-EP, MAM





% add season labels
annotation('textarrow',[.04 .86],[.86 .86],'String','SON',...
    'HeadStyle','none','LineStyle','none','FontSize',16,'FontWeight','bold','TextRotation',90);

annotation('textarrow',[.04 .57],[.57 .57],'String','DJF',...
    'HeadStyle','none','LineStyle','none','FontSize',16,'FontWeight','bold','TextRotation',90);

annotation('textarrow',[.04 .28],[.28 .28],'String','MAM+',...
    'HeadStyle','none','LineStyle','none','FontSize',16,'FontWeight','bold','TextRotation',90);


% a,b,c
annotation('textbox',[.09 .8 .1 .2],'String','a','FontSize',12,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[.39 .8 .1 .2],'String','b','FontSize',12,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[.69 .8 .1 .2],'String','c','FontSize',12,'FontWeight','bold','EdgeColor','none')

annotation('textbox',[.09 .2 .1 .2],'String','g','FontSize',12,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[.39 .2 .1 .2],'String','h','FontSize',12,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[.69 .2 .1 .2],'String','i','FontSize',12,'FontWeight','bold','EdgeColor','none')


annotation('textbox',[.09 .5 .1 .2],'String','d','FontSize',12,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[.39 .5 .1 .2],'String','e','FontSize',12,'FontWeight','bold','EdgeColor','none')
annotation('textbox',[.69 .5 .1 .2],'String','f','FontSize',12,'FontWeight','bold','EdgeColor','none')





%% Save as png 

print ('Fig_2_SSTa_profiles_La_Nina_CP_EP_diff','-r300', '-dpng')



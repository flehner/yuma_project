% ------------------
% yuma_analysis_for_github.m
% Author: Flavio Lehner, flehner@ucar.edu
% See README.txt for further details
% ------------------

close all
clear all

% -- what to plot --
plot1   = 0; % time series
plot2   = 0; % same as plot1, but zoomed in to see backwater effect
plot3   = 0; % WOND vs precip anomalies
plot4   = 0; % show-case plot for Bulletin, various time series
plot5   = 0; % gain forecasts
plot6   = 0; % Contingency table for loss-gain forecast
plot7   = 1; % Contingency plot for loss-gain forecast
% ------------------


pathin = '~/Dropbox/work/yuma_project/'; % adjust as needed

years = 2003:2017; % years to consider initially

% -- streamflow data --
pathin0   = [pathin '/streamflow_data/'];
filein    = 'CRIT-PVID_Order-Actual_final.csv';
data      = read_mixed_csv([pathin0 filein],',');
time      = data(2:end,1);
time_flow = datetime(time, 'InputFormat', 'yyyy-MM-dd');
time_flow = time_flow(1:end-3);
data      = str2double(data(2:end,:));
% -- cut to common length
idx       = find(time_flow==datetime(years(1),1,1)):find(time_flow==datetime(years(end)+1,1,3));
time_flow = time_flow(idx);
data      = data(idx,:);

% -- order and diversion data Imperial Dam:
tmp               = read_mixed_csv([pathin 'streamflow_data/reclamation_orders_and_diversions_2012-2017.csv'],',');
order_all         = str2double(tmp(2:end,4));
diversion_all     = str2double(tmp(2:end,5));
wond              = order_all-diversion_all;
time_order_raw    = tmp(2:end,1);
time_order_all    = datetime(time_order_raw,'InputFormat','yyyy-MM-dd');


weekdays = {'Mo','Tu','We','Th','Fr','Sa','So'};

release = data(1:end-3,2);
arrival = data(4:end,3);
div_crit  = data(:,5);
ret_crit  = data(:,10);
div_pvid  = data(:,8);
ret_pvid1 = data(:,11);
ret_pvid2 = data(:,12);
below_pvd = data(:,13); % flow just below palo verde dam
loss_crit   = ret_crit-div_crit;
loss_pvid1  = ret_pvid1-div_pvid;
loss_pvid2  = (ret_pvid1+ret_pvid2)-div_pvid;
loss_en_route1 = loss_crit + loss_pvid1;
loss_en_route2 = loss_crit + loss_pvid2;
lg = arrival-release;
lg_corrected = lg-loss_en_route2(2:end-2);

drop                  = NaN(length(release),1);
drop(2:end)           = release(2:end)-release(1:end-1);

% -- create and subtract lg climatology
lg_corrected_noleap = lg_corrected;
leap_years = [2004,2008,2012,2016,2020];
for i = 1:length(leap_years)
  d = datetime(leap_years(i),2,29);
  lg_corrected_noleap(find(time_flow==d)) = [];
end
for d = 1:365
  lg_corrected_noleap_clim(d) = nanmean(lg_corrected_noleap(d:365:end));
end
hold on
x = filt121([lg_corrected_noleap_clim lg_corrected_noleap_clim lg_corrected_noleap_clim],1e3); % apply 1-2-1 filter 1000 times to smooth annual cycle
lg_corrected_noleap_clim  = x(366:730);
lg_corrected_leap_clim    = [lg_corrected_noleap_clim(1:59) (lg_corrected_noleap_clim(59)+lg_corrected_noleap_clim(60))/2 lg_corrected_noleap_clim(60:end)];
% -- subtract clim to make anomalies
lg_corrected_anom = [];
for i = 1:length(years)
  idx = find(time_flow==datetime(years(i),1,1)):find(time_flow==datetime(years(i),12,31));
  if ismember(years(i),leap_years)
    lg_corrected_anom = [lg_corrected_anom; lg_corrected(idx)-lg_corrected_leap_clim'];
  else
    lg_corrected_anom = [lg_corrected_anom; lg_corrected(idx)-lg_corrected_noleap_clim'];
  end
end
lg_corrected_orig = lg_corrected;
lg_corrected      = lg_corrected_anom; % <<<<--------------------------------- !!!!!!

% -- calculate day-to-day change in release
for d = 1:7
  tmp1                  = release(d+1:7:end);
  tmp2                  = release(d:7:end);
  s                     = 782;%min(length(tmp1),length(tmp2));
  drop_wd(d,:)          = tmp1(1:s)-tmp2(1:s);
  tmp3                  = lg_corrected(d+1:7:end);
  lg_corrected_wd(d,:)  = tmp3(1:s);
end


% -- predict loss-gain due to backwater effect --
model = 2; % 1=one model for all days, 2=individual model for each weekday
  % -- one prediction model for all days
  xx = length(lg_corrected);
  xi = 1:round(0.9*xx);
  % -- regression
  [b] = regress(lg_corrected(xi),[ones(size(drop(xi))) drop(xi)])
  lg_predicted = b(2)*drop;
  % -- prediction models for individual days of the week
  if model == 2
    for d = 1:7
      x = drop(d:7:end);
      y = lg_corrected(d:7:end);
      xx = length(x);
      xi = 1:round(0.9*xx);
      % -- without intercept
      [b] = regress(y,[ones(size(x)) x])
      lg_predicted(d:7:end) = b(2)*x;
    end
  end

% -- subtract prediction
lg_corrected2 = lg_corrected-lg_predicted;


% -- precip forecasts (includes leap days) --
pathin0  = [pathin '/precip_forecasts/new/'];
p_fcst_tmp0 = [];
p_fcst_tmp  = [];
for yyyy = 2002:2017;
  filein  = ['downscaled_forecasts_' num2str(yyyy) '.nc']
  data    = ncread([pathin0 filein],'forecasts');
  % -- reshuffling
  data    = permute(data,[1 3 2]);
  p_fcst_tmp0 = [p_fcst_tmp0 data];
  % -- reshuffling and conversion from mm/day to cfs
  data    = data * 8925.87*1e6 * (1/86400) * 0.0353147; % area of HUC8 * seconds per day * liter --> cubic feet per second;;
  p_fcst_tmp = [p_fcst_tmp data];
end
p_fcst_tmp0 = permute(p_fcst_tmp0,[2 3 1]);
p_fcst_tmp  = permute(p_fcst_tmp,[2 3 1]);
p_fcst_mm     = [p_fcst_tmp0; p_fcst_tmp0(end,:,:)];
p_fcst      = [p_fcst_tmp; p_fcst_tmp(end,:,:)];

% -- create a time vector --
time_p_fcst = datetime(2002, 1, 1):datetime(2017, 12, 31);


% -- precip observations HUC8 (6-hourly values, mm?) --
p_6h            = ncread([pathin 'ccpa/precip_06h_CCPA_2p5km_LCB_NaN_eq_zero_huc8.nc'],'PRECIP_HUC8');
time_p_6h_raw   = ncread([pathin 'ccpa/time.nc'],'yyyymmddhh_begin');
time_p_6h       = datenum(num2str(time_p_6h_raw),'yyyymmddhh');
for i = 1:length(p_6h)/4
  p_all(i)  = sum(p_6h((i-1)*4+1:i*4));
end
tmp             = round(time_p_6h_raw(1:4:end)/100);
time_p_all      = datetime(datenum(num2str(tmp),'yyyymmdd'),'ConvertFrom','datenum');
p_all_mm        = p_all;
% -- conversion from mm/day to cfs
p_all = p_all * 8925.87*1e6 * (1/86400) * 0.0353147; % area of HUC8, seconds per day, liter-->cubic feet per second


% -- common period
% time_common = 2003:2017;
time_common = datetime(2003, 1, 1):datetime(2017, 12, 31);


% -- save data to csv file --
% -- flow and precip:
[dummy,idx3] = intersect(time_p_all,time_common);
data    = [yyyymmdd(time_common)' release arrival loss_en_route2(2:end-2) lg_corrected_orig lg_corrected lg_corrected2 p_all(idx3)']; % p_fcst(idx3,:,:)];
header  = {'Time' 'Release [P]arker' 'Arrival [I]mperial (3 days later)' 'Loss [B]lythe' 'Loss-gain Imperial (I-P-B)' 'Loss-gain Imperial anomalies' 'Loss-gain Imperial anomalies minus backwater effect' 'Precipitation'};
header  = [header;repmat({','},1,numel(header))];
header  = header(:)';
ofile = [pathin '/data_' num2str(yyyymmdd(time_common(1))) '-' num2str(yyyymmdd(time_common(end))) '.csv']
fid=fopen(ofile,'w');
fprintf(fid,'%s\n',cell2mat(header));
dlmwrite(ofile,data,'-append','precision',8);
fclose(fid);

% -- precip obs and forecast:
lt = 1; % lead time (1-7)
[dummy,idx3] = intersect(time_p_all,time_common);
data    = [yyyymmdd(time_common)' p_all(idx3)' squeeze(p_fcst(idx3,lt,:))];
header  = {'Time' 'Release [P]arker' 'Arrival [I]mperial (3 days later)' 'Loss [B]lythe' 'Loss-gain Imperial (I-P-B)' 'Loss-gain Imperial anomalies' 'Loss-gain Imperial anomalies minus backwater effect' 'Precipitation'};
header  = [header;repmat({','},1,numel(header))];
header  = header(:)';
ofile = [pathin '/data_' num2str(yyyymmdd(time_common(1))) '-' num2str(yyyymmdd(time_common(end))) '.csv']
fid=fopen(ofile,'w');
fprintf(fid,'%s\n',cell2mat(header));
dlmwrite(ofile,data,'-append','precision',8);
fclose(fid);


% return






%% -- plotting ----------------------------------
close all


if plot3 == 1

  ylim1 = [-1000 5000];
  ylim2 = [0 12e4];

  % -- selection based on WOND events --
  idx       = find(time_common==datetime(2012,1,1)):find(time_common==datetime(2017,12,31));
  lg_corrected2_tmp = lg_corrected2(idx);
  idx       = find(time_p_all==datetime(2012,1,1)):find(time_p_all==datetime(2017,12,31));
  p_all_tmp = p_all(idx);

  wl    = 4;
  xlim  = [1 2*wl+1];
  xtick = xlim(1):1:xlim(2);
  xl    = [-wl:1:wl];

  clear('wond_events','wond_events_p','wond_events_lg')
  th = nanmean(wond) + 2*nanstd(wond);
  j = 1;
  i = wl;
  while i < length(wond)-wl
    if wond(i) > th
      wond_events(j,:) = wond(i-wl:i+wl);
      wond_events_lg(j,:) = lg_corrected2_tmp(i-wl:i+wl);
      wond_events_p(j,:) = p_all_tmp(i-wl:i+wl);
      j = j+1;
      i = i+wl+1;
    else
      i = i+1;
    end
  end

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 20 15])

  subplot(2,2,1)
  title(['WOND events n = ' num2str(j)])
  hold on
  h1 = plot(wond_events','Color',[1 .7 .7])
  h2 = plot(mean(wond_events),'r.-','LineWidth',2)
  plot(xtick,xtick*0+nanmean(wond)+2*nanstd(wond),'k--')
  hline(mean(wond),'k')
  ylabel('WOND (cfs)')
  xlabel('Time (days)')
  text(2,(nanmean(wond)+2*nanstd(wond))*1.15,['2\sigma'],'FontSize',12)
  text(wl+.5,nanmean(wond)*1.5,['Mean'],'FontSize',12)
  legend([h1(1) h2(1)],'Individual events','Ensemble mean','location','NorthWest')
  legend boxoff
  set(gca,'XLim',xlim,'XTick',xtick,'XTickLabel',xl,'YLim',ylim1)
  box on

  subplot(2,2,2)
  title(['WOND events n = ' num2str(j)])
  hold on
  h1 = plot(wond_events_p','Color',[.7 .7 1])
  h2 = plot(mean(wond_events_p),'b.-','LineWidth',2)
  plot(xtick,xtick*0+2*nanstd(p_all_tmp),'k--')
  hline(mean(p_all_tmp),'k')
  ylabel('Precipitation (cfs)')
  xlabel('Time (days)')
  text(2,(2*nanstd(p_all_tmp))*1.55,['2\sigma'],'FontSize',12)
  legend([h1(1) h2(1)],'Individual events','Ensemble mean','location','NorthWest')
  legend boxoff
  set(gca,'XLim',xlim,'XTick',xtick,'XTickLabel',xl,'YLim',ylim2)
  box on

  % -- selection based on precip events --
  clear('wond_events','wond_events_p','wond_events_lg')
  th = nanmean(p_all_tmp) + 2*nanstd(p_all_tmp);
  j = 1;
  i = wl;
  while i < length(p_all_tmp)-wl
    if p_all_tmp(i) > th
      wond_events(j,:) = wond(i-wl:i+wl);
      wond_events_lg(j,:) = lg_corrected2_tmp(i-wl:i+wl);
      wond_events_p(j,:) = p_all_tmp(i-wl:i+wl);
      j = j+1;
      i = i+wl+1;
    else
      i = i+1;
    end
  end

  subplot(2,2,3)
  title(['Precipitation events n = ' num2str(j)])
  hold on
  h1 = plot(wond_events','Color',[1 .7 .7])
  h2 = plot(mean(wond_events),'r.-','LineWidth',2)
  plot(xtick,xtick*0+nanmean(wond)+2*nanstd(wond),'k--')
  hline(mean(wond),'k')
  ylabel('WOND (cfs)')
  xlabel('Time (days)')
  text(2,(nanmean(wond)+2*nanstd(wond))*1.15,['2\sigma'],'FontSize',12)
  text(wl+.5,nanmean(wond)*1.5,['Mean'],'FontSize',12)
  legend([h1(1) h2(1)],'Individual events','Ensemble mean','location','NorthWest')
  legend boxoff
  set(gca,'XLim',xlim,'XTick',xtick,'XTickLabel',xl,'YLim',ylim1)
  box on

  subplot(2,2,4)
  title(['Precipitation events n = ' num2str(j)])
  hold on
  h1 = plot(wond_events_p','Color',[.7 .7 1])
  h2 = plot(mean(wond_events_p),'b.-','LineWidth',2)
  plot(xtick,xtick*0+2*nanstd(p_all_tmp),'k--')
  hline(mean(p_all_tmp),'k')
  ylabel('Precipitation (cfs)')
  xlabel('Time (days)')
  text(2,(2*nanstd(p_all_tmp))*1.55,['2\sigma'],'FontSize',12)
  legend([h1(1) h2(1)],'Individual events','Ensemble mean','location','NorthWest')
  legend boxoff
  set(gca,'XLim',xlim,'XTick',xtick,'XTickLabel',xl,'YLim',ylim2)
  box on

  tightfig
  % return
  set(gcf,'PaperPositionMode','auto');
  fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/fig6_wond_analysis'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return

end






if plot1 == 1
  % ------------------------------------
  % Fig 1: general overview of data sets
  close all

  xlim = [datetime(2003,1,1) datetime(2017,12,31)];
  i = 1;
  for y = 2003:2017
    xtick(i) = [datetime(y,1,1)];
    i = i+1;
  end

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 25 16])

  subplot(2,1,1)
  hold on
  h1 = plot(time_common,release,'b')
  h2 = plot(time_common,arrival,'r')
  h3 = plot(time_common,loss_en_route2(2:end-2),'Color',[.5 .5 .5])
  hline(0,'k')
  set(gca,'XLim',xlim,'XTick',xtick,'YLim',[-5000 26000],'YTickLabel',[-5000:5000:25000])
  ylabel('Water flow (cfs)')
  xlabel('Time')
  legend([h1 h2 h3],'Release Parker Dam','Arrival Imperial Dam (3 day lag)','Loss in the Blythe and Palo Verde area (CRIT+PVID; 1 day lag)','location','northwest')
  legend boxoff
  box on

  subplot(2,1,2)
  hold on
  h1 = plot(time_common,lg_corrected,'k')
  hline(0,'w')
  set(gca,'XLim',xlim,'XTick',xtick,'YTick',[-4000:2000:4000])
  ylabel('Water flow (cfs)')
  xlabel('Time')
  legend([h1],'Loss-gain between Parker and Imperial Dam','location','northwest')
  legend boxoff
  box on

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/fig1_huc8_timeseries'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return

end



if plot2 == 1
  % ------------------------------------
  % Fig 3: zoomed in to show backwater effect
  close all

  xlim = [datetime(2016,6,1) datetime(2016,7,31)];

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 15 20])

  subplot(2,1,1)
  hold on
  h1 = plot(time_common,release,'b')
  h2 = plot(time_common,arrival,'r')
  h3 = plot(time_common,loss_en_route2(2:end-2),'Color',[.5 .5 .5])
  h4 = plot(time_common,lg_corrected,'k')
  hline(0,'k')
  set(gca,'XLim',xlim,'YLim',[-5000 26000],'YTickLabel',[-5000:5000:25000])
  ylabel('Water flow (cfs)')
  xlabel('Time')
  legend([h1 h2 h3 h4],...
  'Release Parker Dam',...
  'Arrival Imperial Dam (3 day lag)',...
  'Loss in the Blythe and Palo Verde area (CRIT+PVID; 1 day lag)',...
  'Loss-gain between Parker and Imperial Dam',...
  'location','best')
  legend boxoff
  box on

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/fig1_huc8_timeseries_zoomed_in'];
  print('-r300','-loose', '-depsc', ['' fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return

end



if plot4 == 1

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 20 10])
  hold on

  [dummy,idx3] = intersect(time_p_all,time_common);
  h0f = plot(time_common,rm(squeeze(p_fcst(idx3,1,:)),3,1)*.4,'Color',[1 .8 .8])
  h0 = plot(time_common,rm(p_all(idx3),3)*.4,'Color',[235, 131, 52]/255,'LineWidth',1.5,'LineStyle','-')

  [dummy,idx1] = intersect(time_flow,time_common);
  h1 = plot(time_common,lg_corrected(idx1),'Color',[.7 .7 .7]) % corrected for Blythe diversions
  h2 = plot(time_common,lg_corrected2(idx1),'Color',[0 0 0],'LineStyle','-') % corrected for backwater effect
  h3 = plot(time_common,release(idx1),'b-')
  h4 = plot(time_common,arrival(idx1),'Color',[0 .6 0])

  set(gca,'YLim',[-5e3 20e3])
  hline(0,'k-')
  set(gca,'XLim',[datetime(2015,2,1) datetime(2015,3,27)])
  hline(0,'k-')
  ylabel('Water flux (cfs)')

  legend([h1 h2 h3 h4 h0 h0f(1)],'Loss-gain (Blythe area losses accounted for)','Loss-gain backwater effect removed','Release Parker Dam','Arrival Imperial Dam','3-day precipitation','Ensemble of precipitation forecasts','location','northwest')
  box on

  tightfig

  set(gcf,'PaperPositionMode','auto');
  fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/backwater_effect0'];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  save2pdf(['' fileo '.pdf'])
  return
end






if plot5 == 1
  close all

  f = 0.05; % reg coeff
  lt = 1; % lead time

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 27 7])

  ylim = [-1000 4000];
  xlim = [datetime(2005,2,1) datetime(2005,2,28)];

  subplot(1,3,1)
  title('Gain forecasted/Gain observed')
  hold on
  [dummy,idx3] = intersect(time_p_all,time_common);
  hf  = plot(time_common,rm(squeeze(p_fcst(idx3,lt,:)),3,1)*f,'Color',[1 .8 .8])
  hfm = plot(time_common,nanmean(rm(squeeze(p_fcst(idx3,lt,:)),3,1)*f,2),'Color',[235, 131, 52]/255,'LineWidth',1.5,'LineStyle','-')
  hg  = plot(time_common,lg_corrected2,'k')
  hline(2*nanstd(lg_corrected2),'k--')
  text(xlim(1)+1,2*nanstd(lg_corrected2)*1.15,['2\sigma'],'FontSize',12)
  set(gca,'XLim',xlim,'YLim',ylim)
  ylabel('Loss-gain (cfs)')
  legend([hg hf(1) hfm],'Loss-gain observed','Ensemble forecast','Ensemble forecast mean','location','northwest')
  legend boxoff
  box on

  xlim = [datetime(2012,12,1) datetime(2012,12,31)];
  subplot(1,3,2)
  title('Gain forecasted/No gain observed')
  hold on
  [dummy,idx3] = intersect(time_p_all,time_common);
  hf  = plot(time_common,rm(squeeze(p_fcst(idx3,lt,:)),3,1)*f,'Color',[1 .8 .8])
  hfm = plot(time_common,nanmean(rm(squeeze(p_fcst(idx3,lt,:)),3,1)*f,2),'Color',[235, 131, 52]/255,'LineWidth',1.5,'LineStyle','-')
  hg  = plot(time_common,lg_corrected2,'k')
  hline(2*nanstd(lg_corrected2),'k--')
  set(gca,'XLim',xlim,'YLim',ylim)
  box on

  xlim = [datetime(2013,4,1) datetime(2013,4,30)];
  subplot(1,3,3)
  title('No gain forecasted/Gain observed')
  hold on
  [dummy,idx3] = intersect(time_p_all,time_common);
  hf  = plot(time_common,rm(squeeze(p_fcst(idx3,lt,:)),3,1)*f,'Color',[1 .8 .8])
  hfm = plot(time_common,nanmean(rm(squeeze(p_fcst(idx3,lt,:)),3,1)*f,2),'Color',[235, 131, 52]/255,'LineWidth',1.5,'LineStyle','-')

  hg  = plot(time_common,lg_corrected2,'k')
  hline(2*nanstd(lg_corrected2),'k--')
  set(gca,'XLim',xlim,'YLim',ylim)
  box on

  tightfig

  set(gcf,'PaperPositionMode','auto');
  fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/gain_forecast'];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  save2pdf(['' fileo '.pdf'])
  return

end






if plot6 ==1
  % -- Contingency table ---------------------------------
  % -- choose event: 2 sigma
  gain_th = 2%2*nanstd(lg_corrected2);

  lt = 1;

  fcst  = rm(squeeze(p_fcst_mm(idx3,lt,:)),3,1);
  fcstm = nanmean(rm(squeeze(p_fcst_mm(idx3,lt,:)),3,1),2);

  % -- probability to exceed certain precip level:
  prob1   = nansum([fcst>1],2)/30; % >1 mm/day
  prob3   = nansum([fcst>3],2)/30; % >3 mm/day
  prob5   = nansum([fcst>5],2)/30; % >5 mm/day
  prob10  = nansum([fcst>10],2)/30; % >10 mm/day
  % idx     = find(prob1>0.9);
  % idx     = find(prob3>0.5);
  idx     = find(prob5>0.1);
  % idx     = find(prob10>0.1);
  % idx     = find(fcstm>2);
  sum(lg_corrected2>0)
  ['no of events = ' num2str(sum(lg_corrected2>gain_th))]
  ['no of forecasted events = ' num2str(length(idx))]

  TP = zeros(length(idx),1);
  FP = zeros(length(idx),1);
  for i = 1:length(idx)
    if lg_corrected2(idx(i)-(lt-1))>gain_th || lg_corrected2(idx(i)+1-(lt-1))>gain_th %|| lg_corrected2(idx(i)+2-(lt-1))>gain_th || lg_corrected2(idx(i)+3-(lt-1))>gain_th
      TP(i) = 1;
    else
      FP(i) = 1;
    end
  end
  TP = sum(TP)/length(idx) % true positive
  FP = sum(FP)/length(idx) % false positive

  idx     = find(prob5<=0.1);
  FN = zeros(length(idx),1);
  TN = zeros(length(idx),1);
  for i = 1:length(idx)-3
    if lg_corrected2(idx(i))<gain_th && lg_corrected2(idx(i)+1)<gain_th && lg_corrected2(idx(i)+2)<gain_th && lg_corrected2(idx(i)+3)<gain_th
      TN(i) = 1;
    else
      FN(i) = 1;
    end
  end
  FN = sum(FN)/length(idx) % false negative
  TN = sum(TN)/length(idx) % true negative

  return
end





if plot7 ==1
  % -- Contingency plot ---------------------------------
  close all


  ylim = [0 100];

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 27 7])

  for th = [0 1 2]
    ['th = ' num2str(th)]
    % -- choose threshold in sigma:
    gain_th = th * nanstd(lg_corrected2); % gain threshold
    for lt = 1:5 % loop over lead time
      ['  lt = ' num2str(lt)]
      fcst  = rm(squeeze(p_fcst_mm(idx3,lt,:)),3,1); % all forecasts
      fcstm = nanmean(rm(squeeze(p_fcst_mm(idx3,lt,:)),3,1),2); % ensemble mean forecast
      pj = 0;
      for p = [1 3 5 10] % loop over different strength precipitation forecasts
        pj = pj+1;
        % -- probability to exceed certain precipitation level:
        prob    = nansum([fcst>p],2)/30; % >[1 3 5 10] mm/day
        idx     = find(prob>0); % find time when precip probability is larger than zero
        TP = zeros(length(idx),1);
        FP = zeros(length(idx),1);
        % -- check if there is actually a gain happening up to 2 days after the forecasted precip event
        for i = 1:length(idx)
          if lg_corrected2(idx(i)-(lt-1))>gain_th || lg_corrected2(idx(i)+1-(lt-1))>gain_th || lg_corrected2(idx(i)+2-(lt-1))>gain_th %|| lg_corrected2(idx(i)+3-(lt-1))>gain_th
            TP(i) = 1;
          else
            FP(i) = 1;
          end
        end
        TP = sum(TP)/length(idx); % true positive
        FP = sum(FP)/length(idx); % false positive
        TP_all(lt,pj) = TP*100;
        FP_all(lt,pj) = FP*100;
      end
    end
    % -- climatology:
    clear('tmp')
    xx=3;
    for i = 1:length(lg_corrected2)-xx
      if lg_corrected2(i)>gain_th || lg_corrected2(i+1)>gain_th || lg_corrected2(i+2)>gain_th
        tmp(i)=1;
      else
        tmp(i)=0;
      end
    end
    tmp0 = sum(tmp)/length(lg_corrected2)
    % -- plot:
    subplot(1,3,th+1)
    hold on
    title(['Gains > ' num2str(th) '\sigma'])
    plot(TP_all,'.-','MarkerSize',12)
    hline(tmp0*100,'k--')
    set(gca,'YLim',ylim)
    box on
    xlabel('Lead time (days)')
    if th == 0
      text(1.2,tmp0*100-4,['Climatology'],'FontSize',10)
      ylabel(['Fraction of correctly',char(10),'forecasted events (%)'])
      legend('P(pr > 0 mm/day) > 0','P(pr > 3 mm/day) > 0','P(pr > 5 mm/day) > 0','P(pr > 10 mm/day) > 0','location','southwest')
      legend boxoff
    end
  end

  tightfig

  % return
  set(gcf,'PaperPositionMode','auto');
  fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/gain_forecast_contingency_plot'];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  save2pdf(['' fileo '.pdf'])
  return
end







% --
close all

xlim = [-5e3 5e3];
ylim = [-1 25];

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 15 20])

subplot(2,1,1)
hold on
title('Loss-gain vs precip (2003-2017)')
% plot(lg_corrected(idx1),p_all(idx3),'o','Color',[.7 .7 .7])
% plot(lg_corrected2(idx1),p_all(idx3),'bo')
% [f0,x0,u] = ksdensity(lg_corrected(idx1));
% h1 = line(x0,f0*1e4,'Color',[.5 .5 .5], 'LineWidth',2);
% [f0,x0,u] = ksdensity(lg_corrected2(idx1));
% h1 = line(x0,f0*1e4,'Color',[0 0 1], 'LineWidth',2);
h1 = plot(lg_corrected(idx1),rm(p_all(idx3),3),'o','Color',[.7 .7 .7])
h2 = plot(lg_corrected2(idx1),rm(p_all(idx3),3),'bo')
% [f0,x0,u] = ksdensity(rm(lg_corrected(idx1),3));
% h1 = line(x0,f0*1e4,'Color',[.5 .5 .5], 'LineWidth',2);
% [f0,x0,u] = ksdensity(rm(lg_corrected2(idx1),3));
% h1 = line(x0,f0*1e4,'Color',[0 0 1], 'LineWidth',2);
set(gca,'YLim',ylim,'XLim',xlim)
vline(0,'k')
box on
xlabel('Loss-gain at Imperial (cfs)')
ylabel('3-day precipitation (mm)')
legend([h1 h2],'Uncorrected','Backwater effect removed','location','northwest')
legend boxoff

subplot(2,1,2)
hold on
title('Binned by precip (0mm, 0-5mm, 5-10mm, 10-15mm)')
x1  = lg_corrected(idx1);
x2  = lg_corrected2(idx1);
% y   = p_all(idx3);
% x1  = rm(lg_corrected(idx1),3);
% x2  = rm(lg_corrected2(idx1),3);
y   = rm(p_all(idx3),3);
pp1 = 5;
pp2 = 95;
pp3 = 25;
pp4 = 75;

i = 0;
os = 0;
id = find(y==i);
length(id)
plot([min(x1(id)) max(x1(id))],[i+os i+os],'Color',[.5 .5 .5],'LineStyle',':','LineWidth',1)
plot([min(x2(id)) max(x2(id))],[i+os+1 i+os+1],'Color',[0 0 1],'LineStyle',':','LineWidth',1)
plot([prctile(x1(id),pp1) prctile(x1(id),pp2)],[i+os i+os],'Color',[.5 .5 .5],'LineWidth',2)
plot([prctile(x2(id),pp1) prctile(x2(id),pp2)],[i+os+1 i+os+1],'Color',[0 0 1],'LineWidth',2)
plot([prctile(x1(id),pp3) prctile(x1(id),pp4)],[i+os i+os],'Color',[.5 .5 .5],'LineWidth',5)
plot([prctile(x2(id),pp3) prctile(x2(id),pp4)],[i+os+1 i+os+1],'Color',[0 0 1],'LineWidth',5)

incr = 5;
os = floor(incr/2)+1;
% for i = 0:incr:25 % bins for precip amount
for i = 0:incr:10 % bins for precip amount
  id = find(y>i & y<=(i+incr));
  length(id)
  plot([min(x1(id)) max(x1(id))],[i+os i+os],'Color',[.5 .5 .5],'LineStyle',':','LineWidth',1)
  plot([min(x2(id)) max(x2(id))],[i+os+1 i+os+1],'Color',[0 0 1],'LineStyle',':','LineWidth',1)
  plot([prctile(x1(id),pp1) prctile(x1(id),pp2)],[i+os i+os],'Color',[.5 .5 .5],'LineWidth',2)
  plot([prctile(x2(id),pp1) prctile(x2(id),pp2)],[i+os+1 i+os+1],'Color',[0 0 1],'LineWidth',2)
  plot([prctile(x1(id),pp3) prctile(x1(id),pp4)],[i+os i+os],'Color',[.5 .5 .5],'LineWidth',5)
  plot([prctile(x2(id),pp3) prctile(x2(id),pp4)],[i+os+1 i+os+1],'Color',[0 0 1],'LineWidth',5)
end
set(gca,'YLim',[-1 ylim(2)],'XLim',xlim)
vline(0,'k')
box on
xlabel('Loss-gain at Imperial (cfs)')
ylabel('3-day precipitation (mm)')

set(gcf,'PaperPositionMode','auto');
fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/backwater_effect2'];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
save2pdf(['' fileo '.pdf'])

return




% figure
% hold on
% plot(release,lg_corrected,'rx')
% plot(release(first_Th:7:end),lg_corrected(first_Th:7:end),'o')

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 15 20])

subplot(2,1,1)
hold on
title('Release from Parker Dam')
h = boxplot(drop_wd','plotstyle','traditional','whisker',1.96,'Colors','b','Symbol','.')
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
hline(0,'k')
set(gca,'XTickLabel',{'We->Th';'Th->Fr';'Fr->Sa';'Sa->So';'So->Mo';'Mo->Tu';'Tu->We'})
ylabel('Day-to-day change in Parker release (cfs)')

xlim = [-10e3 10e3];
ylim = [-5e3 5e3];

subplot(2,1,2)
hold on
title('Backwater effect')
d1 = 1; % We->Th
h1 = plot(drop_wd(d1,:),lg_corrected_wd(d1,:),'bo')
d2 = 2; % Th->Fr
h2 = plot(drop_wd(d2,:),lg_corrected_wd(d2,:),'ro')
for d = [3:7]
  h3 = scatter(drop_wd(d,:),lg_corrected_wd(d,:),'MarkerEdgeColor',[.6 .6 .6])%,'MarkerFaceColor',[.7 .7 .7])
  % h3.MarkerFaceAlpha = .2;
  h3.MarkerEdgeAlpha = .5;
end
set(gca,'XLim',xlim,'YLim',ylim)
ax = ancestor(h3, 'axes')
ax.XAxis.Exponent = 0
xtickformat('%.0f')
hline(0,'k')
vline(0,'k')
% -- correlation of day-today drop and loss-gain, for each weekday separatelt
for d = 1:7
  idx1 = find(~isnan(lg_corrected_wd(d,:)));
  corr_wd(d) = corr(drop_wd(d,idx1)',lg_corrected_wd(d,idx1)');
end
% --
legend([h1 h2 h3],...
['We->Th (r=' num2str(corr_wd(d1),'%.2f') ')'],...
['Th->Fr (r=' num2str(corr_wd(d2),'%.2f') ')'],...
['all other weekdays (r=' num2str(nanmean(corr_wd([1 4:7])),'%.2f') ')'],...
'location','SouthWest')
legend boxoff
xlabel('Day-to-day change in Parker release (cfs)')
ylabel('Gains at Imperial 3 days later (cfs)')
box on

% return
set(gcf,'PaperPositionMode','auto');
fileo = ['~/Dropbox/proposals/reclamation_SnT/yuma_shortterm_forecast/fig/backwater_effect1'];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
save2pdf(['' fileo '.pdf'])
% return



% figure
%
% subplot(2,1,1)
% hold on
% idx1 = find(~isnan(lg_corrected_drop_Th));
% plot(drop_Th(idx1),lg_corrected_drop_Th(idx1),'bo')
% idx2 = find(~isnan(lg_corrected_rebound_Fr));
% plot(rebound_Fr(idx2),lg_corrected_rebound_Fr(idx2),'ro')
% [b] = regress(lg_corrected_drop_Th(idx1),[ones(size(drop_Th(idx1))) drop_Th(idx1)])
% plot(drop_Th,b(1)+b(2)*drop_Th,'k','LineWidth',2)
% [b] = regress(lg_corrected_rebound_Fr(idx2),[ones(size(rebound_Fr(idx2))) rebound_Fr(idx2)])
% plot(rebound_Fr,b(1)+b(2)*rebound_Fr,'Color',[.5 .5 .5],'LineWidth',2)
% plot([xlim(2) xlim(1)],[ylim(1) ylim(2)],'k')
% set(gca,'XLim',xlim,'YLim',ylim)
% hline(0,'k')
% vline(0,'k')
% xlabel('Thursday-Wednesday change in flow (cfs)')
% ylabel('Loss-gain (cfs)')
% box on
% corr(drop_Th(idx1),lg_corrected_drop_Th(idx1))
%
% subplot(2,1,2)
% hold on
% idx = find(~isnan(lg_corrected));
% plot(drop,lg_corrected,'o')
% [b] = regress(lg_corrected,[ones(size(drop)) drop])
% plot(drop,b(1)+b(2)*drop,'k','LineWidth',2)
% plot([xlim(2) xlim(1)],[ylim(1) ylim(2)],'k')
% set(gca,'XLim',xlim,'YLim',ylim)
% hline(0,'k')
% vline(0,'k')
% xlabel('Change in flow from previous day (cfs)')
% ylabel('Loss-gain (cfs)')
% box on
% corr(drop(idx),lg_corrected(idx))


figure
subplot(2,1,1)
hold on
plot(lg_corrected,'k')
plot(lg_predicted,'r')
subplot(2,1,2)
hold on
plot(lg_corrected(3:7:end),lg_predicted(3:7:end),'o')




% return



figure
hold on
plot(dt,release,'b')
plot(dt,arrival,'r')
plot(dt,loss_en_route2(2:end-2),'m')
plot(dt,lg,'Color',[.4 .4 .4])
plot(dt,lg_corrected,'k','LineWidth',1.5)
hline(0,'k')
% set(gca,'XLim',[datetime(20160301,'ConvertFrom','yyyymmdd') datetime(20160631,'ConvertFrom','yyyymmdd')],'YLim',[-7e3 14e3])
% legend('Release Parker','Arrival Imperial','Loss-gain Parker-Imperial','loss CRIT+PVID2','loss CRIT+PVID1+PVID2','Loss-gain Parker-Imperial after accounting for Blythe area losses','location','southwest')
legend('Release Parker','Arrival Imperial (3 days lagged)','loss CRIT+PVID','Loss-gain Parker-Imperial','Loss-gain Parker-Imperial after accounting for CRIT+PVID losses','location','southwest')
box on
xlabel('Time')
ylabel('Flow (cfs)')

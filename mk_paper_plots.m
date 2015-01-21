close all, clear all
addpath 'C:\Users\Collin Rennie\Documents\MATLAB\work\general'
%addpath 'E:\matlab_code_uottawa\work\general'
addpath 'tools/'
addpath 'good_data/'
set(0,'defaulttextinterpreter','tex')

dn = 'C:\Users\Stephanie\Documents\my publications\uncertainty_paper\figures\';




%% Figure 1 : stage time series for the two sites
doneF = 1;
if doneF == 0
    % from file 05MJ001.xlsx
    dateW = ['03/06/2008 00:00';
        '03/06/2008 00:00';
        '03/06/2008 01:00';
        '03/06/2008 02:00';
        '03/06/2008 03:00';
        '03/06/2008 03:00';
        '03/06/2008 04:00';
        '03/06/2008 04:55';
        '03/06/2008 05:00';
        '03/06/2008 06:00';
        '03/06/2008 06:00';
        '03/06/2008 07:00';
        '03/06/2008 08:00';
        '03/06/2008 08:50';
        '03/06/2008 09:00';
        '03/06/2008 09:00';
        '03/06/2008 10:00';
        '03/06/2008 10:05';
        '03/06/2008 11:00';
        '03/06/2008 12:00';
        '03/06/2008 12:00';
        '03/06/2008 12:50';
        '03/06/2008 13:00';
        '03/06/2008 14:00';
        '03/06/2008 15:00';
        '03/06/2008 15:00';
        '03/06/2008 15:55';
        '03/06/2008 16:00';
        '03/06/2008 17:00';
        '03/06/2008 18:00';
        '03/06/2008 18:00';
        '03/06/2008 19:00';
        '03/06/2008 19:20';
        '03/06/2008 20:00';
        '03/06/2008 21:00';
        '03/06/2008 21:00';
        '03/06/2008 22:00';
        '03/06/2008 23:00';
        '03/06/2008 23:10';
        '04/06/2008 00:00';
        '04/06/2008 00:00';
        '04/06/2008 01:00';
        '04/06/2008 01:10';
        '04/06/2008 02:00';
        '04/06/2008 03:00';
        '04/06/2008 03:00';
        '04/06/2008 04:00';
        '04/06/2008 05:00';
        '04/06/2008 06:00';
        '04/06/2008 06:00';
        '04/06/2008 07:00';
        '04/06/2008 07:55';
        '04/06/2008 08:00';
        '04/06/2008 09:00';
        '04/06/2008 09:00';
        '04/06/2008 10:00';
        '04/06/2008 10:25';
        '04/06/2008 11:00';
        '04/06/2008 12:00';
        '04/06/2008 12:00';
        '04/06/2008 13:00';
        '04/06/2008 13:25';
        '04/06/2008 14:00';
        '04/06/2008 15:00';
        '04/06/2008 15:00';
        '04/06/2008 16:00';
        '04/06/2008 17:00';
        '04/06/2008 18:00';
        '04/06/2008 18:00';
        '04/06/2008 18:15';
        '04/06/2008 19:00';
        '04/06/2008 20:00';
        '04/06/2008 21:00';
        '04/06/2008 21:00';
        '04/06/2008 21:30';
        '04/06/2008 22:00';
        '04/06/2008 23:00';
        '05/06/2008 00:00';
        
        ];
    
    
    clear ii
    dayW = NaN*ones(1, length(dateW));
    monthW = NaN*ones(1, length(dateW));
    yearW = NaN*ones(1, length(dateW));
    hourW = NaN*ones(1, length(dateW));
    minW = NaN*ones(1, length(dateW));
    for ii = 1:length(dateW)
        dayW(ii) = str2num(dateW(ii, 1:2));
        monthW(ii) = str2num(dateW(ii,4:5));
        yearW(ii) = str2num(dateW(ii,7:10));
        hourW(ii) = str2num(dateW(ii,12:13));
        minW(ii) = str2num(dateW(ii,15:16));
        
    end
    secW = zeros(size(minW));
    WDateNum = datenum(yearW, monthW, dayW, hourW, minW, secW);
    
    
    stage = [231.2248476;
        231.2248476;
        231.2248523;
        231.2248571;
        231.2248618;
        231.2248619;
        231.2248666;
        231.2238709;
        231.2238713;
        231.2238761;
        231.2238761;
        231.2238808;
        231.2238856;
        231.2228895;
        231.2228903;
        231.2228903;
        231.2228951;
        231.2238955;
        231.2238998;
        231.2239046;
        231.2239046;
        231.2249085;
        231.2249093;
        231.2249141;
        231.2249188;
        231.2249188;
        231.2259232;
        231.2259236;
        231.2259283;
        231.2259331;
        231.2259331;
        231.2259378;
        231.2269394;
        231.2269426;
        231.2269473;
        231.2269473;
        231.2269521;
        231.2269568;
        231.2279576;
        231.2279616;
        231.2279616;
        231.2279663;
        231.2289671;
        231.228971;
        231.2289758;
        231.2289758;
        231.2289805;
        231.2289853;
        231.22899;
        231.2289901;
        231.2289948;
        231.2279991;
        231.2279995;
        231.2279968;
        231.2279968;
        231.2279932;
        231.2269917;
        231.2269896;
        231.2269859;
        231.2269859;
        231.2269823;
        231.2259808;
        231.2259787;
        231.2259751;
        231.2259751;
        231.2259715;
        231.2259679;
        231.2259643;
        231.2259642;
        231.2269633;
        231.2269606;
        231.226957;
        231.2269534;
        231.2269534;
        231.2279516;
        231.2279498;
        231.2279462;
        231.2279426;
        ];
    
    WDateNum_05MJ001 = WDateNum;
    stage_05MJ001 = stage;
    clear *W stage
    dateW = [
        '03/06/2008 00:00';
        '03/06/2008 01:00';
        '03/06/2008 01:05';
        '03/06/2008 02:00';
        '03/06/2008 02:30';
        '03/06/2008 03:00';
        '03/06/2008 03:30';
        '03/06/2008 04:00';
        '03/06/2008 05:00';
        '03/06/2008 05:55';
        '03/06/2008 06:00';
        '03/06/2008 06:15';
        '03/06/2008 07:00';
        '03/06/2008 08:00';
        '03/06/2008 08:15';
        '03/06/2008 09:00';
        '03/06/2008 10:00';
        '03/06/2008 11:00';
        '03/06/2008 11:15';
        '03/06/2008 11:35';
        '03/06/2008 12:00';
        '03/06/2008 12:05';
        '03/06/2008 13:00';
        '03/06/2008 14:00';
        '03/06/2008 14:25';
        '03/06/2008 15:00';
        '03/06/2008 15:30';
        '03/06/2008 15:35';
        '03/06/2008 16:00';
        '03/06/2008 17:00';
        '03/06/2008 18:00';
        '03/06/2008 18:25';
        '03/06/2008 19:00';
        '03/06/2008 20:00';
        '03/06/2008 20:15';
        '03/06/2008 21:00';
        '03/06/2008 22:00';
        '03/06/2008 23:00';
        '03/06/2008 23:30';
        '03/06/2008 23:45';
        '04/06/2008 00:00';
        '04/06/2008 00:10';
        '04/06/2008 01:00';
        '04/06/2008 02:00';
        '04/06/2008 02:40';
        '04/06/2008 03:00';
        '04/06/2008 04:00';
        '04/06/2008 05:00';
        '04/06/2008 05:15';
        '04/06/2008 06:00';
        '04/06/2008 07:00';
        '04/06/2008 08:00';
        '04/06/2008 08:40';
        '04/06/2008 09:00';
        '04/06/2008 10:00';
        '04/06/2008 11:00';
        '04/06/2008 11:50';
        '04/06/2008 12:00';
        '04/06/2008 13:00';
        '04/06/2008 13:55';
        '04/06/2008 14:00';
        '04/06/2008 14:20';
        '04/06/2008 15:00';
        '04/06/2008 15:55';
        '04/06/2008 16:00';
        '04/06/2008 17:00';
        '04/06/2008 18:00';
        '04/06/2008 19:00';
        '04/06/2008 20:00';
        '04/06/2008 20:30';
        '04/06/2008 20:45';
        '04/06/2008 21:00';
        '04/06/2008 21:30';
        '04/06/2008 22:00';
        '04/06/2008 23:00';
        '04/06/2008 23:50';
        '05/06/2008 00:00';
        ];
    
    
    clear ii
    for ii = 1:length(dateW)
        dayW(ii) = str2num(dateW(ii, 1:2));
        monthW(ii) = str2num(dateW(ii,4:5));
        yearW(ii) = str2num(dateW(ii,7:10));
        hourW(ii) = str2num(dateW(ii,12:13));
        minW(ii) = str2num(dateW(ii,15:16));
        
    end
    secW = zeros(size(minW));
    WDateNum = datenum(yearW, monthW, dayW, hourW, minW, secW);
    
    
    stage = [268.7444477;
        268.7464029;
        NaN;
        268.7393581;
        268.7343357;
        268.7373132;
        NaN;
        268.7352684;
        268.7312236;
        NaN;
        268.7351788;
        NaN;
        268.7331339;
        268.7310891;
        NaN;
        268.7240443;
        268.733;
        268.7330132;
        268.7170165;
        NaN;
        268.7180265;
        268.7240276;
        268.7160398;
        268.7120531;
        268.7070586;
        268.7070664;
        NaN;
        268.7040741;
        268.7120797;
        268.707093;
        268.7061063;
        268.7121118;
        268.7091195;
        268.7001328;
        NaN;
        268.6981461;
        268.6951594;
        268.6961727;
        268.6901793;
        NaN;
        268.691186;
        NaN;
        268.6901993;
        268.6872126;
        NaN;
        268.6912259;
        268.6862391;
        268.6862524;
        NaN;
        268.6842657;
        268.680279;
        268.6762923;
        NaN;
        268.6803056;
        268.6763189;
        268.6733322;
        268.6693432;
        268.6743455;
        268.6723588;
        NaN;
        268.672372;
        268.6773765;
        268.6693853;
        NaN;
        268.6803986;
        268.6734119;
        268.6724252;
        268.6654385;
        268.6694518;
        268.6644584;
        NaN;
        268.6604651;
        268.6624717;
        268.6594784;
        268.6564916;
        268.6525027;
        268.6525049;
        ];
    
    WDateNum_05PH003 = WDateNum;
    stage_05PH003 = stage;
    
    figure(1), clf(1)
    subplot 211
    plot(WDateNum_05MJ001, stage_05MJ001, 'k', 'linewidth', 1.5)
    hold on
    datetick('x', 'keeplimits')
    ylabel('Stage [m]')
    %xlabel('Date in 2008 [mm/dd]')
    title('Site 05MJ001')
    transStartTimes = [datenum(2008, 06, 04, 11, 59, 22); datenum(2008, 06, 04, 12, 05, 36); datenum(2008, 06, 04, 12, 12, 50); datenum(2008, 06, 04, 12, 19, 23);...
                       datenum(2008, 06, 04, 9, 10, 48); datenum(2008, 06, 04, 9, 15, 41); datenum(2008, 06, 04, 9, 20, 34); datenum(2008, 06, 04, 9, 43, 23); ...
                       datenum(2008, 06, 04, 13, 13, 27); datenum(2008, 06, 04, 13, 23, 09); datenum(2008, 06, 04, 13,29, 30); datenum(2008, 06, 04, 13, 35, 25); ...
                       datenum(2008, 06, 04, 10, 19, 41); datenum(2008, 06, 04, 10, 25, 57); datenum(2008, 06, 04, 10, 31, 12); datenum(2008, 06, 04, 10, 37, 14);]
    trans = ['A-T1'; 'A-T2'; 'A-T3'; 'A-T4'; 'B-T1'; 'B-T2'; 'B-T3'; 'B-T4'; 'C-T1'; 'C-T2';'C-T3';'C-T4'; 'D-T1'; 'D-T2';'D-T3';'D-T4']
   ylims = [231.22 231.228];
    h1(1) = plot([transStartTimes(1) transStartTimes(1)], ylims, '-k', 'color', [.5 .5 .5])
    h1(2) = plot([transStartTimes(2) transStartTimes(2)], ylims, '-k', 'color', [.5 .5 .5])
    h1(3) = plot([transStartTimes(3) transStartTimes(3)], ylims, '-k', 'color', [.5 .5 .5])
    h1(4) =  plot([transStartTimes(4) transStartTimes(4)], ylims, '-k', 'color', [.5 .5 .5])
    text(transStartTimes(1)*.99999998, 231.229, 'A', 'fontsize', 12)
    h1(5) =  plot([transStartTimes(5) transStartTimes(5)], ylims, 'k', 'color', [.5 .5 .5])
    h1(6) =  plot([transStartTimes(6) transStartTimes(6)], ylims, 'k', 'color', [.5 .5 .5])
    h1(7) =  plot([transStartTimes(7) transStartTimes(7)], ylims, 'k', 'color', [.5 .5 .5])
    h1(8) = plot([transStartTimes(8) transStartTimes(8)], ylims, 'k', 'color', [.5 .5 .5])
    text(transStartTimes(5)*.99999998, 231.229, 'B', 'fontsize', 12)
    h1(9) =  plot([transStartTimes(9) transStartTimes(9)], ylims, 'k', 'color', [.5 .5 .5])
    h1(10) =  plot([transStartTimes(10) transStartTimes(10)], ylims, 'k', 'color', [.5 .5 .5])
    h1(11) =  plot([transStartTimes(11) transStartTimes(11)], ylims, 'k', 'color', [.5 .5 .5])
    h1(12) = plot([transStartTimes(12) transStartTimes(12)], ylims, 'k', 'color', [.5 .5 .5])
    text(transStartTimes(9), 231.229, 'C', 'fontsize', 12)
    h1(13) =  plot([transStartTimes(13) transStartTimes(13)], ylims, 'k', 'color', [.5 .5 .5])
    h1(14) =  plot([transStartTimes(14) transStartTimes(14)], ylims, 'k', 'color', [.5 .5 .5])
    h1(15) =  plot([transStartTimes(15) transStartTimes(15)], ylims, 'k', 'color', [.5 .5 .5])
    h1(16) = plot([transStartTimes(16) transStartTimes(16)], ylims, 'k', 'color', [.5 .5 .5])
    text(transStartTimes(13), 231.229, 'D', 'fontsize', 12)
    plot(WDateNum_05MJ001, stage_05MJ001, 'k', 'linewidth', 1.5)
     set(gca,'XGrid','off', 'YGrid', 'on')
    
    subplot 212
    plot(WDateNum_05PH003(stage_05PH003> 0), stage_05PH003(stage_05PH003> 0), 'k', 'linewidth', 1.5)
    hold on
    datetick('x', 'keeplimits')
    ylabel('Stage [m]')
    xlabel('Date in 2008 [mm/dd]')
    title('Site 05PH003')
    
    transStartTimes = [datenum(2008, 06, 03, 13, 15, 12); datenum(2008, 06, 03, 13, 19, 15); datenum(2008, 06, 03, 13, 22, 18); datenum(2008, 06, 03, 13, 25, 22);...
                       datenum(2008, 06, 03, 12, 14, 50); datenum(2008, 06, 03, 12, 18, 0); datenum(2008, 06, 03, 12, 22, 28); datenum(2008, 06, 03, 12, 25, 07); ...
                       datenum(2008, 06, 03, 14, 10, 44); datenum(2008, 06, 03, 14, 14, 27); datenum(2008, 06, 03, 14, 16, 35); datenum(2008, 06, 03, 14, 18, 24); ...
                       datenum(2008, 06, 03, 16, 57, 19); datenum(2008, 06, 03, 17, 0, 48); datenum(2008, 06, 03, 17, 04, 47); datenum(2008, 06, 03, 17,08,25);]
   ylims = [268.65 268.76]; 
                   h1(1) = plot([transStartTimes(1) transStartTimes(1)], ylims, '-k', 'color', [.5 .5 .5])
    h1(2) = plot([transStartTimes(2) transStartTimes(2)], ylims, '-k', 'color', [.5 .5 .5])
    h1(3) = plot([transStartTimes(3) transStartTimes(3)], ylims, '-k', 'color', [.5 .5 .5])
    h1(4) =  plot([transStartTimes(4) transStartTimes(4)], ylims, '-k', 'color', [.5 .5 .5])
    text(transStartTimes(1)*.99999999, 268.77, 'A', 'fontsize', 12)
    h1(5) =  plot([transStartTimes(5) transStartTimes(5)], ylims, 'k', 'color', [.5 .5 .5])
    h1(6) =  plot([transStartTimes(6) transStartTimes(6)], ylims, 'k', 'color', [.5 .5 .5])
    h1(7) =  plot([transStartTimes(7) transStartTimes(7)], ylims, 'k', 'color', [.5 .5 .5])
    h1(8) = plot([transStartTimes(8) transStartTimes(8)], ylims, 'k', 'color', [.5 .5 .5])
    text(transStartTimes(5)*.99999998,  268.77, 'B', 'fontsize', 12)
    h1(9) =  plot([transStartTimes(9) transStartTimes(9)], ylims, 'k', 'color', [.5 .5 .5])
    h1(10) =  plot([transStartTimes(10) transStartTimes(10)], ylims, 'k', 'color', [.5 .5 .5])
    h1(11) =  plot([transStartTimes(11) transStartTimes(11)], ylims, 'k', 'color', [.5 .5 .5])
    h1(12) = plot([transStartTimes(12) transStartTimes(12)], ylims, 'k', 'color', [.5 .5 .5])
    text(transStartTimes(9), 268.77, 'C', 'fontsize', 12)
    h1(13) =  plot([transStartTimes(13) transStartTimes(13)], ylims, 'k', 'color', [.5 .5 .5])
    h1(14) =  plot([transStartTimes(14) transStartTimes(14)], ylims, 'k', 'color', [.5 .5 .5])
    h1(15) =  plot([transStartTimes(15) transStartTimes(15)], ylims, 'k', 'color', [.5 .5 .5])
    h1(16) = plot([transStartTimes(16) transStartTimes(16)], ylims, 'k', 'color', [.5 .5 .5])
    text(transStartTimes(13)*.99999999, 268.77, 'D', 'fontsize', 12) 
    plot(WDateNum_05PH003, stage_05PH003, 'k', 'linewidth', 1.5)
    ylim([ 268.65 268.78])
   set(gca,'XGrid','off', 'YGrid', 'on')
    filename = 'stage_sites_05MJ001_05PH003'
    print ('-depsc', [dn filename])
    
end

 
%low uncertainty and high uncertainty
doneF = 1
if doneF == 0
    figure(1), clf(1), hold on
    % low uncertainty: 05OJ010 (Red River near Erwood, May 4, 2004)
    load good_data/05OJ010_20040504.mat
    Qbreakdown{1:4}
    
    
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    figure(1)
    subplot 211
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    % xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
    
    set(gca, 'box', 'on')
    %axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    set(gca, 'linewidth', 1.5)
    ylabel('Depth [m]', 'fontsize', 16)
    text(-0.1, -5, 'Depth [m]', 'fontsize', 14, 'rotation', 90, 'units', 'normalized')
    text(0, 1.1, ['(a) 05OJ010, Q_{meas}/Q_{tot} = 75%, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized', 'fontsize', 12)
    text(1.04, 1.12, '[m/s]', 'units', 'normalized', 'fontsize', 12)
    
    oldPos = get(gca, 'position')
    % set(gca, 'linewidth', 1.5, 'position', [oldPos(1) 1.2*oldPos(2) oldPos(3:4)])
    

    % high uncertainty
    % Canadian St. Marys Canal near Spring Coulee,  high uncertainty ~5%
    load good_data/05AE026_20100719_O_Connor_AQ1.mat
    Qbreakdown{1:4}
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    figure(1)
    subplot 212
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel(['Distance from left bank (ref: bottom) [m]'], 'fontsize', 14)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    
    %axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    oldPos = get(gca, 'position')
    set(gca, 'linewidth', 1.5, 'position', [oldPos(1) 1.25*oldPos(2) oldPos(3:4)])
    text(0, 1.1, ['(b) 05AE026, Q_{meas}/Q_{tot} = 32%, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized', 'fontsize', 12)
    text(1.04, 1.12, '[m/s]', 'units', 'normalized', 'fontsize', 12)
    
   
    
    filename = 'low_and_high'
    print ('-depsc', [dn filename])
end


% Gatineau River heading including one and two cycle errors
% what I call self dependent error
doneF = 0
if doneF == 0
    load good_data_old/gatineau_20031203_gps_vary_heading_only.mat
    
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    ymax = 0;
    xmin = 10000;
    xmax = 0;
    
    for mm = 1:2
        
        length_reGGA =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel_reGGA(1,:).^2 + (parameters.data{1}.btVel_reGGA(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx_reGGA.^2 + parameters.data{1}.wVely_reGGA.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_reGGA + parameters.data{mm}.lDist;
        
        figure(1)
        subplot (2,1,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Distance from left bank (ref: bottom) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
        % title('Velocity magnitude [m/s]', 'fontsize', 10)
       % axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        caxis([0 2])
        %x = caxis;
        %x = round(100*x)/100;
        
        
        set(gca, 'linewidth', 1.5)
        % set(hc, 'YTick',x)
        
        if mm == 1
            filename = 'gatineau_20031203_xsec_reGGA_v2'
            print ('-depsc', [dn filename])
            print ('-dtiff', [dn filename])
        end
        
        
        % plot pdfs
        figure(2)
        subplot (2,2,mm)
        nr = length(parameters.Q_reGGA{mm,8})
        
        dischargeGGA = parameters.Q_reGGA{mm,8}; % only if only heading is of interest for heading
        
        [n, xout] = hist(parameters.Q_reGGA{mm,8},20); % for heading variation
        pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reGGA)*1.1 > ymax;
            ymax = max(pdfQtest_reGGA)*1.1;
        end
        
        
        % hl(1) = plot([parameters.meanQ_reGGA(mm,8) parameters.meanQ_reGGA(mm,8)], [0 max(pdfQtest_reGGA)*1.1], '-.k');
        hl(1) = plot([parameters.meanQ_reGGA(mm,8) parameters.meanQ_reGGA(mm,8)], [0 1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reGGA, 'b.', 'markersize', 15);
        xplot = linspace(0.99*mean(dischargeGGA), 1.01*mean(dischargeGGA));%linspace(min(discharge), max(discharge));
        % xplot = linspace(min(dischargeGGA), max(dischargeGGA));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeGGA), std(dischargeGGA)), '--r');
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        set(hl, 'linewidth', 1.5)
        if mm == 1
            legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
            legend boxoff
        end
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(dischargeGGA)/mean(dischargeGGA),3)])
        % title('Q GGA')
        
        
        
    end
    xmin = 780; xmax = 820;
    figure(2)
    subplot 221,  xlim([xmin xmax]), ylim([0 1]), text(0, 1.05,'(a) transect from left to right bank', 'units', 'normalized')
    subplot 222,  xlim([xmin xmax]), ylim([0 1]), text(0, 1.05,'(b) transect from right to left bank', 'units', 'normalized')
    %  subplot 223,  xlim([xmin xmax]), ylim([0 1]), text(0, 1.05,'(c)', 'units', 'normalized')
    %  subplot 224,  xlim([xmin xmax]), ylim([0 1]), text(0, 1.05,'(d)', 'units', 'normalized')
    
    filename = 'gatineau_gga_vary_heading'
    print ('-depsc', [dn filename])
    print ('-dtiff', [dn filename])
    
    
end

return
%%

% Gatineau River heading including one and two cycle errors
% what I call self dependent error
doneF = 1
if doneF == 0
    load good_data/gatineau_20031203_gps_vary_heading_only.mat
    Qbreakdown{1:4}
    
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    ymax = 0;
    xmin = 10000;
    xmax = 0;
    
    for mm = 1:4
        
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        %subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 10)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        x = caxis;
        x = round(100*x)/100;
        
        
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
        set(hc, 'YTick',x)
        
        if mm == 1
            filename = 'gatineau_20031203_xsec'
            print ('-depsc', [dn filename])
            print ('-dtiff', [dn filename])
        end
        
        
        % plot pdfs
        figure(2)
        subplot (2,2,mm)
        nr = length(parameters.Q_reGGA{mm,8})
        
        dischargeGGA = parameters.Q_reGGA{mm,8}; % only if only heading is of interest for heading
        
        [n, xout] = hist(parameters.Q_reGGA{mm,8},20); % for heading variation
        pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reGGA)*1.1 > ymax;
            ymax = max(pdfQtest_reGGA)*1.1;
        end
        
        
        hl(1) = plot([parameters.meanQ_reGGA(mm,8) parameters.meanQ_reGGA(mm,8)], [0 max(pdfQtest_reGGA)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reGGA, 'b.', 'markersize', 15);
        xplot = linspace(0.99*mean(dischargeGGA), 1.01*mean(dischargeGGA));%linspace(min(discharge), max(discharge));
        % xplot = linspace(min(dischargeGGA), max(dischargeGGA));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeGGA), std(dischargeGGA)), '--r');
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        set(hl, 'linewidth', 1.5)
        if mm == 1
            legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
            legend boxoff
        end
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(dischargeGGA)/mean(dischargeGGA),3)])
        title('Q GGA')
        
        dischargeVTG = parameters.Q_reVTG{mm,8}; % only if only heading is of interest for heading
        
        [n, xout] = hist(parameters.Q_reVTG{mm,8},20); % for heading variation
        pdfQtest_reVTG = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reVTG)*1.1 > ymax;
            ymax = max(pdfQtest_reVTG)*1.1;
        end
        
        
        
        figure(3)
        subplot (2,2,mm)
        hl(1) = plot([parameters.meanQ_reVTG(mm,8) parameters.meanQ_reVTG(mm,8)], [0 max(pdfQtest_reVTG)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reVTG, 'b.', 'markersize', 15);
        xplot = linspace(0.99*mean(dischargeVTG), 1.01*mean(dischargeVTG));% xplot = linspace(min(dischargeVTG), max(dischargeVTG));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeVTG), std(dischargeVTG)), '--r');
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        
        set(hl, 'linewidth', 1.5)
        if mm == 1
            legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
            legend boxoff
        end
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(dischargeVTG)/mean(dischargeVTG),3)])
        title('Q VTG')
        
    end
    
    figure(2)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_gga_vary_heading'
    print ('-depsc', [dn filename])
    print ('-dtiff', [dn filename])
    
    figure(3)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_vtg_vary_heading'
  %  print ('-depsc', [dn filename])
  %  print ('-dtiff', [dn filename])
    
    
end


% Gatineau River heading normal error
doneF = 1
if doneF == 0
    
    load good_data/gatineau_20031203.mat
    Qbreakdown{1:4}
    
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    ymax = 0;
    xmin = 10000;
    xmax = 0;
    
    for mm = 1:4
        
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        %subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 10)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        x = caxis;
        x = round(100*x)/100;
        
        
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
        set(hc, 'YTick',x)
        
        
        
        % plot pdfs
        figure(2)
        subplot (2,2,mm)
        nr = length(parameters.Q_reGGA{mm,8})
        
        dischargeGGA = parameters.Q_reGGA{mm,8}; % only if only heading is of interest for heading
        
        [n, xout] = hist(parameters.Q_reGGA{mm,8},20); % for heading variation
        pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reGGA)*1.1 > ymax;
            ymax = max(pdfQtest_reGGA)*1.1;
        end
        
        
        hl(1) = plot([parameters.meanQ_reGGA(mm,8) parameters.meanQ_reGGA(mm,8)], [0 max(pdfQtest_reGGA)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reGGA, 'b.', 'markersize', 15);
        xplot = linspace(0.99*mean(dischargeGGA), 1.01*mean(dischargeGGA));
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeGGA), std(dischargeGGA)), '--r');
        
        
        set(hl, 'linewidth', 1.5)
        if mm == 1
            legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
            legend boxoff
        end
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(dischargeGGA)/mean(dischargeGGA),3)])
        title('Q GGA')
        
        dischargeVTG = parameters.Q_reVTG{mm,8}; % only if only heading is of interest for heading
        
        [n, xout] = hist(parameters.Q_reVTG{mm,8},20); % for heading variation
        pdfQtest_reVTG = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reVTG)*1.1 > ymax;
            ymax = max(pdfQtest_reVTG)*1.1;
        end
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        
        figure(3)
        subplot (2,2,mm)
        hl(1) = plot([parameters.meanQ_reVTG(mm,8) parameters.meanQ_reVTG(mm,8)], [0 max(pdfQtest_reVTG)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reVTG, 'b.', 'markersize', 15);
        xplot = linspace(0.99*mean(dischargeVTG), 1.01*mean(dischargeVTG));
        %xplot = linspace(min(dischargeVTG), max(dischargeVTG));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeVTG), std(dischargeVTG)), '--r');
        
        set(hl, 'linewidth', 1.5)
        if mm == 1
            legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
            legend boxoff
        end
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(dischargeVTG)/mean(dischargeVTG),3)])
        title('Q VTG')
        
        
        figure(4)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reBot)*1.1 > ymax;
            ymax = max(pdfQtest_reBot)*1.1;
        end
        
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        %xplot = linspace(min(discharge), max(discharge));
        xplot = linspace(0.99*mean(discharge), 1.01*mean(discharge));
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        title('Q bottom')
        
    end
    
    figure(2)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_gga_normal_heading'
    print ('-depsc', [dn filename])
    print ('-dtiff', [dn filename])
    
    
    figure(3)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_vtg_normal_heading'
    print ('-depsc', [dn filename])
    print ('-dtiff', [dn filename])
    
    figure(4)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_bot_normal_heading'
    print ('-depsc', [dn filename])
    print ('-dtiff', [dn filename])
    
end


% Gatineau River heading normal error
doneF = 1
if doneF == 0
    
    load good_data/gatineau_20031203.mat
    Qbreakdown{1:4}
    
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    ymax = 0;
    xmin = 10000;
    xmax = 0;
    
    for mm = 1:4
        
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        %subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 10)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        x = caxis;
        x = round(100*x)/100;
        
        
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
        set(hc, 'YTick',x)
        
        
        
        % plot pdfs
        figure(2)
        subplot (2,2,mm)
        nr = length(parameters.Q_reGGA{mm,1})
        
        dischargeGGA = parameters.Q_reGGA{mm,1}; % only if only heading is of interest for heading
        
        [n, xout] = hist(parameters.Q_reGGA{mm,1},20); % for heading variation
        pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reGGA)*1.1 > ymax;
            ymax = max(pdfQtest_reGGA)*1.1;
        end
        
        
        hl(1) = plot([parameters.meanQ_reGGA(mm,1) parameters.meanQ_reGGA(mm,1)], [0 max(pdfQtest_reGGA)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reGGA, 'b.', 'markersize', 15);
        xplot = linspace(0.99*mean(dischargeGGA), 1.01*mean(dischargeGGA));
        % xplot = linspace(min(dischargeGGA), max(dischargeGGA));
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeGGA), std(dischargeGGA)), '--r');
        
        
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(dischargeGGA)/mean(dischargeGGA),3)])
        title('Q GGA')
        
        
    end
    
    figure(2)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_gga_all_error'
    % print ('-depsc', [dn filename])
    print ('-dtiff', [dn filename])
    
end


% ALL TRANSECTS FOR MJ001
doneF = 1
if doneF == 0
    figure(1), clf(1), hold on
    mm = 1
    
    load good_data/05MJ001_20080604_angus_colin.mat
    % trans 1 goes L to R
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    subplot (4,1,1)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    % xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    %title('Velocity magnitude [m/s]', 'fontsize', 10)
    % axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1.2])
    %x = caxis;
    %x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    %   set(hc, 'YTick',x)
    text(0.01, 1.3, ['A - T1: used mode 11, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),1),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    text(1.0, 1.25, 'Velocity [m/s]', 'units', 'normalized', 'fontsize', 10)
    
    load good_data/05MJ001_20080604_hardy_karen.mat
    % trans 1 goes L to R
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (4,1,2)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    %xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    %title('Velocity magnitude [m/s]', 'fontsize', 10)
    % axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1.2])
    %x = caxis;
    %x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    % set(hc, 'YTick',x)
    text(0.01, 1.12, ['B - T1: used mode 11, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),1),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    load good_data/05MJ001_20080604_hood_don.mat
    % trans 1 goes L to R
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (4,1,3)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    %xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    % title('Velocity magnitude [m/s]', 'fontsize', 10)
    % axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1.2])
    %x = caxis;
    %x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    %   set(hc, 'YTick',x)
    text(0.01, 1.12, ['C - T1: used mode 11, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),1),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    load good_data/05MJ001_20080604_selinger_dan.mat
    % trans 1 goes R to L
    % because of this you had to change the calculation for length_x and
    % also the xlim parameters
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT(end) - length_refBT + parameters.data{mm}.lDist;
    
    subplot (4,1,4)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel('Distance from left bank (ref: bottom) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    %  title('Velocity magnitude [m/s]', 'fontsize', 10)
    % axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(1) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1.2])
    % x = caxis;
    % x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    %   set(hc, 'YTick',x)
    text(0.01, 1.12, ['D - T1: used mode 5, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),1),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    filename = '05MJ001_20080604_all_xsec'
    print ('-depsc', [dn filename])
    
end





% all transects for Whitemouth at Whitemouth
doneF = 1
if doneF == 0
    figure(1), clf(1), hold on
    mm = 1
    
    load good_data/05PH003_20080603_angus_colin.mat
    % trans 1 goes L to R
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    subplot (4,1,1)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    % xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    % title('Velocity magnitude [m/s]', 'fontsize', 10)
    %  axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1])
    % x = caxis;
    % x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    %set(hc, 'YTick',x)
    text(0.01, 1.2, ['A - T1: used mode 12, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    load good_data/05PH003_20080603_hardy_karen.mat
    % trans 1 goes L to R
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (4,1,2)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    %xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    %title('Velocity magnitude [m/s]', 'fontsize', 10)
    %axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1])
    %x = caxis;
    %x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    %set(hc, 'YTick',x)
    text(0.01, 1.12, ['B - T1: used mode 12, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    load good_data/05PH003_20080603_hood_don.mat
    % trans 1 goes L to R
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (4,1,3)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    %xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    % title('Velocity magnitude [m/s]', 'fontsize', 10)
    % axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    hc = colorbar;
    caxis([0 1])
    % x = caxis;
    % x = round(100*x)/100;
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    % set(hc, 'YTick',x)
    text(0.01, 1.12, ['C - T1: used mode 5, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    load good_data/05PH003_20080603_selinger_dan.mat
    % trans 1 goes R to L, must account for this
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    %   length_x = length_refBT + parameters.data{mm}.lDist;
    length_x = length_refBT(end) - length_refBT + parameters.data{mm}.lDist; % bcs it goes R to L unlike the others
    subplot (4,1,4)
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel('Distance from left bank (ref: bottom) [m]', 'fontsize', 10)
    ylabel('Depth [m]', 'fontsize', 10)
    set(gca, 'box', 'on')
    %  title('Velocity magnitude [m/s]', 'fontsize', 10)
    %  axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    % xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    xlim([0 (length_x(1) + parameters.data{mm}.rDist)])
    hc = colorbar;
    % x = caxis;
    % x = round(100*x)/100;
    caxis([0 1])
    set(gca, 'linewidth', 1.5, 'YTick', [0 1])
    % set(hc, 'YTick',x)
    text(0.01, 1.12, ['D - T1: used mode 12, \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),3), ' m^3/s'], 'units', 'normalized')
    
    
    filename = '05PH003_20080603_all_xsec'
    print ('-depsc', [dn filename])
    
end
return




% Dan Sellinger
% 05PH003_20080603 Whitemouth River near Whitemouth
doneF = 1
if doneF == 0
    
    load good_data/05PH003_20080603_selinger_dan.mat
    Qbreakdown{1:4}
    
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    ymax = 0;
    xmin = 10000;
    xmax = 0;
    
    for mm = 1:4
        
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        %subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 10)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        x = caxis;
        x = round(100*x)/100;
        
        
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
        set(hc, 'YTick',x)
        
        if mm == 1
            filename = '05PH003_20080603_selinger_xsec'
            print ('-depsc', [dn filename])
        end
        
        figure(2)
        subplot (2,2,mm)
        hp = pie(real(parameters.frac(mm,:)));
        [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
        set(ax, 'position', [0.08 0.3 .01 .01], 'units', 'normalized') % cant put it any further left
        text(-0.1, 1.1, ['T', num2str(mm), ': \sigma =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
        legend boxoff
        set(gca, 'fontsize', 8.5)
        
        
        figure(3)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        if max(pdfQtest_reBot)*1.1 > ymax;
            ymax = max(pdfQtest_reBot)*1.1;
        end
        
        
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        % xplot = linspace(min(discharge), max(discharge));
        xplot = linspace(parameters.meanQ_reBot(mm)*0.9, 1.1*parameters.meanQ_reBot(mm));
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        %xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlim([xmin xmax]), ylim([0 ymax])%yli
        xlabel('Discharge [m^3/s]', 'fontsize', 12)
        ylabel('Probability density', 'fontsize', 12)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        
    end
    
    figure(2)
    filename = '05PH003_20080603_selinger_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
    
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05PH003_20080603_selinger_distns'
    print ('-depsc', [dn filename])
    
end



% Colin Angus
% 05PH003_20080603 Whitemouth near Whitemouth
doneF = 0
if doneF == 0
    
    load good_data/05PH003_20080603_angus_colin.mat
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    figure(4), clf(4), hold on
    
    
    for mm = 1:4
        
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
        ylabel('Depth [m]', 'fontsize', 16)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 16)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        colorbar
        set(gca, 'linewidth', 1.5)
        
        %     dn = 'C:\Users\Stephanie\Documents\environment canada\may2014\';
        %     filename = '05MJ001_data'
        %     print ('-dpng', '-r600', [dn filename])
        %
        %
        %     %axis normal
        %     filename = '05MJ001_data_normal'
        %     %print ('-dpng', '-r600', [dn filename])
        %
        
        
        figure(2)
        subplot (2,2,mm)
        hp = pie(real(parameters.frac(mm,:)));
        [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
        set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
        text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
        legend boxoff
        filename = '05MJ001_pie'
        %print ('-dpng', '-r600',[dn filename])
        
        
        figure(3)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        xplot = linspace(min(discharge), max(discharge));
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlabel('Discharge [m^3/s]', 'fontsize', 16)
        ylabel('Probability density', 'fontsize', 16)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        
        filename = '05MJ001_mc'
        %print ('-dpng', '-r600', [dn filename])
        % pause
    end
end


% Karen Hardy
% 05PH003_20080603 Whitemouth near Whitemouth
doneF = 0
if doneF == 0
    
    load good_data/05PH003_20080603_hardy_karen.mat
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    figure(4), clf(4), hold on
    
    
    for mm = 1:4
        
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
        ylabel('Depth [m]', 'fontsize', 16)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 16)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        colorbar
        set(gca, 'linewidth', 1.5)
        
        figure(2)
        subplot (2,2,mm)
        hp = pie(real(parameters.frac(mm,:)));
        [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
        set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
        text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
        legend boxoff
        filename = '05MJ001_pie'
        %print ('-dpng', '-r600',[dn filename])
        
        
        figure(3)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        xplot = linspace(min(discharge), max(discharge));
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlabel('Discharge [m^3/s]', 'fontsize', 16)
        ylabel('Probability density', 'fontsize', 16)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        
    end
end



% Dan Hood
% 05PH003_20080603 Whitemouth near Whitemouth
doneF = 0
if doneF == 0
    
    load good_data/05PH003_20080603_hood_don.mat
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    figure(4), clf(4), hold on
    
    
    for mm = 1:4
        
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
        ylabel('Depth [m]', 'fontsize', 16)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 16)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        colorbar
        set(gca, 'linewidth', 1.5)
        
        figure(2)
        subplot (2,2,mm)
        hp = pie(real(parameters.frac(mm,:)));
        [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
        set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
        text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
        legend boxoff
        filename = '05MJ001_pie'
        %print ('-dpng', '-r600',[dn filename])
        
        
        figure(3)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        xplot = linspace(min(discharge), max(discharge));
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlabel('Discharge [m^3/s]', 'fontsize', 16)
        ylabel('Probability density', 'fontsize', 16)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        
    end
end

% Dan Sellinger
% 05PH003_20080603 Whitemouth near Whitemouth
doneF = 0
if doneF == 0
    
    load good_data/05PH003_20080603_selinger_dan.mat
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    figure(4), clf(4), hold on
    
    
    for mm = 1:4
        
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
        ylabel('Depth [m]', 'fontsize', 16)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 16)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        colorbar
        set(gca, 'linewidth', 1.5)
        
        figure(2)
        subplot (2,2,mm)
        hp = pie(real(parameters.frac(mm,:)));
        [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
        set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
        text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
        legend boxoff
        filename = '05MJ001_pie'
        %print ('-dpng', '-r600',[dn filename])
        
        
        figure(3)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        xplot = linspace(min(discharge), max(discharge));
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlabel('Discharge [m^3/s]', 'fontsize', 16)
        ylabel('Probability density', 'fontsize', 16)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        
    end
end




% Canadian St. Marys Canal near Spring Coulee
doneF = 0
if doneF == 0
    load good_data/05AE026_20100719_O_Connor_AQ1.mat
    
    figure(1), clf(1), hold on
    figure(2), clf(2), hold on
    figure(3), clf(3), hold on
    figure(4), clf(4), hold on
    
    
    for mm = 1:4
        
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
        
        figure(1)
        subplot (2,2,mm)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
        ylabel('Depth [m]', 'fontsize', 16)
        set(gca, 'box', 'on')
        title('Velocity magnitude [m/s]', 'fontsize', 16)
        axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        colorbar
        set(gca, 'linewidth', 1.5)
        
        figure(2)
        subplot (2,2,mm)
        hp = pie(real(parameters.frac(mm,:)));
        [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
        set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
        text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
        legend boxoff
        filename = '05AE_pie'
        %print ('-dpng', '-r600',[dn filename])
        
        
        figure(3)
        subplot (2,2,mm)
        nr = length(parameters.Q_reBot{mm,1})
        discharge = parameters.Q_reBot{mm,1};
        
        [n, xout] = hist(parameters.Q_reBot{mm,1},20);
        pdfQtest_reBot = (n/nr)/diff(xout(1:2));
        
        hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
        hold on
        hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
        xplot = linspace(min(discharge), max(discharge));
        hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
        xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
        xlabel('Discharge [m^3/s]', 'fontsize', 16)
        ylabel('Probability density', 'fontsize', 16)
        display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
        
    end
end

return






% WEIR
doneF = 1 % DONE on 2013/10/15
if doneF == 0
    % load Weir data
    % these data were obtained with a 1200 kHz Rio Grande
    %
    load Weir_8_20120503.mat
    
    % measurement transect index
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    meanQ = nanmean(parameters.meanQ(:,1)) % the mean of the discharge from each transect (no uncertainty included)
    100*nanstd(parameters.meanQ(:,1))/meanQ % percent deviation between transects
    sqrt(nansum(parameters.perErr(:,1).^2))/ sqrt(sum(parameters.transno>0))
    parameters.perErr(:,1)
    return
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    figure
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel('Dist from Left Bank (Ref: BT) [m]', 'fontsize', 16)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    title('Velocity magnitude [m/s]', 'fontsize', 16)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    set(gca, 'linewidth', 1.5)
    
    dn = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\report\figures\';
    filename = 'weir8_welland_canal_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    
    
    figure
    % if labelling pie
    %pie(real(parameters.frac(mm,:)), parameters.name(2:end)')
    % if adding legend
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    %text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ{mm},4), ' m^3/s'], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ(mm),4), ' m^3/s'], 'units', 'normalized')
    
    legend boxoff
    filename = 'weir8_welland_canal_pie'
    print ('-dpng', '-r600',[dn filename])
    %title(['Transect ', num2str{mm}])
    
    
    figure
    nr = length(parameters.Q{mm,1})
    discharge = parameters.Q{mm,1};
    
    [n, xout] = hist(parameters.Q{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    %hl(1) = plot([parameters.meanQ{mm} parameters.meanQ{mm}], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hl(1) = plot([parameters.meanQ(mm) parameters.meanQ(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(105, 175);%(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    %title(['Discharge referenced to bottom track ', parameter, ' is varied'])
    %ylim([0 max(pdfQtest_reBot)])
    %xlim([.9*parameters.meanQ{mm} 1.1*parameters.meanQ{mm}])
    xlim([105 175]), ylim([0 max(pdfQtest_reBot)*1.1])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'weir8_welland_canal_mc'
    print ('-dpng', '-r600', [dn filename])
end




%  Gatineau data
doneF = 1 % not done yet
if doneF == 0
    
    load good_data/gatineau_20031203.mat
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr_reGGA(mm,2:end)./sum(parameters.perErr_reGGA(mm,2:end))
    meanQ_reGGA = nanmean(parameters.meanQ_reGGA(:,1)) % the mean of the discharge from each transect (no uncertainty included)
    100*nanstd(parameters.meanQ_reGGA(:,1))/meanQ_reGGA % percent deviation between transects
    sqrt(nansum(parameters.perErr_reGGA(:,1).^2))/ sqrt(sum(parameters.transno>0))
    parameters.perErr_reGGA(:,1)
    
    length_refGGA =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel_reGGA(1,:).^2 + (parameters.data{1}.btVel_reGGA(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx_reGGA.^2 + parameters.data{1}.wVely_reGGA.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refGGA + parameters.data{mm}.lDist;
    
    figure
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    title('Velocity magnitude [m/s]', 'fontsize', 16)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    set(gca, 'linewidth', 1.5)
    
    dn = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\report\figures\';
    filename = 'gatineau_data';
    % print ('-dpng', '-r600', [dn filename])
    
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reGGA(mm,1)),2),'% of ', num2str(parameters.meanQ_reGGA(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'gatineau_pie';
    % print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q_reGGA{mm,1})
    discharge = parameters.Q_reGGA{mm,1};
    
    [n, xout] = hist(parameters.Q_reGGA{mm,1},20);
    pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ_reGGA(mm) parameters.meanQ_reGGA(mm)], [0 max(pdfQtest_reGGA)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reGGA, 'b.', 'markersize', 15);
    xplot = linspace(mean(discharge)*.95, mean(discharge)*1.05);%linspace(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]), ylim([0 max(pdfQtest_reGGA)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'gatineau_mc';
    % print ('-dpng', '-r600', [dn filename])
    
    load good_data/gatineau_20031203_gps_vary_heading_only
    figure
    nr = length(parameters.Q_reGGA{mm,8})
    discharge = parameters.Q_reGGA{mm,8};
    
    [n, xout] = hist(parameters.Q_reGGA{mm,8},20);
    pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ_reGGA(mm) parameters.meanQ_reGGA(mm)], [0 max(pdfQtest_reGGA)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reGGA, 'b.', 'markersize', 15);
    xplot = linspace(mean(discharge)*.99, mean(discharge)*1.01);%linspace(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]), ylim([0 max(pdfQtest_reGGA)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
end


%  Mueller_M18_streampro
doneF = 1
if doneF == 0
    
    load good_data/Mueller_M18_streampro.mat
    
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    figure
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    title('Velocity magnitude [m/s]', 'fontsize', 16)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    set(gca, 'linewidth', 1.5)
    
    dn = 'C:\Users\Stephanie\Documents\environment canada\may2014\';
    filename = 'mueller_m18_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    %axis normal
    %filename = 'mueller_m18_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'mueller_m18_pie'
    print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q_reBot{mm,1})
    discharge = parameters.Q_reBot{mm,1};
    
    [n, xout] = hist(parameters.Q_reBot{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'mueller_m18_mc'
    print ('-dpng', '-r600', [dn filename])
end

%Mueller RG1200
doneF = 1
if doneF == 0
    
    load good_data/Mueller_RG1200.mat
    
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    figure
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    title('Velocity magnitude [m/s]', 'fontsize', 16)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    set(gca, 'linewidth', 1.5)
    
    dn = 'C:\Users\Stephanie\Documents\environment canada\may2014\';
    filename = 'mueller_RG1200_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    axis normal
    filename = 'mueller_RG1200_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'mueller_RG1200_pie'
    print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q_reBot{mm,1})
    discharge = parameters.Q_reBot{mm,1};
    
    [n, xout] = hist(parameters.Q_reBot{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'mueller_RG1200_mc'
    print ('-dpng', '-r600', [dn filename])
end

% Mueller 478
doneF = 1
if doneF == 0
    
    load good_data/Mueller_478_20111130_streampro.mat
    
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % In case there is a bad ensemble and there is no depth measurement
    botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    figure
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    xlabel(['Dist from L Bank (Ref: BT) [m]'], 'fontsize', 16)
    ylabel('Depth [m]', 'fontsize', 16)
    set(gca, 'box', 'on')
    title('Velocity magnitude [m/s]', 'fontsize', 16)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    colorbar
    set(gca, 'linewidth', 1.5)
    
    dn = 'C:\Users\Stephanie\Documents\environment canada\may2014\';
    filename = 'mueller_478_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    axis normal
    filename = 'mueller_478_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'mueller_478_pie'
    print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q_reBot{mm,1})
    discharge = parameters.Q_reBot{mm,1};
    
    [n, xout] = hist(parameters.Q_reBot{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ_reBot(mm) parameters.meanQ_reBot(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'mueller_478_mc'
    print ('-dpng', '-r600', [dn filename])
end


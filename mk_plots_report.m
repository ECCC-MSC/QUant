close all, clear all
addpath 'C:\Users\Collin Rennie\Documents\MATLAB\work\general'
%addpath 'E:\matlab_code_uottawa\work\general'
addpath 'tools/'
addpath 'good_data/'
set(0,'defaulttextinterpreter','tex')

dn = 'C:\Users\Collin Rennie\Documents\Moore\my publications\uncertainty_paper\figures\';


% Gatineau River heading including one and two cycle errors
doneF = 0
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
            %    filename = '05MJ001_20080604_angus_xsec'
            %      print ('-depsc', [dn filename])
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
        xplot = linspace(min(dischargeGGA), max(dischargeGGA));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeGGA), std(dischargeGGA)), '--r');
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
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
        xplot = linspace(min(dischargeVTG), max(dischargeVTG));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeVTG), std(dischargeVTG)), '--r');
        
        if xplot(1) < xmin
            xmin = xplot(1);
        end
        if xplot(end) > xmax
            xmax = xplot(end);
        end
        set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
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
    
    
    figure(3)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_vtg_vary_heading'
    print ('-depsc', [dn filename])
    
    %     figure(2)
    %     subplot 221,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %     subplot 222,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %     subplot 223,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %     subplot 224,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %
    %     figure(3)
    %     subplot 221,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %     subplot 222,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %     subplot 223,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])
    %     subplot 224,  xlim([0.99*xmin 1.01*xmax]), ylim([0 ymax])

end
return

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
        
        if mm == 1
        %    filename = '05MJ001_20080604_angus_xsec'
      %      print ('-depsc', [dn filename])
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
        xplot = linspace(min(dischargeGGA), max(dischargeGGA));
        
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
        xplot = linspace(min(dischargeVTG), max(dischargeVTG));
        hl(3) = plot(xplot, normpdf(xplot,mean(dischargeVTG), std(dischargeVTG)), '--r');
 
       set(hl, 'linewidth', 1.5)
        legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
        legend boxoff
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
        xplot = linspace(min(discharge), max(discharge));
        
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
    
    
    figure(3)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_vtg_normal_heading'
    print ('-depsc', [dn filename])
    
        figure(4)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = 'gatineau_bot_normal_heading'
    print ('-depsc', [dn filename])
    
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
        
        if mm == 1
        %    filename = '05MJ001_20080604_angus_xsec'
      %      print ('-depsc', [dn filename])
        end
        
      
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
        xplot = linspace(min(dischargeGGA), max(dischargeGGA));
        
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
    print ('-depsc', [dn filename])

end

% ALL TRANSECTS FOR MJ001
doneF = 0
if doneF == 0
    figure(1), clf(1), hold on
  mm = 1
    
    load good_data/05MJ001_20080604_angus_colin.mat
     
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
        x = caxis;
        x = round(100*x)/100;
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
     %   set(hc, 'YTick',x)
        text(0.01, 1.3, 'A: used mode 11', 'units', 'normalized')
        
        
        load good_data/05MJ001_20080604_hardy_karen.mat
     
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
        x = caxis;
        x = round(100*x)/100;  
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
       % set(hc, 'YTick',x)
        text(0.01, 1.3, 'B: used mode 11', 'units', 'normalized')
        
        
        load good_data/05MJ001_20080604_hood_don.mat
     
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
        x = caxis;
        x = round(100*x)/100;
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
     %   set(hc, 'YTick',x)
       text(0.01, 1.3, 'C: used mode 11', 'units', 'normalized')
        
        
        load good_data/05MJ001_20080604_selinger_dan.mat
     
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
     
        subplot (4,1,4)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
      %  title('Velocity magnitude [m/s]', 'fontsize', 10)
       % axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        x = caxis;
        x = round(100*x)/100;
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
     %   set(hc, 'YTick',x)
       text(0.01, 1.3, 'D: used mode 5', 'units', 'normalized')
        
        
            filename = '05MJ001_20080604_all_xsec'
            print ('-depsc', [dn filename])

end
return

% Colin Angus
% 05MJ001_20080604 Assiniboine River near Headingly
doneF = 1
if doneF == 0
    
    load good_data/05MJ001_20080604_angus_colin.mat
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
            filename = '05MJ001_20080604_angus_xsec'
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
        

        hl(1) = plot([parameters.meanQ_reBot(mm,1) parameters.meanQ_reBot(mm,1)], [0 max(pdfQtest_reBot)*1.1], '-.k');
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
    filename = '05MJ001_20080604_angus_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
    subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    
    filename = '05MJ001_20080604_angus_distns'
    print ('-depsc', [dn filename])
    
end

doneF = 1
if doneF == 0
    
    load good_data/05MJ001_20080604_hardy_karen.mat
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
            filename = '05MJ001_20080604_hardy_xsec'
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
    filename = '05MJ001_20080604_hardy_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
  subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05MJ001_20080604_hardy_distns'
    print ('-depsc', [dn filename])
    
end


% Don Hood
% 05MJ001_20080604 Assiniboine River near Headingly
doneF = 1
if doneF == 0
    
    load good_data/05MJ001_20080604_hood_don.mat
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
            filename = '05MJ001_20080604_hood_xsec'
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
    filename = '05MJ001_20080604_hood_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
 subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05MJ001_20080604_hood_distns'
    print ('-depsc', [dn filename])
    
end

% Dan Sellinger
% 05MJ001_20080604 Assiniboine River near Headingly
doneF = 1
if doneF == 0
    
    load good_data/05MJ001_20080604_selinger_dan.mat
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
            filename = '05MJ001_20080604_selinger_xsec'
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
      %  xplot = linspace(min(discharge), max(discharge));
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
    filename = '05MJ001_20080604_selinger_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
 subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05MJ001_20080604_selinger_distns'
    print ('-depsc', [dn filename])
    
end

% all transects for Whitemouth at Whitemouth
doneF = 1
if doneF == 0
    figure(1), clf(1), hold on
  mm = 1
    
    load good_data/05PH003_20080603_angus_colin.mat
     
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
        x = caxis;
        x = round(100*x)/100;
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
        %set(hc, 'YTick',x)
        text(0.01, 1.3, 'A: used mode 12', 'units', 'normalized')
        
        
        load good_data/05PH003_20080603_hardy_karen.mat
     
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
        x = caxis;
        x = round(100*x)/100;  
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
        %set(hc, 'YTick',x)
        text(0.01, 1.3, 'B: used mode 12', 'units', 'normalized')
        
        
        load good_data/05PH003_20080603_hood_don.mat
     
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
        x = caxis;
        x = round(100*x)/100;
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
       % set(hc, 'YTick',x)
        text(0.01, 1.3, 'C: used mode 5', 'units', 'normalized')
        
        
        load good_data/05PH003_20080603_selinger_dan.mat
     
        parameters.frac(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
        
        length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
        vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
        cellDepth = parameters.data{1}.cellDepth;
        botDepth = parameters.data{1}.depthEns;
        % In case there is a bad ensemble and there is no depth measurement
        botDepth(botDepth ==  parameters.data{mm}.ddraft.mean) = NaN;
        length_x = length_refBT + parameters.data{mm}.lDist;
     
        subplot (4,1,4)
        hp = pcolor(length_x, cellDepth, vel);
        hold on
        set(gca, 'ydir', 'rev')
        set(hp,'edgecolor','none')
        shading flat
        plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
        xlabel('Dist from L Bank (Ref: BT) [m]', 'fontsize', 10)
        ylabel('Depth [m]', 'fontsize', 10)
        set(gca, 'box', 'on')
      %  title('Velocity magnitude [m/s]', 'fontsize', 10)
      %  axis equal
        ylim([0 (max(botDepth) + .1*max(botDepth))])
        xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
        hc = colorbar;
        x = caxis;
        x = round(100*x)/100;
        set(gca, 'linewidth', 1.5, 'YTick', [0 1])
       % set(hc, 'YTick',x)
        text(0.01, 1.3, 'D: used mode 12', 'units', 'normalized')
        
        
            filename = '05PH003_20080603_all_xsec'
            print ('-depsc', [dn filename])

end

return
% Colin Angus
% 05PH003_20080603_angus_colin Whitemouth River near Whitemouth
doneF = 1
if doneF == 0
    
    load good_data/05PH003_20080603_angus_colin.mat
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
            filename = '05PH003_20080603_angus_xsec'
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
    filename = '05PH003_20080603_angus_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
 subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05PH003_20080603_angus_distns'
    print ('-depsc', [dn filename])
    
end

% 05PH003_20080603 Karen Hardy Whitemouth River near Whitemouth
doneF = 1
if doneF == 0
    
    load good_data/05PH003_20080603_hardy_karen.mat
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
            filename = '05PH003_20080603_hardy_xsec'
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
        %xplot = linspace(min(discharge), max(discharge));
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
    filename = '05PH003_20080603_hardy_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
 subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05PH003_20080603_hardy_distns'
    print ('-depsc', [dn filename])
    
end

% Don Hood
% 05PH003_20080603 Whitemouth River near Whitemouth
doneF = 1
if doneF == 0
    
    load good_data/05PH003_20080603_hood_don.mat
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
            filename = '05PH003_20080603_hood_xsec'
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
        %xplot = linspace(min(discharge), max(discharge));
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
    filename = '05PH003_20080603_hood_pies'
    print ('-depsc', [dn filename])
    
    figure(3)
 subplot 221,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 222,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 223,  xlim([xmin xmax]), ylim([0 ymax])
    subplot 224,  xlim([xmin xmax]), ylim([0 ymax])
    filename = '05PH003_20080603_hood_distns'
    print ('-depsc', [dn filename])
    
end

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

return

% Colin Angus
% 05PH003_20080603 Whitemouth near Whitemouth
doneF = 1
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
doneF = 1
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
doneF = 1
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
doneF = 1
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


% Gatineau River

doneF = 1
if doneF == 0
    
   load good_data/gatineau_20031203.mat
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
  
return



% Graham Creek Data
doneF = 1;
if doneF == 0
    mm = 1;
    
    measErr = NaN*ones(1,7);
    
    % load Graham Creek Data
    
    load good_data/graham_creek_20130723_streampro.mat
    measErr(1) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    
    load good_data/Weir_8_20120503.mat
    measErr(2) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    
    load good_data/assiniboine_20130515.mat
    measErr(3) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    
    load good_data/assiniboine_05MH005_20110517.mat
    measErr(4) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    
    load good_data/montelimar_20120301.mat
    measErr(5) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    
    load good_data/assiniboine_05MH005_20130516.mat
    measErr(6) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    
    load good_data/gatineau_20031203.mat
    measErr(7) = real(parameters.perErr(mm,1));
    clearvars -except mm measErr
    

    cmin = 0;
    cmax = 4;
    mfactorL = 2;
    mfactorR = 1.5;
    

    figure(1), clf(1), hold on
    set(gcf, 'units', 'normalized')
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % PLOT MONTELIMAR
    load good_data/montelimar_20120301.mat
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
        % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (7,2,1)
    clear hp
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    %caxis([cmin cmax]), 
    caxis([0 max(max(vel))])
    colorbar
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1:2)  mfactorL*oldPos(3:4)])
    %text(1.15,1,['Montelimar: ', num2str(parameters.meanQ{mm},4), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
  title(['Montelimar: ', num2str(parameters.meanQ(mm),4), ' m ^3/s'],'fontsize', 12, 'units', 'normalized')
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot (7,2,2)
    hp5 = pie(real(parameters.frac(mm,:)));
    
    hText = findobj(hp5,'Type','text');
    set(hText, 'String', {''})
        set(gca, 'units', 'normalized')
    
    pieAxis = get(hp5(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');
    
    newRadius = mfactorR*measErr(5)/max(measErr); % 1 for the largest
    newPieAxisPosition = pieAxisPosition .* [1 1 newRadius newRadius];
    
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = 1.2*(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
         
    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.2, .5, ['\sigma = ', num2str(real(parameters.perErr(mm,1)),2),'%'], 'units', 'normalized')

    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%    
    % PLOT GATINEAU
    load good_data/gatineau_20031203.mat
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (7,2,3)
    clear hp
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
  
    set(gca, 'box', 'on')
    %title('Velocity magnitude [m/s]', 'fontsize', 16)
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    %caxis([cmin cmax]), 
   caxis([0 max(max(vel))])
    colorbar
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1:2) mfactorL*oldPos(3:4)])
    %text(1.15,1,['Gatineau: ', num2str(parameters.meanQ{mm},3), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
    title(['Gatineau: ', num2str(parameters.meanQ(mm),3), ' m ^3/s'],'fontsize', 12, 'units', 'normalized')
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot (7,2,4)
    hp7 = pie(real(parameters.frac(mm,:)));
    
    hText = findobj(hp7,'Type','text');
    set(hText, 'String', {''})
    set(gca, 'units', 'normalized')
    
    pieAxis = get(hp7(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');

    
    newRadius = mfactorR*measErr(7)/max(measErr); % 1 for the largest
    newPieAxisPosition = pieAxisPosition .* [1 1 newRadius newRadius];
    
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = 1.2*(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
     
    
    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.2, .5, ['\sigma = ', num2str(real(parameters.perErr(mm,1)),2),'%'], 'units', 'normalized')

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%   
    % PLOT ASSINIBOINE Station  05MH013 20130515
    load good_data/assiniboine_20130515
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;

    
    subplot (7,2,5)
    clear hp
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    set(gca, 'box', 'on')
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    %caxis([cmin cmax]), 
    caxis([0 max(max(vel))])
    colorbar
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1:2) mfactorL*oldPos(3:4)])
    %text(1.15,1,['Assiniboine (05MH013, 2013/05/15) ', num2str(parameters.meanQ{mm},3), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
    title(['Assiniboine 05MH013 ', num2str(parameters.meanQ(mm),3), ' m ^3/s'],'fontsize', 12, 'units', 'normalized')
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot (7,2,6)
    hp3 = pie(real(parameters.frac(mm,:)));
    hText = findobj(hp3,'Type','text');
    set(hText, 'String', {''})
    set(gca, 'units', 'normalized')
    
    pieAxis = get(hp3(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');
    newRadius = mfactorR*measErr(3)/max(measErr); % 1 for the largest
    newPieAxisPosition = pieAxisPosition.*[1 1 newRadius newRadius];
    
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = 1.2*(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
 

    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.2, .5, ['\sigma =  ', num2str(real(parameters.perErr(mm,1)),2),'%'], 'units', 'normalized')


  
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT ASSINIBOINE Station 05MH005 20110517
    load good_data/assiniboine_05MH005_20110517.mat
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
        % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (7,2,7)
    clear hp
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    set(gca, 'box', 'on')
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    %caxis([cmin cmax]), 
    caxis([0 max(max(vel))])
    colorbar
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1:2) mfactorL*oldPos(3:4)])
    %text(1.15,1,['Assiniboine (05MH005, 2011/05/17) ', num2str(parameters.meanQ{mm},4), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
    title(['Assiniboine 05MH005 ', num2str(parameters.meanQ(mm),4), ' m ^3/s'],'fontsize', 12, 'units', 'normalized')
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot (7,2,8)
     toplot =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    hp4 = pie(real(toplot));
    
    hText = findobj(hp4,'Type','text');
    set(hText, 'String', {''})

    pieAxis = get(hp4(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');
    
    newRadius = mfactorR*measErr(4)/max(measErr); % 1 for the largest
    newPieAxisPosition = pieAxisPosition .* [1 1 newRadius newRadius];
    
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = 1.2*(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
         
    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.1, 0.5, ['\sigma =  ', num2str(real(parameters.perErr(mm,1)),2),'% '], 'units', 'normalized')

    
   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT WELLAND CANAL
    load good_data/Weir_8_20120503.mat
    shiftDown = 0.02;
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
        % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;
   

    subplot (7,2,9)
   
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    set(gca, 'box', 'on')
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    %caxis([cmin cmax]), 
    caxis([0 max(max(vel))])
    colorbar
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1) (-shiftDown + oldPos(2)) mfactorL*oldPos(3:4)])
    %text(1.15,1,['Welland Canal ', num2str(parameters.meanQ{mm},3), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
    title(['Welland Canal: ', num2str(parameters.meanQ(mm),3), ' m  ^3/s'],'fontsize', 12, 'units', 'normalized')
     
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    subplot (7,2,10)
    hp2 = pie(real(parameters.frac(mm,:)));
    hText = findobj(hp2,'Type','text');
    set(hText, 'String', {''})
    set(gca, 'units', 'normalized')
    
    pieAxis = get(hp2(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');

    newRadius = mfactorR*measErr(2)/max(measErr); % 1 for the largest
    newPieAxisPosition = pieAxisPosition .* [1 1 newRadius newRadius];
    
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = 1.2*(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
        
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
   
    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.1, 0.5, ['\sigma = ', num2str(real(parameters.perErr(mm,1)),2),'%'], 'units', 'normalized')
  

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % PLOT ASSINIBOINE STATION 05MH005 DATE 20130516
    load good_data/assiniboine_05MH005_20130516.mat
    shiftDown = 0.03;
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
        % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    
    subplot (7,2,11)
    clear hp
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev', 'units', 'normalized')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    set(gca, 'box', 'on')
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    %caxis([cmin cmax]),
    caxis([0 max(max(vel))])
    colorbar
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1) (oldPos(2) - shiftDown)  mfactorL*oldPos(3:4)])
    %text(1.15,1,['Assiniboine (05MH005, 2013/05/16): ', num2str(parameters.meanQ{mm},4), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
    title(['Assiniboine 05MH005: ', num2str(parameters.meanQ(mm),4), ' m  ^3/s'],'fontsize', 12, 'units', 'normalized')

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot (7,2,12)
    hp6 = pie(real(parameters.frac(mm,:)));
    
    hText = findobj(hp6,'Type','text');
    set(hText, 'String', {''})
    set(gca, 'units', 'normalized')
    
    pieAxis = get(hp6(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');
    
    newRadius = mfactorR*measErr(6)/max(measErr); % 1 for the largest
    newPieAxisPosition = pieAxisPosition .* [1 1 newRadius newRadius];
    
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = +1.2*(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
     
    
    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.1, 0.5, ['\sigma = ', num2str(real(parameters.perErr(mm,1)),2),'%'], 'units', 'normalized')

    
   %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    % PLOT GRAHAM CREEK
    load good_data/graham_creek_20130723_streampro.mat
    shiftDown = 0.05;
    length_refBT =  nancumsum(parameters.data{1}.ensDeltaTime.*(sqrt(parameters.data{1}.btVel(1,:).^2 + (parameters.data{1}.btVel(2,:).^2))));
    vel = sqrt(parameters.data{1}.wVelx.^2 + parameters.data{1}.wVely.^2);
    cellDepth = parameters.data{1}.cellDepth;
    botDepth = parameters.data{1}.depthEns;
    
    % include distance from left bank
    length_x = length_refBT + parameters.data{mm}.lDist;
    
    subplot (7,2,13)
    clear hp
    hp = pcolor(length_x, cellDepth, vel);
    hold on
    set(gca, 'ydir', 'rev')
    set(hp,'edgecolor','none')
    shading flat
    plot(length_x, botDepth, 'k-.', 'linewidth', 1.5)
    set(gca, 'box', 'on')
    
    axis equal
    ylim([0 (max(botDepth) + .1*max(botDepth))])
    xlim([0 (length_x(end) + parameters.data{mm}.rDist)])
    %caxis([cmin cmax])
    caxis([0 max(max(vel))])
    colorbar
    set(gca, 'linewidth', 1.5, 'units', 'normalized')
    oldPos = get(gca, 'position');
    set(gca, 'position', [oldPos(1) (oldPos(2)-shiftDown) mfactorL*oldPos(3:4)])
    %text(1.15,1,['Graham Creek: ', num2str(parameters.meanQ{mm},3), ' m^3/s'],'fontsize', 12, 'units', 'normalized')
    title(['Graham Creek: ', num2str(parameters.meanQ(mm),3), ' m ^3/s'],'fontsize', 12, 'units', 'normalized')
    xlabel(['Length (Ref: BT) [m]'], 'fontsize', 14)
    text(-.1, 7, 'Depth [m]', 'fontsize', 14,'rotation', 90, 'units', 'normalized')

    
    %%%%%%%%%%%%%%%%%%%%
    subplot (7,2,14)
    hp1 = pie(real(parameters.frac(mm,:)));
    hold on
    hText = findobj(hp1,'Type','text');
    set(hText, 'String', {''})
    set(gca, 'units', 'normalized')
    
    pieAxis = get(hp1(1), 'Parent');
    pieAxisPosition = get(pieAxis, 'Position');
    
    newRadius = mfactorR*measErr(1)/max(measErr); % 1 for the largest
    newPieAxisPosition = (pieAxisPosition).* [1 1 newRadius newRadius];
   
    deltaXpos = (pieAxisPosition(3) - newPieAxisPosition(3))/2; % Move axis position left or right
    deltaYpos = -shiftDown -(pieAxisPosition(4) - newPieAxisPosition(4)); % Move axis position up or down
    newPieAxisPosition = (pieAxisPosition + [deltaXpos deltaYpos 0 0]) .* [1 1 newRadius newRadius];
       
   
    set(pieAxis, 'Position', newPieAxisPosition); % Set pieAxis to new position
    text(1.1, 0.5, ['\sigma =  ', num2str(real(parameters.perErr(mm,1)),2),'%'], 'units', 'normalized')

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%


    [ax,leg] = legend(hp1(1:2:end), parameters.name(2:end)); 
    set(ax, 'position', [.94 0.8 .001 .005], 'units', 'normalized', 'fontsize', 5)
    legend boxoff

    dn = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\report\figures\';
    filename = 'all_data_compare'
    print ('-dpng', '-r600',[dn filename])
    
    
end


    % Graham Creek Data
    % very shallow
    % 2400 kHz Streampro
doneF = 1
if doneF == 0
    
    load good_data/graham_creek_20130723_streampro.mat
        mm = 1;
    parameters.frac(mm,:) =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    meanQ = nanmean(parameters.meanQ(:,1)) % the mean of the discharge from each transect (no uncertainty included)
    100*nanstd(parameters.meanQ(:,1))/meanQ % percent deviation between transects
    sqrt(nansum(parameters.perErr(:,1).^2))/ sqrt(sum(parameters.transno>0))
    parameters.perErr(:,1)
 
    
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
    filename = 'graham_creek_trans1_data'
    print ('-dpng', '-r600', [dn filename])
    
  %  axis normal
  %  filename = 'graham_creek_trans1_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    
    figure
    % if labelling pie
    %pie(real(parameters.frac(mm,:)), parameters.name(2:end)')
    % if adding legend
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ(mm),3), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'graham_creek_trans1_pie'
    print ('-dpng', '-r600',[dn filename])
    %title(['Transect ', num2str{mm}])
    
    
    figure
    nr = length(parameters.Q{mm,1})
    discharge = parameters.Q{mm,1};
    
    [n, xout] = hist(parameters.Q{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ(mm) parameters.meanQ(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(mean(discharge)*.45, mean(discharge)*1.55);%linspace(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    %title(['Discharge referenced to bottom track ', parameter, ' is varied'])
    %ylim([0 max(pdfQtest_reBot)])
    %xlim([.9*parameters.meanQ{mm} 1.1*parameters.meanQ{mm}])
    xlim([-inf inf]), ylim([0 max(pdfQtest_reBot)*1.1])% ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'graham_creek_trans1_mc'
    print ('-dpng', '-r600', [dn filename])
    
end



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



% Assiniboine: Stn 05MH013 date 2013_05_15
doneF = 1 % DONE as of 2013/10/15
if doneF == 0
    
    load good_data/assiniboine_20130515.mat
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    meanQ = nanmean(parameters.meanQ(:,1)) % the mean of the discharge from each transect (no uncertainty included)
    100*nanstd(parameters.meanQ(:,1))/meanQ % percent deviation between transects
    sqrt(nansum(parameters.perErr(:,1).^2))/ sqrt(sum(parameters.transno>0))
    parameters.perErr(:,1)
    parameters.meanQ(:,1)
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
    filename = 'ass_05MH013_2013_05_15_data'
    print ('-dpng', '-r600', [dn filename])
    
  
    
    figure
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'ass_05MH013_2013_05_15_pie'
    print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q{mm,1})
    discharge = parameters.Q{mm,1};
    
    [n, xout] = hist(parameters.Q{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ(mm) parameters.meanQ(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(210,270);%(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([210 270]), ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'ass_05MH013_2013_05_15_mc'
    print ('-dpng', '-r600', [dn filename])
end


% Assiniboine near Holland Stn 05MH005 date 2011 05 17
% 600 kHz Rio grande
doneF = 1 % DONE 
if doneF == 0
    
    load good_data/assiniboine_05MH005_20110517.mat
    
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    meanQ = nanmean(parameters.meanQ(:,1)) % the mean of the discharge from each transect (no uncertainty included)
    100*nanstd(parameters.meanQ(:,1))/meanQ % percent deviation between transects
    sqrt(nansum(parameters.perErr(:,1).^2))/ sqrt(sum(parameters.transno>0))
    parameters.perErr(:,1)
    
    
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
    
    dn = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\report\figures\';
    filename = 'ass_05MH005_2011_05_17_data'
    print ('-dpng', '-r600', [dn filename])
    
   
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'ass_05MH005_2011_05_17_pie'
    %print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q{mm,1})
    discharge = parameters.Q{mm,1};
    
    [n, xout] = hist(parameters.Q{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ(mm) parameters.meanQ(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(parameters.meanQ(mm) - 180, parameters.meanQ(mm)+180);%(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]),  ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'ass_05MH005_2011_05_17_mc'
    %print ('-dpng', '-r600', [dn filename])
end


% Assiniboine near Holland Stn 05MH005 date 20130516
% 600 kHz Rio grande
doneF = 1 % DONE 
if doneF == 0
    
    load good_data/assiniboine_05MH005_20130516.mat
    
    mm = 1;
    parameters.frac(mm,:) =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    meanQ = nanmean(parameters.meanQ(:,1)) % the mean of the discharge from each transect (no uncertainty included)
    100*nanstd(parameters.meanQ(:,1))/meanQ % percent deviation between transects
    sqrt(nansum(parameters.perErr(:,1).^2))/ sqrt(sum(parameters.transno>0))
    parameters.perErr(:,1)
    
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
    
    dn = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\report\figures\';
    filename = 'ass_05MH005_20130516_data';
    print ('-dpng', '-r600', [dn filename])
    
   
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'ass_05MH005_20130516_pie';
    print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q{mm,1})
    discharge = parameters.Q{mm,1};
    
    [n, xout] = hist(parameters.Q{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ(mm) parameters.meanQ(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
    hold on
    hl(2) = plot(xout, pdfQtest_reBot, 'b.', 'markersize', 15);
    xplot = linspace(parameters.meanQ(mm) - 180, parameters.meanQ(mm)+180);%(min(discharge), max(discharge));
    hl(3) = plot(xplot, normpdf(xplot,mean(discharge), std(discharge)), '--r');
    set(hl, 'linewidth', 1.5)
    legend(hl, 'Winriver II', 'Monte Carlo (MC)', 'Gaussian from MC')
    legend boxoff
    xlim([-inf inf]),  ylim([0 max(pdfQtest_reBot)*1.1])%ylim([-inf inf])
    xlabel('Discharge [m^3/s]', 'fontsize', 16)
    ylabel('Probability density', 'fontsize', 16)
    display(['the uncertainty in percent is ', num2str(100*std(discharge)/mean(discharge),3)])
    
    filename = 'ass_05MH005_20130516_mc';
    print ('-dpng', '-r600', [dn filename])
end


% load Montelimar
% trapezoidal canal
% 1200 kHz rio grande
doneF = 1
if doneF == 0
    
    load good_data/montelimar_20120301.mat
    
 mm = 1;
    parameters.frac(mm,:) =  parameters.perErr(mm,2:end)./sum(parameters.perErr(mm,2:end))
    
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
    
    dn = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\report\figures\';
    filename = 'montelimar_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    axis normal
    filename = 'montelimar_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr(mm,1)),2),'% of ', num2str(parameters.meanQ(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = 'montelimar_pie'
    print ('-dpng', '-r600',[dn filename])
    
    
    figure
    nr = length(parameters.Q{mm,1})
    discharge = parameters.Q{mm,1};
    
    [n, xout] = hist(parameters.Q{mm,1},20);
    pdfQtest_reBot = (n/nr)/diff(xout(1:2));
    
    hl(1) = plot([parameters.meanQ(mm) parameters.meanQ(mm)], [0 max(pdfQtest_reBot)*1.1], '-.k');
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
    
    filename = 'montelimar_mc'
    print ('-dpng', '-r600', [dn filename])
end


%  Gatineau data
doneF = 0 % not done yet
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
return

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

% 05AE026_20100719_O_Connor_AQ1
doneF = 1
if doneF == 0
    
    load good_data/05AE026_20100719_O_Connor_AQ1.mat
    
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
    filename = '05AE026_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    %axis normal
    filename = '05AE026_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05AE026_pie'
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
    
    filename = '05AE026_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 05BN012_20100720
doneF = 1
if doneF == 0
    
    load good_data/05BN012_20100720.mat
    
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
    filename = '05BN012_20100720'
    print ('-dpng', '-r600', [dn filename])
    
    
    %axis normal
    filename = '05BN012_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05BN012_pie'
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
    
    filename = '05BN012_mc'
    print ('-dpng', '-r600', [dn filename])
end


%05DF001_20080918_CK
doneF = 1
if doneF == 0
    
    load good_data/05DF001_20080918_CK.mat
    
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
    filename = '05DF001_data'
    print ('-dpng', '-r600', [dn filename])
    
    
   % axis normal
    filename = '05DF001_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05DF001_pie'
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
    
    filename = '05DF001_mc'
    print ('-dpng', '-r600', [dn filename])
end

%05GG001_20070724_AQ1 
doneF = 1
if doneF == 0
    
    load good_data/05GG001_20070724_AQ1.mat
    
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
    filename = '05GG001_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    axis normal
    filename = '05GG001_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05GG001_pie'
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
    
    filename = '05GG001_mc'
    print ('-dpng', '-r600', [dn filename])
end

% 05JK002_20031021_MQ1
doneF = 1
if doneF == 0
    
    load good_data/05JK002_20031021_MQ1.mat
    
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
    filename = '05JK002_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    axis normal
    filename = '05JK002_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05JK002_pie'
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
    
    filename = '05JK002_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 05JK007_20070723
doneF = 1
if doneF == 0
    
    load good_data/05JK007_20070723.mat
    
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
    filename = '05JK007_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    %axis normal
    filename = '05JK007_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05JK007_pie'
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
    
    filename = '05JK007_mc'
    print ('-dpng', '-r600', [dn filename])
end

% 05LC001_20070726
doneF = 1
if doneF == 0
    
    load good_data/05LC001_20070726.mat
    
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
    filename = '05LC001_data'
    print ('-dpng', '-r600', [dn filename])
    
    
   % axis normal
    filename = '05LC001_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05LC001_pie'
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
    
    filename = '05LC001_mc'
    print ('-dpng', '-r600', [dn filename])
end



% 05MJ001_20080604
doneF = 1
if doneF == 0
    
    load good_data/05MJ001_20080604.mat
    
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
    filename = '05MJ001_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    %axis normal
    filename = '05MJ001_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05MJ001_pie'
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
    
    filename = '05MJ001_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 05OJ010_20040504
% possibly a non-normal distribution
% it looks positively skewed
doneF = 1
if doneF == 0
    
    load good_data/05OJ010_20040504.mat
    
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
    filename = '05OJ010_data'
    print ('-dpng', '-r600', [dn filename])
    
    
  %  axis normal
    filename = '05OJ010_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '05OJ010_pie'
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
    
    filename = '05OJ010_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 07BE001_20030918
doneF = 1
if doneF == 0
    
    load good_data/07BE001_20030918.mat
    
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
    filename = '07BE001_data'
    print ('-dpng', '-r600', [dn filename])
    
    
   % axis normal
    filename = '07BE001_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '07BE001_pie'
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
    
    filename = '07BE001_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 07BK001_0
doneF = 1
if doneF == 0
    
    load good_data/07BK001_0.mat
    
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
    filename = '07BK001_data'
    print ('-dpng', '-r600', [dn filename])
    
    
    %axis normal
    filename = '07BK001_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '07BK001_pie'
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
    
    filename = '07BK001_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 10FB005_20050506
doneF = 0
if doneF == 0
    
    load good_data/10FB005_20050506.mat
    
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
    filename = '10FB005_data'
    print ('-dpng', '-r600', [dn filename])
    
    
  %  axis normal
    filename = '10FB005_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '10FB005_pie'
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
    
    filename = '10FB005_mc'
    print ('-dpng', '-r600', [dn filename])
end


% 10GC003_20050507
doneF = 0
if doneF == 0
    
    load good_data/10GC003_20050507.mat
    
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
    filename = '10GC003_data'
    print ('-dpng', '-r600', [dn filename])
    
    
  %  axis normal
    filename = '10GC003_data_normal'
    %print ('-dpng', '-r600', [dn filename])
    
    figure
    
    hp = pie(real(parameters.frac(mm,:)));
    [ax,leg] = legend(hp(1:2:end), parameters.name(2:end)');
    set(ax, 'position', [0.1 0.3 .01 .01], 'units', 'normalized')
    text(-0.3, 1, ['Total error (1 \sigma) =  ', num2str(real(parameters.perErr_reBot(mm,1)),2),'% of ', num2str(parameters.meanQ_reBot(mm),4), ' m^3/s'], 'units', 'normalized')
    legend boxoff
    filename = '10GC003_pie'
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
    
    filename = '10GC003_mc'
    print ('-dpng', '-r600', [dn filename])
end
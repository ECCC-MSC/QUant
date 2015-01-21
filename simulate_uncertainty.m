clear all, close all

% Last modified 2014/09/17 by SAM changed transData(mm).velSysErrPct = 1.1
% percent  to  [nrows,ncols]  = size(transData(mm).wVelerr);
%      x = reshape(transData(mm).wVelerr, nrows*ncols,1);
%      transData(mm).velSysErr = nanstd(x)
%%%%%%%%%%%%
% 14 July 2014:
% - fixed a mistake with assigning "top" "bot" and "exponent",
% there was an error in the index assignments
% - added "Qbreakdown" to the list of variables that are saved
% - line 436 changed for nn = nStart to ...

% Add paths to folders that contain useful tools
% ensure that these are the correct paths to the folders

% EXAMPLE OF FOLDER PATH:
addpath 'general/'
addpath 'tools/'
addpath 'version_3d4/'

% specify the file you want to process
file = 11;
% If you don't want to process all the transects that were checked in the
% measurement file, then alter the for loop on line 447

% If you want to alter the number of runs, alter line 1088
% START OF SWITCH LOOP WHERE NAMES OF FILES MUST BE ENTERED
switch file
    
    case 1 % Assiniboine Stn 05MH005 date 20110517 high flow % good, done July 8 2014
        % used in EC report
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\05MH005_20110517_aq1\';
        filename = '05MH005_20110517.mmt';
    case 2 % Montelimar good, done on July 9 2014
        % used in EC report
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\montelimar_fr\AmontUsineMO010312JL12WH2_0\';
        filename = 'amontusinemo010312jl12wh2_0_first_4trans.mmt'; %filename = 'AmontUsineMO010312JL12WH2_0.mmt';
    case 3 % good, done 2014-07-08, gatineau - the largest files so the slowest to run includes GPS data
        % used in EC report
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\Gatineau ADCP Data\';
        %filename = 'gatineau_20031203_ej.mmt' ; % there were about 20
        %transects that were checked off so I made a new mmt file
        filename = 'gatineau_20031203_sm_4_trans.mmt' ;
        
    case 4 % Assiniboine stn 05MH005 low flow % good done as of 10 july 2014
        % used in EC report
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\05MH005_RG600_20130516.AQ1\';
        filename = '05MH005_20130516.AQ1.mmt';
    case 5 % assiniboine % good done as of 10 july 2014
        % used in EC report
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\05MH013_RG1200_20130515.AQ1\' ;
        filename = '05MH013_20130515.AQ1_1.mmt';
    case 6 % graham creek % done on July 13 2014
        % used in EC report
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\Graham_Creek\';
        filename = '02kf015_20130723.mmt';
    case 7 % weir on the welland canal % done 13 juillet 2014
        % used in EC report
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\Weir8 2012 rect channel\weir 8_20120503\';
        filename = 'Weir 8_20120503_0.mmt';
    case 8 % Data from D. Mueller done 14 juillet 2014
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\USGS Data\478.20111130.streampro_adcp\';
        filename = '04097540_478_edited.mmt';
    case 9 % Data from D. Mueller
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\USGS Data\M 18 StreamPro\';
        filename = '04233255_18.mmt';
    case 10 %  Data from D. Mueller
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\USGS Data\RG_1200\';
        filename = '02237700_020711_239.mmt';
    case 11 % test of streampro data for S&T, Greg Bickerton
        pathname = 'M:\My Network Documents\MATLAB\ADCP Uncertainty Analysis\QUant Matlab Code\code\data_in\MacKay_River_005_20140925_0\';
        filename = 'MacKay River 005_20140925_0.mmt';
    case 12 % used in paper submitted to JHE
        pathname = 'M:\My Network Documents\MATLAB\ADCP Uncertainty Analysis\QUant Matlab Code\code\data_in\05AE026_20100719_O_Connor_AQ1\';
        filename = '05ae026_20100719.mmt';
    case 13 % DIDN'T work I get the following error: : An invalid XML character (Unicode: 0x8) was found in the element content of the document.
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05BM002_20090818_O_Connor_AQ1\';
        filename = '05BM002_20090818_AQ1.mmt';
        
    case 14 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05BN012_20100720_AQ1\';
        filename = '05BN012_20100720_0.mmt';
    case 15 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05DF001_20080918_AQ1\';
        filename = '05DF001_20080918_CK.mmt';
    case 16 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05GG001_20070724_AQ1\';
        filename = '05gg001_20070724.mmt';
    case 17 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05JK002_20031021_MQ1\';
        filename = '05JK001_20031021.mmt';
        
    case 18 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05JK007_20070723_AQ1\';
        filename = '05jk007_20070723.mmt';
    case 19 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05LC001_20070726_AQ1\';
        filename = '05lc001_20070726.mmt';
    case 20 % good, done, will use for paper
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05MJ001_20080604\Angus_Colin\';
        filename = '05MJ001_20080604.mmt';
    case 21 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05OC001_20040428\';
        filename = '05oc001_20040428.mmt';
    case 22 % used in paper submitted to JHE
        pathname = 'C:\Users\Stephanie\Documents\environment canada\data\ADCP measurements for Analysis\accreditation\05OJ010_20040504\';
        filename = '05oj010_20040504.mmt';
    case 23 %%% old comment: discharge seems to be >5 cms lower than from winriver
        % about 1/10th of the ensembles are lost, maybe that's why
        % otherwise no hiccups, data processed and saved
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\07be001_20030918_AthabascaAthabasca\';
        filename = '07be001_20030918.mmt';
    case 24 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\07BK001_20100616_CM_0_aq1\';
        filename = '07BK001_0.mmt';
    case 25 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\10FB005_GDL_20050506_AQ1\';
        filename = '10fb005_20050506.mmt';
    case 26 %
        pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\10GC003_20050507_AQ1\';
        filename = '10GC003_20050507.mmt';
        
    case 27 % Assiniboine near Headingly MB Angus_Colin % good done as of 10 july 2014
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05MJ001_20080604\Angus_Colin\';
        %pathname = 'G:\environment canada\data\ADCP measurements for Analysis\accreditation\05MJ001_20080604\Angus_Colin\';
        filename = '05MJ001_20080604.mmt';
    case 28 % Assiniboine near Headingly MB Hardy_Karen % good done as of 10 july 2014
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05MJ001_20080604\Hardy_Karen\';
        filename = '05MJ001_20080604.mmt';
    case 29 % Assiniboine near Headingly MB Hood_Dan
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05MJ001_20080604\Hood_Don\';
        filename = 'Assiniboine River near Headingly.mmt';
    case 30 % Assiniboine near Headingly MB Selinger_Dan
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05MJ001_20080604\Selinger_Dan\';
        filename = '05MJ001_20080604.mmt';
    case 31 %Whitemouth near Whitemouth Colin Angus
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05PH003_20080603\Angus_Colin\';
        filename = '05PH003_20080603.mmt';
    case 32
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05PH003_20080603\Hardy_Karen\';
        filename = '05PH003_20080603.mmt';
    case 33
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05PH003_20080603\Hood_Don\';
        filename = '05PH003_20080603.mmt';
    case 34
        pathname = 'C:\Users\Collin Rennie\Documents\Moore\environment canada\data\ADCP measurements for Analysis\accreditation\05PH003_20080603\Selinger_Dan\';
        filename = '05PH003_20080603.mmt';
        
end

%%

% mmt2mat.m is a file written by D. Mueller that "accepts an
% input file name of an MMT file created by TRDI WinRiver II and parses the
% data from the MMT file's xml structure and stores the data in Matlab
% structures."  It outputs the discharge summary and the input parameters
% used during the measurement and processing.


infile = [pathname filename];

[MMT,MMT_Site_Info,MMT_Transects,MMT_Field_Config, ...
    MMT_Active_Config,MMT_Summary_None,MMT_Summary_BT, ...
    MMT_Summary_GGA,MMT_Summary_VTG,MMT_QAQC,MMT_MB_Transects, ...
    MMT_MB_Field_Config,MMT_MB_Active_Config] = mmt2mat(infile);


%%%%%%%%%%%%%%%%%%%%
% the following code gets the value for WP (number of pings/ensemble), WM
% (water mode), WS (cell size) and WO (number of subpings for mode 12) from
% the mmt file use the following order of priority: user commands, wizard
% commands then fixed commands
%%%%%%%%%%%%%%%%%%%%

% Find WP

clear indUse use
use = [];
skip = 0;
try MMT_Active_Config.User_Commands;
catch em
    if strcmp(em.identifier, 'MATLAB:nonExistentField')
        indUse = [];
        skip = 1; % dont try to use the non existent field
    end
end
clear em
if skip == 0
    indUse = find(strncmp('WP', MMT_Active_Config.User_Commands,2));
end


if  ~isempty(indUse)
    use = char(MMT_Active_Config.User_Commands(indUse(1))); % use the first instance only
    
elseif isempty(use)
    skip = 0;
    try MMT_Active_Config.Wizard_Commands;
    catch em
        if strcmp(em.identifier, 'MATLAB:nonExistentField')
            indUse = [];
            skip = 1;
        end
    end
    clear em
    
    if skip == 0
        indUse = find(strncmp('WP', MMT_Active_Config.Wizard_Commands,2));
    end
    
    
    if ~isempty(indUse)
        use = char(MMT_Active_Config.Wizard_Commands(indUse(1)));
    else
        use = MMT_Active_Config.Fixed_Commands(1,:);
        if ~isempty(use)
            [indUse] = find(strncmp('WP', use,2));
            use = char(use(indUse(1)));
        end
    end
end
numPings = str2double(use(3:end));



% find water mode
clear indUse use
use = [];
skip = 0;

try MMT_Active_Config.User_Commands;
catch em
    if strcmp(em.identifier, 'MATLAB:nonExistentField')
        indUse = [];
        skip = 1;
    end
end
clear em

if skip == 0
    indUse = find(strncmp('WM', MMT_Active_Config.User_Commands,2));
end


if  ~isempty(indUse)
    use = char(MMT_Active_Config.User_Commands(indUse(1))); % use the first instance only
    
elseif isempty(use)
    skip = 0;
    try MMT_Active_Config.Wizard_Commands;
        
    catch em
        if strcmp(em.identifier, 'MATLAB:nonExistentField')
            indUse = [];
            skip = 1;
        end
    end
    clear em
    
    if skip == 0
        
        indUse = find(strncmp('WM', MMT_Active_Config.Wizard_Commands,2));
    end
    
    
    if ~isempty(indUse)
        use = char(MMT_Active_Config.Wizard_Commands(indUse(1)));
    else
        use = MMT_Active_Config.Fixed_Commands(1,:);
        if ~isempty(use)
            [indUse] = find(strncmp('WM', use,2));
            use = char(use(indUse(1)));
        end
    end
end

wMode = str2double(use(3:end));

%find WO command values
subPings = 1;
if wMode == 12
    
    clear indUse use
    use = [];
    skip = 0;
    
    try MMT_Active_Config.User_Commands;
    catch em
        if strcmp(em.identifier, 'MATLAB:nonExistentField')
            indUse = [];
            skip = 1;
        end
    end
    clear em
    if skip == 0
        indUse = find(strncmp('WO', MMT_Active_Config.User_Commands,2));
    end
    
    
    if  ~isempty(indUse)
        use = char(MMT_Active_Config.User_Commands(indUse(1))); % use the first instance only
        
    elseif isempty(use)
        skip = 0;
        try MMT_Active_Config.Wizard_Commands;
        catch em
            if strcmp(em.identifier, 'MATLAB:nonExistentField')
                indUse = [];
                skip = 1;
            end
        end
        clear em
        
        if skip == 0
            
            indUse = find(strncmp('WO', MMT_Active_Config.Wizard_Commands,2));
            
        end
        
        
        if ~isempty(indUse)
            use = char(MMT_Active_Config.Wizard_Commands(indUse(1)));
        else
            use = MMT_Active_Config.Fixed_Commands(1,:);
            if ~isempty(use)
                [indUse] = find(strncmp('WO', use,2));
                if ~isempty(indUse)
                    use = char(use(indUse(1)));
                else
                    use = [];
                    % otherwise, the default (WO1,4) was used, that is only
                    % 1 subping
                end
            end
        end
    end
    if ~isempty(use)
        subPings =  str2double(use(3));
    end
end


% find cell size
clear indUse use
use = [];
skip = 0;

try MMT_Active_Config.User_Commands;
catch em
    if strcmp(em.identifier, 'MATLAB:nonExistentField')
        indUse = [];
        skip = 1;
    end
end
clear em

if skip == 0 % if there was no error message
    indUse = find(strncmp('WS', MMT_Active_Config.User_Commands,2));
end


if  ~isempty(indUse)
    use = char(MMT_Active_Config.User_Commands(indUse(1))); % use the first instance only
    
elseif isempty(use)
    skip = 0;
    try MMT_Active_Config.Wizard_Commands;
    catch em
        if strcmp(em.identifier, 'MATLAB:nonExistentField')
            indUse = [];
            skip = 1;
        end
    end
    clear em
    
    if skip == 0
        indUse = find(strncmp('WS', MMT_Active_Config.Wizard_Commands,2));
        
    end
    
    
    if ~isempty(indUse)
        use = char(MMT_Active_Config.Wizard_Commands(indUse(1)));
    else
        use = MMT_Active_Config.Fixed_Commands(1,:);
        if ~isempty(use)
            [indUse] = find(strncmp('WS', use,2));
            use = char(use(indUse(1)));
        end
    end
end

cellSize = str2double(use(3:end));
% %%%%%%%%%%%%%%%%%%%%%%
%   delete extrap_summary.mat %this deletes old summary information if it exists.
% 
%   uiwait(msgbox('The GUI entitiled "extrap - Beta 3.4" is about to open. Click on "Open .mmt file" and select the file you are currently analyzing. Once the code has run, click on "Save Summary". Next, close the GUI so that discharge calculation can continue','Instructions','modal'));
% 
%   handle1 = extrap3;
% 
%   uiwait(handle1) %once the window

%until I have the curve fitting toolbox, input extrap results manually
load extrap_summary.mat
% 
% extrapDev = abs(str2double(QStats{2,3})); % percentage deviation from the default (Power-Power 1/6th exponent)
% 
extrapDev = 0.64; %this is the value for MacKay_River_005_20140925_0

%% possibleParameters is the list of parameters that will be tested

possibleParameters = {'all', 'bottom depth', 'water velocity', 'bottom velocity', ...
    'temperature', 'salinity', 'draft', 'heading', ...
    'magnetic declination', 'left distance', 'right distance', ...
    'left edge vel', 'right edge vel', 'extrap method top and bot', ...
    'left coef', 'right coef', 'missing ensembles'}; % the list of parameters that will be tested as they affect the uncertainty on the discharge

% parameters is the structure file to which the output of the Monte Carlo
% simulations will be saved
parameters.name = {};
parameters.perErr_reBot = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters));
parameters.perErr_reVTG = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters));
parameters.perErr_reGGA = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters));

parameters.meanQ_reBot = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters));
parameters.meanQ_reVTG = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters));
parameters.meanQ_reGGA = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters));

parameters.frac = NaN*ones(sum(MMT_Transects.Checked),length(possibleParameters)-1);
parameters.Q_reBot = {};
parameters.Q_reGGA = {};
parameters.Q_reVTG = {};
parameters.transno = NaN*ones(1,sum(MMT_Transects.Checked));
parameters.data = {};


doneDisplay = 0;
mm = 0;

nTransects = sum(MMT_Transects.Checked == 1);
% nTransects is the number of transects you want to process, it doesn't
% have to be the total number of transects in the measurement file



nStart = find(MMT_Transects.Checked,1, 'first'); % if you want to start at the first checked transect, should use this as the default
Qbreakdown = cell(length(MMT_Transects.Checked));
for nn = nStart%:length(MMT_Transects.Checked) % cycle through the transects to test
    
    
    % bcs you generally only want to use transects in a mmt file that have
    % checkmarks beside them
    if MMT_Transects.Checked(nn) > 0
        mm = mm+1; % index of the structure parameters
        clearvars -except Qbreakdown FitStats file mm doneDisplay MMT* numPings wMode subPings cellSize nn pathname transData parameters uncert possibleParameters nTransects extrapDev % clear all variables except those listed
        
        %%%%%%%%%%%%%%%%%%%%%
        % If you aren't using extrap3.m to determine the best top and
        % bottom extrapolation methods
        
        %         if MMT_Active_Config.Q_Top_Method(nn) == 0;
        %             top = 'Power';
        %         elseif MMT_Active_Config.Q_Top_Method(nn) == 1;
        %             top = 'Constant';
        %         elseif MMT_Active_Config.Q_Top_Method(nn) == 2;
        %             top = '3-point';
        %         end
        %
        %         if MMT_Active_Config.Q_Bottom_Method(nn) == 0;
        %             bot = 'Power';
        %         elseif MMT_Active_Config.Q_Bottom_Method(nn) == 2;
        %             bot = 'no-slip';
        %         end
        %%%%%%%%%%%%%%%%
        % before July 14 2014 was:
        %         top = FitStats{nn+1,2};
        %         bot = FitStats{nn+1,3};
        %         exponent = FitStats{nn+1,4};
        % after July 14 was:
        top = FitStats{mm+1,2};
        bot = FitStats{mm+1,3};
        exponent = FitStats{mm+1,4};
        
        filename = MMT_Transects.Files{nn} % assign the transect to analyse
        transData(mm) = OriginData_sm(filename, pathname, nn, MMT, MMT_Active_Config, top, bot, exponent); % extract the data using OriginData_sm.m
        
        % So that you plot the data first in order to be able to better choose
        % the uncertainy parameters
        
        % list of arguements: input file, boat velocity reference, (1 or 0, use 1 to
        % create a figure of the velocity data, use 0 to not show a figure
        plot_trans_data(transData(mm), 'bot',1)
        
        if isfinite(max(max(transData(mm).btVel_reVTG))) % if there is gps data this value will be non-zero
            gpsData = 1;
            plot_trans_data(transData(mm), 'gga',1)
            plot_trans_data(transData(mm), 'vtg',1)
        else
            gpsData = 0;
        end
        
        
        
        %%%%%%%%
        % Make and display a table of the values read from the mmt file
        %%%%%%%%
        % the position of the figure is given in pixels of your screen
        
        
        % if you want to change the uncertainty with each transect, comment
        % out the following 'if' and its corresponding 'end'
        
        if doneDisplay == 0 % display the user inputs and default uncertainty values
            set(0,'Units','pixels') ;
            scnsize = get(0,'ScreenSize');
            
            % determine the distance to the right and left banks
            beginDist = MMT_Active_Config.Edge_Begin_Shore_Distance(nn);
            endDist = MMT_Active_Config.Edge_End_Shore_Distance(nn);
            beginL = MMT_Active_Config.Edge_Begin_Left_Bank(nn);
            
            if beginL == 0;
                beginR = 1;
            elseif beginL == 1;
                beginR = 0;
            end
            lDist = beginDist*abs(beginL);
            if lDist == 0
                lDist = endDist;
            end
            rDist = beginDist*beginR;
            if rDist == 0
                rDist = endDist;
            end
            
            
            % Display the draft, distance to banks, etc that were used by
            % the ADCP operator
            f = figure('Position', [0.2*scnsize(3) 0.2*scnsize(4) 0.8*scnsize(3) 0.7*scnsize(4)]);
            cnames = {'Values'};
            rnames = {'draft [m]', 'left dist [m]', 'right dist [m]', 'left edge coefficient', 'right edge coefficient', 'top method', 'bottom method', 'exponent'};
            dat = {MMT_Active_Config.Offsets_Transducer_Depth(nn); ...
                lDist; ...
                rDist; ...
                MMT_Active_Config.Q_Left_Edge_Coeff(nn); ...
                MMT_Active_Config.Q_Right_Edge_Coeff(nn); ...
                top; ...
                bot; ...
                MMT_Active_Config.Q_Power_Curve_Coeff(nn);
                };
            columneditable =  false;
            
            clear ll rr begin*
            
            % table size
            t1x = 0.05;
            t1by = 0.25;
            t1width = 0.35;
            t1height = 0.45;
            t1ty = t1by + t1height;
            
            
            tUserInputs =  uitable('Units', 'normalized',  'Position', [t1x t1by t1width t1height], 'Data', dat, ...  %can start out with % 'Parent', fh,
                'ColumnName', cnames, 'RowName', rnames, 'ColumnEditable', columneditable);
            uicontrol('Parent', f, 'Style', 'text', 'units', 'normalized', 'Position', [t1x t1ty t1width 0.05], 'String', 'Values from MMT file', 'fontsize', 20)
            
            %%%%
            % A table that displays the parameters for the normal distributions
            % of uncertainty that you will add to all parameters
            % THE USER WILL SEE THIS TABLE AND BE PROMPTED TO CHANGE THESE VALUES
            %%%%
            cnames = {'Mean noise', 'stdev of noise'};
            rnames = {'draft [m]',  'left edge distance', 'right edge distance', 'heading [deg]', ...
                'magnetic declination  [deg]', 'temperature (deg C)', 'salinity [ppt]', ...
                'extrap top and bottom [% error]', 'left coeff [% error]', 'right coeff [% error]' };
            % suggestions
            dat = { ...
                0 0.05; ... % draft % don't change this to much smaller, we think it may go up and down 5 cm each way due to choppiness (see EJ's comments on deliverable 1)
                0 0.1*lDist; ... % left dist 10% error
                0 0.1*rDist; ... % right dist 10% error
                0 2; ... % heading
                0 2; ... % magnetic declination
                0 2; ... % temperature
                0 2; ...    % salinity
                0 extrapDev; ... % percent deviation in top and bottom extrapolated discharge from extrap3.m
                0 10; ...   % extrap left coeff percent error
                0 10; ...   % extrap right coeff percent error
                % the value of 10% was obtained by analysing results of Dramais 2011 (Engineering degree thesis) and from the excel worksheet "Error Analysis Spreadsheet Modified 2" that I was given with the project documentation
                
                };
            columneditable =  [ true true];
            columnformat = {'numeric', 'numeric'};
            
            % position of the table
            t2x = 0.45;
            t2by = 0.25;
            t2width = 0.5;
            t2height = 0.5;
            t2ty = t2by + t2height;
            
            uicontrol('Style', 'text', 'units', 'normalized', 'Position', [t2x t2ty t2width 0.05], 'String', 'Uncertainty to test: pls modify values then click continue', 'fontsize', 22)
            
            
            
            tUncert =  uitable('Units', 'normalized',  'Position', [t2x t2by t2width t2height], 'Data', dat, ...
                'ColumnName', cnames, 'RowName', rnames, 'ColumnEditable', columneditable);
            
            uicontrol('Style', 'text', 'units', 'normalized', 'Position', [t2x t2ty t2width 0.05], 'String', 'Uncertainty to test', 'fontsize', 22)
            %display('Please modify the uncertainty parameters and then type return. Do not close the GUI')
            %keyboard
            uicontrol('Position', [20 20 200 40], 'String', 'Continue', ...
                'Callback', 'uiresume(gcbf)');
            uiwait(gcf);
            uncert =  get(tUncert, 'Data');
            doneDisplay = 1;
        end
        
        
        %%
        % the next chunk of code assigns mean and standard deviations for the
        % uncertainty on certain parameters
        transData(mm).ddraft.mean = transData(mm).draft;
        transData(mm).ddraft.meanErr = uncert{1,1};
        transData(mm).ddraft.stdErr = uncert{1,2};
        
        transData(mm).llDist.mean = transData(mm).lDist;
        transData(mm).llDist.meanErr = uncert{2,1};
        transData(mm).llDist.stdErr = uncert{2,2};
        
        transData(mm).rrDist.mean = transData(mm).rDist;
        transData(mm).rrDist.meanErr = uncert{3,1};
        transData(mm).rrDist.stdErr = uncert{3,2};
        
        transData(mm).hheading.mean = transData(mm).heading;
        transData(mm).hheading.meanErr = uncert{4,1};
        transData(mm).hheading.stdErr = uncert{4,2};
        
        transData(mm).mmagDec.mean = transData(mm).magDec;
        transData(mm).mmagDec.meanErr = uncert{5,1};
        transData(mm).mmagDec.stdErr = uncert{5,2};
        
        transData(mm).ttemperature.mean = transData(mm).Sensor.temperature_degc';
        transData(mm).ttemperature.meanErr = uncert{6,1};
        transData(mm).ttemperature.stdErr = uncert{6,2};
        
        transData(mm).ssalinity.mean = transData(mm).salinty;
        transData(mm).ssalinity.meanErr = uncert{7,1};
        transData(mm).ssalinity.stdErr = uncert{7,2};
        
        transData(mm).QextrapTBErrPct = uncert{8,2};
        
        transData(mm).lleftCoef.mean = transData(mm).leftCoef;
        transData(mm).lleftCoef.meanErr = uncert{9,1};
        transData(mm).lleftCoef.stdErr = uncert{9,2};
        
        transData(mm).rrightCoef.mean = transData(mm).rightCoef;
        transData(mm).rrightCoef.meanErr = uncert{10,1};
        transData(mm).rrightCoef.stdErr = uncert{10,2};
        
        %%
        %%%%%%%%%%%%%%
        % velocity uncertainty in m/s taken from Table 5, p. 97 of the
        % Winriver II User Guide (P/N 957-6231-00, February 2007) this is
        % the single ping standard deviation, except for mode 12 for which
        % it is the uncertainty for WO20,4 (20 subpings) these are the
        % values for the recommended cell size 600 kHz: WM1=50 cm;
        % WM12,11,5,8 = 25cm 1200 kHz: WM1 = 25cm, WM12 = 10 cm; WM11 = 5
        % cm, WM5 and WM8 = 10 cm WM12 WO20,4
        %%%%%%%%%%%%%%
        switch wMode
            case 1
                vErr = [0.1362; 0.1364];
                
            case 12
                vErr = [0.0624; 0.0695];
                % for workhorse monitor = [0.3703; 0.1882;] for streampro
                % it is apparently 2 mm/s (see line 379)
                
            case 11
                vErr = [0.0074; 0.0134];
                
            case 5
                vErr = [0.0033; 0.0044];
                
            case 8
                vErr = [0.0334; 0.0515];
        end
        
        display ('The single ping velocity uncertainty from table 5 of the Winriver II user guide (Feb 2007) is used')
        
                
        % divide the single ping uncertainty by the sqrt of the number of
        % pings
        if transData(mm).freq == 600
            transData(mm).velErr = vErr(1)/sqrt(numPings);
            if subPings > 1
                transData(mm).velErr = (vErr(1)/sqrt(subPings))/sqrt(numPings);
            end
            
        elseif transData(mm).freq == 1200
            transData(mm).velErr = vErr(2)/sqrt(numPings);
            if subPings > 1
                transData(mm).velErr = (vErr(2)/sqrt(subPings))/sqrt(numPings);
            end
        elseif transData(mm).freq == 2400 % streampro
            transData(mm).velErr = 2E-3; %m/s
            
            % according to p 84 of the May 2013 version of the streampro
            % manual, the velocity accuracy is 1%, +/- 2mm I assume that
            % this is for the default settings, which for the streampro are
            % WM12, WP6 note, these are the commands that were used for the
            % Graham Creek data
        end
        % prior to 8 aout 2014 was:
        % systematic error that does not depend on number of pings
        %        transData(mm).velSysErrPct = 1.1; % from Oberg 2007
        %%%%%%%%%%%%%%
        [nrows,ncols]  = size(transData(mm).wVelerr);
        x = reshape(transData(mm).wVelerr, nrows*ncols,1);
        transData(mm).velSysErr = nanstd(x);
        
        
        % prior to 8 aout 2014 was:  transData(mm).velBtErrPct = 0.5; % from Oberg 2007
        transData(mm).velBtErr = nanstd(transData(mm).btVel(4,:));
        %transData(mm).velBtErrPct = 0.5; % from Oberg 2007
        transData(mm).depthErrPct = 2; % from Simpson 2001
        
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%
        %compute the discharge using the 3 references: bottom, GGA and VTG
        %%%%%%%%%%%%%%%%%%%%%
        % list of arguements: input file, boat velocity reference, (1 or 0)
        % create a figure of the velocity data (1) or not (0)
        Qnodist_reBot = Discharge_sm(transData(mm), 'bot',1);
        if isfinite(abs(Qnodist_reBot.qIntEns -  Qnodist_reBot.qIntEns2))
            forQintEns_reBot = abs(Qnodist_reBot.qIntEns -  Qnodist_reBot.qIntEns2)/ abs(Qnodist_reBot.qIntEns);
        else
            forQintEns_reBot = 0;
        end
        
        if abs(max(max(transData(mm).btVel_reVTG)))>0 % if there is gps data this value will be non-zero
            gpsData = 1;
            Qnodist_reGGA = Discharge_sm(transData(mm), 'gga',1)
            if isfinite(abs(Qnodist_reGGA.qIntEns -  Qnodist_reGGA.qIntEns2))
                forQintEns_reGGA = abs((Qnodist_reGGA.qIntEns -  Qnodist_reGGA.qIntEns2)/Qnodist_reGGA.qIntEns);
            else
                forQintEns_reGGA = 0;
            end
            
            Qnodist_reVTG = Discharge_sm(transData(mm), 'vtg',1)
            if isfinite(abs(Qnodist_reVTG.qIntEns -  Qnodist_reVTG.qIntEns2))
                forQintEns_reVTG = (abs(Qnodist_reVTG.qIntEns -  Qnodist_reVTG.qIntEns2))/Qnodist_reVTG.qIntEns;
            else
                forQintEns_reVTG = 0;
            end
            
        else
            gpsData = 0;
        end
        
        % first do the simulations including uncertainty on all variables
        % and then do the simulations including only uncertainty on one
        % parameter
        
        % possibleParameters:
        % all, bot depth, wat vel, bot vel, temp, salinity, draft, heading, ...
        % mag dec, L dist, R dist L edge vel, R edge vel, extrap method top and bot, ...
        % left coef, right coef, missing ensembles;
        
        pp = 1;
        while pp <= length(possibleParameters)
              
%             pp = 1; % if only want to assess uncertainty of all parameters
%             while pp == 1
            
            % pp = 8; % if only heading ;
            % while pp == 8
            
            clear Qtest*
            
            parameter = possibleParameters{pp};
            
            switch parameter
                
                case 'all'
                    varyBD = 1; % include error in bottom depth
                    varyWV = 1; % include error that is inherent in velocity measurement
                    varyBV = 1; % include error in bottom velocity estimate
                    varyT = 1; % temperature
                    varyS = 1; % salinity
                    varyD = 1; % draft
                    varyH = 1; % heading
                    varyMG = 1; % magnetic declination
                    varyRD = 1; % left distance
                    varyLD = 1; % right distance
                    varyREV = 1; %variability in velocity used for right edge
                    varyLEV = 1; %variability in velocity used for left edge
                    varyQEXTPTB = 1; % uncertainty due to the extrapolation method at top and bottom
                    varyLEC = 1; % uncertainty due to left edge coefficient
                    varyREC = 1; % uncertainty due to right edge coefficient
                    varyQINTENS = 1; % add uncertainty for the interpolated discharge of missing ensembles
                    
                case 'bottom depth'
                    varyBD = 1;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'water velocity'
                    varyBD = 0;
                    varyWV = 1;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'bottom velocity'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 1;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'temperature'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 1;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'salinity'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 1;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'draft'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 1;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'heading'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 1;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'magnetic declination'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 1;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'right distance'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 1;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'left distance'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 1;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'right edge vel'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 1;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'left edge vel'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 1;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'extrap method top and bot'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 1;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'left coeff'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 1;
                    varyREC = 0;
                    varyQINTENS = 0;
                    
                case 'right coeff'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 1;
                    varyQINTENS = 0;
                case 'missing ensembles'
                    varyBD = 0;
                    varyWV = 0;
                    varyBV = 0;
                    varyT = 0;
                    varyS = 0;
                    varyD = 0;
                    varyH = 0;
                    varyMG = 0;
                    varyRD = 0;
                    varyLD = 0;
                    varyREV = 0;
                    varyLEV = 0;
                    varyQEXTPTB = 0;
                    varyLEC = 0;
                    varyREC = 0;
                    varyQINTENS = 1;
                    
            end
            
            nr = 1000; % number of Monte Carlo realizations - recommend using 1000 for test runs
            display(['number or Monte Carlo realizations = ', num2str(nr),3])
            
            Qtest_reBot = NaN*ones(1,nr);
            if gpsData
                Qtest_reGGA = NaN*ones(1,nr);
                Qtest_reVTG = NaN*ones(1,nr);
            end
            for rr = 1:nr
                clear dataOut Q_re*
                
                dataOut = MC_Data(filename, transData(mm), varyT, varyS, varyD, varyH, varyMG, varyRD, varyLD, varyREC, varyLEC, varyWV, varyBV, varyBD);
                
                Q_reBot = Discharge_sm(dataOut, 'bot', 0, varyLEV, varyREV, varyQEXTPTB, varyQINTENS, forQintEns_reBot);
                Qtest_reBot(rr) = real(Q_reBot.qTot); % havent taken the time to figure out why, but there are sometimes complex
                
                if gpsData % if there is gps data
                    Q_reGGA = Discharge_sm(dataOut, 'gga', 0, varyLEV, varyREV, varyQEXTPTB, varyQINTENS, forQintEns_reGGA);
                    Qtest_reGGA(rr) = real(Q_reGGA.qTot);
                    Q_reVTG = Discharge_sm(dataOut, 'vtg',0, varyLEV, varyREV, varyQEXTPTB, varyQINTENS, forQintEns_reVTG);
                    Qtest_reVTG(rr) = real(Q_reVTG.qTot);
                end
                
            end
            
            
            figure
            [n, xout] = hist(Qtest_reBot,20); % make a 20 binned histogram of the discharge values
            pdfQtest_reBot = (n/nr)/diff(xout(1:2)); % convert the histogram data to a pdf
            % plot the pdf of the discharge
            clf(gcf)
            hl(1) = plot([Qnodist_reBot.qTot Qnodist_reBot.qTot], [0 max(pdfQtest_reBot)], 'k'); % the value when no uncertainty is included
            hold on
            hl(2) = plot(xout, pdfQtest_reBot, 'b.');
            xplot = linspace(min(Qtest_reBot), max(Qtest_reBot));
            hl(3) = plot(xplot, normpdf(xplot,mean(Qtest_reBot), std(Qtest_reBot)), '--r'); % plot a normal distribution for comparison
            legend(hl, '"actual" value', 'Monte Carlo', 'Normal distribution')
            title(['Discharge referenced to bottom track ', parameter, ' is varied'])
            xlabel(['Discharge'])
            ylabel(['Probability Density'])
            ylim([0 max(pdfQtest_reBot)])
            % plot out to +/- 5 std
            xlim([Qnodist_reBot.qTot - 5*std(Qtest_reBot) Qnodist_reBot.qTot + 5*std(Qtest_reBot)])
            
            display(['the uncertainty when ', parameter, ' is varied (in %) = ', num2str(100*std(Qtest_reBot)/mean(Qtest_reBot),3)])
           
            if gpsData
                display(['the uncertainty for gga ref data in percent is ', num2str(100*std(Qtest_reGGA)/mean(Qtest_reGGA),3)])
                
                figure(3), clf(3)
                [n, xout] = hist(Qtest_reVTG,20);
                pdfQtest_reVTG = (n/nr)/diff(xout(1:2));
                clf(3)
                hl(1) = plot([Qnodist_reVTG.qTot Qnodist_reVTG.qTot], [0 1], 'k');
                hold on
                hl(2) = plot(xout, pdfQtest_reVTG);
                xplot = linspace(min(Qtest_reVTG), max(Qtest_reVTG));
                hl(3) = plot(xplot, normpdf(xplot,mean(Qtest_reVTG), std(Qtest_reVTG)), '--r');
                legend(hl, '"actual" value', 'Monte Carlo', 'Normal distribution')
                title(['Discharge referenced to VTG ', parameter, ' is varied'])
                
                figure(4), clf(4)
                [n, xout] = hist(Qtest_reGGA,20);
                pdfQtest_reGGA = (n/nr)/diff(xout(1:2));
                clf(4)
                hl(1) =plot([Qnodist_reGGA.qTot Qnodist_reGGA.qTot], [0 1], 'k');
                hold on
                hl(2) = plot(xout,pdfQtest_reGGA); xplot = linspace(min(Qtest_reGGA), max(Qtest_reGGA));
                hl(3) = plot(xplot, normpdf(xplot,mean(Qtest_reGGA), std(Qtest_reGGA)), '--r');
                legend(hl, '"actual" value', 'Monte Carlo', 'Normal distribution')
                title(['Discharge referenced to GGA ', parameter, ' is varied'])
                
                parameters.meanQ_reVTG(mm,pp) = Qnodist_reVTG.qTot; % equivalent to discharge computed with Winriver II
                parameters.Q_reVTG{mm,pp} = Qtest_reVTG; % probability density function of discharge
                parameters.perErr_reVTG(mm,pp) = 100*std(Qtest_reVTG)/mean(Qtest_reVTG); % relative standard deviation
                
                parameters.meanQ_reGGA(mm,pp) = Qnodist_reGGA.qTot; % equivalent to discharge computed with Winriver II
                parameters.Q_reGGA{mm,pp} = Qtest_reGGA; % probability density function of discharge
                parameters.perErr_reGGA(mm,pp) = 100*std(Qtest_reGGA)/mean(Qtest_reGGA); % relative standard deviation
                
            end
            
            parameters.perErr_reBot(mm,pp) = 100*std(Qtest_reBot)/mean(Qtest_reBot); % relative standard deviation
            parameters.meanQ_reBot(mm,pp) = Qnodist_reBot.qTot; % equivalent to discharge computed with Winriver II
            parameters.Q_reBot{mm,pp} = Qtest_reBot; % probability density function of discharge
            
            
            parameters.name(pp) = {parameter}; % parameters that were tested
            %  pause
            pp = pp + 1;
        end
        parameters.transno(mm) = nn; % transect number
        parameters.data{mm} = transData(mm); % input data
    end
    
    %     parameters.frac_reBot(mm,:) =  parameters.perErr_reBot(mm,2:end)./sum(parameters.perErr_reBot(mm,2:end))
    %     if gpsData
    %         parameters.frac_reVTG(mm,:) =  parameters.perErr_reVTG(mm,2:end)./sum(parameters.perErr_reVTG(mm,2:end))
    %         parameters.frac_reGGA(mm,:) =  parameters.perErr_reGGA(mm,2:end)./sum(parameters.perErr_reGGA(mm,2:end))
    %         figure
    %         pie(real(parameters.perErr_reGGA(mm,2:end)), parameters.name(2:end)')
    %         text(-0.05, -0.05, ['One standard deviation is ', num2str(real(parameters.perErr_reGGA(mm,1)),2),'%'], 'units', 'normalized')
    %         title(['Transect ', num2str(mm)])
    %
    %     end
    % moved by sam on 14 july 2014
    %  parameters.transno(mm) = nn; % transect number
    %  parameters.data{mm} = transData(mm); % input data
    
    Qbreakdown{mm} = Qnodist_reBot;
    
    %   pause
end

% to see how the different factors contribute
% meanQ_reBot = nanmean(parameters.meanQ_reBot(:,1)); % the mean of the discharge from each transect (no uncertainty included)
% 
% 100*nanstd(parameters.meanQ_reBot(:,1))/meanQ_reBot; % percent deviation between transects
% 
% sqrt(nansum(parameters.perErr_reBot(:,1).^2))/ sqrt(sum(parameters.transno>0)); % total error of the measurement from monte carlo analysis

% START OF SWITCH LOOP PERTAINING TO SAVING THE OUTPUT.

% Save format: "save" "filename" "parameters to write to that file"

return
switch file
    
    case 1 % Assiniboine Stn 05MH005 date 20110517 high flow
        
        save assiniboine_05MH005_20110517 parameters % file 1  DONE
    case 2
        save ('G:/work/enviro_can/testing_sm/good_data/montelimar_20120301.mat', 'parameters')        % file 2    DONE
    case 3 % done on July 15 2014 good
        % save ('C:\Users\Collin Rennie\Documents\MATLAB\work/enviro_can/testing_sm/good_data/gatineau_20031203.mat', 'parameters', 'Qbreakdown')
        save ('C:\Users\Collin Rennie\Documents\MATLAB\work/enviro_can/testing_sm/good_data/gatineau_20031203_gps_vary_heading_only.mat', 'parameters', 'Qbreakdown')     % file 3    DONE
    case 4
        save ('G:/work/enviro_can/testing_sm/good_data/assiniboine_05MH005_20130516.mat', 'parameters')    % file 4     DONE
    case 5
        save ('G:/work/enviro_can/testing_sm/good_data/assiniboine_20130515.mat', 'parameters')    % file 5    DONE
    case 6
        save ('C:\Users\Collin Rennie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\graham_creek_20130723_streampro.mat', 'parameters') % file 6
    case 7
        save ('C:\Users\Collin Rennie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\Weir_8_20120503.mat', 'parameters')            % file 7    DONE
    case 8
        save ('C:\Users\Collin Rennie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\Mueller_478_20111130_streampro.mat', 'parameters')
    case 9
        save ('C:\Users\Collin Rennie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\Mueller_M18_streampro.mat', 'parameters')
    case 10
        save ('C:\Users\Collin Rennie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\Mueller_RG1200.mat', 'parameters')
    case 11
        
    case 12 % for the sensititivy analysis for JHE paper, do 30%, 10% and 50% for the edge distance
      %  save ('C:\Users\Stephanie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\05AE026_20100719_edge_err_50pct.mat', 'parameters', 'Qbreakdown')
        % save ('C:\Users\Stephanie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\05AE026_20100719_edge_err_10pct.mat', 'parameters', 'Qbreakdown')
        %save ('C:\Users\Collin Rennie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\05AE026_20100719_O_Connor_AQ1.mat', 'parameters', 'Qbreakdown')
        %save 05BM002_20090818_O_Connor_AQ1 parameters
        save ('M:\My Network Documents\MATLAB\ADCP Uncertainty Analysis\QUant Matlab Code\code\data_out\Testing\05AE026_20100719_O_Connor_AQ1.mat', 'parameters', 'Qbreakdown')
        
    case 14
        save 05BN012_20100720 parameters
    case 15
        save 05DF001_20080918_CK parameters
    case 16
        save 05GG001_20070724_AQ1 parameters
        
    case 17
        save 05JK002_20031021_MQ1 parameters
    case 18
        save 05JK007_20070723 parameters
    case 19
        save 05LC001_20070726 parameters
    case 20
        save 05MJ001_20080604_angus_colin parameters
    case 21
        save 05OC001_20040428 parameters
    case 22
       %save ('C:\Users\Stephanie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\05OJ010_20040504_edge_err_10pct.mat', 'parameters', 'Qbreakdown')
        %save ('C:\Users\Stephanie\Documents\MATLAB\work\enviro_can\testing_sm\good_data\05OJ010_20040504_edge_err_50pct.mat', 'parameters', 'Qbreakdown')
        %save ('C:\Users\Collin Rennie\Documents\MATLAB\work/enviro_can/testing_sm/good_data/05OJ010_20040504.mat', 'parameters', 'Qbreakdown')
    case 23
        save 07be001_20030918 parameters
    case 24
        save 07BK001_0 parameters
    case 25
        save 10FB005_20050506 parameters
    case 26
        
    case 27
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05MJ001_20080604_angus_colin.mat', 'parameters')
    case 28 % modified at the end of the day 14 july 2014
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05MJ001_20080604_hardy_karen.mat', 'parameters', 'Qbreakdown')
    case 29
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05MJ001_20080604_hood_don.mat', 'parameters')
    case 30
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05MJ001_20080604_selinger_dan.mat', 'parameters')
    case 31
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05PH003_20080603_angus_colin.mat', 'parameters')
    case 32
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05PH003_20080603_hardy_karen.mat', 'parameters')
    case 33
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05PH003_20080603_hood_don.mat', 'parameters')
    case 34
        save('C:/Users/Collin Rennie/Documents/MATLAB/work/enviro_can/testing_sm/good_data/05PH003_20080603_selinger_dan.mat', 'parameters')
        
end


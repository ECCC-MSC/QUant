classdef OriginData_sm
    %
    % Class definition of original data. Constructor method reads *r.000,
    % *.pd0, or *.mat files from TRDI and SonTek/YSI ADCPs and stores the
    % data with the following properties. The constructor method is the
    % only method for this class.
    % This code is based on David S. Mueller's OriginData.m
    
    % Last modificaitons and validation 2014/08/08 by SAM, changed
    % velSysErrPct to velSysErr and velBtErrPct to velBtErr
    
    properties
        filename        % name of data file
        pathname        % path to data file
        cellDepth       % depth of each cell, may be computed herein
        cellSize        % size of each cell
        depthEns        % mean depth for each ensemble
        beamDepths      % individual beam depths
        rDist
        lDist
        heading
        pitch
        roll
        startBank
        leftCoef        % coefficient for edge estimate
        rightCoef       % coefficient for edge estimate
        leftNumEns2Avg  % number of cells to use for edge estimates
        rightNumEns2Avg % number of cells to use for edge estimates
        
        topMethod
        botMethod
        exponent
        beamAngle
        magDec
        freq
        
        maxCells        % maximum number of cells in a profile
        numCells        % number of cells above side lobe cutoff
        numEns          % number of ensembles
        idxInvalidCells % index to invalid cells in valid ensembles
        idxInvalidEns   % index to invalid ensembles
        perInvalidCells % percentage of cells in valid ensembles that are invalid
        perInvalidEns   % percentage of ensembles that are invalid
        wVelx       % valid velocities in east or x direction
        wVely       % valid velocities in north or y direction
        wVelz       % valid velocities in the vertical direction
        wVelerr       % valid error velocities
        
        wVelx_reGGA % water velocity in east or x direction relative to gga added by sm
        wVely_reGGA % water velocity in north or y direction relative to gga added by sm
        wVelx_reVTG % water velocity in east or x direction relative to vtg added by sm
        wVely_reVTG % water velocity in north or y direction relative to vtg added by sm
        btVel_reGGA
        btVel_reVTG
        
        btVel           % bottom track velocities (2 x ens, 1-east, 2-north) %modified by EJ to keep all 4 ens
        wtVel           % water track velocities (added by EJ)
        draft           % depth of transducers below water surface
        draftUnits      % units for draft
        validData       % flag to indicate that these are valid data
        ensDeltaTime    % time for each ensemble
        
        %added by EJ for theshold/filter sensitivity analysis
        beamsolBT      % specify 3 or 4 beam solution for BT by indicating 3 or 4 (default WinRiver value is 3)
        beamsolWT      % specify 3 or 4 beam soltuion for WT by indicating 3 or 4 (default WinRiver value is 4)
        cellsAboveSL   % matrix for cells above side lobe
        
        Cfg
        Bt
        Wt
        Sensor
        Surface
        temperature
        salinty
        
        %%%%%%%%%%%%%%%%
        % the next 12 variables are not assigned here, but rather in the
        % code simulate_uncertainty.m
        % all parameters except QextrapErrPct are used to add noise in
        % MC_Data.m
        % QextrapErrPct is used in Discharge_sm
        %%%%%%%%%%%%%%%%%
        ddraft
        rrDist
        llDist
        hheading
        mmagDec
        ttemperature
        ssalinity
        velErr % initialized by not used here
        velBtErr % initialized but not used here
        depthErr % initialized but not used here
        QextrapTBErrPct  % initialized but not used here
        lleftCoef % initialized but not used here
        rrightCoef % initialized but not used here
        
        
    end
    methods
        %==================================================================
        function obj=OriginData_sm(filename,pathname,nn, MMT, MMT_Active_Config, top, bot, exponent);
            
            
            obj.magDec = MMT_Active_Config.Offsets_Magnetic_Variation(nn);
            
            %
            % Constructor method
            %==================================================================
            %
            % If no arguments just create object
            % ----------------------------------
            if nargin>0
                
                %--------------------------------------------------------------
                % Read and prepare data from PD0 file
                %--------------------------------------------------------------
                
                fullName=strcat(pathname,filename);
                
                [Hdr, Inst, Cfg, Sensor, Gps, Wt, Bt, Nmea, Gps2, Surface, AutoMode]=readpd0rrss(fullName);
                
                
                obj.freq = Inst.freq(1); % frequency does not change
                obj.Cfg = Cfg;
                obj.Bt = Bt;
                obj.Wt = Wt;
                obj.Sensor = Sensor;
                obj.Surface = Surface;
                obj.temperature = Sensor.temperature_degc';
                obj.salinty = Sensor.salinity_ppt';
                
                obj.draft = MMT_Active_Config.Offsets_Transducer_Depth(nn);
                
                beginL = MMT_Active_Config.Edge_Begin_Left_Bank(nn);
                beginR = 1 - beginL;
                if  MMT_Active_Config.Edge_Begin_Left_Bank(nn) == 1
                    obj.startBank = 'left';
                else
                    obj.startBank = 'right';
                end
                
                
                
                beginDist = MMT_Active_Config.Edge_Begin_Shore_Distance(nn);
                endDist = MMT_Active_Config.Edge_End_Shore_Distance(nn);
                
                beginL = MMT_Active_Config.Edge_Begin_Left_Bank(nn);
                
                if beginL == 0;
                    beginR = 1;
                elseif beginL == -1;
                    beginL = 0;
                    beginR = 1;
                end
                obj.lDist = beginDist*abs(beginL);
                if obj.lDist == 0
                    obj.lDist = endDist;
                end
                obj.rDist = beginDist*beginR;
                if obj.rDist == 0
                    obj.rDist = endDist;
                end
                
                
                
                obj.leftCoef = MMT_Active_Config.Q_Left_Edge_Coeff(nn);
                obj.rightCoef = MMT_Active_Config.Q_Right_Edge_Coeff(nn);
                obj.leftNumEns2Avg = MMT_Active_Config.Q_Shore_Pings_Avg(nn);
                obj.rightNumEns2Avg = MMT_Active_Config.Q_Shore_Pings_Avg(nn);
                obj.beamAngle = Inst.beamAng(1)'; %until march 27 2014 there was an error here, it was always zero
                % if using extrap3.m to determine top and bot coefficients
                obj.topMethod = top;
                obj.botMethod = bot;
                obj.exponent = exponent;
                % if not using extrap3
                %                 obj.topMethod = MMT_Active_Config.Q_Top_Method(nn);
                %                 obj.botMethod = MMT_Active_Config.Q_Bottom_Method(nn);
                %                 obj.exponent = MMT_Active_Config.Q_Power_Curve_Coeff(nn);
                %
                % Check file validity
                % -------------------
                if isstruct(Hdr)
                    obj.heading=Sensor.heading_deg' + obj.magDec;
                    obj.pitch = Sensor.pitch_deg';
                    obj.roll = Sensor.roll_deg';
                    %
                    % Retrieve data to compute sidelobe
                    % ---------------------------------
                    lag=Cfg.lag_cm'./100;
                    pulseLen=Cfg.xmitPulse_cm'./100;
                    regCellSize=Cfg.ws_cm'./100;
                    regCellSize(regCellSize==0)=nan;
                    obj.beamDepths=Bt.depth_m;
                    obj.beamDepths(obj.beamDepths<0.01)=nan; % Screen bad depths reported as zero
                    cell1Dist=Cfg.distBin1_cm'./100;
                    obj.numEns=size(Wt.vel_mps(:,:,1),2);
                    %
                    % Beam angle is used to ID RiverRay data with variable modes and
                    % lags
                    % --------------------------------------------------------------
                    if Inst.beamAng(1)>21
                        lag(AutoMode.Beam1.mode<2 | AutoMode.Beam1.mode>4)=0;
                    end
                    %
                    % surf* data are to accomodate RiverRay. readpd0rr2 sets these
                    % values to nan when reading Rio Grande or StreamPro data
                    % -----------------------------------------------------------
                    surfCells_idx=find(isnan(Surface.no_cells));
                    noSurfCells=Surface.no_cells';
                    noSurfCells(surfCells_idx)=0;
                    maxSurfCells=nanmax(noSurfCells);
                    surfCellSize=Surface.cell_size_cm'./100;
                    surfCell1Dist=Surface.dist_bin1_cm'./100;
                    numRegCells=size(Wt.vel_mps,1);
                    obj.maxCells=maxSurfCells+numRegCells;
                    %
                    % Beam angle is used to ID RiverRay data with variable modes and
                    % lags
                    % --------------------------------------------------------------
                    %if Inst.beamAng(1)<21
                    %
                    % Compute side lobe interference limit
                    % ------------------------------------
                    lagEffect_m=(lag+pulseLen+regCellSize)./2;
                    depthmin=nanmin(obj.beamDepths);
                    lastcell=depthmin.*cosd(Inst.beamAng(1)')-(lagEffect_m);
                    
                    obj.numCells=max([floor(((lastcell-cell1Dist)./regCellSize)+1); zeros(size(lastcell))],[],1);
                    obj.numCells(obj.numCells>numRegCells)=numRegCells;
                    if nanmax(noSurfCells)>0
                        obj.numCells=obj.numCells+noSurfCells;
                    end
                    
                    %
                    % Create matrix with only cells above sidelobe
                    % ---------------------------------------------
                    obj.cellsAboveSL=ones(obj.maxCells,obj.numEns);
                    for jj=1:obj.numEns
                        for i=obj.numCells(jj)+1:obj.maxCells
                            obj.cellsAboveSL(i,jj)=nan;
                        end
                    end
                    %
                    % Compute bottom and water track for beam velocity data
                    % =====================================================
                    
                    %
                    % Create variables for transformation matrix (instrument to earth)

                    % ----------------------------------------------------------------
                         % see correspondence from D. Mueller on 2015-08-12 on why pitch
            % and roll and not included when going from ship to earth
            % coordinates for data acquired by a StreamPro or a Rio Grande
            % ADCP with Winriver II
                    CH=cosd(obj.heading);
                    SH=sind(obj.heading);
                    P=atand(tand(0).*cosd(0));
                    CP=cosd(P);
                    SP=sind(P);
                    CR=cosd(0);
                    SR=sind(0);
                 
    
                    % Convert ship data to earth data
                    % -------------------------------
                    obj.btVel=nan(4,obj.numEns);
                    for ii=1:obj.numEns
                        trans_ie=[((CH(ii).*CR)+(SH(ii).*SP.*SR)) (SH(ii).*CP) ((CH(ii).*SR)-(SH(ii).*SP.*CR));...
                            ((-1.*SH(ii).*CR)+(CH(ii).*SP.*SR)) (CH(ii).*CP) ((-1.*SH(ii).*SR)-(CH(ii).*SP.*CR));...
                            (-1.*CP.*SR) (SP) (CP.*CR)];
                        
                        obj.btVel(1:3,ii)=trans_ie*Bt.vel_mps(1:3,ii);
                        obj.btVel(4,ii)=Bt.vel_mps(4,ii);
                    end
                    obj.wtVel=nan(obj.maxCells,obj.numEns,4);
                    for ii=1:obj.numEns
                        vel_temp=squeeze(Wt.vel_mps(:,ii,:))';
                        trans_ie=[((CH(ii).*CR)+(SH(ii).*SP.*SR)) (SH(ii).*CP) ((CH(ii).*SR)-(SH(ii).*SP.*CR));...
                            ((-1.*SH(ii).*CR)+(CH(ii).*SP.*SR)) (CH(ii).*CP) ((-1.*SH(ii).*SR)-(CH(ii).*SP.*CR));...
                            (-1.*CP.*SR) (SP) (CP.*CR)];
                        
                        vel_earth_temp=trans_ie*vel_temp(1:3,:);
                        vel_earth_temp(4,:)=vel_temp(4,:);
                        obj.wtVel(:,ii,:)=vel_earth_temp';
                    end

                    obj.btVel_reGGA(1,:) = -Gps.ggaVelE_mps;
                    obj.btVel_reGGA(2,:) = -Gps.ggaVelN_mps;
                    obj.btVel_reVTG(1,:) = -Gps.vtgVelE_mps;
                    obj.btVel_reVTG(2,:) = -Gps.vtgVelN_mps;
                    
                    %
                    % Convert water track velocity to water velocity
                    % ----------------------------------------------
                    wVel(:,:,1)=obj.wtVel(:,:,1)-repmat(obj.btVel(1,:),size(obj.wtVel,1),1);
                    wVel(:,:,2)=obj.wtVel(:,:,2)-repmat(obj.btVel(2,:),size(obj.wtVel,1),1);
                    wVel(:,:,3)=obj.wtVel(:,:,3)-repmat(obj.btVel(3,:),size(obj.wtVel,1),1); %added by EJ
                    wVel(:,:,4)=obj.wtVel(:,:,4); %added by EJ
                    
                    obj.wVelx(:,:)=wVel(1:obj.maxCells,:,1).*obj.cellsAboveSL;
                    obj.wVely(:,:)=wVel(1:obj.maxCells,:,2).*obj.cellsAboveSL;
                    obj.wVelz(:,:)=wVel(1:obj.maxCells,:,3).*obj.cellsAboveSL;
                    obj.wVelerr(:,:)=wVel(1:obj.maxCells,:,4).*obj.cellsAboveSL;
                    
                    %====================
                    % ADDED by SAM 2013-07-17
                    %====================
                    % if using gps as the bottom reference
                    wVel_reGGA(:,:,1)=obj.wtVel(:,:,1)-repmat(obj.btVel_reGGA(1,:),size(obj.wtVel,1),1);
                    wVel_reGGA(:,:,2)=obj.wtVel(:,:,2)-repmat(obj.btVel_reGGA(2,:),size(obj.wtVel,1),1);
                    
                    obj.wVelx_reGGA(:,:)=wVel_reGGA(1:obj.maxCells,:,1).*obj.cellsAboveSL;
                    obj.wVely_reGGA(:,:)=wVel_reGGA(1:obj.maxCells,:,2).*obj.cellsAboveSL;
                    
                    wVel_reVTG(:,:,1)=obj.wtVel(:,:,1)-repmat(obj.btVel_reVTG(1,:),size(obj.wtVel,1),1);
                    wVel_reVTG(:,:,2)=obj.wtVel(:,:,2)-repmat(obj.btVel_reVTG(2,:),size(obj.wtVel,1),1);
                    
                    obj.wVelx_reVTG(:,:)=wVel_reVTG(1:obj.maxCells,:,1).*obj.cellsAboveSL;
                    obj.wVely_reVTG(:,:)=wVel_reVTG(1:obj.maxCells,:,2).*obj.cellsAboveSL;
                    
                    
                    
                    %
                    % Develop distance to center of top cell for each ensemble
                    % --------------------------------------------------------
                    depth=obj.beamDepths';
                    if nanmax(noSurfCells)>0
                        dist_cell1_m=surfCell1Dist;
                        dist_cell1_m(surfCells_idx)=cell1Dist(surfCells_idx);
                    else
                        dist_cell1_m=cell1Dist;
                    end
                    %
                    % Combine cell size and cell range from transducer for both
                    % surface and regular cells
                    % ---------------------------------------------------
                    obj.cellDepth=nan(obj.maxCells,obj.numEns);
                    cellSizeAll=nan(obj.maxCells,obj.numEns);
                    for ii=1:obj.numEns
                        if nanmax(noSurfCells)>0
                            numRegCells=obj.maxCells-noSurfCells(ii);
                        else
                            numRegCells=obj.maxCells;
                        end
                        %
                        % Surface cell are present
                        % ------------------------
                        if nanmax(noSurfCells)>0 && noSurfCells(ii)>0
                            obj.cellDepth(1:noSurfCells(ii),ii)=dist_cell1_m(ii)+(0:surfCellSize(ii):(noSurfCells(ii)-1)*surfCellSize(ii))';
                            obj.cellDepth(noSurfCells(ii)+1:end,ii)=obj.cellDepth(noSurfCells(ii),ii)+0.5*surfCellSize(ii)+0.5.*regCellSize(ii)+(0:regCellSize(ii):(numRegCells-1)*regCellSize(ii));
                            cellSizeAll(1:noSurfCells(ii),ii)=repmat(surfCellSize(ii),noSurfCells(ii),1);
                            cellSizeAll(noSurfCells(ii)+1:end,ii)=repmat(regCellSize(ii),numRegCells,1);
                            %
                            % No surface cells
                            % ----------------
                        else
                            %ii
                            obj.cellDepth(1:numRegCells,ii)=dist_cell1_m(ii)+[0:regCellSize(ii):(numRegCells-1)*regCellSize(ii)];
                            cellSizeAll(1:end,ii)=repmat(regCellSize(ii),numRegCells,1);
                        end
                    end
                    % Compute cell depths from water surface
                    % --------------------------------------
                    obj.cellDepth=obj.cellDepth+obj.draft;
                    %
                    % Compute weighted mean depth
                    % ---------------------------
                    w=1-depth./repmat(nansum(depth,2),1,4);
                    obj.depthEns=nansum((depth.*w)./repmat(nansum(w,2),1,4),2)';
                    obj.depthEns=obj.depthEns+obj.draft;
                    %
                    % Compute time for each ensemble
                    % -------------------------------
                    ensTimeSec=Sensor.time(:,1).*3600+Sensor.time(:,2).*60+Sensor.time(:,3)+Sensor.time(:,4)./100;
                    obj.ensDeltaTime=nan(size(ensTimeSec));
                    idxtime=find(~isnan(ensTimeSec));
                    obj.ensDeltaTime(idxtime(2:end))=nandiff(ensTimeSec(idxtime));
                    idx24hr=find(obj.ensDeltaTime<0);
                    obj.ensDeltaTime(idx24hr)=24.*3600+obj.ensDeltaTime(idx24hr);
                    obj.ensDeltaTime=obj.ensDeltaTime';
                    %
                    % Identify invalid ensembles and invalid cells with valid
                    % ensembles
                    % --------------------------------------------------------
                    obj.idxInvalidEns=[];
                    obj.idxInvalidCells=[];
                    for jj=1:obj.numEns
                        test=find(~isnan(obj.wVelx(:,jj)), 1);
                        %
                        % If entire ensemble has no valid data the ensemble is
                        % marked invalid.
                        % -----------------------------------------------------
                        if ~isempty(test)
                            %
                            % Identify invalid cells above sidelobe
                            % -------------------------------------
                            junk=find(isnan(obj.wVelx(1:obj.numCells(jj),jj)))+(jj-1).*obj.maxCells;
                            obj.idxInvalidCells=[obj.idxInvalidCells; junk];
                        else
                            obj.idxInvalidEns=[obj.idxInvalidEns; jj];
                            
                        end
                    end
                    obj.perInvalidEns=(length(obj.idxInvalidEns)./length(obj.depthEns)).*100;
                    obj.perInvalidCells=(length(obj.idxInvalidCells)./(length(find(~isnan(obj.wVelx)))+length(obj.idxInvalidCells))).*100;
                end
                
                if obj.numEns>0
                    obj.validData=1;
                    obj.cellSize=cellSizeAll;
                    obj.filename=filename;
                end
            end
        end
        
    end
    %==================================================================
    methods
        
        function obj=beamsolSensitivity(obj)
            %
            % defines user specified filtering for data
            %==================================================================
            
            obj.beamsolBT=str2double(inputdlg(['specify 3 or 4 beam solution for BT '],'BT beamsol input',1,{num2str(obj.beamsolBT,4)}))
            obj.beamsolWT=str2double(inputdlg(['specify 3 or 4 beam solution for WT '],'WT beamsol input',1,{num2str(obj.beamsolBT,4)}))
            
            if obj.beamsolBT == 4;
                % Remove all three-beam solutions from bottom track
                % -------------------------------------------------
                idx_3b_bt=find(isnan(obj.btVel(4,:)));
                obj.btVel(1,idx_3b_bt)=nan;
                obj.btVel(2,idx_3b_bt)=nan;
                obj.btVel(3,idx_3b_bt)=nan;
            else
            end
            
            if obj.beamsolWT == 4;
                
                for jj=1:obj.maxCells
                    idx_3b_wt=find(isnan(obj.wVelerr(jj,:)));
                    obj.wVelx(jj,idx_3b_wt)=nan;
                    obj.wVely(jj,idx_3b_wt)=nan;
                    obj.wVelz(jj,idx_3b_wt)=nan;
                end
                
                
            end
            %            % Convert water track velocity to water velocity
            %                     % ----------------------------------------------
            wVel(:,:,1)=obj.wtVel(:,:,1)-repmat(obj.btVel(1,:),size(obj.wtVel,1),1);
            wVel(:,:,2)=obj.wtVel(:,:,2)-repmat(obj.btVel(2,:),size(obj.wtVel,1),1);
            wVel(:,:,3)=obj.wtVel(:,:,3)-repmat(obj.btVel(3,:),size(obj.wtVel,1),1); %added by EJ
            wVel(:,:,4)=obj.wtVel(:,:,4); %added by EJ
            
            obj.wVelx(:,:)=wVel(1:obj.maxCells,:,1).*obj.cellsAboveSL;
            obj.wVely(:,:)=wVel(1:obj.maxCells,:,2).*obj.cellsAboveSL;
            obj.wVelz(:,:)=wVel(1:obj.maxCells,:,3).*obj.cellsAboveSL;
            obj.wVelerr(:,:)=wVel(1:obj.maxCells,:,4).*obj.cellsAboveSL;
            
        end
        
        
    end
    
end

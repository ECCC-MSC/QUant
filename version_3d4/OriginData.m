classdef OriginData
    %
    % Class definition of original data. Constructor method reads *r.000,
    % *.pd0, or *.mat files from TRDI and SonTek/YSI ADCPs and stores the
    % data with the following properties. The constructor method is the
    % only method for this class.
    % David S. Mueller, 2/18/2011
    %
    % Modified 6/17/2011
    % 1) Various changes to handle data. May still be a problem with
    % RiverRay data.
    % 2) Cleaned up variable names and added comments.
    %
    %
    % Last modificaitons / validation 5/15/2012
    
    properties
        filename        % name of data file
        pathname        % path to data file
        cellDepth       % depth of each cell, may be computed herein
        cellSize        % size of each cell
        depthEns        % mean depth for each ensemble
        beamDepths      % individual beam depths 
        beginDist
        endDist
        heading
        startBank
        maxCells        % maximum number of cells in a profile
        numCells        % number of cells above side lobe cutoff
        numEns          % number of ensembles
        idxInvalidCells % index to invalid cells in valid ensembles
        idxInvalidEns   % index to invalid ensembles
        perInvalidCells % percentage of cells in valid ensembles that are invalid
        perInvalidEns   % percentage of ensembles that are invalid
        wVelx       % valid velocities in east or x direction
        wVely       % valid velocities in north or y direction
        btVel           % bottom track velocities (2 x ens, 1-east, 2-north)
        draft           % depth of transducers below water surface
        draftUnits      % units for draft
        validData       % flag to indicate that these are valid data
        ensDeltaTime    % time for each ensemble
    end
    methods
        %==================================================================
        function obj=OriginData(filename,pathname,varargin)
        % 
        % Constructor method
        %==================================================================
                    %
        % If no arguments just create object
        % ----------------------------------
        if nargin>0
            %
            % Determine file type
            % -------------------
            fullName=strcat(pathname,filename);
            if strcmp('000',fullName(end-2:end)) % Old WinRiver raw data
                filetype=1;  
            elseif strcmp('mat',fullName(end-2:end)) % RiverSurveyor Live
                filetype=3;
            elseif strcmp('PD0',fullName(end-2:end))||...
                    strcmp('pd0',fullName(end-2:end))% WinRiver II raw data
                filetype=1;
            end
            %
            % Read data file
            % ==============
            %--------------------------------------------------------------
            % Read and prepare data from RiverSurveyor *.mat file
            %--------------------------------------------------------------
            if filetype==3
                load (fullName)
                idx=find(System.Step==3); % Identify data in transect not edges
                if ~isempty(idx)
                    deltaTime=diff(System.Time);
                    obj.ensDeltaTime=deltaTime(idx)';
                    obj.heading=System.Heading(idx)';
                    if Setup.startEdge>0
                        obj.startBank='Right';
                    else
                        obj.startBank='Left';
                    end
                    %
                    % Determine units used in data file
                    % RiverSurveyor prior to version 1.5 had units in the variable name    
                    % -----------------------------------------------------------------
                    if isfield(Summary,'Depth_ft') % English units
                        %
                        % Convert data to meters
                        % ----------------------
                        unitconv=1./3.281;
                        %
                        % Use only data in transect, not edges
                        % ------------------------------------

                        obj.depthEns=Summary.Depth_ft(idx)'.*unitconv;
                        obj.beamDepths=BottomTrack.BT_Beam_Depth_ft(idx)'.*unitconv;
                        obj.beamDepths(obj.beamDepths<0.01)=nan;
                        obj.btVel=-1.*BottomTrack.Bt_Vel_ft_s(idx,1:4)'.*unitconv;
                        obj.wVelx=squeeze(WaterTrack.Velocity_ft_s(:,1,idx)).*unitconv;
                        obj.wVely=squeeze(WaterTrack.Velocity_ft_s(:,2,idx)).*unitconv;
                        obj.beginDist=Setup.Edges_0__DistanceToBank_ft.*unitconv;
                        obj.endDist=Setup.Edges_1__DistanceToBank_ft.*unitconv;
                        %
                        % Develop array of cell size and centers
                        % --------------------------------------
                        obj.maxCells=size(obj.wVelx,1);  
                        cell_size=System.Cell_Size_ft(idx)'.*unitconv;
                        cellSizeAll=repmat(cell_size,obj.maxCells,1);
                        top_of_cells=System.Cell_Start_ft(idx).*unitconv;
                        obj.cellDepth=(repmat((1:obj.maxCells)',1,size(obj.wVelx,2))-0.5).*cellSizeAll+repmat(top_of_cells',obj.maxCells,1);
                        obj.cellDepth(isnan(obj.wVelx))=nan; % Remove depth values for invalid data
                        obj.draft=Setup.sensorDepth;
                    %
                    % SI Units
                    % --------
                    elseif isfield(Summary,'Depth_m') 
                        %
                        % Data in meters
                        % --------------
                        unitconv=1;
                        %
                        % Use only data in transect, not edges
                        % ------------------------------------
                        obj.depthEns=Summary.Depth_m(idx)'.*unitconv;
                        obj.beamDepths=BottomTrack.BT_Beam_Depth_m(idx)'.*unitconv;
                        obj.beamDepths(obj.beamDepths<0.01)=nan;                    
                        obj.btVel=-1.*BottomTrack.Bt_Vel_m_s(idx,1:4)'.*unitconv;
                        obj.wVelx=squeeze(WaterTrack.Velocity_m_s(:,1,idx)).*unitconv;
                        obj.wVely=squeeze(WaterTrack.Velocity_m_s(:,2,idx)).*unitconv;
                        obj.beginDist=Setup.Edges_0__DistanceToBank_m.*unitconv;
                        obj.endDist=Setup.Edges_1__DistanceToBank_m.*unitconv;
                        %
                        % Develop array of cell size and centers
                        % --------------------------------------
                        obj.maxCells=size(obj.wVelx,1);        
                        cell_size=System.Cell_Size_m(idx)'.*unitconv;
                        cellSizeAll=repmat(cell_size,obj.maxCells,1);
                        top_of_cells=System.Cell_Start_m(idx).*unitconv;
                        obj.cellDepth=(repmat((1:obj.maxCells)',1,size(obj.wVelx,2))-0.5).*cellSizeAll+repmat(top_of_cells',obj.maxCells,1);
                        obj.cellDepth(isnan(obj.wVelx))=nan; % Remove depths for invalid data
                        obj.draft=Setup.sensorDepth;
                    %
                    % For RiverSurveyor Version 1.5 or later data
                    % -------------------------------------------
                    elseif isfield(Summary.Units,'Depth') 
                        %
                        % Setup units conversion
                        % ----------------------
                        if strcmp(Summary.Units.Depth,'ft')
                            unitconv=1./3.281;                         
                        else
                            unitconv=1; 
                        end
                        %
                        % Use only data in transect, not edges
                        % ------------------------------------        
                        obj.depthEns=Summary.Depth(idx)'.*unitconv;
                        obj.beamDepths=BottomTrack.BT_Beam_Depth(idx)'.*unitconv;
                        obj.beamDepths(obj.beamDepths<0.01)=nan;
                        obj.btVel=-1.*BottomTrack.BT_Vel(idx,1:4)'.*unitconv;
                        obj.wVelx=squeeze(WaterTrack.Velocity(:,1,idx)).*unitconv;
                        obj.wVely=squeeze(WaterTrack.Velocity(:,2,idx)).*unitconv;
                        obj.beginDist=Setup.Edges_0__DistanceToBank.*unitconv;
                        obj.endDist=Setup.Edges_1__DistanceToBank.*unitconv;
                        %
                        % Develop array of cell size and centers
                        % --------------------------------------        
                        obj.maxCells=size(obj.wVelx,1);         
                        cell_size=System.Cell_Size(idx)'.*unitconv;
                        cellSizeAll=repmat(cell_size,obj.maxCells,1);
                        top_of_cells=System.Cell_Start(idx).*unitconv;
                        obj.cellDepth=(repmat((1:obj.maxCells)',1,size(obj.wVelx,2))-0.5).*cellSizeAll+repmat(top_of_cells',obj.maxCells,1);
                        obj.cellDepth(isnan(obj.wVelx))=nan; % Remove depths for invalid data
                        obj.draft=Setup.sensorDepth;        
                    end
                    obj.numEns=length(idx);
                else
                    obj.numEns=0;
                end
            %--------------------------------------------------------------
            % Read and prepare data from PD0 file
            %--------------------------------------------------------------
            elseif filetype==1 
                if ~isempty(varargin)
                    obj.draft=varargin{1};
                    obj.beginDist=varargin{2};
                    obj.endDist=varargin{3};
                    obj.startBank=varargin{4};
                end

                %
                % Read file
                % ---------
                fullName=strcat(pathname,filename);
                 [Hdr, Inst, Cfg, Sensor,~, Wt, Bt, ~, ~, Surface, AutoMode]=readpd0rrss(fullName);

                %
                % Check file validity
                % -------------------
                if isstruct(Hdr)
                    obj.heading=Sensor.heading_deg;
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
                        
                   % else
%lagrr=[AutoMode.Beam1.lag_length_cm'; AutoMode.Beam2.lag_length_cm'; AutoMode.Beam3.lag_length_cm'; AutoMode.Beam4.lag_length_cm']./100;
%pulseLenrr=[AutoMode.Beam1.trans_length_cm'; AutoMode.Beam2.trans_length_cm'; AutoMode.Beam3.trans_length_cm'; AutoMode.Beam4.trans_length_cm']./100;
                   % end
                    
                    %
                    % Create matrix with only cells above sidelobe
                    % ---------------------------------------------
                   cellsAboveSL=ones(obj.maxCells,obj.numEns);
                    for j=1:obj.numEns
                        for i=obj.numCells(j)+1:obj.maxCells
                            cellsAboveSL(i,j)=nan;
                        end
                    end
                    %
                    % Compute bottom and water track for beam velocity data
                    % =====================================================
                    % All RiverRay data are collected in beam coordinates
                    % ---------------------------------------------------
                    if strcmp(Cfg.coordSys(1,1:4),'Beam')
                        %
                        % Combine Surface Velocity Bins and Regular Velocity Bins into 
                        % one matrix
                        % -------------------------------------------------------------
                        vel_mps_beam=nan(obj.maxCells,obj.numEns,4);
                        if maxSurfCells>0
                            vel_mps_beam(1:maxSurfCells,:,:)=Surface.vel_mps(1:maxSurfCells,:,:);
                        end
                        for iCell=1:obj.numEns
                            vel_mps_beam(noSurfCells(iCell)+1:noSurfCells(iCell)+numRegCells,iCell,:)=Wt.vel_mps(1:numRegCells,iCell,:);
                        end
                        %
                        % Create transformation matrix (beam to instrument)
                        % -------------------------------------------------
                        a=1./(2.*sind(Inst.beamAng(1)));
                        b=1./(4.*cosd(Inst.beamAng(1)));
                        d=a./sqrt(2);
                        if strcmp(Inst.pat,'Convex')
                            c=1;
                        else
                            c=-1;
                        end
                        trans_bi=[-1.*c.*a c.*a 0 0; 0 0 c.*a -1.*c.*a; b b b b; d d -1.*d -1.*d];
                        %
                        % Create variables for transformation matrix (instrument to earth)
                        % ----------------------------------------------------------------
                        CH=cosd(Sensor.heading_deg);
                        SH=sind(Sensor.heading_deg);
                        P=atand(tand(Sensor.pitch_deg).*cosd(Sensor.roll_deg));
                        CP=cosd(P);
                        SP=sind(P);
                        CR=cosd(Sensor.roll_deg);
                        SR=sind(Sensor.roll_deg);
                        %
                        % Convert beam water data to earth water data
                        % Only 4 beam solutions are valid
                        % NEED TO MOVE THIS TO ANOTHER METHOD
                        % -------------------------------------------
                        wtVel=nan(obj.maxCells,obj.numEns,4);
                        for ii=1:obj.numEns  
                            vel_temp=trans_bi*squeeze(vel_mps_beam(:,ii,:))';   
                            trans_ie=[((CH(ii).*CR(ii))+(SH(ii).*SP(ii).*SR(ii))) (SH(ii).*CP(ii)) ((CH(ii).*SR(ii))-(SH(ii).*SP(ii).*CR(ii)));...
                                    ((-1.*SH(ii).*CR(ii))+(CH(ii).*SP(ii).*SR(ii))) (CH(ii).*CP(ii)) ((-1.*SH(ii).*SR(ii))-(CH(ii).*SP(ii).*CR(ii)));...
                                    (-1.*CP(ii).*SR(ii)) (SP(ii)) (CP(ii).*CR(ii))];
                            vel_earth_temp=trans_ie*vel_temp(1:3,:);
                            vel_earth_temp(4,:)=vel_temp(4,:);
                            wtVel(:,ii,:)=vel_earth_temp';
                        end
                        %
                        % Convert beam bottom track to earth bottom track
                        % Only 4 beam solutions are valid
                        % NEED TO MOVE THIS TO ANOTHER METHOD
                        % -----------------------------------------------
                        obj.btVel=nan(4,obj.numEns);
                        for ii=1:obj.numEns
                            vel_temp=trans_bi*Bt.vel_mps(:,ii);   
                            trans_ie=[((CH(ii).*CR(ii))+(SH(ii).*SP(ii).*SR(ii))) (SH(ii).*CP(ii)) ((CH(ii).*SR(ii))-(SH(ii).*SP(ii).*CR(ii)));...
                                    ((-1.*SH(ii).*CR(ii))+(CH(ii).*SP(ii).*SR(ii))) (CH(ii).*CP(ii)) ((-1.*SH(ii).*SR(ii))-(CH(ii).*SP(ii).*CR(ii)));...
                                    (-1.*CP(ii).*SR(ii)) (SP(ii)) (CP(ii).*CR(ii))];                  
                            vel_earth_temp=trans_ie*vel_temp(1:3);
                            vel_earth_temp(4)=vel_temp(4);
                            obj.btVel(:,ii)=vel_earth_temp';
                        end
                    %
                    % If not beam coordinates (typically Rio or StreamPro
                    % just apply heading.
                    % ---------------------------------------------------
                    else
                        
                        %
                        % Create variables for transformation matrix (instrument to earth)
                        % ----------------------------------------------------------------
                        CH=cosd(obj.heading);
                        SH=sind(obj.heading);
                        P=atand(tand(0).*cosd(0));
                        CP=cosd(P);
                        SP=sind(P);
                        CR=cosd(0);
                        SR=sind(0);
                        %
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
                        wtVel=nan(obj.maxCells,obj.numEns,4);
                        for ii=1:obj.numEns  
                            vel_temp=squeeze(Wt.vel_mps(:,ii,:))';
                            trans_ie=[((CH(ii).*CR)+(SH(ii).*SP.*SR)) (SH(ii).*CP) ((CH(ii).*SR)-(SH(ii).*SP.*CR));...
                                    ((-1.*SH(ii).*CR)+(CH(ii).*SP.*SR)) (CH(ii).*CP) ((-1.*SH(ii).*SR)-(CH(ii).*SP.*CR));...
                                    (-1.*CP.*SR) (SP) (CP.*CR)];
                                                 
                            
                           
                            vel_earth_temp=trans_ie*vel_temp(1:3,:);
                            vel_earth_temp(4,:)=vel_temp(4,:);
                            wtVel(:,ii,:)=vel_earth_temp';
                        end
                    end
                    %
                    % Remove all three-beam solutions from bottom track
                    % -------------------------------------------------
                    idx_3b_bt=find(isnan(obj.btVel(4,:)));
                    obj.btVel(1,idx_3b_bt)=nan;
                    obj.btVel(2,idx_3b_bt)=nan;
                    obj.btVel(3,idx_3b_bt)=nan;
                    %
                    % Remove all three-beam solutions from water track
                    % ------------------------------------------------
                     for jj=1:obj.maxCells
                         idx_3b_wt=find(isnan(wtVel(jj,:,4)));
                         wtVel(jj,idx_3b_wt,1)=nan;
                         wtVel(jj,idx_3b_wt,2)=nan;
                         wtVel(jj,idx_3b_wt,3)=nan;            
                     end
                    %
                    % Convert water track velocity to water velocity
                    % ----------------------------------------------
                    wVel(:,:,1)=wtVel(:,:,1)-repmat(obj.btVel(1,:),size(wtVel,1),1);
                    wVel(:,:,2)=wtVel(:,:,2)-repmat(obj.btVel(2,:),size(wtVel,1),1);
                    obj.wVelx(:,:)=wVel(1:obj.maxCells,:,1).*cellsAboveSL;
                    obj.wVely(:,:)=wVel(1:obj.maxCells,:,2).*cellsAboveSL;
                    %
                    % Correct bottom track vector directions
                    % Raw BT vectors are reference to a stationary
                    % ADCP and must be flipped to reflect boat motion
                    % ------------------------------------------------
                    obj.btVel=-1.*obj.btVel(1:2,:);        
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
                    for j=1:obj.numEns
                        test=find(~isnan(obj.wVelx(:,j)), 1);
                        %
                        % If entire ensemble has no valid data the ensemble is
                        % marked invalid.
                        % -----------------------------------------------------
                        if ~isempty(test)
                           %
                           % Identify invalid cells above sidelobe
                           % -------------------------------------
                           temp=find(isnan(obj.wVelx(1:obj.numCells(j),j)))+(j-1).*obj.maxCells;
                           obj.idxInvalidCells=[obj.idxInvalidCells; temp];
                        else
                            obj.idxInvalidEns=[obj.idxInvalidEns; j];
                        end
                    end
                    obj.perInvalidEns=(length(obj.idxInvalidEns)./length(obj.depthEns)).*100;
                    obj.perInvalidCells=(length(obj.idxInvalidCells)./(length(find(~isnan(obj.wVelx)))+length(obj.idxInvalidCells))).*100;
                end
            end 
            if obj.numEns>0
                obj.validData=1;
                obj.cellSize=cellSizeAll;
                obj.filename=filename;
            end
        end
        end
    end
    
end
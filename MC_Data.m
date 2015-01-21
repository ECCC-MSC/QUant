classdef MC_Data
    % Code is based on Dave Mueller's code OriginData.m (2012/06/05)
    %
    % It takes the output of the class OriginData_sm as input.
    %
    % Modified 2013/09/18 by S.A. Moore in order to add uncertainty to the
    % following variables: adcp draft, heading, magnetic declination, temperature,
    % salinity, distance to start bank, distance to end bank, begin velocity,
    % end velocity
    %
    % Uncertainty related to the extrapolation method has not been added
    % because the user can simply use Dave Mueller's code extrap3.m to see how
    % discharge estimates change if the extrapolation method changes
    % (I used version Beta 3.4 (2012/05/15))
    %
    % uncertainty is added to the various parameters by adding noise with a
    % user specified probability distribution to each parameter.
    
    % The discharge is then calculated as normal, using the code
    % Discharge_sm.m for each realization of each variable.
    %
    %
    % Last modificaitons / validation 2014/08/08 no longer use a systematic
    % velocity error of 1.1%, but do sqrt((dataIn.velErr./Wt.vel_mps).^2 + (dataIn.velSysErr/Wt.vel_mps).^2)
    
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
        startBank
        leftCoef        % coefficient for edge estimate
        rightCoef       % coefficient for edge estimate
        leftNumEns2Avg  % number of cells to use for edge estimates
        rightNumEns2Avg % number of cells to use for edge estimates
        
        topMethod
        botMethod
        exponent
        
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
        T              % temperature
        S              % salinity
        QextrapTBErrPct %percent error in the top and bottom extrapolation, obtained with extrap3.m
        
    end
    methods
        %==================================================================
        function obj=MC_Data(filename, dataIn, varyT, varyS, varyD, varyH, varyMG, varyRD, varyLD, varyREC, varyLEC, varyWV, varyBV, varyBD)
            
            obj.topMethod = dataIn.topMethod;
            obj.botMethod = dataIn.botMethod;
            obj.exponent = dataIn.exponent;
            
            oldDraft = dataIn.ddraft.mean; % the original draft before addition of error, important, do not delete
            obj.draft =  dataIn.ddraft.mean; % initialize, you'll add variation
            obj.numEns = dataIn.numEns;
            oldCellDepth = dataIn.cellDepth;
            obj.depthEns = dataIn.depthEns;
            obj.beamDepths = dataIn.beamDepths;
            obj.ensDeltaTime = dataIn.ensDeltaTime;
            
            obj.startBank = dataIn.startBank;
            obj.lDist = dataIn.llDist.mean;
            obj.rDist = dataIn.rrDist.mean;
            
            obj.leftCoef = dataIn.lleftCoef.mean; % as of Feb 18 I add error to these values
            obj.rightCoef = dataIn.rrightCoef.mean;
            
            obj.heading = dataIn.hheading.mean;
            
            
            obj.leftNumEns2Avg  = dataIn.leftNumEns2Avg;
            obj.rightNumEns2Avg = dataIn.rightNumEns2Avg;
            
            beamAngle = dataIn.beamAngle;
            
            obj.btVel_reGGA = dataIn.btVel_reGGA;
            obj.btVel_reVTG = dataIn.btVel_reVTG;
            
            Cfg = dataIn.Cfg;
            Bt = dataIn.Bt;
            Wt = dataIn.Wt;
            Sensor = dataIn.Sensor;
            
            obj.QextrapTBErrPct = dataIn.QextrapTBErrPct;
            
            
            
            lag = Cfg.lag_cm'./100;
            pulseLen = Cfg.xmitPulse_cm'./100;
            regCellSize = Cfg.ws_cm'./100;
            regCellSize(regCellSize==0)=nan;
            cell1Dist = Cfg.distBin1_cm'./100;
            numRegCells=size(Wt.vel_mps,1);
            obj.maxCells=numRegCells;
            
         
            
            
            surfCells_idx=find(isnan(dataIn.Surface.no_cells));
            noSurfCells=dataIn.Surface.no_cells';
            noSurfCells(surfCells_idx)=0;
            maxSurfCells=nanmax(noSurfCells);
            surfCellSize=dataIn.Surface.cell_size_cm'./100;
            surfCell1Dist=dataIn.Surface.dist_bin1_cm'./100;
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
            lastcell=depthmin.*cosd(beamAngle)-(lagEffect_m);
          
            obj.numCells=max([floor(((lastcell-cell1Dist)./regCellSize)+1); zeros(size(lastcell))],[],1);
            obj.numCells(obj.numCells>numRegCells)=numRegCells;
            if nanmax(noSurfCells)>0
                obj.numCells=obj.numCells+noSurfCells;
            end
            
            
            %++++++++++++++++
            % now modify the parameters
            % take a random sample for each ensemble
            %++++++++++++++++
            
            obj.draft = dataIn.ddraft.mean + varyD*random('norm', dataIn.ddraft.meanErr, dataIn.ddraft.stdErr,1, obj.numEns);
            
            % if error is independent of heading
            obj.heading = dataIn.hheading.mean + varyH*random('norm', dataIn.hheading.meanErr, dataIn.hheading.stdErr,1, obj.numEns) + varyMG*random('norm', dataIn.mmagDec.meanErr, dataIn.mmagDec.stdErr,1, obj.numEns);
            
            %%%%%%%%%%%%%%%%%%%%
            % SAM April 10 2014
            % the uncertainty on the compass heading changes with heading
            % (recall Marsden proceedings form Utah)
            
            % Personal communication from R. Marsden sent April 14, 2014:
            %%%%%%%%%%%%%%%%%
            %             The heading output by a magnetic compass can be described by
            %  
            % Houtput = Htrue + A1*cos(Htrue + delta1) + A2*cos(2*Htrue +delta2)
            %  
            % where 
            %  
            % Htrue is the actual heading the would be reported if there were no errors, 
            %  
            % A1 and A2 are the magnitude of the 'one cycle' and 'two cycle' errors in degrees
            %  
            % delta1 and delta2 are the phases of the two error types. 
            %  
            % For a well calibrated compass, A1 and A2 are on the order of 1 degree or less and the overall compass error is 1degree RMS or less.
            % compass calibration should remove the one cycle and two cycle
            % errors
            % but some error remains (see figure3 of Marsden - HMEM2012-000235).
            % assume that A1 = A2 = .25
            % assume that delta 1 = -15 deg and delta 2 = -25 deg (phase of
            % one and two cycle 
           
            % last used by SAMoore 15 juillet 2014
            % UNCOMMENT ME IF YOU WANT TO USE A HEADING DEPENDENT ERROR
            %obj.heading = dataIn.hheading.mean + varyH*random('norm', 1,1, 1,length(dataIn.hheading.mean)).*cosd(dataIn.hheading.mean - 15) + varyH*random('norm', 1,1, 1,length(dataIn.hheading.mean)).*cosd(2*dataIn.hheading.mean - 25) + varyMG*random('norm', dataIn.mmagDec.meanErr, dataIn.mmagDec.stdErr,1, obj.numEns);
           
   
          %  phis = linspace(0,360,100);
          %  noise =  1.*cosd(phis - 15) + 1*cosd(2*phis - 25); with a
          %  noise of 1 degree on both components
        
            %%%%%%%%%%%%%%%%%%%%%
            
            % recall line 127 of OriginData_sm:
            %  obj.heading=Sensor.heading_deg + obj.magDec;, therefore you
            %  must add the noise from both heading and magnetic
            %  declination
            
            
            c = Sensor.sos_mps'; % this is the value of the sos used to calculate beam depth, size and all velocities
            obj.T = Sensor.temperature_degc';
            obj.S = Sensor.salinity_ppt';
            
            obj.T = obj.T + varyT*random('norm', dataIn.ttemperature.meanErr, dataIn.ttemperature.stdErr,1,obj.numEns);
            obj.S = obj.S + varyS*random('norm', dataIn.ssalinity.meanErr, dataIn.ssalinity.stdErr,1,obj.numEns);
            % recalculate the sound speed for the new temperature and
            % salinity
            c2 = round(1449.2 + 4.6*obj.T - 0.055*obj.T.^2 + 0.00029*obj.T.^3 + (1.34 - 0.01*obj.T).*(obj.S - 35));
            % If no uncertainty, c2 should be equivalent to Sensor.sos_mps
            % NB: Winriver rounds to the nearest one
            
            
            [nbins,nens,nbeams] = size(Wt.vel_mps);
            
            
            obj.depthEns = (c2./c).*(obj.depthEns - oldDraft + obj.draft);
            cell1Dist = (c2./c).*cell1Dist;
            
            mfactor = repmat(c2./c, nbins,1);
            Wt.vel_mps(:,:,1) = mfactor.*Wt.vel_mps(:,:,1);
            Wt.vel_mps(:,:,2) = mfactor.*Wt.vel_mps(:,:,2);
            Wt.vel_mps(:,:,3) = mfactor.*Wt.vel_mps(:,:,3);
            Wt.vel_mps(:,:,4) = mfactor.*Wt.vel_mps(:,:,4);
            
            Bt.vel_mps(1,:) =(c2./c).*Bt.vel_mps(1,:);
            Bt.vel_mps(2,:) = (c2./c).*Bt.vel_mps(2,:);
            Bt.vel_mps(3,:) = (c2./c).*Bt.vel_mps(3,:);
            Bt.vel_mps(4,:) = (c2./c).*Bt.vel_mps(4,:);
            
            
            
            % ADD RANDOM UNCERTAINTY RESULTING FROM ACOUSTIC METHOD
            % error is given in m/s
            % add random error (in cm/s, function of number of pings) and
            % systematic error
            %%%%%%%%%%
            % the error is the root mean square of the relative errors
            % (systematic and ping dependent)
            % last modified 2014/08/08 by SAM
          
            
            Wt.vel_mps = Wt.vel_mps + varyWV*Wt.vel_mps.*random('norm', 0, sqrt((dataIn.velErr./Wt.vel_mps).^2 + (dataIn.velSysErr./Wt.vel_mps).^2), nbins, nens, nbeams) ;% last three numbers give the size     
            Bt.vel_mps = Bt.vel_mps + varyBV.*random('norm', 0, dataIn.velBtErr, nbeams, nens);

            %Wt.vel_mps = Wt.vel_mps + varyWV*Wt.vel_mps.*random('norm', 0, sqrt((dataIn.velErr./Wt.vel_mps).^2 + (dataIn.velSysErrPct/100).^2), nbins, nens, nbeams) ;% last three numbers give the size
           % Bt.vel_mps = Bt.vel_mps + varyBV*Bt.vel_mps.*random('norm', 0, dataIn.velBtErrPct, nbeams, nens);
            
            
            
            % ADD RANDOM UNCERTAINTY TO THE DEPTH,(SEE P 116 OF SIMPSON 2001)
            
            obj.depthEns = obj.depthEns + varyBD*obj.depthEns.*random('norm', 0, dataIn.depthErrPct/100,1,obj.numEns);
            
            
            obj.lDist = obj.lDist + varyLD*random('norm', dataIn.llDist.meanErr, dataIn.llDist.stdErr,1);
            obj.rDist = obj.rDist + varyRD*random('norm', dataIn.rrDist.meanErr, dataIn.rrDist.stdErr,1);
            
            obj.leftCoef = obj.leftCoef + varyLEC*random('norm', obj.leftCoef*dataIn.lleftCoef.meanErr/100, obj.leftCoef*dataIn.lleftCoef.stdErr/100,1);
            obj.rightCoef = obj.rightCoef + varyREC*random('norm', obj.rightCoef*dataIn.rrightCoef.meanErr/100, obj.rightCoef*dataIn.rrightCoef.stdErr/100,1);
            
            
            for ind = 1:obj.numEns
                obj.cellDepth(:,ind) = oldCellDepth(:,ind) - oldDraft + obj.draft(ind);
                obj.beamDepths(:,ind) = obj.beamDepths(:,ind) -oldDraft + obj.draft(ind);
                cell1Dist(ind) = cell1Dist(ind) - oldDraft + obj.draft(ind);
            end
            
            
            %
            % Create matrix with only cells above sidelobe
            % ---------------------------------------------
            obj.cellsAboveSL=ones(obj.maxCells,obj.numEns);
            for j=1:obj.numEns
                for i=obj.numCells(j)+1:obj.maxCells
                    obj.cellsAboveSL(i,j)=nan;
                end
            end
            
            
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
            
            
            
            wVel_reGGA(:,:,1)=obj.wtVel(:,:,1)-repmat(obj.btVel_reGGA(1,:),size(obj.wtVel,1),1);
            wVel_reGGA(:,:,2)=obj.wtVel(:,:,2)-repmat(obj.btVel_reGGA(2,:),size(obj.wtVel,1),1);
            
            obj.wVelx_reGGA(:,:)=wVel_reGGA(1:obj.maxCells,:,1).*obj.cellsAboveSL;
            obj.wVely_reGGA(:,:)=wVel_reGGA(1:obj.maxCells,:,2).*obj.cellsAboveSL;
            
            wVel_reVTG(:,:,1)=obj.wtVel(:,:,1)-repmat(obj.btVel_reVTG(1,:),size(obj.wtVel,1),1);
            wVel_reVTG(:,:,2)=obj.wtVel(:,:,2)-repmat(obj.btVel_reVTG(2,:),size(obj.wtVel,1),1);
            
            obj.wVelx_reVTG(:,:)=wVel_reVTG(1:obj.maxCells,:,1).*obj.cellsAboveSL;
            obj.wVely_reVTG(:,:)=wVel_reVTG(1:obj.maxCells,:,2).*obj.cellsAboveSL;
            
            
            %
            % Combine cell size and cell range from transducer for both
            % surface and regular cells
            % ---------------------------------------------------
            cellSizeAll=nan(obj.maxCells,obj.numEns);
            for ii=1:obj.numEns
                
                numRegCells=obj.maxCells;
                cellSizeAll(1:end,ii)=repmat(regCellSize(ii),numRegCells,1);
                
            end
            
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
            
            
            if obj.numEns>0
                obj.validData=1;
                obj.cellSize=cellSizeAll;
                obj.filename=filename;
            end
        end
        
    end
    
end
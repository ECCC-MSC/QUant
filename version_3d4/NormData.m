classdef NormData
    %
    % Class definition for normalized data. The input data are normalized
    % for relative depth and relative unit discharge or velocity. The
    % constuctor method allows an object to be formed without any data or
    % creates a new object from OriginData or completes the properties of
    % an object that is already of class NormData. The constuctor method
    % also allows only a portion of the data to be used in the
    % normalization process by specifying the data extent.
    % David S. Mueller, 2/18/2011
    %
    % Modified 6/18/2011, dsm
    % 1) Added DisplayNos
    % 2) Cleaned up variable names and added comments
    %
    % Modified 10/17/2011, dsm
    % 3) Added validData index to identify data meeting threshold criteria.
    % Changed threshold from a number to a percent of the median number
    % of valid cells in a normalized cell. This is applied to each transect
    % and to the composite so that the composite will need proportionaly
    % more in each normalized cell since it is compiling all transects.
    %
    % Last modifications / validations 5/15/2012 dsm
    
    properties
        fileName             % name of transect file
        cellDepthNormalized  % normalized depth of cell
        unitNormalized       % normalized discharge or velocity for all depth cells
        unitNormalizedMed    % median of normalized data within 5% partitions
        unitNormalizedNo     % number of data points in each median
        unitNormalizedz      % relative depth for each median (5% increments)
        unitNormalized25     % value for which 25% of normalized values are smaller
        unitNormalized75     % value for which 25% of normalized values are larger
        dataType             % type of data (v, q, V, or Q)
        dataExtent
        validData            % index of median values with point count greater than threshold cutoff
    end
    methods
        %==================================================================
        function obj=NormData(datain,datatype,threshold,varargin)
        %
        % Constructor method
        % If no arguments are provided the object is formed with no data in
        % the properties. If datain is of class OriginData the normalized
        % data are computed from the original data. If datain is of class
        % NormData then the unitNormalized, cellDepthNormalized, and
        % dataType properties must be previously defined and the remaining
        % properties will be computed.
        %==================================================================
            %
            % If no arguments just create object
            % ----------------------------------
            if nargin>0
                %
                % If the data extent is not defined set dataextent to zero
                % to trigger all data to be used.
                % --------------------------------------------------------
                if ~isempty(varargin)
                    dataextent=varargin{1};
                else
                    dataextent=[0 100];
                end
                %
                % Check data class for datain. If class is OriginData
                % process the original data.
                % ---------------------------------------------------
                if isa(datain,'OriginData')
                    %
                    % Assign variables from datain
                    % ----------------------------
                    filename=datain.filename;
                    cellDepth=datain.cellDepth;
                    depthEns=datain.depthEns;
                    maxCells=datain.maxCells;
                    wVelx=datain.wVelx;
                    wVely=datain.wVely;
                    btVel=datain.btVel;
                    %
                    % Compute normalized cell depth by average depth in each ensemble
                    % --------------------------------------------------------------
                    normcellDepth=(cellDepth(1:maxCells,:))./repmat(depthEns,maxCells,1);
                    %
                    % If datatype is discharge compute unit discharge for
                    % each cell
                    % ---------------------------------------------------
                    if strcmp(datatype,'q')
                        %
                        % Compute the cross product for each cell
                        % ---------------------------------------
                        unit=(wVelx.*repmat(btVel(2,:),maxCells,1)-...
                               wVely.*repmat(btVel(1,:),maxCells,1));
                    else
                        %
                        % Compute mean velocity components in each ensemble
                        % -------------------------------------------------
                        wVelMean1=nanmean(wVelx);
                        wVelMean2=nanmean(wVely);
                        %
                        % Compute a unit vector in the mean flow direction for each ensemble
                        % ------------------------------------------------------------------
                        [dir ~]=cart2pol(wVelMean1,wVelMean2);
                        [unitVec(:,1) unitVec(:,2)]=pol2cart(dir,1);
                        %
                        % Compute the velocity magnitude in the direction of the mean velocity
                        % of each ensemble using the dot product
                        % --------------------------------------------------------------------
                        unit=nan(size(wVelx,1),size(btVel,2));
                        for i=1:size(wVelx,1)
                            unit(i,:)=dot([wVelx(i,:); wVely(i,:)],unitVec');
                        end
                    end
                    %
                    % Compute total 
                    % -------------
                    meas=nansum(nansum(unit));
                    %
                    % Adjust to postive value
                    % -------------------------
                    if meas<0
                        unit=unit.*-1;
                    end
                    %
                    % Compute normalize unit values
                    % -----------------------------
                    unitNorm=unit./repmat(abs(nanmean((unit))),size(unit,1),1);
                    %
                    % Apply extents if they have been specified
                    % -----------------------------------------
                    if dataextent(1)~=0 || dataextent(2)~=100
                        %
                        % unit discharge is computed here because the
                        % unitNorm could be based on velocity
                        % -------------------------------------------
                        unit=(wVelx.*repmat(btVel(2,:),maxCells,1)-...
                               wVely.*repmat(btVel(1,:),maxCells,1));
                        unitens=nansum(unit);
                        unittotal=nancumsum(unitens);
                        %
                        % Adjust so total discharge is positive
                        % -------------------------------------
                        if unittotal(end)<0
                            unittotal=unittotal.*-1;
                        end
                        %
                        % Apply extents
                        % -------------
                        unitlower=unittotal(end).*dataextent(1)./100;
                        unitupper=unittotal(end).*dataextent(2)./100;
                        idxextent=find(unittotal>unitlower & unittotal<unitupper);
                        unitNorm=unitNorm(:,idxextent);
                        normcellDepth=normcellDepth(:,idxextent);
                    end
                    %
                    % If whole profile is negative make positive
                    % ------------------------------------------
                    idxNeg1=nan(size(unitNorm,2),1);
                    idxNeg2=nan(size(unitNorm,2),1);
                    for c=1:size(unitNorm,2)
                        idxNeg1(c)=length(find(unitNorm(:,c)<0));
                        idxNeg2(c)=length(find(~isnan(unitNorm(:,c))));
                    end
                    idxNeg=idxNeg1==idxNeg2;
                    unitNorm(:,idxNeg)=unitNorm(:,idxNeg).*-1;  
                %
                % If datain is of class NormData assign variables
                % -----------------------------------------------
                else
                    unitNorm=datain.unitNormalized;
                    normcellDepth=datain.cellDepthNormalized;
                    datatype=datain.dataType;
                    filename=datain.fileName;
                end          
                %
                % Compute median values for profile
                % ---------------------------------
                avgInterval=0:0.05:1;
                unitNormMed=nan(1,length(avgInterval)-1);
                unitNormMedNo=nan(1,length(avgInterval)-1);
                unit25=nan(1,length(avgInterval)-1);
                unit75=nan(1,length(avgInterval)-1);
                avgz=nan(1,length(avgInterval)-1);
                for i=1:length(avgInterval)-1
                    idx=find(normcellDepth>avgInterval(i) & normcellDepth<=avgInterval(i+1));
                    unitNormMed(i)=nanmedian((unitNorm(idx)));
                    unitNormMedNo(i)=sum(~isnan(unitNorm(idx)));
                    unit25(i)=prctile(unitNorm(idx),25);
                    unit75(i)=prctile(unitNorm(idx),75);
                    %avgz(i)=1-nanmean(avgInterval(i:i+1));  
                    avgz(i)=1-nanmean(normcellDepth(idx));
                end
                cutoff=nanmedian(unitNormMedNo(unitNormMedNo>0)).*(threshold./100);
                valid=find(unitNormMedNo>cutoff);
                %
                % Data Stored
                % -----------
                obj.fileName=filename;
                obj.cellDepthNormalized=normcellDepth;
                obj.unitNormalized=unitNorm;
                obj.unitNormalizedMed=unitNormMed;
                obj.unitNormalizedNo=unitNormMedNo;
                obj.unitNormalized25=unit25;
                obj.unitNormalized75=unit75;
                obj.dataType=datatype;
                obj.unitNormalizedz=avgz;
                obj.dataExtent=dataextent;
                obj.validData=valid;
            end
        end 
        %==================================================================
        function plotRaw (normData,h)
        % 
        % plots the unitNormalized data for each depth cell
        %==================================================================
            %
            % Plot data
            % ---------  
            hold (h,'on')
            for i=1:length(normData)
                plot(h,normData(i).unitNormalized,1-normData(i).cellDepthNormalized,'.','MarkerFaceColor',...
                [0.8275 0.8275 0.8275],'MarkerEdgeColor',[0.8275 0.8275 0.8275]);   
            end
            hold off
        end
        %==================================================================
        function plotSurface (normData,h)
        % 
        % plots the unitNormalized data for each depth cell
        %==================================================================
            %
            % Plot data
            % ---------  
            hold (h,'on')
            for i=1:length(normData)
                plot(h,normData(i).unitNormalized(1,:),1-normData(i).cellDepthNormalized(1,:),'og');   
            end
            hold off
        end
        %==================================================================
        function plotMed (normData,h)
        %
        % plots the median values as squares and horizontal lines between
        % the 25 and 75 percentiles for each transect in normData.
        %==================================================================
            %
            % Data used
            % ---------           

            %
            % Process each transect in normData
            % ---------------------------------
            for i=1:length(normData)
                %
                % Assign variables from normData
                % ------------------------------
                unitNormMed=normData(i).unitNormalizedMed;
                avgz=normData(i).unitNormalizedz;
                unitNormNo=normData(i).unitNormalizedNo;
                unit25=normData(i).unitNormalized25;
                unit75=normData(i).unitNormalized75;
                datatype=normData(i).dataType;
                validData=normData(i).validData;
                %
                % If data are from a composite profile use different
                % plotting parameters.
                % --------------------------------------------------
                if strcmp(datatype,'Q') || strcmp(datatype,'V')
                    %
                    % Plot all median values in red
                    % --------------------
                    hold (h,'on')
                    plot(h,unitNormMed,avgz,'sr','MarkerFaceColor','r')
                    %
                    % Plot innerquartile range bars around medians in red
                    % ---------------------------------------------------
                    for j=1:20
                        hold (h,'on')
                        plot(h,[unit25(j) unit75(j)],[avgz(j) avgz(j)],'-r')
                    end
                    %
                    % Plot combined median values that meet in black
                    % -----------------------------------------------
                    hold (h,'on')
                    plot(h,unitNormMed(validData),avgz(validData),'sk','MarkerFaceColor','k')
                    %
                    % Plot innerquartile range bars around median values that meet criteria
                    % in blue
                    % ---------------------------------------------------------------------
                    for k=1:length(validData)
                        hold (h,'on')
                        plot(h,[unit25(validData(k)) unit75(validData(k))],[avgz(validData(k)) avgz(validData(k))],'-k','LineWidth',2)
                    end
                %
                % If data are from individual transects use blue line color
                % ---------------------------------------------------------
                else
                    %
                    % Plot all median values in red
                    % --------------------
                    hold (h,'on')
                    plot(h,unitNormMed,avgz,'sr')
                    %
                    % Plot innerquartile range bars around medians in red
                    % ---------------------------------------------------
                    idxnan=find(~isnan(unit25));
                    for j=1:length(idxnan)
                        hold (h,'on')
                        plot(h,[unit25(idxnan(j)) unit75(idxnan(j))],[avgz(idxnan(j)) avgz(idxnan(j))],'-r')
                    end
                    %
                    % Plot median values that meet thresholds in blue
                    % -----------------------------------------------
                    hold (h,'on')
                    plot(h,unitNormMed(validData),avgz(validData),'sb')
                    %
                    % Plot innerquartile range bars around median values that meet criteria
                    % in blue
                    % ---------------------------------------------------------------------
                    for j=1:length(validData)
                        hold (h,'on')
                        plot(h,[unit25(validData(j)) unit75(validData(j))],[avgz(validData(j)) avgz(validData(j))],'-b')
                    end
                end
            end
            hold off
        end
        %==================================================================
        function plotScale (normData,h)
        %
        % Sets the scale for the plot based on the data plotted.
        %==================================================================
            %
            % Compute the overall minimum value of the 25 percentile and
            % the overall maximum value of the 75 percentile.
            % ----------------------------------------------------------
            minavgall=nan(length(normData),1);
            maxavgall=nan(length(normData),1);
            for i=1:length(normData)
                minavgall(i)=nanmin(normData(i).unitNormalized25(normData(i).validData));
                maxavgall(i)=nanmax(normData(i).unitNormalized75(normData(i).validData));
            end
            %
            % Use percentiles to set axes limits
            % ----------------------------------
            minavg=min(minavgall);
            maxavg=max(maxavgall);
            if minavg>0 && maxavg>0
                minavg=0;
                lower=0;
            else
                lower=minavg.*1.2;
            end
            if maxavg<0 && minavg<0
                upper=0;
            else
                upper=maxavg.*1.2;
            end                
            %
            % Scale axes
            % ----------
            ylim(h,[0 1]);
            xlim(h,[lower upper]);
            box (h,'on')
        end
        %==================================================================
        function plotLabel (normData,h)
        %
        % Label axes based on type of data plotted.
        %==================================================================
            %
            % Set datatype to either dischage (q) or velocity (v)
            % ---------------------------------------------------
            datatype= normData(1).dataType;
            %
            % Label axes
            % ----------
            if strcmpi(datatype,'q')
                xlabel(h,'Normalized Unit Q')
            else
                xlabel(h,'Normalized Velocity')
            end
            ylabel(h,'Normalized Distance from Streambed');
                
        end
        %==================================================================
        function DisplayNos (normData,h)
        %
        % Display number of cells in each median value in text boxes
        % defined by h vector.
        %==================================================================
        tableData=[normData.unitNormalizedz',normData.unitNormalizedNo'];
        set(h,'Data',tableData);
        end
    end
end
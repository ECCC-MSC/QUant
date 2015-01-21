classdef SelectFit < FitData
    %
    % Class definition for data class that contains all of the automated
    % extrapolation method selection information. This inherits all the
    % properties from FitData and addes the automated selection data to it.
    % David S. Mueller 6/17/2011
    %
    %
    % Last modificaitons / validation 5/15/2012
    
    properties
        fitMethod    % User selected method Automatic or Manual
        botMethodAuto % Selected extrapolation for top
        topMethodAuto % Selected extrapolation for bottom
        exponentAuto  % Selected exponent
        topfitr2     % Top fit custom r^2
        topmaxdiff   % Maximum difference between power and 3-pt at top
        botdiff      % Difference between power and no slip at z=0.1
        botrsqr      % Bottom fit r^2
        fitrsqr      % Selected fit of selected power/no slip fit
        nsexponent   % No slip optimized exponent
        ppexponent   % Power power optimized exponent
        topr2
        

    end
    methods
        %==================================================================        
        function obj=SelectFit(normalized,fitMethod,varargin)
        % 
        % Constructor method
        % Requires a single normalized data set from NormData.
        % Accesses FitData with top Power, bottom Power, and Optimized.
        % Includes logic to make automatic selection of best extrapolation
        % method.
        %==================================================================
        obj=obj@FitData();
        if nargin>0
            % Retrieve variables
            % ------------------
            validData=normalized.validData;
            obj=varargin{1};
            ppobj=FitData(normalized,'Power','Power','optimize');
            obj.ppexponent=ppobj.exponent;
            obj.residuals=ppobj.residuals;
            obj.rsqr=ppobj.rsqr;
            obj.exponent95confint=ppobj.exponent95confint;
            obj.fitMethod=fitMethod;
            obj.fileName=normalized.fileName;
            obj.dataType=normalized.dataType;
            % Inherit from FitData
            % --------------------
            if strcmpi(fitMethod,'Automatic')

                % More the 6 cells are required to compute an optimized fit. For
                % fewer than 7 cells the default power/power fit is selected due to
                % lack of sufficient data for a good analysis.
                % -----------------------------------------------------------------
                if length(obj.residuals)> 6

                    % Default fit set to power/power
                    % ------------------------------
                    obj.topMethodAuto='Power';
                    obj.botMethodAuto='Power';

                    % Evaluate difference in data and power fit at water surface
                    % using a linear fit through the top 4 median cells.
                    % ----------------------------------------------------------
                    linfit=regstats(normalized.unitNormalizedMed(validData(1:4)),normalized.unitNormalizedz(validData(1:4)),'linear',{'beta','r','rsquare'});
                    dsmfitr2=1-(sum(linfit.r.^2)./mean(abs(linfit.r)));
                    obj.topfitr2=dsmfitr2;
                    obj.topr2=linfit.rsquare;

                    % Select the best power fit
                    % -------------------------
                    if abs(obj.topfitr2)<0.8 || obj.topr2<0.9 ||...
                       0.1667>obj.exponent95confint(1) && 0.1667<obj.exponent95confint(2)||...
                       ppobj.rsqr<0.8
                        ppobj=FitData(normalized,'Power','Power','manual',0.1667);
                    end

                    % Set default selected exponent to optimized exponent
                    % ---------------------------------------------------
                    obj.exponentAuto=ppobj.exponent;
                    obj.fitrsqr=ppobj.rsqr;                
                    obj.topmaxdiff=ppobj.u(end)-(sum(linfit.beta));

                    % Evaluate the difference at the bottom between power using
                    % the whole profile and power using only the bottom third.
                    % ---------------------------------------------------------
                    nsobj=FitData(normalized,'Constant','No Slip','optimize');
                    obj.nsexponent=nsobj.exponent;
                    obj.botrsqr=nsobj.rsqr;
                    obj.botdiff=ppobj.u(ppobj.z==0.1)-nsobj.u(nsobj.z==0.1);

                    % Begin automatic selection logic
                    % ===============================
                    if ((abs(obj.topfitr2)>0.8 || obj.topr2>0.9) && abs(obj.topmaxdiff)>0.1) ...                        
                         && (obj.topmaxdiff>0 || normalized.unitNormalizedMed(validData(1))-ppobj.u(end)>0.05)...
                         || ((abs(obj.botdiff)>0.1)&& obj.botrsqr>0.6)...
                         || (sign(normalized.unitNormalizedMed(validData(1)))~=sign(normalized.unitNormalizedMed(validData(end))))

                        % No Slip at bottom
                        % -----------------
                        obj.botMethodAuto='No Slip';
                        if nsobj.rsqr>0.8
                            obj.exponentAuto=nsobj.exponent;
                            obj.fitrsqr=nsobj.rsqr;
                        else
                            obj.exponentAuto=0.1667;
                            obj.fitrsqr=nan;
                        end
                        if ~isempty(nsobj.exponent95confint)
                            obj.exponent95confint(1)=nsobj.exponent95confint(1);
                            obj.exponent95confint(2)=nsobj.exponent95confint(2);
                        else
                            obj.exponent95confint(1)=nan;
                            obj.exponent95confint(2)=nan;
                        end
                        obj.topMethodAuto='Constant';
                    else

                        % Power/Power fit
                        % --------------- 
                        if ~isempty(ppobj.exponent95confint)
                            if (0.1667>ppobj.exponent95confint(1) && 0.1667<ppobj.exponent95confint(2)) || ppobj.rsqr<0.8
                                obj.exponentAuto=0.1667;
                            end
                        end
                    end
                else

                    % Insufficient data for anything other than the default
                    % -----------------------------------------------------
                   obj.topMethodAuto='Power';
                   obj.botMethodAuto='Power';
                   obj.exponentAuto=0.1667;
                end
                update=FitData(normalized,obj.topMethodAuto,obj.botMethodAuto,'manual',obj.exponentAuto);            
  
            else
                
                % Manual fit computations and assignments
                % ---------------------------------------
                update=FitData(normalized,obj.topMethod,obj.botMethod,obj.expMethod,obj.exponent); 
                nsobj=FitData(normalized,'Constant','No Slip','optimize');
                obj.nsexponent=nsobj.exponent;

            end
            obj.topMethod=update.topMethod;
            obj.botMethod=update.botMethod;
            obj.exponent=update.exponent;
            obj.coef=update.coef;
            obj.u=update.u;
            obj.z=update.z;
            obj.expMethod=update.expMethod;
            obj.residuals=update.residuals;
        end
        end
        %==================================================================
        function Display(selectData)
        %
        % Display automated selection statistics in a separate window.
        %==================================================================
            f=open ('DisplayFitData.fig');
            h=guihandles(f);
            nrows=size(selectData,2);
            tableData=cell(nrows,8);
            for j=1:nrows
                tableData(j,1)={selectData(j).fileName};
                tableData(j,2)={selectData(j).exponent};
                tableData(j,3)={selectData(j).rsqr};
                tableData{j,4}=selectData(j).exponent95confint(1);
                tableData{j,5}=selectData(j).exponent95confint(2);
                tableData{j,6}=selectData(j).sseTop;
                tableData{j,7}=selectData(j).topmaxdiff;
                tableData{j,8}=selectData(j).toprsquared;
            end
            set(h.uitable1,'Data',tableData)
            drawnow;
        end
        %==================================================================
        function DisplayTable(selectData,h)
        %
        % Displays automated selection statistics in predefined table with 
        % h as handle.
        %==================================================================
            nrows=size(selectData,2);
            tableData=cell(nrows,8);
            for j=1:nrows
                tableData(j,1)={selectData(j).fileName};
                tableData(j,2)={num2str(selectData(j).topMethod,'%8s')};
                tableData{j,3}=num2str(selectData(j).botMethod,'%8sf');
                tableData(j,4)={num2str(selectData(j).exponentAuto,'%6.4f')};
                tableData{j,5}=num2str(selectData(j).fitrsqr,'%6.4f');
                tableData{j,6}=num2str(selectData(j).topmaxdiff,'%6.4f');
                if selectData(j).topr2>0.9
                    tableData{j,7}=num2str(selectData(j).topr2,'%6.4f');  
                else
                    tableData{j,7}=num2str(selectData(j).topfitr2,'%6.4f');
                end
                tableData{j,8}=num2str(selectData(j).botdiff,'%6.4f');
                tableData{j,9}=num2str(selectData(j).botrsqr,'%6.4f');
            end
            set(h,'Data',tableData)
            drawnow;
        end        
    end
end
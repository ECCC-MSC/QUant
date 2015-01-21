classdef FitData
    %
    % Class definition for profile extrapolation fit properties. This class
    % consists of a constructor method and a method to plot the fit as a
    % solid line. 
    % Data required for the constructor method include data of class
    % NormData, threshold for the minimum number of points for a valid
    % median, top extrapolation method, bottom extrapolation method, type
    % of fit, and if a manual fit, the exponent.
    % David S. Mueller, 2/18/2011
    %
    % Modified 6/17/2011, dsm
    % 1) Added fit statistics
    %
    % Modified 10/17/2011, dsm
    % 2) Moved application of threshold criteria to NormData as property
    % validData.
    %
    % Last modificaitons / validation 5/15/2012
    
    properties
        fileName           % name of transect file
        topMethod          % top extrapolation method
        botMethod          % bottom extrapolation method
        coef               % power fit coefficient
        exponent           % power fit exponent
        u                  % fit values of the variable
        z                  % distance from the streambed for fit variable
        expMethod          % method to determine exponent (default, optimize, or manual)
        dataType           % type of data (velocity or unit discharge)
        exponent95confint  % 95% confidence intervals for optimized exponent
        residuals          % residuals from fit
        rsqr               % adjusted r^2 for optimized exponent
    end
    methods
        %==================================================================
        function obj=FitData(normData,top,bot,method,varargin)
        %
        % Constructor method
        % normData but be of call NormData
        % threshold is the minimum number of points in a median value
        % required for the value to be used in the fit process.
        % top is the top extrapolation method (constant or power)
        % bot is the bottom extrapolation method (power or noslip)
        % method is the method used to define the exponent (default, optimize, or
        % manual). default is 1/6.
        % varargin defines the user supplied exponent if the method is
        % manual.
        %==================================================================
        %
        % If no arguments just create object
        % ----------------------------------
        if nargin>0
            %
            % Assign data from normData
            % -------------------------
            unitNormNo=normData.unitNormalizedNo;
            avgz=normData.unitNormalizedz;
            y=normData.unitNormalizedMed;
            idxz=normData.validData;
            zc=nan;
            %
            % Initialize fit boundaries
            % -------------------------
            lowerbnd=[-Inf 0.01];
            upperbnd=[Inf   1];
            %
            % Apply threshold to determine median values that exceed threshold values
            % -----------------------------------------------------------------------
            if ~isempty(idxz)
                idxpower=idxz;
                %
                % Select median values to use in extrapolation based on extrapolation
                % methods selected and create fit output data arrays
                % --------------------------------------------------------------------
                fitcombo=[top bot];
                switch fitcombo
                    case ('PowerPower') 
                        obj.z=0:0.01:1;
                        obj.z=obj.z';
                        zc=nan;
                        uc=nan;
                    case ('ConstantPower') 
                        obj.z=0:0.01:max(avgz(idxz));
                        obj.z=[obj.z' ; nan];
                        zc=max(avgz(idxz))+0.01:0.01:1;
                        zc=zc';
                        uc=repmat(y(idxz(1)),size(zc));
                    case ('ConstantNo Slip') 
                        %
                        % Optimize Constant / No Slip if sufficient cells
                        % are available.
                        % -----------------------------------------------
                        if strcmpi(method,'optimize')
                            idx=idxz(1+end-floor(length(avgz(idxz))/3):end);
                            if length(idx)<4
                                method='default';
                            end
                        %
                        % Compute Constant / No Slip using WinRiver II and
                        % RiverSurveyor Live default cells
                        % ------------------------------------------------
                        else
                            idx=find(avgz(idxz)<=0.2); 
                            if isempty(idx)
                                idx=idxz(end);
                            else
                                idx=idxz(idx);
                            end
                        end
                        %
                        % Configure u and z arrays
                        % ------------------------
                        idxns=idx;
                        obj.z=0:0.01:avgz(idxns(1));
                        obj.z=[obj.z' ; nan];
                        idxpower=idx;
                        zc=max(avgz(idxz))+0.01:0.01:1;
                        zc=zc';
                        uc=repmat(y(idxz(1)),size(zc));
                end
                %
                % Determine exponent
                % ------------------                        
                zfit=avgz(idxpower);
                yfit=y(idxpower);
                % 
                % Check data validity
                % -------------------
                ok_ = isfinite(zfit) & isfinite(yfit);
                if ~all( ok_ )
                    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
                        'Ignoring NaNs and Infs in data' );
                end
                obj.exponent=nan;
                obj.exponent95confint=[nan nan];
                obj.rsqr=nan;
                switch lower(method)
                    case ('manual')
                        obj.exponent=varargin{1};
                        model=['x.^' num2str(obj.exponent)];
                        ft_=fittype({model},'coefficients',{'a1'});
                        fo_ = fitoptions('method','LinearLeastSquares');
                    case ('default')
                        obj.exponent=1./6;
                        model=['x.^' num2str(obj.exponent)];
                        ft_=fittype({model},'coefficients',{'a1'});
                        fo_ = fitoptions('method','LinearLeastSquares');
                    case ('optimize')

                        % Set fit options
                        % ---------------
                        fo_ = fitoptions('method','NonlinearLeastSquares','Lower',lowerbnd,'Upper',upperbnd);                        
                        ft_ = fittype('power1');                        
                        %
                        % Set fit data
                        % ------------
                        strt=yfit(ok_);
                        st_ = [strt(end) 1./6 ];
                        set(fo_,'Startpoint',st_);
                end
                        %
                        % Fit data
                        % --------
                        if length(ok_)>1
                            [cf, gof, ~] = fit(zfit(ok_)',yfit(ok_)',ft_,fo_);
                            %
                            % Extract exponent and confidence intervals from fit
                            % --------------------------------------------------
                            if strcmpi(method,'optimize')
                                obj.exponent=cf.b;
                                if obj.exponent<0.05 
                                    obj.exponent=0.05;
                                end
                            end
                            if ~isempty(confint(cf))
                                exponent95ci=confint(cf);
                                if strcmpi(method,'optimize')
                                    exponent95ci=exponent95ci(:,2);
                                end
                                obj.exponent95confint=exponent95ci;
                                obj.rsqr=gof.rsquare; 
                            else
                                exponent95ci=nan;
                                exponent95ci=nan;
                                obj.exponent95confint=nan;
                                obj.rsqr=nan;  
                            end
                        end
                
                %
                % Fit power curve to appropriate data
                % -----------------------------------
                obj.coef=((obj.exponent+1).*0.05.*nansum(y(idxpower)))./...
                    nansum(((avgz(idxpower)+0.5.*0.05).^(obj.exponent+1))-((avgz(idxpower)-0.5.*0.05).^(obj.exponent+1)));
                %
                % Compute residuals
                % -----------------
                obj.residuals=y(idxpower)-obj.coef.*avgz(idxpower).^(obj.exponent);            
                %
                % Compute values (velocity or discharge) based on exponent and compute
                % coefficient
                % --------------------------------------------------------------------
                obj.u=obj.coef.*obj.z.^(obj.exponent);
                if ~isnan(zc)
                    obj.u=[obj.u ; uc];
                    obj.z=[obj.z ; zc];
                end

                %
                % Assign variables to object properties
                % -------------------------------------
                obj.fileName=normData.fileName;
                obj.topMethod=top;
                obj.botMethod=bot;
                obj.expMethod=method;
                obj.dataType=normData.dataType;
            else
                obj.exponent=nan;
                obj.exponent95confint=[nan nan];
                obj.rsqr=nan;
                obj.fileName=normData.fileName;
                obj.topMethod=top;
                obj.botMethod=bot;
                obj.expMethod=method;
                obj.dataType=normData.dataType;
            end
        end
        end
        %==================================================================
        function plotfit(fitData,h,varargin)
        %
        % If input to the plot method is from a single transect 
        % (dataType=q or v) the line is blue. If the input to the plot method
        % is a composite of multiple transects (dataType=Q or V) then the line
        % is a heavy black line.
        %==================================================================
            %
            % Plot specified extraplation
            % ---------------------------
            
            for i=1:length(fitData)
                if strcmp(fitData(i).dataType,'Q') || strcmp(fitData(i).dataType,'V')
                    hold (h,'on')
                    plot(h,fitData(i).u,fitData(i).z,'-k','LineWidth',2); 
                else
                    hold (h,'on')
                    if ~isempty(varargin)
                        if strcmp(varargin{i},'Left')
                            plot(h,fitData(i).u,fitData(i).z,'-m');
                        else
                            plot(h,fitData(i).u,fitData(i).z,'-b');
                        end
                    else
                        plot(h,fitData(i).u,fitData(i).z,'-b');
                    end
                end
            end
            hold off
        end
    end
end



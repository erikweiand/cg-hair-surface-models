% Calculate water contact angles with errors
% E. Weiand - 01/2022

clear
clc

%% parameters

% file name to be read (.mat from readTraj.m)
fileName = 'tabTraj_nvt_pt1to2.mat';

% binning range and bin size
dAmin = 95;     % [A^2], cylinder surface area increase
Nbins = 500;    % number of radial bins
Nz = 310;       % number of surface-normal bins
dz = 0.5;       % [A], vertical bin size
zMin = 5;     	% [A], minimum surface normal threshold for binning
zMax = 125;     % [A], maximum point for density fitting

% domain dimensions
lx = 357.336;   % [A]
ly = 309.462;   % [A]
lz = 250.0;     % [A]

% timestep control
nMin = 200;         % first timestep
nTimesteps = 100;   % number of timesteps for averaging

% contact angle binning
binSize = 5;        % no. of timesteps per average

% surface-normal point for contact angle measurement
yIntercept = 23.28; % [A]

% switch contact angle measurement (> 90° or < 90°, for values close to
% 90°, visualize the fits and also use calcContactAngle.m first)
isHydrophobic = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computation (do not modify)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data file
load(fileName)

% Avogadro number
Na = 6.02214086e23;     % [1/mol]

% single water unit molar weight
Mw = 0.072;             % [kg/mol]

% radial binning vectors
i = 1:Nbins;
edgesR = sqrt(i.*dAmin/pi);
r = (edgesR(1:end-1) + edgesR(2:end))/2;
lenR = length(r);

% z direction binning
i = 1:Nz;
edgesZ = zMin + dz.*i;
z = (edgesZ(1:end-1) + edgesZ(2:end))/2;
lenZ = length(z);

% bin volume (constant in r direction per definition)
dV = dz .* dAmin .* (edgesR(2) - edgesR(1));

% global binning matrix
beadCount = zeros(lenZ,lenR);
posVar = zeros(nTimesteps,3);

%% cloud shifting and binning
nMax = nMin+nTimesteps-1;

ljCutoff = 9;      % offset to be excluded from yIntercept for binning
zMinFit = yIntercept + ljCutoff; %zMin;
kMin = find(z>zMinFit,1);
zeFitFull = ones(lenZ,nTimesteps).*(-1);
ang = linspace(0,pi,50)';
dxPlot = 10;        % line length for plotting
nBins = floor(nTimesteps/binSize);
contactAngle = zeros(nBins,1);
contactAngleMin = zeros(nBins,1);
contactAngleMax = zeros(nBins,1);
currBin = 1;

beadCountCurr = zeros(lenZ,lenR);

for n=nMin:nMax
    currData = sortedData{n};
    currData = currData(currData(:,3) > 4 & ...
        currData(:,3) < 7,:);   % exclude graphene
    currPos = currData(:,4:6);
    
    xo = 0;
    yo = 0;

    % use configuration that is not separated by period boundaries
    shiftedPosXp = mod(currPos(:,1) + lx/2,lx);
    shiftedPosXm = mod(currPos(:,1) - lx/2,lx);
    varX = [var(shiftedPosXm,1,1) var(currPos(:,1),1,1) ...
        var(shiftedPosXp,1,1)];
    [minVarX, minIndX] = min(varX);
    if minIndX == 1
        currPos(:,1) = shiftedPosXm;
    elseif minIndX ==3
        currPos(:,1) = shiftedPosXp;
    end
    shiftedPosYp = mod(currPos(:,2) + ly/2,ly);
    shiftedPosYm = mod(currPos(:,2) - ly/2,ly);
    varY = [var(shiftedPosYm,1,1) var(currPos(:,2),1,1) ...
        var(shiftedPosYp,1,1)];
    [minVarY, minIndY] = min(varY);
    if minIndY == 1
        currPos(:,2) = shiftedPosYm;
    elseif minIndY ==3
        currPos(:,2) = shiftedPosYp;
    end

    % shift variance-maximized droplet by l/2 to obtain minimum variance
    xMin = [currPos(:,1) currPos(:,2) currPos(:,3)];
    
    % compute position variance in Cartesian directions
    posVar(n,:) = var(xMin,1,1);
    
    % calculate center of mass
    centMass = mean(xMin,1);

    % transform into cylindrical coordinates
    R = sqrt((xMin(:,1)-centMass(1)).^2 + (xMin(:,2)-centMass(2)).^2);

    % perform z and radial wise binning
    for k=1:lenZ
        inPos = R(xMin(:,3) >= edgesZ(k) & xMin(:,3) < edgesZ(k+1));

        if ~isempty(inPos)
            [N,currEdges] = histcounts(inPos,edgesR);

            beadCount(k,:) = beadCount(k,:) + N;
            beadCountCurr(k,:) = beadCountCurr(k,:) + N;
        end
    end
    
    % perform ensemble averaging
    if mod(n-nMin+1,binSize) == 0
        zeFitCurr = ones(lenZ,1).*(-1);
        
        for k=1:lenZ
            if k > kMin
                densCurr = (beadCountCurr(k,:) ./ binSize ./ Na .* Mw) ...
                    ./ (dV*(1e-10)^3);   % mass density

                [xData, yData] = prepareCurveData(r, densCurr);

                % Set up fittype and options.
                ft = fittype( '0.5*a*(1-tanh(2*(x-b)/10))', ...
                    'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.StartPoint = [0.278498218867048 0.546881519204984];

                % Fit model to data.
                [fitresult, gof] = fit( xData, yData, ft, opts );

                coeffs = coeffvalues(fitresult);
                densFitCurr = coeffs(1);
                zeFitCurr(k) = coeffs(2);
            end
        end
            
        ind = zeFitCurr(:) > 10 & (z < zMax)';
        
        [xData, yData] = prepareCurveData(z(ind)', zeFitCurr(ind));
        ft = fittype( 'sqrt(abs(b-(x-a))*(b+(x-a)))', ...
            'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares', ...
            'Lower', [-Inf 60] );
        opts.Display = 'Off';
        opts.StartPoint = [0 200];
        opts.Robust = 'LAR';

        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        coeffs = coeffvalues(fitresult);
        xhat = coeffs(1);
        yhat = 0;
        rhat = coeffs(2);

        % 95 percent conf bounds
        confBounds = confint(fitresult);
        xhatMin = confBounds(1,1);
        xhatMax = confBounds(2,1);
        rhatMin = confBounds(1,2);
        rhatMax = confBounds(2,2);

        % current contact angle
        xIntercept = sqrt(rhat^2 - (xhat-yIntercept)^2);
        xInterceptMin = sqrt(rhatMin^2 - (xhatMin-yIntercept)^2);
        xInterceptMax = sqrt(rhatMax^2 - (xhatMax-yIntercept)^2);
        if isHydrophobic == true
            contactAngle(currBin) = 90 + acos(xIntercept/rhat)*180/pi;
            contactAngleMin(currBin) = 90 + ...
                acos(xInterceptMin/rhatMin)*180/pi;
            contactAngleMax(currBin) = 90 + ...
                acos(xInterceptMax/rhatMax)*180/pi;
        else
            contactAngle(currBin) = 90 - acos(xIntercept/rhat)*180/pi;
            contactAngleMin(currBin) = 90 - ...
                acos(xInterceptMin/rhatMin)*180/pi;
            contactAngleMax(currBin) = 90 - ...
                acos(xInterceptMax/rhatMax)*180/pi;
        end
        currBin = currBin + 1;
        
        beadCountCurr = zeros(lenZ,lenR);
    end
end

% calculate water density profile in surface-normal direction
densWater = (beadCount ./ Na .* Mw) ./ nTimesteps ./ (dV*(1e-10)^3);

%% mean properties
meanZeFit = mean(zeFitFull,2);
meanContactAngle_avCircles = mean(contactAngle);
minContactAngle = min(contactAngleMin);
maxContactAngle = max(contactAngleMax);
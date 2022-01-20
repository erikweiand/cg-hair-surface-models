% Calculate a single, average liquid contact angle
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

% surface-normal point for contact angle measurement
yIntercept = 23.28; % [A]

% switch contact angle measurement (> 90° or < 90°, use visualizaton to
% determine the flag value)
isHydrophobic = true;

% visual center of mass debugging; create a "CheckMin/" subdir if true
visualPosDebug = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computation (do not modify)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data file
load(fileName)

% Avogadro number
Na = 6.02214086e23;     % [1/mol]

% single water bead molar weight
Mw = 0.024;             % [kg/mol]

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
for n=nMin:nMax
    currData = sortedData{n};
    currData = currData(currData(:,3) > 4 & ...
        currData(:,3) < 7,:);   % exclude graphene
    currPos = currData(:,4:6);

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
    
    % plot variance-minimized bead map for counter-checking
    if visualPosDebug == true
        f = figure('visible','off');
        scatter(xMin(:,1),xMin(:,2),1,'.')
        hold on
        scatter(centMass(1),centMass(2),5,'v','MarkerFaceColor','r')
        axis equal
        saveas(f,['CheckMin/t_' num2str(n) '_checkMin'],'png')
        close(f)
    end

    % transform into cylindrical coordinates
    R = sqrt((xMin(:,1)-centMass(1)).^2 + (xMin(:,2)-centMass(2)).^2);

    % perform z and radial wise binning
    for k=1:lenZ
        inPos = R(xMin(:,3) >= edgesZ(k) & xMin(:,3) < edgesZ(k+1));

        if ~isempty(inPos)
            [N,currEdges] = histcounts(inPos,edgesR);

            beadCount(k,:) = beadCount(k,:) + N;
        end
    end
end

densWater = (beadCount ./ Na .* Mw) ./ nTimesteps ./ (dV*(1e-10)^3);

%% perform z wise equimolar distance points from sigmoid fitting
ljCutoff = 9;      % offset from surface for density fitting
zMinFit = yIntercept + ljCutoff; %zMin;
kMin = find(z>zMinFit,1);
nBeadsZ = sum(beadCount,2);
densFit = zeros(lenZ,1);
zeFit = ones(lenZ,1).*(-1);
for k=kMin:lenZ
    if nBeadsZ(k) > 100
        [xData, yData] = prepareCurveData(r, densWater(k,:));

        % Set up fittype and options.
        ft = fittype( '0.5*a*(1-tanh(2*(x-b)/10))', ...
            'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.278498218867048 0.546881519204984];

        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        
        coeffs = coeffvalues(fitresult);
        densFit(k) = coeffs(1);
        zeFit(k) = coeffs(2);
    end
end

% circle fit
ang = linspace(0,pi,50)';
ind = zeFit > 10 & (z < zMax)';
c = [z(ind)' zeFit(ind) ones(length(z(ind)'), 1)] ...
    \-(z(ind)'.^2+zeFit(ind).^2);   %least squares fit
xhat = -c(1)/2;
yhat = -c(2)/2;
rhat = sqrt(xhat^2+yhat^2-c(3));

%% plotting

% plot density map
f=figure(1);
f.Color = "white";
set(f,'position',[0,0,900,600])
x = linspace(0,r(end),300);
y = z;
[X,Y] = meshgrid(x,y);
[X2,Y2] = meshgrid(r,z);
Vq = interp2(X2,Y2,densWater,X,Y);
imagesc(gca,[0 r(end)],[zMin z(end)],Vq)
set(gca,'YDir','normal')
set(gca,'FontSize',14)
hold on
axis equal
colormap(flipud(bone(256)))
h = colorbar;
ylabel(h, '\rho [kg/m^3]','FontSize',18)

% plot scatter and circular fit
plot(rhat*sin(ang)+yhat,rhat*cos(ang)+xhat,'r','linewidth',2)
scatter(zeFit(ind),z(ind),10,'o','MarkerFaceColor','b')
grid on
box on
Ang = char(197);
xlabel(['r [' Ang ']'],'FontSize',18)
ylabel(['z [' Ang ']'],'FontSize',18)
axis equal
xlim([0 140])
ylim([0 150])

% calculate contact angle
dxPlot = 10;        % line length for plotting
xIntercept = sqrt(rhat^2 - (xhat-yIntercept)^2);
if isHydrophobic == true
    contactAngle = 90 + acos(xIntercept/rhat)*180/pi;
else
    contactAngle = 90 - acos(xIntercept/rhat)*180/pi;
end
dyPlot = tan((180-contactAngle)/180*pi);
plot([xIntercept+yhat-dxPlot xIntercept+yhat+dxPlot], ...
    [-dxPlot*dyPlot+yIntercept dxPlot*dyPlot+yIntercept],'b','LineWidth',2)

% add label
labelText = ['$\varphi=' num2str(contactAngle,'%10.2f') '^{\circ}$'];
annotation(f,'textbox',...
    [0.5 0.2 0.25 0.05],...
    'String',labelText,...
    'FitBoxToText','off','EdgeColor','none', ...
    'Interpreter','latex','FontSize',18);

%% print to figure
% print(f,'contactAngle_0p65nm_water_0p25dam_C1','-dpng','-r300')
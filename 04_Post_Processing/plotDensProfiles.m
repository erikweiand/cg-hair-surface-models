% Plot monolayer and water density profiles
% E. Weiand - 01/2022

clear
clc

%% parameters

% file name to be read (.mat from readTraj.m)
fileName = 'tabTraj_nvt_pt1to2.mat';

% binning range and bin size
Nz = 600;     % number of surface-normal bins
dz = 0.25;    % surface-normal bin size
zMin = 0;     % [A], minimum surface-normal coordinate for binning

% domain dimensions
lx = 357.336;   % [A]
ly = 309.462;   % [A]
lz = 250.0;     % [A]

% timestep control
nMin = 200;         % first timestep
nTimesteps = 100;   % number of timesteps for averaging

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computation (do not modify)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data file
load(fileName)

% Avogadro number
Na = 6.02214086e23;     % [1/mol]

% single hex bead molar weight
Mw = 0.072;             % [kg/mol]

% single fatty acid bead molar weight
Mc = 0.072;             % [kg/mol]
Mt = 0.090;             % [kg/mol]


%% cloud shifting and binning
nMax = nMin+nTimesteps-1;

% z direction binning
i = 1:Nz;
edgesZ = zMin + dz.*i;
z = (edgesZ(1:end-1) + edgesZ(2:end))/2;
lenZ = length(z);

% bin volume (constant in r direction per definition)
dV = dz .* lx .* ly;

% global binning matrix
beadCountWater = zeros(lenZ,1);
beadCountFA = zeros(lenZ,1);
beadCountFAterm = zeros(lenZ,1);

for n=nMin:nMax
    currData = sortedData{n};
    currDataWater = currData(currData(:,3) > 4 & currData(:,3) < 7,:);
    currDataFA = currData(currData(:,3) < 4 | currData(:,3) == 11,:);
    currDataFAterm = currData(currData(:,3) == 10,:);
    currPosWater = currDataWater(:,4:6);
    currPosFA = currDataFA(:,4:6);
    currPosFAterm = currDataFAterm(:,4:6);
    
    for k=1:lenZ
        inPosWater = currPosWater(:,3) >= edgesZ(k) & ...
            currPosWater(:,3) < edgesZ(k+1);
        inPosFA = currPosFA(:,3) >= edgesZ(k) & ...
            currPosFA(:,3) < edgesZ(k+1);
        inPosFAterm = currPosFAterm(:,3) >= edgesZ(k) & ...
            currPosFAterm(:,3) < edgesZ(k+1);
        beadCountWater(k,:) = beadCountWater(k,:) + sum(inPosWater);
        beadCountFA(k,:) = beadCountFA(k,:) + sum(inPosFA);
        beadCountFAterm(k,:) = beadCountFAterm(k,:) + sum(inPosFAterm);
    end
end

% mass densities
densWater = (beadCountWater ./ Na .* Mw) ./ nTimesteps ./ (dV*(1e-10)^3);
densFA = (beadCountFA ./ Na .* Mc + beadCountFAterm ./ Na .* Mt) ...
    ./ nTimesteps ./ (dV * (1e-10)^3);

% cumulative densities
densFAcum = cumtrapz(z,densFA);
densWaterCum = cumtrapz(z,densWater);

%% plotting

f=figure(1);
f.Color = "white";
set(f,'position',[0,0,600,450])

yyaxis left
plot(z/10,densFA,'LineWidth',1.5)
xlim([0 14])
ylim([0 1600])
ylabel('\rho [kg/m^3]','FontSize',18)

yyaxis right
plot(z/10,densWater,'LineWidth',1.5)
ylim([0 400])

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

grid on
box on
xlabel(['z [nm]'],'FontSize',18)
% ylabel('\rho [kg/m^3]','FontSize',18)

legend({'Monolayer','Water'},'FontSize',12,'Location','northeast')

% Create arrow
annotation(f,'arrow',[0.350000000000001 0.450000000000001],...
    [0.746666666666667 0.746666666666667],...
    'Color',[0.525098039215686 0.865882352941177 0],...
    'LineWidth',1.5);

% Create arrow
annotation(f,'arrow',[0.363333333333338 0.263333333333338],...
    [0.144444444444444 0.144444444444444],...
    'Color',[0 0.27921568627451 0.474509803921569],...
    'LineWidth',1.5);

title('\chi_N = 0.25','FontWeight','normal','FontSize',16)

% cumulative lipid density (useful for determining surface coordinate, by
% definition integral(rho)/integral(rho_inf) = 0.99
f2 = figure(2);
plot(z,densFAcum./densFAcum(end))
yline(0.99)

%% print to figure
% print(f,'densityProfiles_0p65nm_0p25dam_water_C1','-dpng','-r300')
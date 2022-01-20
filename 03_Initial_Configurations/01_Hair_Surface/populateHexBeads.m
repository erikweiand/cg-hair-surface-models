% Graft fatty acids in a triangular structure inside a periodic box,
% E. Weiand - 01/2022,
% Based on: https://github.com/JeroenMulkers/lammps_graphene/blob/master/graphene.m

% General notes:
%   - The random damage seeds are NOT repeatable (i.e. launching this
%   script twice will give you two different damage patterns. Consider
%   using MATLAB's "rng" function if this is of interest for you.

clear
clc

%% Input parameters

% interatomic distance
a  = 2.836;     % [Angstroem]

% number damage ratio (0 <= chiN <= 1)
damageRatio = 0.5;

% periodic box size
lx = 357.336;
ly = 309.462;
lz = 150;

% number of repetitions in the x direction (must be
% sufficiently large to compensate for rotation)
nx = 75;        

% number of repetitions in the y direction (must be
% sufficiently large to compensate for rotation)
ny = 60;

% output file name
filename = 'fatty_acids.data';

%% Calculation

% angle for rotation along surface-normal axis
% to obey periodicity
phi = atan(sqrt(3) / 2 / 4.5);

% size of the unit cell
A = sqrt(84) / 4 * a;
B = sqrt(3) * A;

% fatty acid bead separation distance in z
dF = 4.0;

% minimum bead position
zMin = 0.0;

% number of fatty acid chains in basic building block
nFA = 2;

% Coordinates of the fatty acid beads (single chain)
nBeads = 8;
fattyAcidPos = zeros(nBeads, 3);
fattyAcidPos(:,3) = linspace(zMin, dF * nBeads, nBeads);

% Coordinates of the 4 atoms in the unit cell
base = [ fattyAcidPos;
         fattyAcidPos + [A/2 , B/2 , 0.0]];

% Total number of atoms
N = size(base, 1) * nx * ny;

% Calculate the coordinates of the atoms in the layer
coords = zeros(N, 3);
id = 0;
for ix=0:(nx - 1)
    for iy=0:(ny - 1)
        for iatom=1:size(base, 1)
            id = id + 1;
            coords(id,:) = base(iatom, :)+[ix * A, iy * B, 0];
        end
    end
end

% translate coordinates prior to rotation
coords = coords - [2 * A, 12 * B, 0];

% rotate points along the surface-normal axis
coords = [coords(:, 1) .* cos(phi) - coords(:, 2) .* sin(phi), coords(:, 1) .* sin(phi) + coords(:, 2) .* cos(phi), coords(:, 3)];

% remove points outside periodic bounds
coords(coords(:,1) < 1e-3, :)= [];
coords(coords(:,2) < 1e-3, :)= [];
coords(coords(:,1) > lx - 1e-3, :)= [];
coords(coords(:,2) > ly - 1e-3, :)= [];

% total number of FA chains
N = length(coords) / nBeads;

% number of fatty acid chains to be removed
Nrem = floor(N * damageRatio);
remPerm = randperm(N);
indSulf = remPerm(1:Nrem);
indSulf = sort(indSulf)';

% bead types in fatty acid
fattyAcidTypes = [9 3 2 1 1 1 1 10]';   % 1: C1 bead
                                        % 2: C5 bead
                                        % 3: P5 bead
                                        % 9: ghost graphene atoms
                                        % 10: terminal group (C1 or C2)
sulfTypes = [9 3 11]';                  % 11: SO4 bead (Qa)

% number of beads per chain
nBeads = length(fattyAcidTypes);
nBeadsSulf = length(sulfTypes);

% replace coordinates of removed FA chains
indSulfCurr = 1;
coordsNew = [];
for i=1:N
    if (damageRatio > 0) && (i == indSulf(indSulfCurr))
        currCoords = coords((i - 1) * nBeads + 1:(i - 1) * nBeads + nBeadsSulf, :);
        coordsNew = [coordsNew; currCoords];
        if indSulfCurr < Nrem
            indSulfCurr = indSulfCurr + 1;
        end
    else
        currCoords = coords((i - 1) * nBeads + 1:i * nBeads,:);
        coordsNew = [coordsNew; currCoords];
    end
end
coords = coordsNew;

% merge types
types = [];
indSulfCurr = 1;
for i=1:N
    if (damageRatio > 0) && (i == indSulf(indSulfCurr))
        types = [types; sulfTypes];
        if indSulfCurr < Nrem
            indSulfCurr = indSulfCurr + 1;
        end
    else
        types = [types; fattyAcidTypes];
    end
end
nBeadsTot = length(types);

fattyAcidCharges = zeros(nBeads, 1);
sulfCharges = [0; 0; -1.0];

% merge charges
charges = [];
indSulfCurr = 1;
for i=1:N
    if (damageRatio > 0) && (i == indSulf(indSulfCurr))
        charges = [charges; sulfCharges];
        if indSulfCurr < Nrem
            indSulfCurr = indSulfCurr + 1;
        end
    else
        charges = [charges; fattyAcidCharges];
    end
end

% fatty acid bonds
fattyAcidBonds = [1 2;
    2 3;
    3 4;
    4 5;
    5 6;
    6 7;
    7 8];

% sulfonate group bond
sulfBonds = [1 2;
    2 3];

% fatty acid angles
fattyAcidAngles = [1 2 3;
    2 3 4;
    3 4 5;
    4 5 6;
    5 6 7;
    6 7 8];

sulfAngles = [1 2 3];

nBondsPerChain = nBeads;
nAnglesPerChain = nBeads - 1;

nBondsTot = length(fattyAcidBonds);
nBondsTotSulf = 2;
nAnglesTot = length(fattyAcidAngles);
nAnglesTotSulf = 1;

bonds = [];
angles = [];
indSulfCurr = 1;
beadOffset = 0;
for i=1:N
    if (damageRatio > 0) && (i == indSulf(indSulfCurr))
        currBonds = sulfBonds + beadOffset;
        currAngles = sulfAngles + beadOffset;
        bonds = [bonds; currBonds];
        angles = [angles; currAngles];
        if indSulfCurr < Nrem
            indSulfCurr = indSulfCurr + 1;
        end
        beadOffset = beadOffset + nBeadsSulf;
    else
        currBonds = fattyAcidBonds + beadOffset;
        currAngles = fattyAcidAngles + beadOffset;
        bonds = [bonds; currBonds];
        angles = [angles; currAngles];
        beadOffset = beadOffset + nBeads;
    end
end

% write results
fid = fopen(filename,'w');
fprintf(fid,'Written in MATLAB. Graphene a=%g\n\n',a);
fprintf(fid,'%g atoms\n',nBeadsTot);
fprintf(fid,'11 atom types\n');
fprintf(fid,'%g bonds\n',length(bonds));
fprintf(fid,'3 bond types\n');
fprintf(fid,'%g angles\n',length(angles));
fprintf(fid,'3 angle types\n\n');
fprintf(fid,'0 %g xlo xhi\n',lx);
fprintf(fid,'0 %g ylo yhi\n',ly);
fprintf(fid,'0 %g zlo zhi\n\n',lz);
fprintf(fid,'Masses\n\n');
fprintf(fid,'1 72.0\n');
fprintf(fid,'2 72.0\n');
fprintf(fid,'3 72.0\n');
fprintf(fid,'4 72.0\n');
fprintf(fid,'5 24.0\n');
fprintf(fid,'6 24.0\n');
fprintf(fid,'7 72.0\n');
fprintf(fid,'8 72.0\n');
fprintf(fid,'9 24.0\n');
fprintf(fid,'10 90.0\n');
fprintf(fid,'11 72.0\n\n');

% write Atoms
fprintf(fid,'Atoms # full\n\n');
for i=1:nBeadsTot
    fprintf(fid,'%g %g %g %g %g %g %g 0 0 0\n',i,i,types(i),charges(i),coords(i,:));
end

% write Bonds
fprintf(fid,'\nBonds\n\n');
for i=1:length(bonds)
    fprintf(fid,'%g 1 %g %g\n',i,bonds(i,1), bonds(i,2));
end

% write Angles
fprintf(fid,'\nAngles\n\n');
for i=1:length(angles)
    fprintf(fid,'%g 1 %g %g %g\n',i,angles(i,:));
end
fclose(fid);
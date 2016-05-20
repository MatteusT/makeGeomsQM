function allZmats = getG09geoms(dir,startGeom,increments)
% Parsing the output file or the log file from gaussian
% startGeom can be inputed as 'last' in whichcase only the last geometry
% will be outputed
fid1 = fopen(dir,'r');
t1 = textscan(fid1,'%s');
text = t1{1};
if nargin < 2
startGeom = 30;
end
if nargin < 3
increments = 10;
end

phrase = {'NAtoms='};
loc = findText(text,phrase);
natom = str2double(text{loc(1)+1});

phrase = {'Input','orientation:'};
loc = findText(text,phrase);
atomOrder = zeros(1,natom);
icurr = loc(1)+16;
for iatom = 1:natom
    atomOrder(iatom) = str2double(text{icurr});
    icurr = icurr + 6;
end
%

phrase = {'Symbolic','Z-matrix:'};
loc = findText(text,phrase);
bondC = zeros(natom,1);
angleC = zeros(natom,1);
dihedralC = zeros(natom,1);
zmat = ZMatrix;
zmat.make_atom( atomOrder(1), 0, 0, 0);
icurr = loc+10;
bondC(2) = str2double(text{icurr});
zmat.make_atom( atomOrder(2), bondC(2), 0, 0);
icurr = icurr+ 3;
bondC(3) = str2double(text{icurr});
icurr = icurr +2;
angleC(3)= str2double(text{icurr});
zmat.make_atom( atomOrder(3), bondC(3), angleC(3), 0);

for iatom = 4:natom
    icurr = icurr + 3;
    bondC(iatom) = str2double(text{icurr});
    angleC(iatom)= str2double(text{icurr+2});
    dihedralC(iatom) = str2double(text{icurr+4});
    icurr = icurr+5;
    zmat.make_atom( atomOrder(iatom), bondC(iatom), angleC(iatom), dihedralC(iatom));
end

phrase = {'Input','orientation:'};
loc = findText(text,phrase);
rcart = zeros(3,natom);
ic = 1;
if strcmp(startGeom,'last')
startGeom = length(loc);
increments = 1;
end
for i = startGeom:increments:length(loc)
    icurr = loc(i)+18;
    for iatom = 1:natom
        rcart(1,iatom) = str2double(text{icurr});
        rcart(2,iatom) = str2double(text{icurr+1});
        rcart(3,iatom) = str2double(text{icurr+2});
        icurr= icurr + 6;
    end
    
    [bonds,angles,dihedrals] = cart2Zmat(rcart,bondC,angleC,dihedralC);
    zmat.pars.bond_pars = bonds;
    zmat.pars.ang_pars = angles;
    zmat.pars.di_pars = dihedrals;
    
    allZmats{ic} = saveZmat(zmat);
    ic = ic + 1;
    % Convert from Bohr radii to Angstroms
    % rcart = rcart / 1.889726124565062;
end
end

function zmat = saveZmat(zmat)
save('tempZ.mat','zmat')
clear zmat;
load('tempZ.mat');
end
%
%

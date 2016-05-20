clear classes
addpath('C:\Users\Matteus\master_MSQC') % Has dependencies to MSQC
dir = 'C:/Users/Matteus/quamboML/makeGeomsQM/hydroCarbons/'; % Directory where you have the molecular structures
mol = {'C2H6','C3H8','C4H10','tC4H10'}; % Name of the files containing the molecular structures
ic = 1;
for imol = 1:length(mol)
        molDir = [dir,mol{imol}];
        zmatOpt(imol) = getZmat(molDir);
        ic = ic + 1;
end
%%
bdev = 0.2; % bond deviation
adev = 10;  % angle deviation
ddev = 40;  % dihedral deviation
ndist = 30; % number of deviations per molecule

zmats = cell(length(zmatOpt),ndist+1);
zmats(:,1) =  zmatOpt';

for iopt = 1:length(zmatOpt)
for i = 1:ndist
  dGeom = distortGeom(zmatOpt{iopt});
    dGeom.distortM1(bdev,adev,ddev);
   zmats{iopt,i+1} = dGeom.zmat;
end
end

%%
save([dir,'hydroCarbons_zmats.mat'],'zmats')
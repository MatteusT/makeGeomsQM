classdef distortGeom
    %DISTORTGEOM: random disortion of any geometry inputed as a zmat class
    %   Input: zmat, bondDev, angDev and diDev
    %          For now the only rotation is along the first diheral and you
    %          can set the diDistortion to anything between -180 and 180
    
    properties
        zmat
        zparsOrig
        bondList
        angleList
        dihedralList
    end
    
    methods
        function obj = distortGeom(objIn)
            obj.zmat = objIn.deepCopy;
            obj.zparsOrig = obj.zmat.pars;
            isBonded = obj.zmat.isBonded;
            natoms = size(isBonded,1);
            obj.bondList = [];
            ic = 0;
            for i =1:natoms
                for j = (i+1):natoms
                    if (isBonded(i,j))
                        ic = ic+1;
                        obj.bondList(1,ic) = i;
                        obj.bondList(2,ic) = j;
                    end
                end
            end
            obj.angleList = [];
            nbonds = size(obj.bondList,2);
            ic = 0;
            for i =1:nbonds
                for j = (i+1):nbonds
                    uniqAtoms = unique([obj.bondList(:,i);obj.bondList(:,j)]);
                    if (length(uniqAtoms) == 3)
                        ic = ic + 1;
                        central = intersect(obj.bondList(:,i),obj.bondList(:,j));
                        outside = sort(setdiff(uniqAtoms,central));
                        obj.angleList(1,ic) = outside(1);
                        obj.angleList(2,ic) = central;
                        obj.angleList(3,ic) = outside(2);
                    end
                end
            end
            
        end
        
        function res = NB1dist(obj)
            rcart = obj.zmat.toCart;
            natoms = size(rcart, 2);
            %           res = zeros(natoms*(natoms-1)/2, 1);
            ic = 0;
            res = NaN;
            for i =1:natoms
                for j = 1:natoms
                    if ((~obj.zmat.isBonded(i,j))&& i~=j)
                        ic = ic + 1;
                        res(ic) = norm(rcart(:, i) - rcart(:, j));
                    end
                end
            end
        end
        
        function res = NB2dist(obj)
            rcart = obj.zmat.toCart;
            natoms = size(rcart, 2);
            %           res = zeros(natoms*(natoms-1)/2, 1);
            ic = 0;
            res = Inf;
            b2 = obj.zmat.isBonded*obj.zmat.isBonded;
            for i =1:natoms
                for j = 1:natoms
                    if ((~obj.zmat.isBonded(i,j))&& i~=j && ~b2(i,j))
                        ic = ic + 1;
                        res(ic) = norm(rcart(:, i) - rcart(:, j));
                    end
                end
            end
        end
        
        function res = bondAngles(obj)
            rcart = obj.zmat.toCart;
            nang = size(obj.angleList,2);
            res = zeros(nang,1);
            for ic = 1:nang
                bond1 = rcart(:, obj.angleList(1,ic)) - ...
                    rcart(:, obj.angleList(2,ic));
                bond2 = rcart(:, obj.angleList(3,ic)) - ...
                    rcart(:, obj.angleList(2,ic));
                res(ic) = acosd(dot(bond1,bond2)/norm(bond1)/norm(bond2));
            end
        end
        function res = bondLengths(obj)
            rcart = obj.zmat.toCart;
            nbonds = size(obj.bondList,2);
            res = zeros(nbonds,1);
            for ic = 1:nbonds
                i = obj.bondList(1,ic);
                j = obj.bondList(2,ic);
                res(ic) = norm(rcart(:,i)-rcart(:,j));
            end
        end
        function distortM1(obj,bondDev,angDev,diDev) 
            bls = obj.zparsOrig.bond_pars;
            angs = obj.zparsOrig.ang_pars;
            origAllAngs = obj.bondAngles;
            bldev = bondDev * 2 * (rand(size(bls)) - 0.5);
            obj.zmat.pars.bond_pars = bls + bldev;
            cont = 1;
            while cont == 1
            adev = angDev * 2 * (rand(size(angs)) -0.5);
            ddev = diDev * 2 * (rand(1) -0.5);
            obj.zmat.pars.ang_pars = angs + adev;
            obj.zmat.pars.di_pars(1) = obj.zmat.pars.di_pars(1) + ddev;
            if (max(abs(origAllAngs-obj.bondAngles))<angDev+0.01 && max(abs(obj.NB2dist))> 3)
            cont = 0;
%           disp([max(abs(origAllAngs-obj.bondAngles)) min(abs(origAllAngs-obj.bondAngles))])
            end
            end
        end
    end
    
end


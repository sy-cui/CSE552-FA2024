% Assume defined: TRACK_ELEM_IDX, TRACK_QUAD_INDEX

if nen == 3
    nel = 3;
elseif nen == 4
    if ix(TRACK_ELEM_IDX,nen) == 0
        nel = 3;
    else
        nel = 4;
    end
elseif nen == 6
    nel = 6;
else
    if ix(TRACK_ELEM_IDX,nen) == 0
        nel = 6;
    else
        nel = 9;
    end
end

% Extract patch nodal coordinates
[xr, ~] = extract_elem_var(NodeTable,TRACK_ELEM_IDX,ndm,nel,ix);
[xc, ~] = extract_elem_var(NodeCurr,TRACK_ELEM_IDX,ndm,nel,ix);

u = reshape(xc-xr,[],1);

% Extract patch material properties
ma = ix(TRACK_ELEM_IDX,nen+1);
mateprop = MateT(ma,:);

E = mateprop(1);
nu= mateprop(2);

% Material moduli tensor
if PSPS == 's' % Plane Stress
    Dmat = E/(1-nu^2)*[1  nu 0;
                       nu 1  0;
                       0  0  (1-nu)/2];

else % Plane Strain
    Dmat = E/(1+nu)/(1-2*nu)*[1-nu nu    0;
                              nu   1-nu  0;
                              0    0     (1-2*nu)/2];     
end

% Load Guass Integration Points
lint = 4;

if nel == 3
    [Wgt,r,s] =  intpntt(TRACK_QUAD_IDX,lint,0);
else
    [Wgt,r,s] =  intpntq(TRACK_QUAD_IDX,lint,0);
end

[Qxy, Jfe] = shpg_2d(shpl_2d(r,s,nel),xr,nel);

if nel == 3
    Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        
            0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3)
            Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3)];
else
    Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        Qxy(1,4) 0 
            0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3) 0        Qxy(2,4)
            Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3) Qxy(2,4) Qxy(1,4)];
end

strain = Bmat * u;
stress = Dmat * strain;

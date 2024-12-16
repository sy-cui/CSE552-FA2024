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
[xr, ElemFlag] = extract_elem_var(NodeTable,TRACK_ELEM_IDX,ndm,nel,ix);
[xc, ~]        = extract_elem_var(NodeCurr,TRACK_ELEM_IDX,ndm,nel,ix);

% Extract patch material properties
ma = ix(TRACK_ELEM_IDX,nen+1);
mateprop = MateT(ma,:);

% Load Guass Integration Points
lint = 4;

if nel == 3
    [Wgt,r,s] =  intpntt(TRACK_QUAD_IDX,lint,0);
else
    [Wgt,r,s] =  intpntq(TRACK_QUAD_IDX,lint,0);
end

% Deformation gradient
[F, J] = def_grad_2d(r,s,xc,xr);

% Cauchy stress
sigma = cauchy(F,J,mateprop);

% Green-Lagrange strain
GLS = 0.5*(F'*F - eye(2));

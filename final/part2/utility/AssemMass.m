% Assemble element mass into model mass
%
% Copyright (C) Arif Masud and Tim Truster
%
% 10/2009
% UIUC

locind1 = 0;
for ie1 = 1:nel
    for l1 = 1:ndf
        locind1 = locind1 + 1;
        grow = EDOFT(locind1);
        if grow <= neq
            locind2 = 0;
            for ie2 = 1:nel
                for l2 = 1:ndf
                    locind2 = locind2 + 1;
                    gcol = EDOFT(locind2);
                    if gcol <= neq
                        Mdd(grow, gcol) = Mdd(grow, gcol) + ElemM(locind1,locind2); %#ok<SAGROW>
                    else %gcol > neq
                        gcol = gcol - neq;
                        Mdf(grow, gcol) = Mdf(grow, gcol) + ElemM(locind1,locind2); %#ok<SAGROW>
                    end %if gcol
                end %l2
            end %ie2
        else %grow > neq
            grow = grow - neq;
            locind2 = 0;
            for ie2 = 1:nel
                for l2 = 1:ndf
                    locind2 = locind2 + 1;
                    gcol = EDOFT(locind2);
                    if gcol <= neq
                        Mfd(grow, gcol) = Mfd(grow, gcol) + ElemM(locind1,locind2); %#ok<SAGROW>
                    else %gcol > neq
                        gcol = gcol - neq;
                        Mff(grow, gcol) = Mff(grow, gcol) + ElemM(locind1,locind2); %#ok<SAGROW>
                    end %if gcol
                end %l2
            end %ie2
        end %if grow
    end %l1
end %ie1
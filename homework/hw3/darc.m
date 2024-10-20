function [dfDu, dfDl] = darc(Du, Dl, dK, b, c)
    f = arc(Du,Dl,dK,b,c);
    dfDu = c * dK * Du ./ f;
    dfDl = b * Dl / f;
end
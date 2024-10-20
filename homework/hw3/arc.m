function f = arc(Du, Dl, dK, b, c)
    f = sqrt(c*Du'*dK*Du + b*Dl^2);
end
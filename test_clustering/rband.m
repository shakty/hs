function ry = rband(ri, band, inner)
    if (inner)
        ry = sqrt((ri^2*pi - band)/pi);
    else
        ry = sqrt((ri^2*pi + band)/pi);
    end
end
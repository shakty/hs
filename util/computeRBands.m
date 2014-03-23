function band = computeRBands(startFrom, BAND, limit)

    band = [startFrom;startFrom];
    i = 1;

    while (band(1,i) < limit)
        band(2,i) = rband(band(1,i), BAND, 0);        
        i = i + 1;        
        band(1,i) = band(2,i-1);
    end
    band = band(:,1:end-1);

end


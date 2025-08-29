function [wind_v] = wind(alt)

alt_data =        [0 330 600 750 900 1500 2000 3000 4200 5500 7000 9000 10000 11700 13500 30000];
wind_data = [8 10   10  10  9   6    4    9    15   19   17   18   22     25    18    19];

wind_v = interp1(alt_data, wind_data, alt, 'linear', 'extrap');
end


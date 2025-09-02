function [wind_v, wind_dir] = wind(alt)

alt_data =  [0 100 600 750 900 1500 2000 3000 4200 5500 7000 9000 10000 11700 13500 30000]; % m/s
wind_data = [8 10   10  10  9   6    4    9    15   19   17   18   22     25    18    19];  % m
dir_data =  [cos(10) sin(10); cos(0) sin(0); cos(20) sin(20); cos(50) sin(50); cos(110) sin(110); cos(100) sin(100); cos(120) sin(120); cos(120) sin(120); cos(150) sin(150); cos(130) sin(130); cos(180) sin(180); cos(200) sin(200); cos(240) sin(240); cos(230) sin(230); cos(230) sin(230); cos(70) sin(70)];  % unit vect

wind_v = interp1(alt_data, wind_data, alt, 'linear', 'extrap');

[~, idx] = min(abs(alt_data - alt));
wind_dir = dir_data(idx,:);
end


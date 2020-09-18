function RunNoiseGenerator(noisecolor, noisenr)
% generates colored noise and saves to csv

keys = {'white','pink','blue'}; values = [1,2,3];
M = containers.Map(keys,values);

noise = [];
for i = 1:noisenr
    noise = [noise; NoiseGenerator(60,1000,M(noisecolor))]; 
end 
csvwrite(strcat(noisecolor, '_noise.csv'), noise)
end

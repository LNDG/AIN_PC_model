function PowFreqCorr=RunNoiseGenerator(noisecolor, range)
% function to run noise generator and saves to csv. Returns the log log 
% correlation of power and frequency
    
    % noisecolor: 'white', 'pink' or 'blue'
    % range: range of frequences in noise [lowerbound, higherbound]

keys = {'white','pink','blue'}; values = [1,2,3];
M = containers.Map(keys,values);
fs = 1000;
noise = [];
PowFreqCorr = [];
for i = 1:4
    noise = [noise; NoiseGenerator(3,fs,M(noisecolor),range(1),range(2),8)];
    [pxx, fx] = pwelch(noise(i,:),hann(100),[],[range(1):1:range(2)],fs);
    plot(fx, pxx)
    lpxx = log(pxx);
    lfx = log(fx);
    PowFreqCorr = [PowFreqCorr, corr(lfx',lpxx')];
end 
csvwrite(strcat(noisecolor, '_noise_', string(range(1)), '-', string(range(2)), '.csv'), noise)
end

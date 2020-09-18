function NoiseVector = NoiseGenerator(stimtime,samplingrate,noisecolour)%,lowcut,highcut,filtorder)
%This function generates noise vectors of different colours of any length and filters them (via butterworth filtering) to any desired frequency range. Douglas Garrett, March 2014.
%    
    %function NoiseVector = NoiseGenerator(stimtime,samplingrate,noisecolour,lowcut,highcut,filtorder)
        % NoiseVector - output row vector of noise samples.
        % stimtime - total desired stimulation time (seconds).
        % samplingrate - sampling rate of stimulator (Hz).
        % noisecolour - enter as 1 (=white), 2 (=pink), or 3(=blue).
        % lowcut - low frequency cutoff (Hz); highpass filter boundary.
        % highcut - high frequency cutoff (Hz); lowpass filter boundary.
        % filtorder - butterworth filter order. Recommended order = 8.
 
%White noise has a flat spectrum, pink noise falls off at 3 dB per octave, and blue noise increases at 3 dB per octave. 

% Define the length of the noise vector.
% Ensure that the M is even.
N = stimtime*samplingrate;
if rem(N,2)
    M = N+1;
else
    M = N;
end

% Generate white noise of length M (by default, randn produces noise of
% mean = 0 and stdev=1). If noisecolour is specified 'pink' (==2) or 'blue' (==3), we will then filter
% white noise to achieve coloured noise.
NoiseVector = randn(1, M);
    
    if noisecolour==2||noisecolour==3
            % FFT
            NoiseVector = fft(NoiseVector);

            % prepare a vector for 1/f multiplication
            NumUniquePts = M/2 + 1;
            n = 1:NumUniquePts;
                    
                if noisecolour==2
                    % multiply the left half of the spectrum so the power spectral density
                    % is inversely proportional to the frequency by factor 1/f, i.e. the
                    % amplitudes are inversely proportional to 1/sqrt(f)
                    NoiseVector(1:NumUniquePts) = NoiseVector(1:NumUniquePts)./sqrt(n);

                elseif noisecolour==3            
                    % multiply the left half of the spectrum so the power spectral density
                    % is proportional to the frequency by factor f, i.e. the
                    % amplitudes are proportional to sqrt(f)
                    NoiseVector(1:NumUniquePts) = NoiseVector(1:NumUniquePts).*sqrt(n);
                
                end
                
            % prepare a right half of the spectrum - a copy of the left one,
            % except the DC component and Nyquist frequency - they are unique
            NoiseVector(NumUniquePts+1:M) = real(NoiseVector(M/2:-1:2)) -1i*imag(NoiseVector(M/2:-1:2));

            % IFFT
            NoiseVector = ifft(NoiseVector);

            % prepare final output vector
            NoiseVector = real(NoiseVector(1, 1:N));

            % normalise to mean=0 and Stdev=1, to equate with white noise.
            NoiseVector = zscore(NoiseVector);
    end

%now filter...use two successive filters here as there seems to be a bug
%with Matlab's bandpass option.
%[B,A] = butter(filtorder,lowcut/(samplingrate/2),'high'); 
%NoiseVector  = filtfilt(B,A,NoiseVector); clear A B;

%[B,A] = butter(filtorder,highcut/(samplingrate/2),'low');
%NoiseVector  = filtfilt(B,A,NoiseVector); clear A B

end

%be very careful when simulating blue noise in very low frequency range
%with sparse data (say 600 sec at 1.55 Hz sampling rate, and .01-.10
%bandpass range. Seems to be major filter artifacts at the boundaries of
%the times series, likely dirven by .01 Hz filter end. Pink and white seem ok though. Look into this further if
%needed? Moving blue noise to .025 Hz greatly reduces the issue.


%functions adapted from: 
    %1. Noise Generation with MATLAB Implementation. Author: M.Sc. Eng. Hristo Zhivomirov 07/30/13. See http://www.mathworks.co.uk/matlabcentral/fileexchange/42919-pink-red-blue-and-violet-noise-generation-with-matlab-implementation/content/pinknoise.m  

    %2. Thomas Grandy's in-house filt_but_lp and _hp functions.
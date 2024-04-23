%------------------------------------------%
%        *** SSSP - HOMEWORK #3 ***        %
%------------------------------------------%
%     Emulation of the Leslie Speaker      %
%------------------------------------------%
% Name: Brusca Alfredo                     %
% Student ID: 10936149                     %
% Name: Pomarico Riccardo                  %
% Student ID: 10661306                     %
%------------------------------------------%

clear; close all; clc;

%% modulation speed
mod_speed = 'chorale';

%% Read the input file
[x, Fs] = audioread('HammondRef.wav');
x=x(:,1);           % take left channel only

%% FX parameters 
switch lower(mod_speed)
    
    case {'chorale'}
        freq=2;
            
    case {'tremolo'}
        freq=6;
        
    otherwise
        error('mod_speed \"%s\" not found.', mod_speed)

end

%% Apply FX
y  = leslie(x, Fs, freq);

%% Avoid any (possible) clipping
y = rescale(y,-1.,1.);

%% Playback
%audiowrite([mod_speed,'.wav'], y, Fs);
soundsc(y, Fs)

%% Read the reference audio file
dir_name = 'Leslie_ref';
addpath(dir_name);
[y_ref, ~] = audioread(fullfile(dir_name, strcat(mod_speed,'.wav')));

%% Display the MSE
MSE = mean(abs(y.'-y_ref).^2);
MSE_str = sprintf('MSE: %g', MSE);
disp(MSE_str)

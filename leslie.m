function [y,y_lpf,y_hpf,y_hp_sdf] = leslie(x, Fs, freq)

%% SASP Homework #3 Group composition:
% Brusca Alfredo 10936149
% Pomarico Riccardo 10661306
% 
% Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

% length of the input signal
N = length(x);

% global modulator parameters
alpha=0.9;
% tremble spectral delay filter parameter 
Ms_t=0.2;
Mb_t=-0.75;
N_sdf_t=4;
% bass spectral delay filter parameter 
Ms_b=0.04;
Mb_b=-0.92;
N_sdf_b=3;

% cross-over network design
fc=800;                 % cutoff frequency

%TODO: compute the coefficients for the two 4th order butterworth filters
%with cutoff frequency fc
[b_lp, a_lp]= butter(4, fc/(Fs/2), 'low'); %LPF design
[b_hp, a_hp]= butter(4, fc/(Fs/2), "high");  %HPF design

% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state=zeros(N_sdf_t+1,1);
hpf.in=zeros(4,1);
% lp filter buffers
lpf.state=zeros(N_sdf_b+1,1);
lpf.in=zeros(4,1);
% treble sdf filter buffers
sdf_h.state=zeros(N_sdf_t,1);
sdf_h.in=zeros(N_sdf_t+1,1);
% bass sdf filter buffers
sdf_b.state=zeros(N_sdf_b,1);
sdf_b.in=zeros(N_sdf_b+1,1);

% modulators
m_b = Ms_b*sin(2*pi*(freq)*(1:N)/Fs)+Mb_b; % bass modulator
m_t = Ms_t*sin(2*pi*(freq+0.1)*(1:N)/Fs)+Mb_t; % tremble modulator
for i=1:N_sdf_b
    Ni_b(i) = nchoosek (N_sdf_b,i);
end
for i=1:N_sdf_t
    Ni_t(i) = nchoosek (N_sdf_t,i);
end

%sample processing
for n=1:N

    % compute crossover network filters outputs
    
    y_lpf = b_lp(1)*x(n) + b_lp(2)*lpf.in(1) + b_lp(3)*lpf.in(2) + b_lp(4)*lpf.in(3) + b_lp(5)*lpf.in(4)- a_lp(2)*lpf.state(1) - a_lp(3)*lpf.state(2) - a_lp(4)*lpf.state(3) - a_lp(5)*lpf.state(4);
    y_hpf = b_hp(1)*x(n) + b_hp(2)*hpf.in(1) + b_hp(3)*hpf.in(2) + b_hp(4)*hpf.in(3) + b_hp(5)*hpf.in(4)- a_hp(2)*hpf.state(1) - a_hp(3)*hpf.state(2) - a_hp(4)*hpf.state(3) - a_hp(5)*hpf.state(4);
    for i = N_sdf_b+1:-1:2
        lpf.state(i) = lpf.state(i-1);
    end
    for i = N_sdf_t+1:-1:2
        hpf.state(i) = hpf.state(i-1);
    end
    for i = 4:-1:2
        hpf.in(i) = hpf.in(i-1);
        lpf.in(i) = lpf.in(i-1);
    end

    lpf.in(1) = x(n);
    hpf.in(1) = x(n);
    lpf.state(1) = y_lpf;
    hpf.state(1) = y_hpf;
    sdf_b.in = lpf.state;
    sdf_h.in = hpf.state;

    % compute bass SDF output
    y_lp_sdf = sdf_b.in(N_sdf_b+1);
    for i = 1 : N_sdf_b
        y_lp_sdf = y_lp_sdf + Ni_b(i)*(m_b(n)^i)*(sdf_b.in(1+N_sdf_b-i)-sdf_b.state(i));
    end


    % compute treble SDF output
    y_hp_sdf = sdf_h.in(N_sdf_t+1);
    for i = 1 : N_sdf_t
        y_hp_sdf = y_hp_sdf + Ni_t(i)*(m_t(n)^i)*(sdf_h.in(1+N_sdf_t-i)-sdf_h.state(i));
    end

    for i = N_sdf_b:-1:2
        sdf_b.state(i) = sdf_b.state(i-1);
    end
    sdf_b.state(1) = y_lp_sdf;

    for i = N_sdf_t:-1:2
        sdf_h.state(i) = sdf_h.state(i-1);
    end
    sdf_h.state(1) = y_hp_sdf;

    % implement AM modulation block*
    y_lp_am=(1+alpha*(m_b(n)))*y_lp_sdf;
    y_hp_am=(1+alpha*(m_t(n)))*y_hp_sdf;

    y(n)= y_lp_am + y_hp_am;
    

end

end


function X=Wforward(P)

options.wavelet_type = 'haar'; % kind of wavelet
options.wavelet_vm = 4; % number of vanishing moments
Jmin = 3; %  minimum scale

X=perform_wavelet_transform(P, Jmin, +1, options);
end

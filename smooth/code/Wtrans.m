function P=Wtrans(X)

options.wavelet_type = 'haar'; % kind of wavelet
options.wavelet_vm = 4; % number of vanishing moments
Jmin = 3; %  minimum scale

P=perform_wavelet_transform(X, Jmin, -1, options);
end


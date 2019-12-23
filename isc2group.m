function [isc_tg, perf] = isc2group(X, G)
% ISC2GROUP computes neural ISC with respect to two groups with separate
% attentional focus. Similar to [1].
%   [ISC_TG] = ISC2GROUP(X, G) computes the most correlated
%   components separately for two groups and computes the correlations of
%   all participants with the responses of these groups. Input (X) contains
%   electroencephalograpic responses, which should be a 3-dimensional array
%   with size (T x D x N), where T is the number of samples in time, D is
%   the number of EEG channels and N is the number of participants. G
%   should be an N-dimensional array containing values 1 and 2,
%   corresponding to the group participants correspond to. This function
%   plots results in a bar plot, showing mean subject-to-group within and
%   between-groups as well as the relation between the separate scores of
%   an individual.
%
%   Dependencies: CorrCA
%
%   Copyright: Lucas C. Parra, adapted by Ivo Stuldreher
%   ivo.stuldreher@tno.nl
%
%   Last saved: 19-06-2019


%   REFERENCES 
%   [1] Ki, J. J., Kelly, S. P., & Parra, L. C. (2016). Attention strongly
%   modulates reliability of neural responses to naturalistic narrative
%   stimuli. Journal of Neuroscience, 36(10), 3092-3101.

% obtain number of channels (D) and number of participants (X)
[~,D,N] = size(X);

% validate G meets the requirements
errmsg = sprintf('G should have length %i and should contain only values 1 and 2', N);
assert(length(G) == N && all(G == 1 | G == 2), errmsg);

% compute cross-covariance between all N participants
fprintf('Computing cross-covariances between all %i partcipants...\n', N)
Rij = permute(reshape(nancov(X(:,:)),[D N  D N]),[1 3 2 4]);
fprintf('Finished\n')

isc_tg = zeros(D,N,2); % pre-assign variable
perf = zeros(1,2); 
for g = 1 : 2
    
    fprintf('Processing group %i/2...\n', g);
        
    % create test set for group g
    Xt = X(:,:,G == g);
    [vt] = corrca(Xt, 'shrinkage', 0.5);
    nbytes = fprintf('processing participant 0 of %d', N);
   
    for p = 1 : N
        while nbytes > 0
            fprintf('\b')
            nbytes = nbytes - 1;
        end
        nbytes = fprintf('processing participant %d of %d', p, N);
        
        % if participant p is group group g, adapt test set to not contain
        % participant p for bias
        if G(p) == g
            X_train = X(:,:,setdiff(find(G == g), p));
            [vi] = corrca(X_train, 'shrinkage', 0.5);
        else
            vi = vt;
        end
        
        % compute isc of participant p with respect to group g
        Rw = 0; 
        Rb = 0;
        N_group = length(setdiff(find(G == g), p));
        for j = setdiff(1 : N, find(G ~= g))
            if p~=j 
                Rw = Rw+1/(N_group-1)*(Rij(:,:,p,p)+Rij(:,:,j,j)); 
                Rb = Rb+1/(N_group-1)*(Rij(:,:,p,j)+Rij(:,:,j,p)); 
            end
        end
        isc_tg(:,p,g) = diag(vi'*Rb*vi)./diag(vi'*Rw*vi);
        
    end
        
end

for g = 1 : 2
    
    perf(g) = (2/N) * length(find(squeeze(sum(isc_tg(1:3,G == g, g))) > squeeze(sum(isc_tg(1:3,G == g, setdiff(1:2,g))))));
    
end
   
end

function Xr = phaserandomized(X)
% Generate phase randomized surrogate data Xr that preserves spatial and
% temporal correlation in X, following Prichard D, Theiler J. Generating 
% surrogate data for time series with several simultaneously measured 
% variables. Physical review letters. 1994 Aug 15;73(7):951.

[T,D,N] = size(X);

Tr = round(T/2)*2; % this code only works if T is even; make it so
for i = 1:N
    Xfft = fft(X(:,:,i),Tr); % will add a zero at the end if uneven length
    Amp = abs  (Xfft(1:Tr/2+1,:)); % original amplitude
    Phi = angle(Xfft(1:Tr/2+1,:)); % orignal phase
    Phir = 4*acos(0)*rand(Tr/2-1,1)-2*acos(0); % random phase to add
    tmp(2:Tr/2,:) = Amp(2:Tr/2,:).*exp(sqrt(-1)*(Phi(2:Tr/2,:)+repmat(Phir,1,D))); % Theiler's magic
    tmp = ifft([Xfft(1,:); tmp(2:Tr/2,:); Xfft(Tr/2+1,:); conj(tmp(Tr/2:-1:2,:))]); % resynthsized keeping it real
    Xr(:,:,i) = tmp(1:T,:,:); % grab only the original length
end

end
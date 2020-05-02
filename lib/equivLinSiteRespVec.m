function [freq_array,af,t_out,accel_on_surface,...
	      new_profile,out_a,out_v,out_d,out_gamma,out_tau,max_avd,max_gt,tf,G_mtrx,D_mtrx]...
         = equivLinSiteRespVec(profile,motion,curve,fig_visible,boundary,output_or_not)
% (c) Wei Li & Jian Shi
%
% [Inputs]
%    profile
%    motion
%    curve
%    fig_visible: 'off' (default) and 'on'
%    boundary: 'elastic' (default) and 'rigid'
%    output_or_not: 'n' (default) and 'y'
%
% [Outputs]
%    freq_array
%    af: equivalent linear amplification function (magnitude)
%    t_out: time array, same as input time array
%    accel_on_surface: only one column
%    new_profile: new (denser) soil profile after re-discretization
%    out_a: acceleration time history of every layer
%    out_v
%    out_d
%    out_gamma
%    out_tau
%    max_avd: maximum acceleration, velocity, and displacement at every
%             interface (including the soil-rock interface)
%    max_gt: maximum strain (gamma) and stress (tau) of every layer
%    tf: equivalent linear transfer function (complex)
%
%
% Please note that the "motion" should be incident motion, regardless of 
% the "boundary" type.
%
% First edition by Wei Li, 3/19/2007
% Last Modified, 12/14/2013 (changed maximum iteration time from 20 to 10)
% 
% * * *  Renamed from tfEQLv3 into the current name on July, 2015
%
% ---------- 4/17/2018: Vectorized version. --------------------

if nargin < 6, output_or_not = 'n'; end
if nargin < 5, boundary = 'elastic'; end
if nargin < 4, fig_visible = 'on'; end

I = 1i;

%% Part 1.1: Data preperation -- useful constants
% % % % boundary = 1; % 1 -> elastic, 2-> rigid, 3 -> borehole
tolerance = 0.075; % allowable tolerance for G and D (damping ratio)
R_gamma = 0.65; % used for computing characteristic strain
nr_iter = 10; % maximum number of iteration, regardless of convergance

%% Part 1.2: Data preperation -- import earthquake motion
if ischar(motion)
    quake = importdata(motion);  % time history -> quake
else
    quake = motion;
end
n = length(quake);

t = quake(:,1); % time array
accel_in = quake(:,2); % acceleration (m/s^2) -> a
dt = t(2)-t(1); % time interval of input motion

if floor(n/2) == n/2 % if n is even
    flag = 0;
    t = [t;t(end)+dt]; % pad one point at the end
    accel_in = [accel_in;accel_in(end)]; % to make sure that n is odd
    n = n+1;  % n is guaranteed to be odd
else  % if n is odd
    flag = 1;
end

ACCEL_IN = fft(accel_in);  % n x 1
N = length(ACCEL_IN); % N is definitely odd. N and n are identical!

f = (1:1:N)'./(N*dt);

%% Part 1.3: Data preperation -- Import soil profile
if ischar(profile)
    site_profile = importdata(profile); % soil propertis ->profile
    [file_dir,~,~] = fileparts(profile);
else
    site_profile = profile;
end
h = site_profile(:,1); % thickness of each layer (m)
vs = site_profile(:,2); % shear wave velocity of each layer (m/s)
D = site_profile(:,3); % initial value of damping ratio (unit: 1)
rho = site_profile(:,4); % mass density of soil of each layer (unit: kg/m3)
material_nr = site_profile(:,5); % material number

% Part 1.3.1 -- Divide initial profile into more sublayers * * *
[h,vs,D,rho,material_nr] = stratify(h,vs,D,rho,material_nr);
new_profile = [h,vs,D,rho,material_nr];
% * * * * * * * * * * * * * * * * * * * * * * * * 

nr_layer = length(h);  % it includes the rock layer
G_max = rho.*vs.^2; % G_max (unit: Pa)
G = G_max; %initial value of G

% Part 1.3.2 -- layer boundary depth & layer midpoint depth arrays
layer_boundary_depth = zeros(nr_layer,1);
layer_midpoint_depth = zeros(nr_layer-1,1);

for tt = 1 : 1 : nr_layer-1
	layer_boundary_depth(tt+1) = layer_boundary_depth(tt) + h(tt);
	layer_midpoint_depth(tt) = (layer_boundary_depth(tt+1)+layer_boundary_depth(tt))/2;
end

%% Part 1.4: Data preperation -- import modulus/damping curves
if ischar(curve)
    degr = importdata(curve); % degr curve -> degr
else
    degr = curve;
end

n_obs = size(degr,1); % number of points in a curve (=10)
strain_G = zeros(n_obs,nr_layer-1);
G_vector = zeros(n_obs,nr_layer-1);
strain_D = zeros(n_obs,nr_layer-1);
D_vector = zeros(n_obs,nr_layer-1);

for k = 1 : 1 : nr_layer-1
    strain_G(:,k) = degr(:,(material_nr(k)-1)*4+1)/100.0; % unit: from "percent" to "1"
    G_vector(:,k) = degr(:,(material_nr(k)-1)*4+2);
    strain_D(:,k) = degr(:,(material_nr(k)-1)*4+3)/100.0; % unit: from "percent" to "1"
    D_vector(:,k) = degr(:,(material_nr(k)-1)*4+4)/100.0; % unit: from "percent" to "1"
end

% Part 1.4.1 -- Adjust the damping curve according to array D
for k = 1 : 1 : nr_layer-1
    D_vector(:,k) = D_vector(:,k) - D_vector(1,k) + D(k);
end

%% Part 2: Start of iteration
G_mtrx = zeros(nr_layer-1,nr_iter+1); % matrix that stores G of all iterations
D_mtrx = zeros(nr_layer-1,nr_iter+1); % matrix that stores D of all iterations
G_mtrx(:,1) = G(1:end-1);
D_mtrx(:,1) = D(1:end-1);

fprintf('\n');
for i_iter = 1 : 1 : nr_iter
    fprintf('Iteration No.%d.',i_iter);
    %% Part 2.1: Linear transfer function
    % 2.1.0: Make a denser frequency array (decrease df), added 4/28/2016
    df = f(2)-f(1);
    max_f = max(f);
    f_oversample_factor = 15; % *** This is to eliminate frequency aliasing! (4/28/2016) *** (It has to be an odd number.)
    df_ = df / f_oversample_factor;  % new delta f (= 0.05 Hz)
    f = (df_ : df_ : max_f)';  % up-sample the frequency array
    N = N * f_oversample_factor;  % scale up N as well
    omega = 2*pi*f; % angular frequency array, N x 1
    
    % 2.1.1: Complex impedance ratio
    alpha = zeros(nr_layer-1,1);
    for j = 1 : 1: nr_layer-1
        alpha(j) = (rho(j)*sqrt(G(j)*(1+2*I*D(j))/rho(j)))/...
            (rho(j+1)*sqrt(G(j+1)*(1+2*I*D(j+1))/rho(j+1))); % complex impedance ratio
    end
    if strcmp(boundary,'rigid')
        alpha(nr_layer-1) = 0; % =1: rigid bedrock, =0: elastic bedrock
    end
    % 2.1.2: Complex shear wave velocity (Kramer, p260)
    vs_star = sqrt(G.*(1+2*I*D)./rho);  % nr_layer x 1
    
    % 2.1.3: Complex wave number (Kramer, p260)
    
    %k_star = zeros(N/2+0.5, nr_layer); % preallocation
    %for k = 1 : 1: nr_layer  % k: layer index
    %    for j = 1 : 1 : N/2+0.5  % j: frequency index
    %        k_star(j,k) = omega(j)/vs_star(k);
    %    end
    %end
    
    vs_star_recip = 1./vs_star;  % element-wise reciprocal, nr_layer x 1
    k_star = omega(1 : N/2+0.5) * transpose(vs_star_recip);  % (N/2+0.5) x nr_layer
    
    % 2.1.3 Calculate A & B (Kramer, p269)
    A = zeros(N/2+0.5,nr_layer);
    B = zeros(N/2+0.5,nr_layer);
    A(:,1) = 1;
    B(:,1) = 1;
    for k = 1 : 1 : nr_layer-1  % layer index
        
        %------------ vectorized version ----------------
        A(:,k+1) = 0.5*A(:,k) * (1+alpha(k)) .* exp(I*k_star(:,k)*h(k))...
                  +0.5*B(:,k) * (1-alpha(k)) .* exp(-I*k_star(:,k)*h(k)); % left half
        B(:,k+1) = 0.5*A(:,k) * (1-alpha(k)) .* exp(I*k_star(:,k)*h(k))...
                  +0.5*B(:,k) * (1+alpha(k)) .* exp(-I*k_star(:,k)*h(k)); % left half
              
        %for j = 1 : 1: N/2+0.5  % frequency index
        %    A(j,k+1) = 0.5*A(j,k)*(1+alpha(k))*exp(I*k_star(k,j)*h(k))...
        %        +0.5*B(j,k)*(1-alpha(k))*exp(-I*k_star(k,j)*h(k)); % left half
        %    B(j,k+1) = 0.5*A(j,k)*(1-alpha(k))*exp(I*k_star(k,j)*h(k))...
        %        +0.5*B(j,k)*(1+alpha(k))*exp(-I*k_star(k,j)*h(k)); % left half
        %end
    end
    % 2.1.4 Calculate linear transfer function
    H_ss = zeros(N/2+0.5,nr_layer); % single-sided transfer function
    H_append = zeros(N/2-0.5,nr_layer); % the other half
    for k = 1 : 1 : nr_layer
        H_ss(:,k) = (A(:,k) + B(:,k)) ./ A(:,nr_layer); % Transfer Function
        
        %for j = 1 : 1 : N/2+0.5
        %    H_ss(j,k) = (A(j,k)+B(j,k))/A(j,nr_layer); % Transfer Function
        %end
               
        H_ss(1,k) = real(H_ss(1,k)); % see note (1) below
        H_append(:,k) = flipud(H_ss(2:end,k));
        H_append(:,k) = conj(H_append(:,k));
    end

    H = [H_ss;H_append];
    
    H = H(1:f_oversample_factor:end,:);  % downsample every layer of H (added 4/28/2016)
    f = f(1:f_oversample_factor:end); % downsample f (added 4/28/2016)
    N = N / f_oversample_factor; % scale down N (added 4/28/2016)

    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    % Notes: 
    % (1) The value of transfer function when freq = 0 should be either 1
    %     or 2 (a real number). Because of discretization errors, tf_ss(1)
    %    is not a real number, but very close.
    %
    % (2) Two examples:
    %    [a] fft([1 2 3 4 5 6 7 8]') = 
    %          36
    %          -4 + 9.6569i
    %          -4 + 4.0000i
    %          -4 + 1.6569i
    %          -4
    %          -4 - 1.6569i
    %          -4 - 4.0000i
    %          -4 - 9.6569i
    %
    %    [b] fft([1 2 3 4 5 6 7]') = 
    %          28.0          
    %          -3.5 + 7.2678i
    %          -3.5 + 2.7912i
    %          -3.5 + 0.7989i
    %          -3.5 - 0.7989i
    %          -3.5 - 2.7912i
    %          -3.5 - 7.2678i
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    %% Part 3: Response acceleration of each layer (from input motion & linear TF)
    %ACCEL_OUT = zeros(N,nr_layer);
    % peak = zeros(nr_layer,1); % max. accel. in each layer (absolute value)
    % accel = zeros(n,nr_layer); 
	%accel_out = zeros(n,nr_layer); % n by nr_layer array, each column for one layer
	%veloc = zeros(n,nr_layer);
    %displ = zeros(n,nr_layer);
	
    ACCEL_OUT = H .* repmat(ACCEL_IN, [1,nr_layer]);  % (N x nr_layer) * (N x nr_layer)
    %ACCEL_OUT = H .* ACCEL_IN;  % no repmat needed for newer MATLAB versions
    accel_out = ifft(ACCEL_OUT);  % column-wise ifft (each column is a layer)
    veloc = cumsum(accel_out) * dt;  % perform cumsum to each column
    displ = cumsum(veloc) * dt;  % perform cumsum to each column
    
    % baseline correction of displacement time history:
    offset = (1:1:n)'/(n+1) * displ(n,:);  % (N x 1) * (1 * nr_layer) = (N x nr_layer)
    displ = displ - offset;  % N x nr_layer
    
    %for k = 1 : 1 : nr_layer
    %    for j = 1 : 1 : N
    %        ACCEL_OUT(j,k) = H(j,k)*ACCEL_IN(j); % Calculate the response in the frequency domain
    %    end
    %    
    %    accel_out(:,k) = ifft(ACCEL_OUT(:,k)); % accel_out is guaranteed to be real
    %    % peak(k) = max(abs(accel(:,k)));
    %    %%% veloc(:,k) = integr(accel_out(:,k),dt);
    %    %%% displ(:,k) = integr(veloc(:,k),dt);
    %    veloc(:,k) = cumsum(accel_out(:,k)) * dt; % numerical integration
    %    displ(:,k) = cumsum(veloc(:,k)) * dt;
    %    % baseline correction of displacement time history
    %    for j = 1 : 1 : n
    %        displ(j,k) = displ(j,k)-j/(n+1)*displ(n,k);
    %    end
    %end
    
    %% Part 4: Strain time history of each layer
    strain = zeros(n,nr_layer-1);
    eff_strain = zeros(nr_layer-1,1); % effective strain (Kramer, pp.271-272)
    for k = 1 : 1: nr_layer-1
        strain(:,k) = (displ(:,k)-displ(:,k+1))/h(k); % strain at midpoint of each layer (unit: 1)
        eff_strain(k) = R_gamma * max(abs(strain(:,k))); % unit of strain: 1 (not percent)
    end

    %% Part 5: New modulus & damping ratio based on effective strain
    G_new = zeros(nr_layer-1,1);
    D_new = zeros(nr_layer-1,1);    
    G_relative_diff = zeros(nr_layer-1,1);
    D_relative_diff = zeros(nr_layer-1,1);
    for k = 1 : 1: nr_layer-1
        strainG_interp = [0;strain_G(:,k);1e2];
        Gvector_interp = [1;G_vector(:,k);min(1e-4,G_vector(end,k))];
        strainD_interp = [0;strain_D(:,k);1e2];
        Dvector_interp = [D_vector(1,k);D_vector(:,k);D_vector(end,k)];
        G_new(k) = G_max(k) * interp1(strainG_interp,Gvector_interp,eff_strain(k),'linear');  % manually set strain range
                %  ^ Interpolation needs to start from G_max(k), otherwise
                %  the shear modulus values would get smaller and smaller
                %  and eventually to 0.
        D_new(k) = interp1(strainD_interp,Dvector_interp,eff_strain(k),'linear'); % and no extrapolation is allowed
        %G_relative_diff(k) = abs(G(k)-G_new(k)) / G_new(k);
        %D_relative_diff(k) = abs(D(k)-D_new(k)) / D_new(k);
    end
    G_relative_diff = abs(G(1:end-1) - G_new) ./ G_new;
    D_relative_diff = abs(D(1:end-1) - D_new) ./ D_new;
    G(1:end-1) = G_new;
    D(1:end-1) = D_new;
    G_mtrx(:,i_iter+1) = G_new;
    D_mtrx(:,i_iter+1) = D_new;
    fprintf('  G_diff = %7.2f%%, D_diff = %7.2f%%\n',max(G_relative_diff)*100,max(D_relative_diff)*100);
    
    %% Part 6: End of iteration check
    if max(max(G_relative_diff),max(D_relative_diff)) < tolerance
        fprintf('---------- Convergence achieved ----------\n');
        break;
    end
end

%% Part 7.1: Calculate stress from strain
%strain_fft = zeros(N,nr_layer-1);  % N is guarenteed to be odd
stress_fft = zeros(N,nr_layer-1);  % N is guarenteed to be odd
%stress = zeros(N,nr_layer-1);

%for j = 1 : 1 : nr_layer-1
%    strain_fft(:,j) = fft(strain(:,j));
%    stress_fft(1:N/2+0.5,j) = strain_fft(1:N/2+0.5,j)*(G(j)*(1+2*I*D(j)));
%    for k = 1:N/2-0.5
%        stress_fft(N+1-k,j) = conj(stress_fft(k+1,j));
%    end
%    stress_fft(1,j) = real(stress_fft(1,j));
%    stress(:,j) = ifft(stress_fft(:,j));
%end

strain_fft = fft(strain);  % column-wise fft, N x (nr_layer-1)
modulus = transpose(G .* (1 + 2*I*D));  % 1 x nr_layer
modulus = modulus(1:end-1);  % 1 x (nr_layer-1)
modulus_repmat = repmat(modulus, [N/2+0.5,1]);
%modulus_repmat = modulus;  % no repmat needed for newer MATLAB versions
stress_fft(1:N/2+0.5,:) = strain_fft(1:N/2+0.5,:) .* modulus_repmat;
stress_fft(N/2+0.5+1:end, :) = flipud(conj(stress_fft(2:N/2+0.5, :)));
stress_fft(1,:) = real(stress_fft(1,:));
stress = ifft(stress_fft);

%% Part 7.2: Odd-even check
if flag == 0 % if originally the input motion length is even
    freq_array = f(1:N/2-0.5);
    tf = H(1:N/2-0.5,1)/2;
    af = abs(H(1:N/2-0.5,1))/2;
    t_out = t(1:end-1);
	accel_out = accel_out(1:end-1,:); % ditch the whole last row
	veloc = veloc(1:end-1,:);
	displ = displ(1:end-1,:);
	strain = strain(1:end-1,:);
	stress = stress(1:end-1,:);
    accel_on_surface = accel_out(:,1); % the first column
else  % flag == 1
    freq_array = f(1:N/2+0.5);
    tf = H(1:N/2+0.5,1)/2;
    af = abs(H(1:N/2+0.5,1))/2;
    t_out = t;
    accel_on_surface = accel_out(:,1);
end

%% Plot equivalent linear transfer function
if strcmpi(fig_visible,'on')
    f0 = findF0( [ freq_array, af ] );
    h_tf = figure('visible',fig_visible);
    semilogx(freq_array,af,'k','linewidth',1.5);
    grid on;
    xlim([0 max(freq_array)]); xlabel('Frequency (Hz)');
    if strcmp(boundary,'elastic')
        title('Equivalent Linear AF (S to Outcrop)');
    end
    if strcmp(boundary,'rigid')
        title('Equivalent Linear TF (S to Borehole)');
    end
    annotation('textbox',[.15 .8 .1 .1],'String',sprintf('f_0 = %.2f Hz',f0),'BackgroundColor','w');
end

%% Plot input and output accelerations together
if strcmpi(fig_visible,'on')
    h_accel = figure('visible',fig_visible);
    plot(t_out,accel_on_surface,'b',t,accel_in,'r');
    legend('Output','Input','location','northeast');%,'BackgroundColor','w');
    grid on;
    xlim([0 max(t)]); xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
end

%% Observe modulus reduction process and damping ratio changing process
% % figure('visible','on');plot((1:1:nr_iter+1)',G_mtrx(:,:),'.-'); xlim([1 i_iter]);
% % figure('visible','on');plot((1:1:nr_iter+1)',D_mtrx(:,:),'.-'); xlim([1 i_iter]);

%% Part 7: Post processing
out_gamma = strain;
out_tau = stress;
out_a = accel_out;
out_v = veloc;
out_d = displ;

% max_v = zeros(nr_layer,1);
% max_d = zeros(nr_layer,1);
% max_gamma = zeros(nr_layer-1,1);
% max_tau = zeros(nr_layer-1,1);
% for j = 1 : nr_layer-1
    % max_tau(j) = max(abs(out_tau(:,j)));
    % max_gamma(j) = max(abs(out_gamma(:,j)));
% end
% for j = 1 : nr_layer
    % max_v(j) = max(abs(out_v(:,j)));
    % max_d(j) = max(abs(out_d(:,j)));
% end

max_a = max(abs(accel_out)).'; % (absolute) maximum value of each column of accel_out
max_v = max(abs(veloc)).';  % nr_layer x 1 (after transposition)
max_d = max(abs(displ)).';  % nr_layer x 1 (after transposition)
max_gamma = max(abs(out_gamma)).'; % (nr_layer-1) x 1 (after transposition)
max_tau = max(abs(out_tau)).'; % (nr_layer-1) x 1 (after transposition)

max_avd = [layer_boundary_depth,max_a,max_v,max_d];
max_gt = [layer_midpoint_depth,max_gamma,max_tau];

%%
if strcmp(output_or_not,'y')
    dlmwrite(fullfile(file_dir,'equiv_linear_tf.dat'),[ f, abs(H(1,:))'/2 ],'delimiter','\t');
    saveas(h_tf,fullfile(file_dir,'equiv_linear_tf.png'));
    
    dlmwrite(fullfile(file_dir,'accel_on_surface.dat'),[ t_out, accel_on_surface ],'delimiter','\t');
    saveas(h_accel,fullfile(file_dir,'accel_on_surface.png'));

%     accel = applyBandpass(out_a(:,1),0.1,15,0.01);
%     accel = taperTukey(accel);
%     surf = [t,accel];
%     dlmwrite(fullfile(file_dir,'surface_accel.dat'),surf,'delimiter','\t');
%     h = figure('visible','off');
%     plot(t,accel); title('Acceleration (m/s^2) time history @ surface  ');
%     saveas(h,fullfile(file_dir,'acc_surface.png'));
    
    strn2 = [t,out_gamma(:,2)];
    strn4 = [t,out_gamma(:,4)];
    strn8 = [t,out_gamma(:,8)];
    strs2 = [t,out_tau(:,2)];
    strs4 = [t,out_tau(:,4)];
    strs8 = [t,out_tau(:,8)];
    dlmwrite(fullfile(file_dir,'strain2.dat'),strn2,'delimiter','\t');
    dlmwrite(fullfile(file_dir,'strain4.dat'),strn4,'delimiter','\t');
    dlmwrite(fullfile(file_dir,'strain8.dat'),strn8,'delimiter','\t');
    dlmwrite(fullfile(file_dir,'stress2.dat'),strs2,'delimiter','\t');
    dlmwrite(fullfile(file_dir,'stress4.dat'),strs4,'delimiter','\t');
    dlmwrite(fullfile(file_dir,'stress8.dat'),strs8,'delimiter','\t');
    
    h=figure('visible','off');
    subplot(311); plot(t,strn2(:,2)*100); title('strain (%) time history @ depth=7.5m  ');
    subplot(312); plot(t,strn4(:,2)*100); title('strain (%) time history @ depth=15.5m  ');
    subplot(313); plot(t,strn8(:,2)*100); title('strain (%) time history @ depth=29.5m  ');
    saveas (h,fullfile(file_dir,'strain_his.png'));
    close(h);
    
    h=figure('Visible','off');
    subplot(311); plot(t,strs2(:,2)/1000); title('stress (kPa) time history @ depth=7.5m  ');
    subplot(312); plot(t,strs4(:,2)/1000); title('stress (kPa) time history @ depth=15.5m  ');
    subplot(313); plot(t,strs8(:,2)/1000); title('stress (kPa) time history @ depth=29.5m  ');
    saveas (h,fullfile(file_dir,'stress_time_history.png'));
    close(h);
    
    h=figure('Visible','off');
    subplot(321); plot(strn2(:,2)*100,strs2(:,2)/1000); title('hysteretic loop [stress(kPa) vs strain (%)] @ depth=7.5m  ');
    subplot(323); plot(strn4(:,2)*100,strs4(:,2)/1000); title('hysteretic loop [stress(kPa) vs strain (%)] @ depth=15.5m  ');
    subplot(325); plot(strn8(:,2)*100,strs8(:,2)/1000); title('hysteretic loop [stress(kPa) vs strain (%)] @ depth=29.5m  ');
    saveas (h,fullfile(file_dir,'hysteresis_loop.png'));
    close(h);
    
    % [n_max,nlayer] = size(out_a);
    % acc_f = zeros(n_max,nlayer);
    % max_a = zeros(1,nlayer);
    % for j = 1:nlayer
        % acc_f(:,j) = applyBandpass(out_a(:,j),0.1,15,0.01); % note: same as Wei's "getfilter.m" function
        % acc_f(:,j) = applyTaper(acc_f(:,j)); % note: same as Wei's "taper.m" function
        % for k = 1:n_max
            % if abs(acc_f(k,j)) > max_a(j)
                % max_a(j) = abs(acc_f(k,j));
            % end
        % end
    % end
    % h=figure('Visible','off');
    % subplot(231); plot(max_gamma*100,0-layer_depth,'-bs'); xlabel('Max. Shear Strain (%)'); ylabel('Depth (m)');
    % subplot(232); plot(max_tau/1000,0-layer_depth,'-bs'); xlabel('Max. Shear Stress (kPa)'); ylabel('Depth (m)');
    % subplot(234); plot(max_a/100,0-node_depth,'-bs'); xlabel('Max. Accel. (cm/s^2)'); ylabel('Depth (m)');
    % subplot(235); plot(max_v*100,0-node_depth,'-bs'); xlabel('Max. Relative Veloc. (cm/s)'); ylabel('Depth (m)');
    % subplot(236); plot(max_d*100,0-node_depth,'-bs'); xlabel('Max. Relative Displ. (cm)'); ylabel('Depth (m)');
    % saveas(h,fullfile(file_dir,'iteration.png'));
    % close(h);
end

end

function y = applyTaper(x)
% Taper through 5% tukeywin
    n = length(x);
    y = x.*tukeywin(n,0.05);
end

function [index,maxValue] = findF0(x) % Local function: find f_0
if isvector(x) == 1
    maxValue = max(x);
    l = length(x);
    for ii = 1 : 1 : l
        if x(ii) >= maxValue
            break
        end
    end
    index = ii;
else
    if min(size(x)) > 2
        tolerance('Required dimensions are 1x$, $x1, 2x$ or $x1.');
    end
    if size(x,1) == 2
        maxValue = max(x(2,:));
        l = length(x(2,:));
        for ii = 1 : 1 : l
            if x(2,ii) >= maxValue
                break
            end
        end
        index = x(1,ii);
    end
    if size(x,2) == 2
        maxValue = max(x(:,2));
        l = length(x(:,2));
        for ii = 1 : 1 : l
            if x(ii,2) >= maxValue
                break
            end
        end
        index = x(ii,1);
    end
end

end

function [acc_f] = applyBandpass(acc,low,high,dt)
% bandpass filter for a time series
%
% function [acc_f]=getfilter(acc,low,high,dt)
%   acc_f: filtered time series
%   acc:   input time series
%   low:   low cut frequency
%   high:  high cut frequency
%   dt:    time step of input time series

nyq = 1/(2.0*dt); % nyquist frequency
n = 4; 
wn = [low,high]/nyq;

[b,a] = butter(n,wn);  % Filter coefficients (numerator/denominator)

hdf1 = dfilt.df2(b,a);
acc_f = filter(hdf1,acc);

end

function y = integr(x,dt)
dt2 = dt/2;
y = zeros(size(x));
y(1) = 0;
for i = 2:length(x)
    y(i) = y(i-1) + (x(i-1)+x(i))*dt2;
end
end



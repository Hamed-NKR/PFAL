clc
clear
close all

% location of previously saved aggregate data
fdir_in = 'D:\Hamed\CND\PhD\Publication\DLCA2\mainscatter_sigmapp10\FLAT';
fname_in = 'FLAT-26NOV24';

varnames = {'pars_out'}; % varaiables to be imported

% load aggregate data
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

% === Load Aggregate Data ===
pp = pars_out(1).pp{503};  % aggregate number 503

d = pp(:,2);  % diameters
x = [pp(:,3), pp(:,4), pp(:,5)];  % positions

figure(1)
UTILS.PLOTPP(x(:,1),x(:,2),x(:,3),d)

% === Normalize Radii and Positions ===
r = d / 2;
rm = mean(r);
r = r / rm;
x = x / rm;
x = x - mean(x,1);  % center the aggregate

% === Collapse Simulation Parameters ===
k = 0.2;
k_decay = 1;
eps = 0.5;
v_decay = 0.98;
jit = 0.0;
jit_decay = 0.98;
Dt = 0.01;
nt = 1000;
nj = 8;

N = size(x,1);
sig = (r + r') / (2^(1/6));
v = zeros(N,3);
m = (4/3) * pi * r.^3;

for jj = 1:nj
    for ii = 1:nt
        jit = jit * jit_decay;
        k = k * k_decay;

        cen = sum(x .* r, 1) ./ sum(r);
        F = -k * (x - cen);

        d0 = squareform(pdist(x));
        d0(d0 == 0) = Inf;

        dvec = permute(x, [1 3 2]) - permute(x, [3 1 2]);
        ds = sqrt(sum(dvec.^2, 3));
        ds(ds == 0) = Inf;

        d_unit = dvec ./ ds;
        rmin = 0.3;
        d0_clamped = max(d0, rmin);

        Fvw_mag = 48 * eps * ((sig ./ d0_clamped).^12 ./ d0_clamped - 0.5 * (sig ./ d0_clamped).^6 ./ d0_clamped);
        Fvw_mag(d0 > 1.5 .* sig) = 0;

        Fvw_mag_exp = repmat(Fvw_mag, 1, 1, 3);
        Fvw = squeeze(sum(Fvw_mag_exp .* d_unit, 2));
        F = F + Fvw + jit * randn(N,3);

        a = F ./ m;
        v = v + a * Dt;
        x = x + v * Dt;
        v = v * v_decay;
    end
end

figure(2)
UTILS.PLOTPP(x(:,1),x(:,2),x(:,3),2*r)
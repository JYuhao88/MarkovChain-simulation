clear;
% simulation parameter
T = 3000;
traj = 80000;
state = zeros(1, traj);
next = state;

%% model parameter
rho = 4e-1;
d = 1e-1;
r = rho/d;
p = 0.2;
q = 1- p;

ptime = zeros(1,traj);
jt = 0;
while min(ptime) <T
    pjump = ptime <T;
    state(pjump) = next(pjump);
    lamda = 1./(rho*(1-p) + d*state(pjump));
    waitt = exprnd(lamda);
    ptime(pjump) = ptime(pjump) + waitt;
    
    dice = rand([1, traj]);
    pback = state .*d ./(rho*(1-p)+state*d);
    bjump = dice < pback;
    fjump = ~bjump;
    dice_ = (dice .*(state.*d + rho*(1-p))-state.*d)/rho + p;

    fstep = geoinv(dice_, p);
    fjump = fjump &pjump;
    bjump = bjump &pjump;
    next(fjump) = state(fjump) + fstep(fjump);
    next(bjump) = state(bjump) - 1;
    next(~pjump) = state(~pjump);
end
sta = tabulate(state);

% plot
close all
figure;
bar(sta(:, 1), sta(:, 3)/100);
hold on;
x = 0:max(state);
% y = binopdf(x, x+r-1, 1-p)*p;
y = gamma(r+x)./(gamma(r).*gamma(x+1)).*(1-p).^x.*p.^r;
plot(x, y, '-o');
title('single-death');


clear;
% simulation parameter
T = 3000;
traj = 80000;
state = zeros(1, traj);
next = state;

%% model parameter
rho = 5e-1;
d = 2e-1;
r = rho/d;

ptime = zeros(1,traj);
while min(ptime) < T
    pjump = (ptime < T);
    state(pjump) = next(pjump);
    lamda = 1./(rho + d*state(pjump));
    waitt = exprnd(lamda);
    ptime(pjump) = ptime(pjump) + waitt;
    
    dice = rand([1, traj]);
    fjump = (dice<= rho./(rho+state*d));
    bjump = ~fjump;
    
    fjump = fjump &pjump;
    bjump = bjump &pjump;
    next(fjump) = state(fjump) + 1;
    next(bjump) = state(bjump) - 1;
    next(~pjump) = state(~pjump);
end

% while min(ptime) < T
%     pjump = (ptime < T);
%     state(pjump) = next(pjump);
%     
%     ftime = exprnd(1./(rho*ones(1, traj)));
%     btime = exprnd(1./(state*d));
%     btime(btime==0) = 1e5;
%     
%     waitt = min([ftime; btime], [], 1);
%     ptime(pjump) = ptime(pjump) + waitt(pjump);
%     
%     fjump = (ftime<btime);
%     bjump = ~fjump;
%     
%     fjump = fjump &pjump;
%     bjump = bjump &pjump;
%     next(fjump) = state(fjump) + 1;
%     next(bjump) = state(bjump) - 1;
%     next(~pjump) = state(~pjump);
% end

sta = tabulate(state);

% plot
close all
figure;
bar(sta(:, 1), sta(:, 3)/100);
hold on;
x = 0:max(state);
%y = r.^x.*exp(-r)./factorial(x);
y = poisspdf(x, rho/d);
plot(x, y, '-o');
title('birth-death');


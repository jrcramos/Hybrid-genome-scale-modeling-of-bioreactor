% this codes simulates 9 synthetic cultivation conditions using a dynamic
% model in "odefun.m"

close all
load pars.mat
pars(end+1:end+21)=0; %set feed as 0 initially in the parameter vector
for ibatch=1:9
    [tspan, x0, f]=get_batch_condition(ibatch); % get time, initial condition and feed 

    pars(end-20:end)=f; % update for each experiment
    odeopts=odeset('AbsTol',1e-8,'RelTol',1e-7,'NonNegative',1:47);
    ofun=@(t,s)odefun(t,s,pars);
    % codegen odefun -args {0,zeros(47,1,'double'),zeros(182,1,'double')}
    % ofun=@(t,s)odefun_mex(t,s,pars);
    [t,y]=ode15s(ofun,tspan,x0,odeopts);
    figure('Name',sprintf('Br%d [extracellular states]',ibatch),'units','normalized','outerposition',[0 0.1 0.5 0.75])
    for i=1:25
        subplot(5,5,i)
        plot(t,y(:,i),'b')
    end
end
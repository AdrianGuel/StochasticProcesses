%Information length from stochastic simulation of the OU process
%      dx = theta * ( mu - x(t) ) dt + sigma dW,   
%      x(0) = x0,
% using an initial gaussian distribution p(x,0)~N(0,sigma)

close all
clear all
clc

set(gcf,'color','w');

ax1 = subplot(2,2,[2 4]);
%set(ax1,'YScale','log')
%set(ax1,'XScale','log')
hold(ax1,'on')
grid(ax1,'on')
xlabel(ax1,'$x\sim \mathcal{N}(\langle x\rangle,\Sigma)$','Interpreter','Latex','FontSize', 14)
ylabel(ax1,'$p(x;t)$','Interpreter','Latex','FontSize', 14)
xlim(ax1,[-0.5 2])

ax2 = subplot(2,2,1);
%set(ax2,'YScale','log')
%set(ax1,'XScale','log')
hold(ax2,'on')
grid(ax2,'on')
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax2,'$\mathcal{E}(t)$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(2,2,3);
%set(ax3,'YScale','log')
%set(ax1,'XScale','log')
hold(ax3,'on')
grid(ax3,'on')
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax3,'$\mathcal{L}(t)$','Interpreter','Latex','FontSize', 14)


%Stochastic simulation
theta=3;
mu=0;
sigma=0.1;
D=(sigma^2);
x0=1;
N=1e5;
tmax=5;
[t,x]=ornstein_uhlenbeck_euler_maruyama ( theta, mu, sigma, x0, tmax, N);

%deterministic solution of mu and sigma
%tspan = [0 5];   

opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) Sigmax(t,y,theta,D), t,  sigma, opts);
[t,M] = ode45(@(t,y) Mux(t,y,theta), t, x0, opts);
%plot ( t, x, 'k-' )


xaux = linspace(-0.5,2,N);
xi=zeros(length(xaux),length(t));
f=zeros(length(xaux),length(t));
z=zeros(length(xaux),length(t));

for k=1:length(t)
    %[f(:,k),xi(:,k)] = ksdensity(x(:,k),xaux);
    pd=fitdist(x(:,k),'normal');
    f(:,k)= pdf(pd,xaux);
    z(:,k) = normpdf(xaux,M(k),sqrt(y(k)));
    if mod(k,100)==0
      %plot(ax1,xi(:,k),f(:,k),'k',xaux,z(:,k),'b');
      plot(ax1,xaux,f(:,k),'k',xaux,z(:,k),'b--');
    end
end

% E=zeros(1,length(t));

Ts=diff(t);
% for k=2:length(t)
%     E(k) =  4*(sum(sqrt(z(:,k))-sqrt(z(:,k-1)))/Ts(1))^2;
% end

%%Computation of IL
[fxz,fyz] = gradient(sqrt(z),Ts(1));
Ez=trapz(xaux,4*fxz.^2,1);
ILz=cumtrapz(t,Ez);

[fx,fy] = gradient(sqrt(f),Ts(1));
E=trapz(xaux,4*fx.^2,1);
IL=cumtrapz(t,E);

Et=((-theta*M).^2)./y+0.5*((-2*theta*y+D)./y).^2;
ILt=cumtrapz(t,Et);
%E=cumtrapz(xi,f)
plot(ax2,t,E,'r',t,Et,'k',t,Ez,'b--')
plot(ax3,t,IL,'r',t,ILt,'k',t,ILz,'b--')

leg1 = legend(ax1,{'Estimated','Theoretical'});
leg2 = legend(ax3,{'Estimated','Theoretical a','Theoretical b'});

function dydt = Sigmax(t,y,theta,D)
   dydt=-2*theta*y+D;
end
 
function dydt = Mux(t,y,theta)
   dydt=-theta*y;
end

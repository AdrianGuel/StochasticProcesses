%Information length from stochastic simulation of Kramers equation
%      dx = y dt + sqrt(D_11) dW,   
%      dy = (-omega²*x-gamma*y) dt + sqrt(D_22) dW,
%      x(0) = x0, y(0) = y0,
% using an initial gaussian distribution p(x,y,0)~N(0,Sigma)

close all
clear all
clc

set(gcf,'color','w');

ax1 = subplot(2,2,[2 4]);
%set(ax1,'YScale','log')
%set(ax1,'XScale','log')
hold(ax1,'on')
grid(ax1,'on')
xlabel(ax1,'$x\sim \mathcal{N}(\mu_x, \mu_y,\Sigma)$','Interpreter','Latex','FontSize', 14)
ylabel(ax1,'$y\sim \mathcal{N}(\mu_x, \mu_y,\Sigma)$','Interpreter','Latex','FontSize', 14)
xlim(ax1,[-1 3])

ax2 = subplot(2,2,1);
%set(ax2,'YScale','log')
%set(ax1,'XScale','log')
hold(ax2,'on')
grid(ax2,'on')
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)
%ylabel(ax2,'$\mathcal{E}(t)$','Interpreter','Latex','FontSize', 14)
ylabel(ax2,'$x(t)$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(2,2,3);
%set(ax3,'YScale','log')
%set(ax1,'XScale','log')
hold(ax3,'on')
grid(ax3,'on')
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)
%ylabel(ax3,'$\mathcal{L}(t)$','Interpreter','Latex','FontSize', 14)
ylabel(ax3,'$y(t)$','Interpreter','Latex','FontSize', 14)

%Stochastic simulation
omega=1;
gamma=2;
n=2;
A=[0,1;-omega^2,-gamma];
sigma0=[0.01 0; 0 0.01];
D=[0.1 0; 0 0.1];
x0=1;
y0=2;
X0=[x0,y0];
N=1e2;
tmax=5;
[t,x,y]=Kramers_euler_maruyama (omega,gamma, D, x0, y0, sigma0, tmax, N);    

% 
% %deterministic solution of mu and sigma
% %tspan = [0 5];   
% 
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t,yx] = ode45(@(t,y) Sigmax(t,y,A,D,n), t,  sigma0, opts);
[t,M] = ode45(@(t,y) Mux(t,y,A), t, X0, opts);


plot(ax1,x,y,'k.',M(:,1),M(:,2),'b')
plot(ax2,t,x,'k')
plot(ax3,t,y,'k')


x1 = linspace(-1,3,100);
x2 = linspace(-1,3,100);
[X1,X2] = meshgrid(x1,x2);
xi = [X1(:) X2(:)]; 
f=zeros(length(x1)*length(x2),length(t));
for k=1:length(t)
    f(:,k)= ksdensity([x(:,k) y(:,k)],xi);
%    if mod(k,100)==0
%        plot3(xi(:,1),xi(:,2),f,'k')
%         X = reshape(ep(:,1),length(x1),length(x2));
%         Y = reshape(ep(:,2),length(x1),length(x2));
%         Z = reshape(f,length(x1),length(x2));
%         contour (ax1,X,Y,Z,1,'k')
end   

% 
% % E=zeros(1,length(t));
% 
% Ts=diff(t);
% % for k=2:length(t)
% %     E(k) =  4*(sum(sqrt(z(:,k))-sqrt(z(:,k-1)))/Ts(1))^2;
% % end
% 
% %%Computation of IL
% [fxz,fyz] = gradient(sqrt(z),Ts(1));
% Ez=trapz(xaux,4*fxz.^2,1);
% ILz=cumtrapz(t,Ez);
% 
% [fx,fy] = gradient(sqrt(f),Ts(1));
% E=trapz(xaux,4*fx.^2,1);
% IL=cumtrapz(t,E);
% 
% Et=((-theta*M).^2)./y+0.5*((-2*theta*y+D)./y).^2;
% ILt=cumtrapz(t,Et);
% %E=cumtrapz(xi,f)
% plot(ax2,t,E,'r',t,Et,'k',t,Ez,'b--')
% plot(ax3,t,IL,'r',t,ILt,'k',t,ILz,'b--')
% 
% leg1 = legend(ax1,{'Estimated','Theoretical'});
% leg2 = legend(ax3,{'Estimated','Theoretical a','Theoretical b'});
% 
function dydt = Sigmax(t,y,A,D,n)
   At=A';
   aux2=reshape(y,2,2);
   dydt=kron(eye(n,n),A)*aux2(:)+kron(eye(n,n),aux2)*At(:)+2*D(:);
end
%  
function dydt = Mux(t,y,A)
   dydt=A*y;
end

%cubic stochastic differential equation under the laplace assumption
%Adrian Guel 2022
%Reference of the Laplace assumption "Population dynamics under the laplace
%assumption", Andre C. Marreiros et al. 2009

clearvars;
close all;
clc;

fig=figure('visible','on');
set(fig, 'Position',  [615,328,800,454])
set(gcf,'color','w');
ax1 = subplot(3,2,[2 4]);
hold(ax1,'on')
grid(ax1,'on')
xlabel(ax1,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax1,'$x$','Interpreter','Latex','FontSize', 14)
zlabel(ax1,'$p(\mathbf{x};t)$','Interpreter','Latex','FontSize', 14)
axis(ax1,'square')
view(ax1,3)

ax2 = subplot(3,2,1);
hold(ax2,'on')
grid(ax2,'on')
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax2,'$\mu(t)$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(3,2,3);
hold(ax3,'on')
grid(ax3,'on')
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax3,'$\Sigma(t)$','Interpreter','Latex','FontSize', 14)

ax4 = subplot(3,2,5);
hold(ax4,'on')
grid(ax4,'on')
xlabel(ax4,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax4,'$\Gamma^2(t)$','Interpreter','Latex','FontSize', 14)

ax5 = subplot(3,2,6);
hold(ax5,'on')
grid(ax5,'on')
xlabel(ax5,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax5,'$\mathcal{L}(t)$','Interpreter','Latex','FontSize', 14)


y0=[0.5;1e-1];
D=1e-3;
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);
tspan=[0 25];
[t,L] = ode45(@(t,y) LaplacianA(t,y,D), tspan, y0, opts);


    plot(ax2,t,L(:,1))
    plot(ax3,t,L(:,2))
    
    xaux = -2:.01:2;
    z=[];l=1;
    for k=1:10:length(t)
    z(:,l) = normpdf(xaux,L(k,1),sqrt(L(k,2)));
    p1=plot3(ax1,t(k)*ones(1,length(xaux)),xaux,z(:,l),'k');
    %p1=plot(xaux,z(:,l),'k');
    p1.Color(4) = 0.05;
    hold on
    l=l+1;
    end
    E=((-L(:,1).^3-3*L(:,1).*L(:,2)).^2)./L(:,2)+0.5*((-6*L(:,2).*L(:,1).^2+2*D)./L(:,2)).^2;
    IL=cumtrapz(t,sqrt(E));
    plot(ax4,t,E)
    plot(ax5,t,IL)

function dydt=LaplacianA(t,y,D)
    dydt = zeros(2,1);
    dydt(1)=-y(1)^3-3*y(1)*y(2);
    dydt(2)=-6*y(2)*(y(1)^2)+2*D;
end

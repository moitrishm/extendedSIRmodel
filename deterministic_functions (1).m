%%%%%%%%%%%21 March 2021%%%%%%%
function deterministic_functions(n)
f1=figure(1);

set(f1,'PaperUnits','centimeters','Units','centimeters','PaperPosition',[0 0 80 50],'Position',[0 0 80 50])

figure('GraphicsSmoothing','on');

set(0, 'defaultaxesfontsize', 12);

set(0, 'defaulttextfontsize', 12);

set(0, 'defaultlinelinewidth', 1.5);

mu_0=0.2;
T1=0.4;    %the width of the first peak in the double square
T2=0.1;    %the width of the second peak in the double square
T=T1+T2;   %the width of the usual square function
first=0.25;    %the starting point of the first peak (in terms of the proportion of year)
second=0.75; %the starting point of the second peak
p=200;      %the shape parameter for the family of functions
a=0.01;
t1=10*365;

ft = linspace(0, t1,t1*10); %pulse simulated for 5 years
switch n
    case 'sinusoidal'
        f = (mu_0)*(1+a*cos(2*pi*ft/365));
    case 'square'
        f = (mu_0)*abs(square((2*pi*(ft))/(365) + (2.5*pi) , T*100))+(mu_0*a)*square(2*pi*(ft)/(365) +(2.5*pi) , T*100);
    case 'double_square'
        f =  (mu_0)*abs(square((2*pi*(ft))/(365)- (first*10*pi), T1*100))+(mu_0*a)*square(2*pi*(ft)/(365) -(first*10*pi) , T1*100)+(mu_0*a)*(square((2*pi*ft/365) - (second*10*pi), T2*100)+1);
    case 'family'
        f = (mu_0) * (1 + a*sign(cos(2*pi*ft/365)).*abs(cos(2*pi*ft/365)).^(p));
    case 'wide_sine'
        f = ((mu_0)*(1-2*a*((sin(pi*ft./365)).^(8))))+0.0009;
    case 'narrow_sine'
        f = ((mu_0)*(1+2*a*((sin(pi*ft./365)).^(8))))-0.0009;
    otherwise
        warning('Unexpected function type')

end



tspan=[0, t1]; % simulate system for 5 years

%%initial conditions vector
N1=5000;
I10=100;
R10=0;
S10=N1-I10-R10;
N2=50000;
I20=1000;
S20=N2-I20;
y0=[S10, I10, R10, S20, I20];



%%solve and plot year-wise: t/365 for better visualisation
[t,y]=ode23s(@(t,y) sir_rhs(t,y, ft,f),tspan,y0);


subplot(3, 1, 1);
plot (t/365, y(:,2)/N1);
ylim([0.28 0.39]);
xlim([0 10]);
xlabel ('time');
ylabel ('I_h');

subplot(3, 1, 2);
plot (t/365,y(:,5)./(y(:,4) + y(:,5)));
ylim([0.16 0.23]);
xlim([0 10]);
xlabel ('time');
ylabel ('I_v');

subplot(3, 1, 3);
plot (ft/365,f);
xlim([0 1]);
xlabel ('time');
ylabel ('mu(t)');

%function to specify system equations
function dydt = sir_rhs(t,y,ft,f)  
f = interp1(ft,f,t);  %%interpolates values of forcing function 'f' at simulation times 't'
dydt=zeros(5,1);

%specify parameters as a vector
beta_hv=0.03;
mu_v=0.2;
beta_vh=0.1;
alpha=0.0011;
mu_h=0.1;
gamma=0.01;
pars=[mu_h, beta_hv, alpha, gamma, beta_vh, mu_v];


% %%system equations 
dydt(1)=pars(1)*(y(1)+y(2)+y(3)) - pars(2)*y(5) * y(1) /(y(1)+y(2)+y(3)) + pars(3) * y(3) - pars(1) * y(1);
dydt(2)= pars(2) * y(5) * y(1) / (y(1)+y(2)+y(3)) - pars(4) * y(2) - pars(1)* y(2);
dydt(3)=pars(4) * y(2) - pars(3) * y(3) -  pars(1)* y(3);
dydt(4)=f.*(y(4)+y(5))- pars(5)* y(2) * y(4) / (y(1)+y(2)+y(3)) - pars(6)* y(4); %'.*' because only one element of 1D array f is multiplied
dydt(5)=pars(5) * y(2) * y(4) / (y(1)+y(2)+y(3)) - pars(6) * y(5); 
end
end
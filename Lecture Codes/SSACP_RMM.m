% H. Reed Ogrosky
% Virginia Commonwealth University
% June 2018

% This script implements SSA-CP on the Real-time Multivariate MJO (RMM) 
% indices.  For equation numbers in comments, refer to Ogrosky et al., 2019

% Reference for SSA-CP:
% Ogrosky, H.R., Stechmann, S.N., Chen, N., and Majda, A.J. (2019)
% "Singular Spectrum Analysis with Conditional Predictions for State
% Estimation and Forecasting". Geophys. Res. Lett. 46. 
% https://doi.org/10.1029/2018GL081100

% Reference for RMM indices
% Wheeler, M.C., and Hendon, H.H. (2004) "An all-season real-time 
% multivariate MJO index: Development of an index for monitoring and
% prediction." Mon. Wea. Rev. 132, 1917-1932.

% The RMM indices can be obtained online at Bureau of Meteorology website
% (http://www.bom.gov.au/climate/mjo/).

clear all
close all

%%%%%% Parameter values %%%%%%
% Embedding window length
M = 51;

% Maximum number of modes to keep?
maxmcount=10;     % Must be less than or equal to M*D, where D is the number of spatial dimensions

% How many modes do you want to sum for your reconstruction plot?
maxmodetoplot=2;  % Must be less than or equal to maxmcount

% Which spatial grid point reconstruction do you want to plot?
dtoplot=1;        % Must be a positive integer less than or equal to D
                  % 1 plots RMM1, 2 plots RMM2

%%%%%% End Parameter values %%%%%%

% Load data here - this script creates two row vectors, RMM1 and RMM2
load_RMM_indices

% Form (D-dimensional) time series x
x=[RMM1; RMM2];

% Remove last 600+2*(M-1) days to match Fig. 1 of Ogrosky et al., 2019
x=x(:,1:end-600-2*(M-1));  

% D=number of spatial dimensions, N=length of time series
[D,N]=size(x);

if maxmcount>M*D error('maxmcount may not exceed M*D'); end
if maxmodetoplot>maxmcount error('maxmodetoplot may not exceed maxmcount'); end
if dtoplot>D error('dtoplot may not exceed D'); end

% Remove mean
for dcount=1:D
    x(dcount,:)=x(dcount,:)-mean(x(dcount,:));
end

t=(-M+2:N+M-1)';

% Create time-delayed embedding of data matrix (eq. 1)
X=zeros(D*M,N-M+1);
for mcount=1:M
    X(1+(mcount-1)*D:mcount*D,:)=x(:,mcount:N-M+mcount);
end

% Calculate covariance matrix
C=X*X'/(N-M+1);

% Find eigenvalues of covariance matrix (eq. 2)
[V,Lambda]=eig(C);
Lambda=diag(Lambda);      
[Lambda,ind]=sort(Lambda,'descend'); 
V=V(:,ind);             

% Calculate principal components (eq. 3)
Phi=X'*V(:,1:maxmcount); 

% Calculate extended principal components (eqs. 7, 8, and related text)
Phitilde=[zeros(M-1,maxmcount); Phi; zeros(M-1,maxmcount)];
for mcount=1:M-1
    y1=X(1+D*mcount:end,end);                       % y1 - known data
    C21=C(end-D*mcount+1:end,1:end-D*mcount);       % (see eq. 8)
    C11=C(1:end-D*mcount,1:end-D*mcount);           % (see eq. 8)
    mu21=C21*(C11\y1);                              % (eq. 7)
    y=[y1; mu21];                                   % column of Xtilde
    Phitilde(end-M+1+mcount,:)=y'*V(:,1:maxmcount); % Phitilde=Xtilde'*V
end

% Calculate the extended reconstructed components (eq. 9)
ztilde=zeros(N,maxmcount,D);
for modecount=1:maxmcount 
    for dcount=1:D
        phivproduct=flipud(Phitilde(:,modecount)*V(dcount:D:end,modecount)'); % invert projection - first channel
        for ncount=1-(M-1):N+(M-1) 
            ztilde(ncount+(M-1),modecount,dcount)=mean(diag(phivproduct,-N+ncount) );
        end
    end
end

% Make SSA-CP portions of figures 1a, 1d from Ogrosky et al., 2019
figfontsize=14;
plotwindowlength=100;  % How far before N should plot start?
yl=2.3;                % Sets limits of y-axis

figure;
set(gcf,'Position',[100,100,400,200]);
hold on;
fill([N-M+1 N N N-M+1],[-yl -yl yl yl]*1.1,[1 0.95 0.95]);  % Creates state-estimation shaded region
fill([N N+M-1 N+M-1 N],[-yl -yl yl yl]*1.1,[1 0.95 1]);     % Creates forecasting shaded region

% Plot reconstructed component. (Red) times prior to N; (magenta) times after N
plot(t,sum(ztilde(:,1:maxmodetoplot,dtoplot),2),'r','Linewidth',1.5);
plot(t(end-M+1:end),sum(ztilde(end-M+1:end,1:maxmodetoplot,dtoplot),2),'m','Linewidth',2);

title(strcat(['RMM',int2str(dtoplot),' (Modes 1-',int2str(maxmodetoplot),')']));
set(gca,'XTick',[N-2*(M-1):(M-1):N+M-1]);
set(gca,'XTickLabel',{'N-2(M-1)','N-M+1','N','N+M-1'});
ylim([-yl yl]);
xlabel('Day');
xlim([N-plotwindowlength N+M-1]);
set(gca,'fontsize',figfontsize);

figureHandle = gcf;
set(gca,'fontsize',figfontsize);
set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
set(gcf,'Units','points');
set(gcf,'PaperUnits','points');
sizepaper=get(gcf,'Position');
sizepaper=sizepaper(3:4);
set(gcf,'PaperSize',sizepaper);
set(gcf,'PaperPosition',[0,0,sizepaper(1),sizepaper(2)]);

% filenamestring=strcat(['RMMindices']);
% saveas(gcf,filenamestring,'epsc');

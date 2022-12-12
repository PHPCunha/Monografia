function modelo = plsmodel2(Xcal,ycal,Xtest,ytest,options)
% Esta funcao tem como principal objetivo facilitar o uso de PLS entre
% leigos em quimiometria, de forma clara e rapida.
%   Detailed explanation goes here

% [model] = plsmodel(Xcal,ycal,Xtest,ytest,options);
% [model] = plsmodel(Xcal,ycal,options);
% [model] = plsmodel(model,Xtest,ytest);
%

%options=[];
%options.Xpretreat      = {'center'}; % Pretratamento da matriz Xcal;;
%options.vene           = 5;      % Tamanho da janela de validacao cruzada;
%options.vl             = 20;     % Numero de variaveis latentes;

%% Preparacao do modelo.
if nargin == 3;
   if isstruct(Xcal);
      error("Nao foi programado, ainda, calculo de Test com modelo.");
   else
      narg = 1;
      %narg = "Calibracao";
      options = Xtest;
   end
elseif nargin == 5;
      narg = 2;
      %narg = "Teste";
else
      error("Adicionado numero incorreto de inputs.");
end

%% Notando erros;

if size(Xcal,1) ~= size(ycal,1);
   error("Amostras de Calibracao com dimensoes incorretas.");
end

%% Desenvolvimento do modelo.
if narg == 1;
   Xcal=pretrat(Xcal,[],options.Xpretreat); % pretreat Xcal
   CV=plscv(Xcal,ycal,options.vl,options.vene,options.Xpretreat);
   plot(1:options.vl,CV.RMSECV,'b'), hold on
   plot(1:options.vl,CV.RMSECV,'ko','MarkerFaceColor','b')
   xlabel('Latent Variable Number','FontSize',14);
   ylabel('RMSECV','FontSize',14);

   modelo.RMSECV = CV.RMSECV;
   modelo.R2cv = CV.Q2;
   modelo.Ycv = [ycal CV.Ypred];
   modelo.options = options;

   return %Finalizando a partir daqui.

%% Desenvolvimento da Previcao criando o modelo do zero.
elseif  narg == 2;
   if size(Xtest,1) ~= size(ytest,1);
   error("Amostras de Teste com dimensoes incorretas.");
   end

   PLS=pls(Xcal,ycal,options.vl,options.Xpretreat);
   CV=plscv(Xcal,ycal,options.vl);

   Xcal2=[Xcal ones(size(Xcal,1),1)];
   ycal2=Xcal2*PLS.regcoef(:,options.vl);

   Xtest2=[Xtest ones(size(Xtest,1),1)];
   ytest2=Xtest2*PLS.regcoef(:,options.vl);

   ycv2 = CV.Ypred(:,end);

   RMSEC=sqrt((sum((ycal-ycal2).^2))./(length(ycal2)-options.vl-1));
   R2c = corrcoef(ycal,ycal2); R2c=R2c(1,2)^2;
   biasc = sum(ycal-ycal2)/length(ycal);

   RMSEP=sqrt((sum((ytest-ytest2).^2))./(length(ytest2)));
   R2p = corrcoef(ytest,ytest2); R2p=R2p(1,2)^2;
   biasp = sum(ytest-ytest2)/length(ytest);

   RMSECV=sqrt((sum((ycal-ycv2).^2))./(length(ycv2)-options.vl-1));
   R2cv = corrcoef(ycal,ycv2); R2cv=R2cv(1,2)^2;
   biascv = sum(ycal-ycv2)/length(ycal);
else
   error("Algo nao esta certo.");
end

%% Finalizando modelo
modelo.data    = date;
modelo.RMSEC   = RMSEC;
modelo.RMSECV  = RMSECV;
modelo.RMSEP   = RMSEP;
modelo.R2c     = R2c;
modelo.R2cv    = R2cv;
modelo.R2p     = R2p;
modelo.Ycal    = [ycal ycal2];
modelo.Ycv     = [ycal ycv2];
modelo.Ytest   = [ytest ytest2];
%modelo.bias.c  = biasc;
%modelo.bias.cv = biascv;
%modelo.bias.p  = biasp;
modelo.options = options;

%% Outras informa��es;
modelo.PLS.loading    = PLS.X_loadings;
modelo.PLS.regcoef    = PLS.regcoef;
modelo.PLS.vip        = PLS.VIP;
modelo.coef.B = PLS.coef.B;
modelo.coef.T = PLS.coef.T;
modelo.coef.P = PLS.coef.P;
modelo.coef.W = PLS.coef.W;

end

%% Sub-fuction

% PLScv %%Enxugar
function CV=plscv(X,y,A,K,method,PROCESS,order)
%+++ K-fold Cross-validation for PLS
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The maximal number of latent variables for cross-validation
%            K: fold. when K=m, it is leave-one-out CV
%       method: pretreatment method. Contains: autoscaling,
%               pareto,minmax,center or none.
%      PROCESS: =1 : print process.
%               =0 : don't print process.
%+++ Order: =0  sorted, default. For CV partition.
%           =1  random.
%           =2  original.
%+++ Output: Structural data: CV
%+++ Hongdong Li, Oct. 16, 2008.
%+++ Tutor: Prof. Yizeng Liang, yizeng_liang@263.net.
%+++ Contact: lhdcsu@gmail.com.

if nargin<7;order=0;end
if nargin<6;PROCESS=1;end
if nargin<5;method='center';end
if nargin<4;K=10;end
if nargin<3;A=3;end

if order==0
  [y,indexyy]=sort(y);
  X=X(indexyy,:);
elseif order==1
  indexyy=randperm(length(y));
  X=X(indexyy,:);
  y=y(indexyy);
elseif order==2
  indexyy=1:length(y);
  X=X(indexyy,:);
  y=y(indexyy);
end

[Mx,Nx]=size(X);
A=min([size(X) A]);
yytest=nan(Mx,1);
YR=nan(Mx,A);

groups = 1+rem(0:Mx-1,K);
for group=1:K

    calk = find(groups~=group);
    testk = find(groups==group);

    Xcal=X(calk,:);
    ycal=y(calk);

    Xtest=X(testk,:);
    ytest=y(testk);

    %   data pretreatment
    [Xs,xpara1,xpara2]=pretreat(Xcal,method);
    [ys,ypara1,ypara2]=pretreat(ycal,'center');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [B,Wstar,T,P,Q]=plsnipals(Xs,ys,A);   % no pretreatment.

    yp=[];
    for j=1:A
        B=Wstar(:,1:j)*Q(1:j);
        %+++ calculate the coefficient linking Xcal and ycal.
        C=ypara2*B./xpara2';
        coef=[C;ypara1-xpara1*C;];
        %+++ predict
        Xteste=[Xtest ones(size(Xtest,1),1)];
        ypred=Xteste*coef;
        yp=[yp ypred];
    end

    YR(testk,:)=[yp];
    yytest(testk,:)=ytest;
    if PROCESS==1; fprintf('The %dth group finished.\n',group); end;
end

%+++ return the original order
YR(indexyy,:)=YR;
y(indexyy)=y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% mean and sd of squared error %%%%%%%%%%%%
  error=YR-repmat(y,1,A);
  error2=error.^2;
  error2_MEAN=sum(error2)/Mx;
  error2_SD= sqrt(sum((error2-repmat(mean(error2),Mx,1)).^2)/(Mx-1)); % unbiased estimator

  %+++ calculate Q2
  cv=sqrt(error2_MEAN);
  [RMSEP,index]=min(cv);index=min(index);
  SST=sumsqr(yytest-mean(y));
  for i=1:A
    SSE=sumsqr(YR(:,i)-y);
    Q2(i)=1-SSE/SST;
  end

  %+++ take standard deviation into account for LV chosen
  indexSD=find(error2_MEAN <= min(error2_MEAN)+error2_SD(index));
  indexSD=min(indexSD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+++ output  %%%%%%%%%%%%%%%%
  CV.method=method;
  CV.Ypred=YR;
  CV.predError=error;
  CV.RMSECV=cv;
  CV.Q2=Q2;
  CV.RMSECV_min=RMSEP;
  CV.Q2_max=Q2(index);
  CV.optLV=index;
  CV.note='*** The following is based on global min MSE + 1SD';
  CV.RMSECV_min_1SD=cv(indexSD);
  CV.Q2_max_1SD=Q2(indexSD);
  CV.optLV_1SD=indexSD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% PLS   %%Enxugar
function PLS=pls(X,y,A,method)
%+++  PLS=pls(x0,y0,A,method);
%+++  Input:
%     X,y: sample data and y-value to predict
%     A: number of PLS components
%     method: pretreat method for X, either "center" or "autoscaling". y is
%     always centered in our libPLS package.
%+++  Ouput : is a structural array which are explained at the end of this code.
%+++  Hongdong Li, June 1,2008.
%+++  Tutor:Yizeng Liang, yizeng_liang@263.net.
%+++  Contact: lhdcsu@gmail.com.

if nargin<4;method='center';end
if nargin<3;A=2;end;

[Mx,Nx]=size(X);

%+++ check effectiveness of A.
A=min([Mx Nx A]);

%+++ data pretreatment
[Xs,xpara1,xpara2]=pretreat(X,method);
[ys,ypara1,ypara2]=pretreat(y,'center');

%+++ Use the pretreated data to build a PLS model
[B,Wstar,T,P,Q,W,~,~]=plsnipals(Xs,ys,A);
% notice that here, B is the regression coefficients linking Xs and ys.

%+++ perform target projection;
%[tpt,tpw,tpp,SR]=tp(Xs,B);

%+++ calculate VIP *********************
VIP=vip(Xs,ys,T,W,Q);

%+++ get regression coefficients that link X and y (original data) ************
coef=zeros(Nx+1,A);
for j=1:A
    Bj=Wstar(:,1:j)*Q(1:j);
    C=ypara2*Bj./xpara2';
    coef(:,j)=[C;ypara1-xpara1*C;];
end

%+++ ********************************************
x_expand=[X ones(Mx,1)];
ypred=x_expand*coef(:,end);
error=ypred-y;
%********************************************
SST=sum((y-mean(y)).^2);
SSE=sum((y-ypred).^2);
R2=1-SSE/SST;

%[fm.RMSEC,fm.R2c,fm.biasc]=rmse(ycal(:,1),PLS.y_fit,nvl);

RMSEC=sqrt((sum((y-ypred).^2))./(length(y)-A-1));
R2c = corrcoef(y,ypred); R2c=R2c(1,2)^2;
biasc = sum(y-ypred)/length(y);

%+++ Output**************************************
PLS.regcoef=coef;
PLS.X_scores=T;
PLS.X_loadings=P;
PLS.VIP=VIP;
PLS.ypred=ypred;
PLS.RMSEC = RMSEC;
PLS.R2c = R2c;
PLS.biasc = biasc;
PLS.coef.B = B;
PLS.coef.T = T;
PLS.coef.P = P;
PLS.coef.W = W;

end
% Paulo R Filgueiras - 14/12/2015
% Edit Pedro H. P. Cunha - 16/09/2022

function VIP=vip(X,y,T,W,Q)

%+++ Calculate the vip for each variable to the response;
%+++ T,W,Q are output from plsnipals.m
%+++ T: socres, which can be obtained by pls_nipals.m
%+++ W: weight, which can be obtained by pls_nipals.m
%+++ Q: Y-loadings
%+++ VIP=sqrt(p*q/s);
%+++      where, p is a constant, denoting the number of variables in x
%                q stands for the explained variance of Y by each variable
%                s represents the total variance explained by A components
%+++ Reference: Tahir Mehmood et al, Chemometrics and Intelligent Laboratory Systems 118 (2012)62?69
%+++ HDLi


s=diag(T'*T*Q*Q');
%initializing
[m,p]=size(X);
[m,h]=size(T);
%+++ calculate VIP;
VIP=[];
for i=1:p
    weight=[];
    for j=1:h
        weight(j,1)= (W(i,j)/norm(W(:,j)))^2;
    end
    q=s'*weight;  % explained variance by variable i
    VIP(i)=sqrt(p*q/sum(s));
end

%+++
end

% NISPALS
function [B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(X,Y,A)
%+++ The NIPALS algorithm for both PLS-1 (a single y) and PLS-2 (multiple Y)
%+++ X: n x p matrix
%+++ Y: n x m matrix
%+++ A: number of latent variables
%+++ Code: Hongdong Li, lhdcsu@gmail.com, Feb, 2014
%+++ reference: Wold, S., M. Sj�str�m, and L. Eriksson, 2001. PLS-regression: a basic tool of chemometrics,
%               Chemometr. Intell. Lab. 58(2001)109-130.

varX=sum(sum(X.^2));
varY=sum(sum(Y.^2));
for i=1:A
    error=1;
    u=Y(:,1);
    niter=0;
    while (error>1e-8 && niter<1000)  % for convergence test
        w=X'*u/(u'*u);
        w=w/norm(w);
        t=X*w;
        q=Y'*t/(t'*t);  % regress Y against t;
        u1=Y*q/(q'*q);
        error=norm(u1-u)/norm(u);
        u=u1;
        niter=niter+1;
    end
    p=X'*t/(t'*t);
    X=X-t*p';
    Y=Y-t*q';

    %+++ store
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;
    Q(:,i)=q;

end

%+++ calculate explained variance
R2X=diag(T'*T*P'*P)/varX;
R2Y=diag(T'*T*Q'*Q)/varY;

Wstar=W*(P'*W)^(-1);
B=Wstar*Q';
Q=Q';

%+++
end

% Funcao do matlab
function [s,n] = sumsqr(x)
%SUMSQR Sum of squared elements of a matrix or matrices.
%
%  [S,N] = sumsqr(X) returns S the sum of all squared finite elements in M,
%  and the number of finite elements N.  M may be a numeric matrix or a
%  cell array of numeric matrices.
%
%  For example:
%
%    s = sumsqr([1 2;3 4])
%    [s,n] = sumsqr({[1 2; NaN 4], [4 5; 2 3]})
%
%  See also MEANSQR, SUMABS, MEANABS.

% Mark Beale, 1-31-92
% Copyright 1992-2010 The MathWorks, Inc.

if nargin < 1,error(message('nnet:Args:NotEnough'));end

if isreal(x)
  notFinite = find(~isfinite(x));
  x(notFinite) = 0;
  s = x.*x;
  s = sum(s(:));
  n = numel(x) - length(notFinite);
elseif iscell(x)
  numx = numel(x);
  s = zeros(1,numx);
  n = zeros(1,numx);
  for i=1:numx
    [s(i),n(i)] = sumsqr(x{i});
  end
  s = sum(s);
  n = sum(n);
else
  error(message('nnet:NNData:NotNumOrCellNum'))
end

end

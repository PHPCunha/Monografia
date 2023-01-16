function modelo = lev_res(modelo,Xcal,ycal,Xtest,ytest)
%Cria grafico de leverage e residuos absolutos do modelo para identificacao
%de outliers
%Exemplo:
% lev_res(modelo,Xcal,ycal,Xtest,ytest);
% modelo = lev_res(modelo,Xcal,ycal,Xtest,ytest);
%
%% Futuras atualizacoes
% Melhor descricao da funcao.
% Explicacao com base na literatura.
%

% Leverage
conf = confidence_interval(Xcal,Xtest,ycal,modelo.options.Xpretreat,modelo.options.vl,0.95,modelo)

lev = conf.leverage_cal;
lev_test = conf.leverage_test;
lev_limite = conf.leverage_limit;

%Calcula os residuos absolutos de calibracao e teste
res_cal   = abs(ycal - modelo.Ycal(:,2));
res_test  = abs(ytest - modelo.Ytest(:,2));
res_cal2  = ycal - modelo.Ycal(:,2);
res_test2 = ytest - modelo.Ytest(:,2);

%Plota os graficos de leverage em funcao de residuo absoluto de calibracao
%e teste
plot(lev,res_cal,'b.','MarkerSize',12), hold on
plot(lev_test,res_test,'r.','MarkerSize',12), hold off

%Nomeia os eixos x e y
xlabel('Leverage','FontSize',16)
ylabel('Residuo absoluto','FontSize',16)
set(gca,'FontSize',18);

%Cria linhas verdes limitanto alto leverage e alto residuo absoluto
vline(lev_limite,'k');
hline(2*modelo.RMSEP,'k');
res_limite = 2*modelo.RMSEP;
legend('Calibracao','Previsao');

% Verificando que amostras podem ser consideradas outlier.
Sample_cal = find( lev_limite < lev & res_limite < res_cal);
Sample_test = find( lev_limite < lev_test & res_limite < res_test);

modelo.lev_res.lev_cal    = lev;
modelo.lev_res.lev_test   = lev_test;
modelo.lev_res.lev_limite = lev_limite;
modelo.lev_res.res_cal    = res_cal;
modelo.lev_res.res_test   = res_test;
modelo.lev_res.res_limite = res_limite;
modelo.lev_res.res_cal2   = res_cal2;
modelo.lev_res.res_test2  = res_test2;
modelo.lev_res.Sample_cal = Sample_cal;
modelo.lev_res.Sample_test = Sample_test;


%% Sub-Rotinas

function conf = confidence_interval(Xcal,Xtest,ycal,Xpretreat,vl,alpha,modelo)
% confidence = confidence_interval(Xcal,Xtest,ycal,Xpretreat,vl,alpha)
%
% Rotina para calcular o intervalo de confianca em modelos de regressao PLS.
% Os calculos sao baseados na norma ASTM E1655-05 (2012).
% NORMA: ASTM E1655-05 (2012) Standard Practices for Infrared Multivariate
%        Quantitative Analysis.
% input
% Xcal  : Matriz X de calibracao.
% Xtest : Matriz X de teste.
% ycal  : vetor y de calibracao.
% Xpretreat : metodo de preprocessamento dos dados X.
%             Xcal e Xtest serao pre-processados conforme construcao do
%             modelo de calibracao. Ex: {'center'}
%  vl  : escalar com o numero de variaveis latentes.
% alpha : nivel de abrangencia do intervalo de confianca. default: 0.95
%
% output
% confidence : estrutura contendo o intervalo de confianca e o leverage das
%              amostras
% saidas:
%   + cal : matriz contendo [valor previsto(PLS) intervalo de confianca] calibracao
%   + test: matriz contendo [valor previsto(PLS) intervalo de confianca]  teste.
%   + leverage_cal  : vetor contendo o leverage das amostras de calibracao.
%   + leverage_test : vetor contendo o leverage das amostras de teste.
%   + leverage_limit : valor limite de leverage segundo a norma ASTM e 1655-05(2012).
%
% Exemplo:
% confidence = confidence_interval(Xcal,Xtest,ycal,{'msc';'center'},3,0.95);
%
% Paulo R. Filgueiras  - 11/09/2014
%

if nargin<6, alpha=0.95; end

% -------------------- C�lculo dos valores previstos --------------------
%[ycal_prev,coef] = submodel_pls(Xcal,ycal,Xcal,Xpretreat,vl,ycal); %previsao das amostras de calibracao.
ycal_prev = modelo.Ycal(:,2);
coef.R = modelo.coef.W; %PLS.W;
coef.P = modelo.coef.P; %PLS.X_loadings;
coef.T = modelo.coef.T; %PLS.X_scores;
ytest_prev = modelo.Ytest(:,2); %= submodel_pls(Xcal,ycal,Xtest,Xpretreat,vl,ones(size(Xtest,1),1)); %previsao das amostras de teste.


% -------------------- Standard Error of Calibration  SEC --------------------
RMSEC = sqrt((sum((ycal-ycal_prev).^2))/(length(ycal)-vl-1));
% norma ASTM E1655-05 (2012) pag.12

% --------------------  leverage  --------------------
[~,Xtest_pretreat]=pretrat(Xcal,Xtest,Xpretreat);

tv=Xtest_pretreat*coef.R(:,1:vl)*inv(coef.P(:,1:vl)'*coef.R(:,1:vl));

h_limit=3*(vl)/size(Xcal,1); % leverage limite, segundo a norma ASTM E1655-05 (2012) pag.15

% Leverage from calibration set
for ki=1:size(Xcal,1)
    h_cal(ki,1)=coef.T(ki,1:vl)*inv(coef.T(:,1:vl)'*coef.T(:,1:vl))*coef.T(ki,1:vl)';
end

% Leverage from test set
for ki=1:size(Xtest,1)
    h_test(ki,1)=tv(ki,1:vl)*inv(coef.T(:,1:vl)'*coef.T(:,1:vl))*tv(ki,1:vl)';
end
% h_limit : leverage limite
%  h_cal  : leverage das amsostras de calibracao.
%  h_test : leverage das amsostras de teste.
% --------------------  end leverage computation  --------------------


% --------------------  Estatastica t-student  --------------------
t_stat=abs(tinv((1-alpha)/2,length(ycal)-1));


% --------------------  Intervalo de confianca  --------------------
% Calculado segundo a norma ASTM E1655-05 (2012) pag.14
conf.cal  =[ycal_prev  t_stat*RMSEC*sqrt(1+h_cal)];
conf.test =[ytest_prev t_stat*RMSEC*sqrt(1+h_test)];
conf.leverage_cal = h_cal;
conf.leverage_test = h_test;
conf.leverage_limit = h_limit;
conf.Ttest = tv;




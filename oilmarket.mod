// A STRUCTURAL MODEL OF THE OIL MARKET
// European Central Bank

// 1. DECLARING VARIABLES

    // 1.1 ENDOGENOUS VARIABLES

    var  Y K C I O L B P R Ch Oh Kh Ih Xh Bh Rh Ct Ot Kt It Xt Bt Rt lambda Gz at ah CA CASA CAP dlogIh dlogCh dloglambda
    dlogIP dlogO Oshare dlogROP dlogOt dlogOh dlogC;

    // 1.2 EXOGENOUS VARIABLES

    varexo eps_Z eps_at eps_ah;

    // 1.3 PARAMETERS
 
    parameters   alpha beta gamma delta eta nu omega Z Zh Zt GzSS GoSS gamma_t psi Bss  Btss Bhss Rss
                 rho_Z rho_at rho_ah stderr_eps_Z stderr_eps_at stderr_eps_ah;  

// 2. PARAMETERS AND STEADY-STATE VALUES

    // 2.1 CALIBRATED PARAMETERS

    alpha   = 0.67;                 % labor share in final goods production
    beta    = 1.01^(-1/12);         % time discount rate
    omega   = 1;                    % labor supply elasticity            
    GzSS    = 1.03^(1/12);          % TFP growth rate oil-importer
    GoSS    = 0.99818;              % productivity growth shortfall in oil production w.r.t TFP growth
    delta   = 1.10^(1/12) - 1;      % capital depreciation in oil production  
    eta     = 4;                    % inverse price elasticity of oil demand: 4 high elasticity VS 21 low elasticty
    gamma_t = 0.4;                  % variable input share in fringe producers production function
    psi     = 0.001;                % elasticity of the interest rates to the change in indebtness   


    // 2.2 DYNAMIC PARAMETERS (don't affect the steady-state)

    
    rho_Z   =  0.001;                  % annual persistence of TFP growth rate shock
    rho_at  =  0;                    % annual persistence of fringe producers productivity shock
    rho_ah  =  0.5;                  % annual persistence of dominant producer preference shock
    stderr_eps_Z  =  0.004;          % standard deviation of the TFP innovation 
    stderr_eps_at =  0.05;           % standard deviation of the innovation to frinde producers productivity    0.05 high elasticty 0.04 low elasticity
    stderr_eps_ah =  0.001;          % standard deviation of preferences of dominant producer

    // 2.3 STEADY-STATE MOMENTS
    
    % SAshare =  0.123;               % Share of dominant producer in production 
    SAshare =  0.4;               % Share of dominant producer in production 
    oilshare = 0.05;                % Share of oil imports / GDP
    Btss     = 0;                   % Debt over GDP of fringe producers
    Bhss     = 0;                   % Debt over GDP of dominant producer
    capacity = 1.5;                 % dominant producer capacity utilization 75%

    // 2.4 STEADY-STATE VALUES
    
    Rss = GzSS / beta;
    Zt  = 1;                        % initial level of fringe productivity (normalization)
    Pss = (Rss-(1-delta))^(1-gamma_t) / (Zt * gamma_t^gamma_t * (1-gamma_t)^(1-gamma_t));
    
    Yss = 1;
    Kss = (1-alpha) * GzSS * Yss / (GzSS/beta - (1-delta));
    Iss = Kss * (1-(1-delta)*GzSS^(-1));
    Oss = oilshare * Yss / Pss;
    
    Otss = (1-SAshare) * Oss;
    Ohss = SAshare * Oss;
    
    Xtss = gamma_t * Pss *Otss;
    Ktss = (1-gamma_t) * GzSS * Pss *Otss / (GzSS/beta - (1-delta));
    Itss = Ktss * (1-(1-delta)*GzSS^(-1));
    Ctss = Pss *Otss - Itss - Xtss - Btss * (1- Rss * GzSS^(-1));
    
    lambdass = -Pss *Ohss / (1/eta * Oss + gamma_t / (1-gamma_t) * Otss);
    gamma = 1 / (1 + (GzSS/beta - (1-delta))*Ktss / capacity / Xtss /GzSS);  % variable input share in dominant producer production function
    
    Xhss = gamma * (Pss + lambdass) *Ohss;
    Khss = (1-gamma) * GzSS * (Pss + lambdass) *Ohss / (GzSS/beta - (1-delta));
    Ihss = Khss * (1-(1-delta)*GzSS^(-1));
    Chss = Pss *Ohss - Ihss - Xhss - Bhss * (1- Rss * GzSS^(-1));
    Zh   = Ohss / Xhss^gamma / (Khss*GzSS^(-1))^(1-gamma);                   % dominant producer productivity advantage
    
    Css = Yss -Iss  - Itss - Xtss - Ctss -  Ihss - Xhss - Chss;
    Lss = (alpha *Yss / Css)^(1/(omega+1));
    Z  = 1/Lss * (Yss / (Kss*GzSS^(-1))^(1-alpha))^(1/alpha);                % initial level of TFP (normalization)
    nu  = Pss * Oss^eta / Css;                                               % efficiency level
    Bss = -Btss -Bhss;   
 

  
// 3. DECLARATION OF THE NONLINEAR DSGE MODEL 

model; 

    // 3.1 CONVERTING PERSISTENCE PARAMETERS FROM ANNUAL TO MONTHLY FREQUENCY

    # rho_Zm =  rho_Z^(1/12);     % monthly TFP growth persistence
    # rho_atm = rho_at^(1/12);    % monthly fringe producers productivity persistence
    # rho_ahm = rho_ah^(1/12);    % monthly dominant producer productivity persistence

    // 3.2 OIL-IMPORTING REGION CONDITIONS

    Y   = (Z* L)^(alpha) * (K(-1)/Gz)^(1-alpha);   %(K(-1)/Gz)^(1-alpha);
    K   = (1-delta)* (K(-1)/Gz) + I;
    C   = P* O^(eta)/nu;
    1   = C* L^(omega+1)/(alpha* Y);
    1   = beta* C/C(+1)* ((1-delta)/Gz(+1)+(1-alpha)* Y(+1)/K);
    Y   = C + I + P * O + B - R(-1)* B(-1)/Gz;
    1   = beta* C/C(+1)/Gz(+1)* R;
 
    // 3.3 FRINGE CONDITIONS (variables with tilde, "t")
   
    Ot  = exp(at)* Zt* Xt^(gamma_t)* (Kt(-1)/Gz)^(1-gamma_t);
    Kt  = (1-delta)* (Kt(-1)/Gz) + It;
    Xt  = gamma_t* P* Ot;
    1   = beta* Ct/Ct(+1)* ((1-delta)/Gz(+1) + (1-gamma_t)* P(+1)* Ot(+1)/Kt);
    1   = beta* Ct/Ct(+1)/Gz(+1)* Rt;
    Ct  = P* Ot -It -Xt -Bt + Rt(-1) * Bt(-1)/Gz;

    // 3.4 DOMINANT PRODUCER CONDITIONS (variables with hat, "h")

    Oh  = Zh* Xh^(gamma)* (Kh(-1)/Gz)^(1-gamma);
	Kh  = (1-delta)* (Kh(-1)/Gz) + Ih;
    Xh  = gamma* Oh * (lambda+P);
    1   = exp(ah(+1)-ah)*beta* Ch/Ch(+1)* ((1-delta)/Gz(+1)+ (1-gamma)* (P(+1)+lambda(+1))* Oh(+1)/Kh);
    1   = exp(ah(+1)-ah)*beta* Ch/Ch(+1)/Gz(+1)* Rh;
    Ch  = P* Oh -Ih -Xh -Bh + Rh(-1) * Bh(-1)/Gz;
    lambda  = -P * Oh / (1/eta *O + gamma_t /(1-gamma_t) * Ot);
 
    // 3.5 MARKET CLEARING CONDITIONS
   
    O   = Oh + Ot;
    B   = 0;
    Bt   = 0;
    Bh   = 0;
   % R   = Rss + psi * (exp(-(B-Bss))-1);
   % Rt  = Rss + psi * (exp(-(Bt-Btss))-1);
   % Rh  = Rss + psi * (exp(-(Bh-Bhss))-1);

    // 3.6 SHOCK PROCESSES

    log(Gz/GzSS) = rho_Zm*log(Gz(-1)/GzSS) + eps_Z;      % TFP growth shock process
    at   = rho_atm* at(-1) + eps_at  ;                   % Fringe producers productivity shock process
    ah   = rho_ahm* ah(-1) + eps_ah  ;                   % Dominant producer productivity shock process

    // 3.7 SOME DATA TRANSFORMATIONS 

    dlogIP    = log(Y/Y(-1)*Gz);
    dlogO     = log(O/O(-1)*Gz*GoSS);
    dlogROP   = log(P/P(-1)/GoSS);
    Oshare    = Oh/O;
    dlogOt    = log(Ot/Ot(-1)*Gz*GoSS);
    dlogOh    = log(Oh/Oh(-1)*Gz*GoSS);
    dlogC     = log(C/C(-1)*Gz);
    CA        = 1 - (C + I)/Y; 
    CASA      = 1 - (Ch + Ih + Xh)/(P* Oh);
    CAP       = (Xh/Kh)^gamma / (Xt/Kt)^gamma_t;
    dlogIh    = log(Ih/Ih(-1)*Gz);
    dlogCh    = log(Ch/Ch(-1)*Gz);
    dloglambda = log(lambda/lambda(-1)/GoSS);

end; 

// 4. INITIAL VALUES FOR STEADY-STATE COMPUTATION

    initval; 

    Y    =  Yss;
    K    =  Kss;
    C    =  Css;
    I    =  Iss;
    O    =  Oss;
    L    =  Lss;
    B    =  Bss;
    P    =  Pss;
    R    =  Rss;
    Ch   =  Chss;
    Oh   =  Ohss;
    Kh   =  Khss;
    Ih   =  Ihss;
    Xh   =  Xhss;
    Bh   =  Bhss;
    Rh   =  Rss;
    Ct   =  Ctss;
    Ot   =  Otss;
    Kt   =  Ktss;
    It   =  Itss;
    Xt   =  Xtss;
    Bt   =  Btss;
    Rt   =  Rss;
    lambda =  lambdass;
    Gz   = GzSS;
    at   = 0;
    ah   = 0;
    dlogIP  = log(GzSS);
    dlogC    = log(GzSS);
    dlogO   = log(GzSS*GoSS);
    Oshare  = Ohss / Oss;
    dlogROP = -log(GoSS);
    dlogOt  = log(GzSS*GoSS);
    dlogOh  = log(GzSS*GoSS);
    CA      = 1 - (Css + Iss)/Yss; 
    CASA    = 1 - (Chss + Ihss + Xhss)/(Pss* Ohss);
    CAP     = (Xhss/Khss)^gamma / (Xtss/Ktss)^gamma_t;
    dlogIh  = log(GzSS);
    dlogCh  = log(GzSS);
    dloglambda = -log(GoSS);

end; 

steady;
%check;

// 5. SHOCKS

shocks;
        
     // 5.1 DETERMINISTIC SHOCKS

     // 5.1.1 SHOCK TO FRINGE PRODUCERS PRODUCTIVITY (SUPPLY SHOCK)

      % var eps_at;
      % periods 1:37   38:96    97:120; 
      % values 0.0025  0.0040 0.0030;
 
     
     // 5.1.2 SHOCK TO DOMINANT SUPPLIER PREFERENCES (SUPPLY SHOCK)

       var eps_ah;
       periods 1:1 2:120;
       values  0  0;   

     // 5.1.3 SHOCK TO TFP
   
      var eps_Z; 
      periods 1:1  2:120; 
      values  -0.03 0 ;
      
      // 5.2 STOCHASTIC SHOCKS

     % var eps_Z    = stderr_eps_Z^2;           % TFP growth innovation
     % var eps_at   = stderr_eps_at^2;          % Fringe producers productivity innovation
     % var eps_ah   = stderr_eps_ah^2;          % Dominant producer productivity innovation

end;

// 6. EXECUTION

simul(periods=120);     % IF DETERMINISTIC SHOCKS
    
     %stoch_simul(order=1, irf=60, periods=120) K  Kh  Kt at ah Gz Y  P  Ot Oh dlogROP dlogIP dlogOt dlogOh Oshare;      % IF STOCHASTIC SHOCKS

time = (2015+1/24):1/12:2030;     
time = time(1:122);
style = 'r';               
line = 1;                  


figure(1)
subplot(3,1,1)
oil_price = 100*(P/Pss-1);
plot(time,oil_price,style,'Linewidth', line)
ylabel('% change')
title('Oil price')
xlim([2015 2018])
hold on

subplot(3,1,2)
%prod = 9*exp(cumsum(dlogOh-log(GzSS*GoSS)))-9;
prod = 100*exp(cumsum(dlogOh-log(GzSS*GoSS)))-100;
plot(time,prod,style,'Linewidth', line)
%ylabel('mbd')
ylabel('%')
title('Saudi Arabia production')
xlim([2015 2018])
hold on

subplot(3,1,3)
prod = 100*exp(cumsum(dlogOt-log(GzSS*GoSS)))-100;
plot(time,prod,style,'Linewidth', line)
ylabel('%')
title('Fringe production')
xlim([2015 2018])
hold on


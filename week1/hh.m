function dudt = hh(t, u, p, IAppFun)
    v = u(1);
    n = u(2);
    m = u(3);
    h = u(4);

    Cm  = p(1);
    gNa = p(2);
    gK  = p(3);
    gL  = p(4);
    ENa = p(5);
    EK  = p(6);
    EL  = p(7);
    phi = p(8);

    INa  = -gNa*m^3*h*(v-ENa);
    IK   = -gK*n^4*(v-EK);
    IL   = -gL*(v-EL);
    IApp = IAppFun(t);

    an =  0.01*(v+55)/(1-exp(-(v+55)/10));
    bn =  0.125*exp(-(v+65)/80);
    am =  0.1*(v+40)/(1-exp(-(v+40)/10));
    bm =  4*exp(-(v+65)/18);
    ah =  0.07*exp(-(v+65)/20);
    bh =  1/(1+exp(-(v+35)/10));
    
    dudt = zeros(4,1);
    dudt(1) = INa + IK + IL + IApp;
    dudt(2) = phi * ( an*(1-n) - bn*n); 
    dudt(3) = phi * ( am*(1-m) - bm*m);
    dudt(4) = phi * ( ah*(1-h) - bh*h);
end

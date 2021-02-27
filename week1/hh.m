function dudt = HodgkinHuxley(p)

    Cm  = p(1);
    gNa = p(2);
    gK  = p(3);
    gL  = p(4);
    ENa = p(5);
    EK  = p(6);
    EL  = p(7);
    phi = p(8);

    dudt = 0;
end

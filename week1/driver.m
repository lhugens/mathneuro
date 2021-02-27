p(1) =     1;           % Cm:   membrane capacitance [microFarads/cm^2]
p(2) =   120;           % gNa:  sodium conductance [milliSiemens/cm^3]
p(3) =    36;           % gK:   potassium conductance [milliSiemens/cm^3]
p(4) =   0.3;           % gL:   leak conductance [milliSiemens/cm^3]
p(5) =    50;           % eNa:  sodium Nernst potential [milliVolts]
p(6) =   -77;           % eK:   potassium Nernst potential [milliVolts] 
p(7) = -54.4;           % eL:   leak reversal potential [milliVolts]
p(8) = 3^((20-6.3)/10);   % phi:  temperature factor, see ET, equation 1.44.


HodgkinHuxley(p);

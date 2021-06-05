function [I_turn V_turn] = turn(Vs, a)
    [r1 r2 i1 i2] = J_eigen(Vs, a);
    idx = find(r2>0);
    V_turn = min(Vs(find(r2>0)));
    I_turn = I_fixed(V_turn, a);
end

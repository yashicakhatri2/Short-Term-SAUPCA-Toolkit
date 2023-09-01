function OUT = substituteFunctions(IN, rfun, TAfun, Eccfun)
    IN = subs(IN, 'r', rfun);
    IN = subs(IN, 'TA', TAfun);
    IN = subs(IN, 'Ecc', Eccfun);
    OUT = IN;
end
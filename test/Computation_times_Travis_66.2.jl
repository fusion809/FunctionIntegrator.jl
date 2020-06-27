# Taken from https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/354648732
L_chebyshev1=[1.228102; 0.380927; 0.343716; 0.088164; 0.036754; 0.088224; 0.038897; 0.041772; 0.101443; 0.090530; 0.130748; 0.092721; 0.052740];
L_chebyshev2=[0.097074; 8.306089; 0.093900; 0.091194; 0.036205; 0.083073; 0.038763; 0.042806; 0.113077; 0.098064; 0.057584; 0.084314; 0.053570];
L_chebyshev3 = [0.167595; 3.375615; 0.087435; 0.086833; 0.042152; 0.083704; 0.044630; 0.047584; 0.111407; 0.086341; 0.062481; 0.083567; 0.064393];
L_chebyshev4 = [0.102224; 4.064065; 0.086987; 0.085582; 0.041838; 0.086670; 0.043695; 0.046258; 0.104307; 0.086315; 0.061660; 0.083707; 0.060925];
L_legendre=[0.099377; 4.215644; 0.102441; 0.062042; 0.025503; 0.064393; 0.030050; 0.029080; 0.064633; 0.063727; 0.038881; 0.055240; 0.038958];
L_lobatto=[0.179238; 19.749167; 0.088166; 0.068954; 0.000061; 0.071718; 0.000102; 0.000317; 0.072961; 0.105261; 0.000111; 0.021398; 0.000097];
L_radau = [0.136686; 20.713057; 0.087354; 0.073406; 0.000068; 0.071638; 0.000109; 0.000463; 0.073130; 0.097385; 0.000084; 0.019310; 0.000060];
L_simpsons=[0.046564; 13.614002; 0.042339; 0.024079; 0.000013; 0.029465; 0.000060; 0.000057; 0.027102; 0.036607; 0.019030; 0.045476; 0.017362];
L_trapezoidal=[0.037420; 13.258997; 0.035951; 0.015603; 0.000708; 0.020648; 0.001492; 0.001060; 0.108480; 0.024120; 0.011874; 0.086951; 0.010772];

N = length(L_simpsons);

# RMS of times
rms_chebyshev1 = sqrt(L_chebyshev1'*L_chebyshev1/N);
rms_chebyshev1_wo_simppen = sqrt((L_chebyshev1[1]^2+L_chebyshev1[3:N]'*L_chebyshev1[3:N])/(N-1));
rms_chebyshev2 = sqrt(L_chebyshev2'*L_chebyshev2/N);
rms_chebyshev2_wo_simppen = sqrt((L_chebyshev2[1]^2+L_chebyshev2[3:N]'*L_chebyshev2[3:N])/(N-1));
rms_chebyshev3 = sqrt(L_chebyshev3'*L_chebyshev3/N);
rms_chebyshev3_wo_simppen = sqrt((L_chebyshev3[1]^2+L_chebyshev3[3:N]'*L_chebyshev3[3:N])/(N-1));
rms_chebyshev4 = sqrt(L_chebyshev4'*L_chebyshev4/N);
rms_chebyshev4_wo_simppen = sqrt((L_chebyshev4[1]^2+L_chebyshev4[3:N]'*L_chebyshev4[3:N])/(N-1));
rms_legendre = sqrt(L_legendre'*L_legendre/N);
rms_legendre_wo_simppen = sqrt((L_legendre[1]^2+L_legendre[3:N]'*L_legendre[3:N])/(N-1));
rms_lobatto_wo_simppen = sqrt((L_lobatto[1]^2+L_lobatto[3:N]'*L_lobatto[3:N])/(N-1));
rms_radau_wo_simppen = sqrt((L_radau[1]^2+L_radau[3:N]'*L_radau[3:N])/(N-1));
rms_simpsons_wo_simppen = sqrt((L_simpsons[1]^2+L_simpsons[3:N]'*L_simpsons[3:N])/(N-1));
rms_trapezoidal_wo_simppen = sqrt((L_trapezoidal[1]^2+L_trapezoidal[3:N]'*L_trapezoidal[3:N])/(N-1));

rms_wo_simppen = [rms_chebyshev1_wo_simppen; rms_chebyshev2_wo_simppen; rms_chebyshev3_wo_simppen; rms_chebyshev4_wo_simppen; rms_legendre_wo_simppen; rms_lobatto_wo_simppen; rms_radau_wo_simppen; rms_simpsons_wo_simppen; rms_trapezoidal_wo_simppen];

println("rms_chebyshev1 is $(rms_chebyshev1)")
println("rms_chebyshev2 is $(rms_chebyshev2)")
println("rms_chebyshev3 is $(rms_chebyshev3)")
println("rms_chebyshev4 is $(rms_chebyshev4)")
println("rms_legendre   is $(rms_legendre)")

println("rms_chebyshev1_wo_simppen is  $(rms_chebyshev1_wo_simppen)")
println("rms_chebyshev2_wo_simppen is  $(rms_chebyshev2_wo_simppen)")
println("rms_chebyshev3_wo_simppen is  $(rms_chebyshev3_wo_simppen)")
println("rms_chebyshev4_wo_simppen is  $(rms_chebyshev4_wo_simppen)")
println("rms_legendre_wo_simppen is    $(rms_legendre_wo_simppen)")
println("rms_lobatto_wo_simppen is     $(rms_lobatto_wo_simppen)")
println("rms_radau_wo_simppen is       $(rms_radau_wo_simppen)")
println("rms_simpsons_wo_simppen is    $(rms_simpsons_wo_simppen)")
println("rms_trapezoidal_wo_simppen is $(rms_trapezoidal_wo_simppen)")

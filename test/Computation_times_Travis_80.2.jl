# Taken from https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/355708079
# 10th entry is simppen
L_chebyshev1=[1.435044; 0.307644; 0.045131; 0.104719; 0.116603; 0.090129; 0.088141; 0.114219; 0.162800; 0.441775; 0.102337; 0.050055; 0.048896];
L_chebyshev2=[0.197077; 0.068784; 0.043831; 0.099813; 0.111966; 0.125740; 0.103428; 0.098729; 0.064575; 8.991385; 0.096214; 0.052295; 0.049748];
L_chebyshev3 = [0.155732; 0.072039; 0.049798; 0.111535; 0.099392; 0.101607; 0.111045; 0.111373; 0.070177; 3.208550; 0.092763; 0.057960; 0.055562];
L_chebyshev4 = [0.146024; 0.070966; 0.049619; 0.092735; 0.098091; 0.100671; 0.099425; 0.099531; 0.070195; 4.146543; 0.094798; 0.056763; 0.053876];
L_jacobi = [21.200547; 0.741378; 0.079454; 0.468545; 0.220772; 2.786482; 0.531896; 0.446746; 0.479568; 22.058539; 0.242638; 1.933591; 0.694156];
L_legendre=[0.105352; 0.076702; 0.030384; 0.071562; 0.078758; 0.078509; 0.088236; 0.090079; 0.043331; 4.383104; 0.059413; 0.036632; 0.033564];
L_lobatto=[0.192820; 0.012460; 0.000103; 0.080504; 0.143451; 0.117507; 0.086438; 0.091315; 0.000142; 21.927058; 0.025977; 0.000409; 0.000141];
L_radau = [0.138636; 0.012143; 0.000040; 0.094062; 0.129707; 0.104701; 0.086289; 0.087083; 0.000057; 22.282400; 0.025089; 0.000345; 0.000109];
L_rectangle_midpoint = [0.073889; 0.000017; 0.000330; 0.019159; 0.022585; 0.030032; 0.018576; 0.021532; 0.000027; 1.863286; 0.035673; 0.000545; 0.000903];
L_simpsons=[0.060178; 0.040657; 0.000012; 0.023218; 0.041971; 0.028887; 0.025965; 0.031452; 0.018607; 13.831230; 0.047495; 0.000064; 0.000078];
L_trapezoidal=[0.161483; 0.034553; 0.000577; 0.017627; 0.032416; 0.022288; 0.019312; 0.023317; 0.013583; 13.361007; 0.077226; 0.001042; 0.001420];

N = length(L_simpsons);

# RMS of times
rms_chebyshev1 = sqrt(L_chebyshev1'*L_chebyshev1/N);
rms_chebyshev1_wo_simppen = sqrt((L_chebyshev1[1:9]'*L_chebyshev1[1:9]+L_chebyshev1[11:N]'*L_chebyshev1[11:N])/(N-1));
rms_chebyshev2 = sqrt(L_chebyshev2'*L_chebyshev2/N);
rms_chebyshev2_wo_simppen = sqrt((L_chebyshev2[1:9]'*L_chebyshev2[1:9]+L_chebyshev2[11:N]'*L_chebyshev2[11:N])/(N-1));
rms_chebyshev3 = sqrt(L_chebyshev3'*L_chebyshev3/N);
rms_chebyshev3_wo_simppen = sqrt((L_chebyshev3[1:9]'*L_chebyshev3[1:9]+L_chebyshev3[11:N]'*L_chebyshev3[11:N])/(N-1));
rms_chebyshev4 = sqrt(L_chebyshev4'*L_chebyshev4/N);
rms_chebyshev4_wo_simppen = sqrt((L_chebyshev4[1:9]'*L_chebyshev4[1:9]+L_chebyshev4[11:N]'*L_chebyshev4[11:N])/(N-1));
rms_jacobi = sqrt(L_jacobi'*L_jacobi/N);
rms_jacobi_wo_simppen = sqrt((L_jacobi[1:9]'*L_jacobi[1:9]+L_jacobi[11:N]'*L_jacobi[11:N])/(N-1));
rms_legendre = sqrt(L_legendre'*L_legendre/N);
rms_legendre_wo_simppen = sqrt((L_legendre[1:9]'*L_legendre[1:9]+L_legendre[11:N]'*L_legendre[11:N])/(N-1));
rms_lobatto_wo_simppen = sqrt((L_lobatto[1:9]'*L_lobatto[1:9]+L_lobatto[11:N]'*L_lobatto[11:N])/(N-1));
rms_radau_wo_simppen = sqrt((L_radau[1:9]'*L_radau[1:9]+L_radau[11:N]'*L_radau[11:N])/(N-1));
rms_rectangle_midpoint = sqrt(L_rectangle_midpoint'*L_rectangle_midpoint/N);
rms_rectangle_midpoint_wo_simppen = sqrt((L_rectangle_midpoint[1:9]'*L_rectangle_midpoint[1:9]+L_rectangle_midpoint[11:N]'*L_rectangle_midpoint[11:N])/(N-1));
rms_simpsons_wo_simppen = sqrt((L_simpsons[1:9]'*L_simpsons[1:9]+L_simpsons[11:N]'*L_simpsons[11:N])/(N-1));
rms_trapezoidal_wo_simppen = sqrt((L_trapezoidal[1:9]'*L_trapezoidal[1:9]+L_trapezoidal[11:N]'*L_trapezoidal[11:N])/(N-1));

rms_wo_simppen = [rms_chebyshev1_wo_simppen; rms_chebyshev2_wo_simppen; rms_chebyshev3_wo_simppen; rms_chebyshev4_wo_simppen; rms_jacobi_wo_simppen; rms_legendre_wo_simppen; rms_lobatto_wo_simppen; rms_radau_wo_simppen; rms_rectangle_midpoint_wo_simppen; rms_simpsons_wo_simppen; rms_trapezoidal_wo_simppen];

println("rms_chebyshev1 is $(rms_chebyshev1)")
println("rms_chebyshev2 is $(rms_chebyshev2)")
println("rms_chebyshev3 is $(rms_chebyshev3)")
println("rms_chebyshev4 is $(rms_chebyshev4)")
println("rms_jacobi     is $(rms_jacobi)")
println("rms_legendre   is $(rms_legendre)")

println("rms_chebyshev1_wo_simppen is  $(rms_chebyshev1_wo_simppen)")
println("rms_chebyshev2_wo_simppen is  $(rms_chebyshev2_wo_simppen)")
println("rms_chebyshev3_wo_simppen is  $(rms_chebyshev3_wo_simppen)")
println("rms_chebyshev4_wo_simppen is  $(rms_chebyshev4_wo_simppen)")
println("rms_jacobi_wo_simppen is      $(rms_jacobi_wo_simppen)")
println("rms_legendre_wo_simppen is    $(rms_legendre_wo_simppen)")
println("rms_lobatto_wo_simppen is     $(rms_lobatto_wo_simppen)")
println("rms_radau_wo_simppen is       $(rms_radau_wo_simppen)")
println("rms_rectangle_midpoint_wo_simppen is $(rms_rectangle_midpoint_wo_simppen)")
println("rms_simpsons_wo_simppen is    $(rms_simpsons_wo_simppen)")
println("rms_trapezoidal_wo_simppen is $(rms_trapezoidal_wo_simppen)")

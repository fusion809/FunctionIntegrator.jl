# Taken from https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/359262615
# 10th entry is simppen
L_adaptive_simpsons = [0.144420; 0.048034; 0.000163; 0.019791; 0.031644; 0.022732; 0.022157; 0.023056; 0.000081; 12.635982; 0.002164; 0.000840; 0.000220];
L_chebyshev1 = [1.374153; 0.295273; 0.042028; 0.097700; 0.085278; 0.091216; 0.085430; 0.086918; 0.136373; 0.405410; 0.088461; 0.048735; 0.045503];
L_chebyshev2 = [0.189887; 0.063260; 0.041871; 0.095741; 0.109335; 0.101715; 0.111041; 0.100269; 0.082595; 8.732129; 0.090618; 0.050086; 0.045789];
L_chebyshev3 = [0.166648; 0.067508; 0.048010; 0.090679; 0.093410; 0.123339; 0.094715; 0.107830; 0.068438; 3.184987; 0.088439; 0.056318; 0.051559];
L_chebyshev4 = [0.133982; 0.067320; 0.048107; 0.093402; 0.105892; 0.101815; 0.095580; 0.095526; 0.067627; 4.097529; 0.090069; 0.054310; 0.050882];
L_jacobi = [20.730626; 0.696972; 0.074401; 0.472530; 0.214534; 2.607443; 0.521480; 0.440099; 0.493098; 21.869202; 0.162223; 1.875650; 0.579464];
L_legendre = [0.104126; 0.075988; 0.029669; 0.069286; 0.074264; 0.075576; 0.072840; 0.071971; 0.043136; 4.255756; 0.058376; 0.034468; 0.030826];
L_lobatto = [0.191745; 0.012332; 0.000122; 0.090188; 0.139901; 0.105984; 0.083651; 0.081196; 0.000138; 21.580829; 0.092406; 0.000301; 0.000141];
L_radau = [0.140783; 0.012576; 0.000053; 0.079066; 0.108806; 0.111685; 0.096150; 0.081501; 0.000065; 21.579555; 0.023179; 0.000421; 0.000113];
L_rectangle_midpoint = [0.074701; 0.000019; 0.000266; 0.018933; 0.020185; 0.018622; 0.018862; 0.018857; 0.000027; 1.564159; 0.029700; 0.000414; 0.000740];
L_simpsons = [0.057260; 0.018206; 0.000013; 0.020936; 0.029926; 0.022812; 0.023652; 0.024098; 0.018203; 13.383746; 0.036649; 0.000066; 0.000052];
L_simpsons38 = [0.056655; 0.040500; 0.000012; 0.023536; 0.031635; 0.026131; 0.024863; 0.025426; 0.017806; 14.101716; 0.040174; 0.000060; 0.000041];
L_trapezoidal = [0.141598; 0.034595; 0.000494; 0.015096; 0.021626; 0.017495; 0.016807; 0.017933; 0.012839; 12.977622; 0.069035; 0.000937; 0.001293];
N = length(L_simpsons);

# RMS of times
rms_adaptive_simpsons_wo_simppen = sqrt((L_adaptive_simpsons[1:9]'*L_adaptive_simpsons[1:9]+L_adaptive_simpsons[11:N]'*L_adaptive_simpsons[11:N])/(N-1));
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
rms_simpsons38_wo_simppen = sqrt((L_simpsons38[1:9]'*L_simpsons38[1:9]+L_simpsons38[11:N]'*L_simpsons38[11:N])/(N-1));
rms_trapezoidal_wo_simppen = sqrt((L_trapezoidal[1:9]'*L_trapezoidal[1:9]+L_trapezoidal[11:N]'*L_trapezoidal[11:N])/(N-1));

rms_wo_simppen = [rms_chebyshev1_wo_simppen; rms_chebyshev2_wo_simppen; rms_chebyshev3_wo_simppen; rms_chebyshev4_wo_simppen; rms_jacobi_wo_simppen; rms_legendre_wo_simppen; rms_lobatto_wo_simppen; rms_radau_wo_simppen; rms_rectangle_midpoint_wo_simppen; rms_simpsons_wo_simppen; rms_trapezoidal_wo_simppen];

println("rms_chebyshev1 is $(rms_chebyshev1)")
println("rms_chebyshev2 is $(rms_chebyshev2)")
println("rms_chebyshev3 is $(rms_chebyshev3)")
println("rms_chebyshev4 is $(rms_chebyshev4)")
println("rms_jacobi     is $(rms_jacobi)")
println("rms_legendre   is $(rms_legendre)")

println("rms_adaptive_simpsons_wo_simppen is    $(rms_adaptive_simpsons_wo_simppen)")
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
println("rms_simpsons38_wo_simppen is    $(rms_simpsons38_wo_simppen)")
println("rms_trapezoidal_wo_simppen is $(rms_trapezoidal_wo_simppen)")

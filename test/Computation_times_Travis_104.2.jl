# Taken from https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/364059665
# 10th entry is simppen
L_adaptive_simpsons = [; ];
L_chebyshev1 = [; ];
L_chebyshev2 = [; ];
L_chebyshev3 = [; ];
L_chebyshev4 = [; ];
L_jacobi = [; ];
L_legendre = [; ];
L_lobatto = [; ];
L_radau = [; ];
L_rectangle_midpoint = [; ];
L_simpsons = [; ];
L_simpsons38 = [; ];
L_trapezoidal = [; ];
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

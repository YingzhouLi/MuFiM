A = getHfd2D(5,5);

AMF = Multifrontal(A);

Aapprox = AMF*eye(size(A));
disp('MF*I')
norm(A-Aapprox)/norm(full(A))

Aapprox = eye(size(A))*AMF;
disp('I*MF')
norm(A-Aapprox)/norm(full(A))

Ainvapprox = eye(size(A))/AMF;
disp('I/MF')
norm(eye(size(A))-Ainvapprox*A)

Ainvapprox = AMF\eye(size(A));
disp('MF\I')
norm(eye(size(A))-Ainvapprox*A)
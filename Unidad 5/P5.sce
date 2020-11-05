// Ejercicio 1


// Resuelve el sistema lineal Ax=b con el método de jacobi
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
function x = jacobi(A, b, x, eps, maxIter)
    D = diag(diag(A))
    Dinv = inv(D)
    LU = A - D
    for i = 1:maxIter
        xNuevo = Dinv * (b - LU * x)
        if norm(xNuevo - x) < eps
            x = xNuevo
            // mostramos las iteraciones
            disp (i)
            break
        end
        x = xNuevo
    end
endfunction


// Resuelve el sistema lineal Ax=b con el método de gauss_seidel
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
function x = gauss_seidel(A, b, x, eps, maxIter)
    n = size(A, 1)
    xNuevo = x
    for contador = 1:maxIter
        xNuevo(1) = 1/A(1, 1) * (b(1) - A(1, 2:n) * x(2:n))
        for i = 2:n-1
            xNuevo(i) = 1/A(i, i) * (b(i) - A(i, 1:i-1) * xNuevo(1 : i-1) - A(i, i+1:n) * x(i+1:n))
        end
        xNuevo(n) = 1/A(n, n) * (b(n) - A(n, 1:n-1) * xNuevo(1:n-1))
        if norm(xNuevo - x) < eps
            x = xNuevo
            // mostramos las iteraciones
            disp (contador)
            break
        end
        x = xNuevo
    end
endfunction


// --> A = [0 2 4; 1 -1 -1; 1 -1 2];

// --> b = [0 .375 0]';


// a) TODO
// b) TODO


// c)

// --> A = [1 -1 0; -1 2 -1; 0 -1 1.1];

// --> b = [0 1 0]'; 

// --> jacobi (A, b, [0 0 0]', .01, 1000)

//    171.
//  ans  =

//    10.789086
//    10.798673
//    9.8082602

// --> gauss_seidel(A, b, [0 0 0]', .01, 1000)

//    97.
//  ans  =

//    10.873565
//    10.879312
//    9.8902836




// Ejercicio 2

// --> A = [10 1 2 3 4; 1 9 -1 2 -3; 2 -1 7 3 -5; 3 2 3 12 -1; 4 -3 -5 -1 15];

// --> b = [12 -27 14 -17 12]';

// --> jacobi(A, b, zeros(5, 1), 1e-6, 1000)

//    67.                                   <-N° de iteraciones
//  ans  =

//    1.0000016
//   -2.0000015
//    2.9999973
//   -1.9999996
//    0.9999981



// --> gauss_seidel(A, b, zeros(5, 1), 1e-6, 1000)

//    38.                                   <-N° de iteraciones
//  ans  =

//    1.0000009
//   -2.0000007
//    2.9999987
//   -1.9999999
//    0.9999992



// Ejercicio 3 TODO



// Ejercicio 4

function [x,a] = gausselimPP(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana con pivoteo parcial.
    
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A b]; // Matriz aumentada
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(a(k,k));  //pivoteo
        for i=k+1:n
            if abs(a(i,k))>amax then
                kpivot = i; amax = a(k,i);
            end;
        end;
        temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
        
        for i=k+1:n
            for j=k+1:n+1
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            end;
            for j=1:k        // no hace falta para calcular la solución x
                a(i,j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end;
    end;
    
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k);
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction


// --> N = 500;
// --> A = 8 * eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1)+ diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
// --> b = ones(N,1);


// --> tic(); gausselimPP(A,b); t = toc()
//  t  = 
//    189.51243



// --> tic(); gauss_seidel(A, b, zeros(N, 1), 1e-6, 10000); t = toc();
//    18.
// --> t
//  t  = 
//    0.179603



// --> tic(); gauss_seidel(A, b, zeros(N, 1), 1e-12, 10000); t = toc();
//    45.
// --> t
//  t  = 
//    0.4146840




// --> N=1000;
// --> A = 8 * eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1)+ diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
// --> b = ones(N,1);

// tic(); gausselimPP(A,b); t = toc()
//  t  = 
//    1561.3197

// -->  tic(); gauss_seidel(A, b, zeros(N, 1), 1e-6, 10000); t = toc();
//    18.
// --> t
//  t  = 
//    0.4799447


// -->  tic(); gauss_seidel(A, b, zeros(N, 1), 1e-12, 10000); t = toc();
//     45.
// --> t
//  t  = 
//    1.1827675


// Ejercicio 5

// Resuelve el sistema lineal Ax=b con el método de sobrerrelajación (SOR)
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
// con w coeficiente de relajación
function x = SOR(A, b, w, x, eps, maxIter)
    n = size(A, 1)
    xNuevo = x
    for contador = 1:maxIter
        xNuevo(1) = w/A(1, 1) * (b(1) - A(1, 2:n) * x(2:n)) + (1 - w) * x(1)
        for i = 2:n-1
            xNuevo(i) = w/A(i, i) * (b(i) - A(i, 1:i-1) * xNuevo(1 : i-1) - A(i, i+1:n) * x(i+1:n)) + (1 - w) * x(i)
        end
        xNuevo(n) = w/A(n, n) * (b(n) - A(n, 1:n-1) * xNuevo(1:n-1)) + (1 - w) * x(n)
        if norm(xNuevo - x) < eps
            x = xNuevo
            // mostramos las iteraciones
            disp (contador)
            break
        end
        x = xNuevo
    end
endfunction


// --> A = [4 3 0; 3 4 -1; 0 -1 4];

// --> b = [24 30 -24]';


// a)

// --> gauss_seidel(A, b, [0 0 0]', 1e-7, 1000)
//    36.
//  ans  =
//    3.0000001
//    3.9999999
//   -5.       


// b)

// TODO CALCULAR W Y UTILIZARLO PARA RESOLVER Ax=b

// --> Dinv = inv(diag(diag(A)));
// --> TJ = (eye(A) - Dinv * A);
// --> rho = (max(abs(spec( TJ ))));
// --> w = 2 / (1 + sqrt(1 - rho ^ 2 ));

// --> SOR(A, b, w, [0 0 0]', 1e-7, 1000)
//    16.
//  ans  =
//    3.
//    4.
//   -5.

// Notemos que tardó 16 iteraciones en vez de 36 con el normal
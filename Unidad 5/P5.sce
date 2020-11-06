// Resuelve el sistema lineal Ax=b con el método iterativo de forma vectorial con matriz N
// tal que Nx^(k+1) = (N-A)x^(k) + b
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
// Además muestra la cantidad de iteraciones realizadas si se llegó a la condición de parada
function x = metodoIterativoVectorial(N, A, b, x, eps, maxIter)
    Ninv = inv(N)
    NmA = N - A
    for i = 1:maxIter
        xNuevo = Ninv *(NmA * x + b)
        if norm(xNuevo - x) < eps 
            x = xNuevo
            // mostramos las iteraciones
            disp (i)
            return
        end
        x = xNuevo
    end
    x = %nan
endfunction




// Ejercicio 1

// Resuelve el sistema lineal Ax=b con el método de jacobi
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
// Además muestra la cantidad de iteraciones realizadas si se llegó a la condición de parada
function x = jacobi(A, b, x, eps, maxIter)
    N = diag(diag(A))
    x = metodoIterativoVectorial(N, A, b, x, eps, maxIter)
endfunction


// Resuelve el sistema lineal Ax=b con el método de gauss_seidel
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
// Además muestra la cantidad de iteraciones realizadas si se llegó a la condición de parada
function x = gauss_seidel(A, b, x, eps, maxIter)
    N = tril(A)
    x = metodoIterativoVectorial(N, A, b, x, eps, maxIter)
endfunction



// a) TODO
// b) TODO


// c)


// --> A = [0 2 4; 1 -1 -1; 1 -1 2];
// --> b = [0 .375 0]';



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



// Ejercicio 3 
// TODO



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
//    0.0086329



// --> tic(); gauss_seidel(A, b, zeros(N, 1), 1e-12, 10000); t = toc();
//    45.
// --> t
//  t  = 
//    0.0101253



// --> N=1000;
// --> A = 8 * eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1)+ diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
// --> b = ones(N,1);

// tic(); gausselimPP(A,b); t = toc()
//  t  = 
//    1561.3197

// --> tic(); gauss_seidel(A, b, zeros(N, 1), 1e-6, 10000); t = toc();
//    18.
// --> t
//  t  = 
//    0.050344

// --> tic(); gauss_seidel(A, b, zeros(N, 1), 1e-12, 10000); t = toc();
//    45.
// --> t
//  t  = 
//    0.0622858



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


function x = SOR(A, b, w, x, eps, maxIter)
    N = tril(A)
    Ninv = inv(N)
    NmA = N - A
    for i = 1:maxIter
        xNuevo = w * Ninv *(NmA * x + b) + (1 - w) * x
        if norm(xNuevo - x) < eps 
            x = xNuevo
            // mostramos las iteraciones
            disp (i)
            return
        end
        x = xNuevo
    end
    x = %nan
endfunction



// --> A = [4 3 0; 3 4 -1; 0 -1 4];
// --> b = [24 30 -24]';


// a)

// --> gauss_seidel(A, b, [0 0 0]', 1e-8, 1000)
//    41.
//  ans  =
//    3.
//    4.
//   -5.       


// b)

// --> Dinv = inv(diag(diag(A)));
// --> TJ = (eye(A) - Dinv * A);
// --> rho = (max(abs(spec( TJ ))));
// --> w = 2 / (1 + sqrt(1 - rho ^ 2 ));

// --> SOR(A, b, w, [0 0 0]', 1e-8, 1000)
//    31.
//  ans  =
//    3.
//    4.
//   -5.

// Notemos que tardó 31 iteraciones en vez de 41 con el normal
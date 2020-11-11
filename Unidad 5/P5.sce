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
    n = size(A, 1)
    for iter = 1:maxIter // Iteramos
        for i = 1:n      // Aplicamos el algoritmo por filas
            xNuevo(i) = 1/A(i, i) * (b(i) - A(i, :) * x + A(i, i) * x(i))
        end
        if norm(xNuevo - x) < eps // Criterio de corte
            x = xNuevo
            disp (iter) // mostramos las iteraciones
            return
        end
        x = xNuevo  // Reasignamos
    end
    x = %nan
endfunction


// Resuelve el sistema lineal Ax=b con el método de gauss_seidel
// comenzando desde x con una tolerancia de eps
// con un máximo de iteraciones maxIter
// Además muestra la cantidad de iteraciones realizadas si se llegó a la condición de parada
function x = gauss_seidel(A, b, x, eps, maxIter)
    n = size(A, 1)
    xNuevo = x
    for iter = 1:maxIter // Iteramos
        // Aplicamos el algoritmo por filas
        xNuevo(1) = 1/A(1, 1) * (b(1) - A(1, 2:n) * x(2:n)) // Primera fila
        for i = 2:n-1    // Filas 2 a n
            xNuevo(i) = 1/A(i, i) * (b(i) - A(i, 1:i-1) * xNuevo(1:i-1) - A(i, i+1:n) * x(i+1:n))
        end
        xNuevo(n) = 1/A(n, n) * (b(n) - A(n, 1:n-1) * xNuevo(1:n-1)) // Última fila

        if norm(xNuevo - x) < eps // Criterio de corte
            x = xNuevo
            disp (iter) // mostramos las iteraciones
            return
        end
        x = xNuevo  // Reasignamos
    end
    x = %nan
endfunction



// a) y b)

// --> A = [0 2 4; 1 -1 -1; 1 -1 2];
// --> b = [0 .375 0]';


// Dado que la diagonal de A no es invertible
// y la triangular inferior asociada a A tampoco, 
// no podemos aplicar los algoritmos
// de jacobi y gauss_seidel
// Entonces redefinimos el sistema como


// --> A = [1 -1 -1; 0 2 4; 1 -1 2];
// --> b = [.375 0 0]';

// A no es diagonal dominante por lo cual no podemos aplicar T3 ni T4

// Jacobi:
// Utilizando C1:
// --> N = diag(diag(A));
// --> m = eye(A) - inv(N) * A;
// --> rho = (max(abs(spec( m ))))
//  rho  =
//    1.3440402
// Que es mayor que 1, por lo tanto no podemos asegurar su convergencia

// Gauss-Seidel:
// Utilizando C1:
// --> N = tril(A);
// --> m = eye(A) - inv(N) * A;
// --> rho = (max(abs(spec( m ))))
//  rho  =
//    2.
// Que es mayor que 1, por lo tanto no podemos asegurar su convergencia




// --> A = [1 -1 0; -1 2 -1; 0 -1 1.1];
// --> b = [0 1 0]'; 

// A no es diagonal dominante por lo cual no podemos aplicar T3 ni T4

// Jacobi: 

// Utilizando el Teorema 1:
// --> N = diag(diag(A));
// --> m = eye(A) - inv(N) * A;
// --> norm(m)
//  ans  =
//    1.3514608
// Que es mayor que 1, por lo tanto no podemos asegurar su convergencia con T1


// Utilizando el C1:
// --> rho = (max(abs(spec( m ))))
//  rho  =
//    0.9770084
// Lo cual es menor que 1 => converge para cualquier vector inicial 


// Gauss-Seidel:
// Utilizando el Teorema 1:
// --> N = tril(A);
// --> m = eye(A) - inv(N) * A;
// --> norm(m)
//  ans  =
//    1.2781759
// Que es mayor que 1, por lo tanto no podemos asegurar su convergencia con T1

// Utilizando el C1:
// --> rho = (max(abs(spec( m ))))
//  rho  =
//    0.9545455
// Lo cual es menor que 1 => converge para cualquier vector inicial 


// c)

// Utilizo el sistema con las filas alternadas:
// --> A = [1 -1 1; 0 2 4; 1 -1 2];
// --> b = [.375 0 0]';

// --> jacobi(A, b, zeros(3, 1), 1e-2, 1000)
//  ans  =
//    Nan

// --> gauss_seidel(A, b, zeros(3, 1), 1e-2, 1000)
//  ans  =
//    Nan

// Vemos que no convergen comenzando desde x0 = 0


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

// A(i, i) = 2      para todo i
// A(i, j) = -1     para i,j / |i-j| = 1
// A(i, j) = 0      para el resto

// Al estar aplicando el método de Gauss-Seidel,
// La matriz N será
// [
//     2  0 0 ... 0
//     -1 2 0 ... 0
//     0 -1 2 ... 0
//     ............
//     0 ... 0 -1 2
// ]

// Y su inversa aplicando eliminacion gaussiana nos queda de la forma
// [
//     2^-1  0    0    ...   0
//     2^-2  2^-1 0    ...   0
//     2^-3  2^-2 2^-1 ...   0
//     .........................
//     2^-n 2^-(n-1)...2^-2  2^-1
// ]

// Por lo que (N^-1) * A queda
// [
//     1  2^-1    0       ...         0
//     0  1-2^-2  2^-1    ...         0
//     0  2^-3    1-2^-2  ...         0
//     .........................
//     0  2^-n    2^-(n-1) ...  2^-3  1 - 2^-2
// ]

// Y al final cuando se lo restamos a la identidad
// nos queda I - (N^-1) * A:
// [
//     0  2^-1  0       ...         0
//     0  2^-2  2^-1    ...         0
//     0  2^-3  2^-2    ...         0
//     .........................
//     0  2^-n  2^-(n-1) ...  2^-3  2^-2
// ]

// Crea la matriz de iteración de A de tamaño nxn
// para el método de Gauss-Seidel 
function M = Ejercicio3(n)
    A = diag(ones(1, n) * 2) + diag(ones(1, n-1) * -1, -1) + diag(ones(1, n-1) * -1, 1)
    N = tril(A)
    M = eye(n, n) - inv(N) * A
endfunction


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
//    17.
//  ans  =
//    3.
//    4.
//   -5.

// Notemos que tardó 17 iteraciones en vez de 41 con el normal

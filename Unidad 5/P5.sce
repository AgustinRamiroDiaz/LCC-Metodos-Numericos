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

//    67.
//  ans  =

//    1.0000016
//   -2.0000015
//    2.9999973
//   -1.9999996
//    0.9999981



// --> gauss_seidel(A, b, zeros(5, 1), 1e-6, 1000)

//    38.
//  ans  =

//    1.0000009
//   -2.0000007
//    2.9999987
//   -1.9999999
//    0.9999992
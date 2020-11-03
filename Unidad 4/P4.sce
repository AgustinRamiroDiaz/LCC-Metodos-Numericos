// 1

function x = resolverTriangularSuperior(A, b)
    [n,m] = size(A) 
    x(n) = b(n) / A(n, n)
    for i = n-1 : -1 : 1
        x(i) = (b(i) - A(i, i + 1 : n) * x(i + 1 : n)) / A(i, i)
    end
endfunction

A = [2 0 3; 0 3 2; 0 0 3]
b = [1 2 3]'

resolverTriangularSuperior(A, b)
// [-1, 0, 1]'


function x = resolverTriangularInferior(A, b)
    [n,m] = size(A) 
    x(1) = b(1) / A(1, 1)
    for i = 2 : n
        x(i) = (b(i) - A(i, 1 : i - 1) * x(1 : i - 1)) / A(i, i)
    end
endfunction

A = [3 0 0; 3 2 0; 2 0 3]
b = [1 2 3]'

resolverTriangularInferior(A, b)
// [0.3333333 0.5 0.7777778]'


// 2
// a
function [x,a] = gausselim(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana sin pivoteo.  
    
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
    
    // Eliminación progresiva
    n = nA;
    for k=1:n-1
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

// b
// i
A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]' 

gausselim(A, b)
//  [-1. 2. 0. 1.]'

// ii
A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]'

gausselim(A, b)
//  [Nan Nan Nan Nan]'


// iii

A = [1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2]
b = [2 1 0 -3]'

gausselim(A, b)
// [-4. 0.6666667 -7. 1.3333333]'


// c
function [x,a,SR,MD] = gausselimCount(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // Además cuenta la cantidad de operaciones realizadas
    // La función implementa el método de Eliminación Gaussiana sin pivoteo.  
    SR = 0
    MD = 0

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
    
    // Eliminación progresiva
    n = nA;
    for k=1:n-1
        for i=k+1:n
            for j=k+1:n+1
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
                SR = SR + 1
                MD = MD + 2
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
            SR = SR + 1
            MD = MD + 1
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
        SR = SR + 1
        MD = MD + 1
    end;
endfunction

gausselimCount([1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2], [2 1 0 -3]')

// x  = 
// -4.       
//  0.6666667
// -7.       
//  1.3333333
// a  = 
//  1.   1.   0.   4.    2.
//  0.  -1.  -1.  -7.   -3.
//  0.   0.   3.   21.   7.
//  0.   0.   0.  -3.   -4.
// SR  = 
//  29.
// MD  = 
//  49.


// d
function [x,a] = gausselimCorta(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana sin pivoteo.  
    
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
    
    // Eliminación progresiva
    n = nA;
    for k=1:n-1
        for i=k+1:n
            a(i,k+1:n+1) = a(i,k+1:n+1) - a(k,k+1:n+1)*a(i,k)/a(k,k);
            a(i,1:k) = 0;  // no hace falta para calcular la solución x
        end;
    end;

    Aprima = a(:, 1:n)
    bprima = a(:, n + 1)

    // Sustitución regresiva
    x(n) = bprima(n) / Aprima(n, n)
    for i = n-1 : -1 : 1
        x(i) = (bprima(i) - Aprima(i, i + 1 : n) * x(i + 1 : n)) / Aprima(i, i)
    end
endfunction

gausselimCorta([1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2], [2 1 0 -3]')

// ans  =

//   -4.       
//    0.6666667
//   -7.       
//    1.3333333


// 4
function d = determinante(A)
    [n,m] = size(A) 
    
    if n<>m then
        error('gausselim - La matriz debe ser cuadrada');
        abort;
    end
        
    a = A
    // Eliminación progresiva
    for k=1:n-1
        for i=k+1:n
            a(i,k+1:n) = a(i,k+1:n) - a(k,k+1:n)*a(i,k)/a(k,k);
            a(i,1:k) = 0; 
        end;
    end;

    d = prod(diag(a))
endfunction

determinante([1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2])
// ans = 9








// 5

// a
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
    
    Aprima = a(:, 1:n)
    bprima = a(:, n + 1)

    // Sustitución regresiva
    x(n) = bprima(n) / Aprima(n, n)
    for i = n-1 : -1 : 1
        x(i) = (bprima(i) - Aprima(i, i + 1 : n) * x(i + 1 : n)) / Aprima(i, i)
    end
endfunction

// b
// i
A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]' 

gausselimPP(A, b)
//  [-1. 2. 0. 1.]'

// ii
A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]'

gausselimPP(A, b)
//  [-7 3 2 2]'


// iii
A = [1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2]
b = [2 1 0 -3]'

gausselimPP(A, b)
// [-4. 0.6666667 -7. 1.3333333]'




// 6
function x = resolverDiagonal(A, b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;

    n = nA

    for k = 1:n
        x(k) = b(k) / A(k, k)
    end

endfunction

function [x, cop] = resolverTridiagonal(A, b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;

    n = nA
    cop = 0     // cantidad de operaciones

    // Borro la diagonal inferior
    for k=2:n
        multiplicador = A(k,k-1) / A(k-1,k-1)
        
        A(k, k) = A(k, k) - A(k-1, k) * multiplicador
        A(k, k-1) = 0
        
        b(k) = b(k) - b(k-1) * multiplicador
        
        cop = cop + 5
    end;
    
    // Borro la diagonal superior
    for k=n-1:-1:1
        multiplicador = A(k,k+1) / A(k+1,k+1)
        
        A(k, k+1) = 0;
        
        b(k) = b(k) - b(k+1) * multiplicador

        cop = cop + 3
    end;

    x = resolverDiagonal(A, b)
endfunction

// --> A = [1 2 0 0 0; 3 4 5 0 0; 0 6 7 8 0; 0 0 9 10 11; 0 0 0 12 13];
// --> b = [1 2 3 4 5]';
// --> [x, cop] = resolverTridiagonal(A,b)
//  x  = 
//    0.122449 
//    0.4387755
//   -0.0244898
//    0.0673469
//    0.322449 
//  cop  = 
//    32.





// Ejercicio 7

// Dada una matriz A obtiene la factorizacion PA=LU 
// a partir de la eliminacion de Gauss con pivoteo parcial
function [L, U, P] = factorizacionPALU(A)
    [n,m] = size(A) 
    
    if n<>m then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    end

    U = A
    L = eye(A)
    P = eye(A)
    for k = 1:m-1
        ipivot = k; umax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k)) > umax then
                ipivot = i; umax = A(k,i);
            end;
        end;
        temp = U(ipivot, k:m); U(ipivot, k:m) = U(k, k:m); U(k, k:m) = temp;
        temp = L(ipivot, 1:k-1); L(ipivot, 1:k-1) = L(k, 1:k-1); L(k, 1:k-1) = temp;
        temp = P(ipivot, :); P(ipivot, :) = P(k, :); P(k, :) = temp;
        
        for j = k+1:m
            L(j, k) = U(j, k) / U(k, k)
            U(j, k:m) = U(j, k:m) - L(j, k) * U(k, k:m)
        end
    end
endfunction
    
// --> A = [2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8];

// --> [L, U, P] = factorizacionPALU(A)
//  L  = 

//    1.          0.          0.          0.
//    1.3333333   1.          0.          0.
//    0.6666667   0.7142857   1.          0.
//    0.3333333   0.5714286   0.3333333   1.
//  U  = 

//    6.   7.          9.          8.       
//    0.  -2.3333333  -3.         -5.6666667
//    0.   0.         -0.8571429  -0.2857143
//    0.   0.          0.          0.6666667
//  P  = 

//    0.   0.   0.   1.
//    0.   0.   1.   0.
//    0.   1.   0.   0.
//    1.   0.   0.   0.





// Ejercicio 8

// a)

// --> A = [1.012 -2.132 3.104; -2.132 4.096 -7.013; 3.104 -7.013 0.014];
// --> [L, U, P] = factorizacionPALU(A)
//  L  = 

//    1.          0.          0.
//   -0.6868557   1.          0.
//    0.3260309  -0.2142473   1.
//  U  = 

//    3.104  -7.013       0.014    
//    0.     -0.7209188  -7.003384 
//    0.      0.          1.5989796
//  P  = 

//    0.   0.   1.
//    0.   1.   0.
//    1.   0.   0.

// --> [L, U] = lu(A)
//  L  = 

//    0.3260309  -0.2142473   1.
//   -0.6868557   1.          0.
//    1.          0.          0.
//  U  = 

//    3.104  -7.013       0.014    
//    0.     -0.7209188  -7.003384 
//    0.      0.          1.5989796



// b)

// --> A = [-2.1756 4.0231 -2.1732 5.1967; -4.0231 6.0000 0 1.1973; -1.0000 5.2107 1.1111 0; 6.0235 7.0000 0 4.1561];

// --> [L, U, P] = factorizacionPALU(A)
//  L  = 

//    1.          0.          0.          0.
//   -0.6679007   1.          0.          0.
//   -0.3611854   0.6136965   1.          0.
//   -0.1660164   0.596968   -0.5112737   1.
//  U  = 

//    6.0235   7.          0.       4.1561   
//    0.       10.675305   0.       3.9731622
//    0.       0.         -2.1732   4.2595067
//    0.       0.          0.       0.4959041
//  P  = 

//    0.   0.   0.   1.
//    0.   1.   0.   0.
//    1.   0.   0.   0.
//    0.   0.   1.   0.

// --> [L, U] = lu(A)
//  L  = 

//   -0.3611854   0.6136965   1.          0.
//   -0.6679007   1.          0.          0.
//   -0.1660164   0.596968   -0.5112737   1.
//    1.          0.          0.          0.
//  U  = 

//    6.0235   7.          0.       4.1561   
//    0.       10.675305   0.       3.9731622
//    0.       0.         -2.1732   4.2595067
//    0.       0.          0.       0.4959041

// La diferencia parece radicar en que la funcion lu de Scilab 
// no utiliza LU = PA sino LU = A
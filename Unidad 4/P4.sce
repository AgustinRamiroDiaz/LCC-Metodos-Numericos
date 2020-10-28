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

function x = resolverTridiagonal(A, b)
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

    for k=2:n
        a(k, :) = a(k, :) - a(k-1, :)*a(k,k-1)/a(k-1,k-1);
    end;
    
    Aprima = a(:, 1:n)
    bprima = a(:, n + 1)

    // Sustitución regresiva
    x(n) = bprima(n) / Aprima(n, n)
    for i = n-1 : -1 : 1
        x(i) = (bprima(i) - Aprima(i, i + 1 : n) * x(i + 1 : n)) / Aprima(i, i)
    end
endfunction

// TODO contar la cantidad de operaciones






// 6

function [L U] = factorizacionPALU(A)
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

endfunction

// TODO testear con la matriz A dada por el ejercicio
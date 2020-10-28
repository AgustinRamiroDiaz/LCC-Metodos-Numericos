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
            SD = SD + 1
            MD = MD + 1
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
        SD = SD + 1
        MD = MD + 1
    end;
endfunction
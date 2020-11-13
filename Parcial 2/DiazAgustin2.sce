// Item a)

function [x,a] = gausselim(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y la matriz b.
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
    m = mA + mb;
    for k=1:n-1              // recorremos las filas
        for i=k+1:n          // cada fila se la restamos a las filas sucesivas
            for j=k+1:m    // recorremos las columnas
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k); // restamos 
            end;
            for j=1:k        // no hace falta para calcular la solución x
                a(i,j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end;
    end;
    
    // Sustitución regresiva
    x(n, :) = a(n, n+1:m)/a(n,n);
    for i = n-1:-1:1
        // acumulamos la suma para poder hacer la sustitución
        sumk = zeros(1, mb)
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k, :);
        end;
        // sustituimos con la suma
        x(i, :) = (a(i, n+1:m)-sumk)/a(i,i);
    end;
endfunction


// Item b)
A = [1 2 3; 3 -2 1; 4 2 -1];
b = [14 9 -2; 2 -5 2; 5 19 12];

x = gausselim(A, b)
disp("Item b)", x)

// Item c)
x = gausselim(A, eye(3, 3))
disp("Item c)", x)

// Resultado de Scilab:

// "Item b)"

// 1.   2.   2.
// 2.   5.   1.
// 3.  -1.  -2.

// "Item c)"

// 0.      0.1428571   0.1428571
// 0.125  -0.2321429   0.1428571
// 0.25    0.1071429  -0.1428571
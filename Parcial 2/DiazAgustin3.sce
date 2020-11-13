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
            return
        end
        x = xNuevo  // Reasignamos
    end
    x = %nan
endfunction



disp("N = 100")
N = 100;
A = 8 * eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1)+ diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
b = ones(N,1);

tic(); gausselimPP(A,b); t = toc()
disp("Gauss", t)  

tic(); gauss_seidel(A, b, zeros(N, 1), 1e-6, 10000); t = toc();
disp("Gauss-Seidel, eps = 1e-6", t)

tic(); gauss_seidel(A, b, zeros(N, 1), 1e-11, 10000); t = toc();
disp("Gauss-Seidel, eps = 1e-11", t)



disp("N = 1000");
A = 8 * eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1)+ diag(ones(N-3,1),3) + diag(ones(N-3,1),-3);
b = ones(N,1);

tic(); gausselimPP(A,b); t = toc()
disp("Gauss", t)  


tic(); gauss_seidel(A, b, zeros(N, 1), 1e-6, 10000); t = toc();
disp("Gauss-Seidel, eps = 1e-6", t)


tic(); gauss_seidel(A, b, zeros(N, 1), 1e-11, 10000); t = toc();
disp("Gauss-Seidel, eps = 1e-11", t)

// Resultado de Scilab:

// "N = 100"

// "Gauss"

//  1.5089352

// "Gauss-Seidel, eps = 1e-6"

//  0.019358

// "Gauss-Seidel, eps = 1e-11"

//  0.0436458

// "N = 1000"

// "Gauss"

//  1.5098143

// "Gauss-Seidel, eps = 1e-6"

//  0.0191995

// "Gauss-Seidel, eps = 1e-11"

//  0.0430633
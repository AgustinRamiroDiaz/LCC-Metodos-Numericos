// Resuelve el sistema Ax=b donde A es una matriz triangular superior
// a través de sustitución regresiva
function x = resolverTriangularSuperior(A, b)
    [n,m] = size(A) 
    x(n) = b(n) / A(n, n)
    for i = n-1 : -1 : 1
        x(i) = (b(i) - A(i, i + 1 : n) * x(i + 1 : n)) / A(i, i)
    end
endfunction

// Resuelve el sistema Ax=b donde A es una matriz triangular superior
// a través de sustitución progresiva
function x = resolverTriangularInferior(A, b)
    [n,m] = size(A) 
    x(1) = b(1) / A(1, 1)
    for i = 2 : n
        x(i) = (b(i) - A(i, 1 : i - 1) * x(1 : i - 1)) / A(i, i)
    end
endfunction


function [U,ind] = cholesky(A)
    // Factorización de Cholesky.
    // U'*U = A
    // Trabaja únicamente con la parte triangular superior.
    //
    // ind = 1  si se obtuvo la factorización de Cholesky.
    //     = 0  si A no es definida positiva
    //
    //******************
    eps = 1.0e-8
    //******************
    
    n = size(A,1)
    U = zeros(n,n)
    
    t = A(1,1)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end

    // Obtenemos el primer elemento aplicando el caso base del algoritmo
    U(1,1) = sqrt(t)
    for j = 2:n
        U(1,j) = A(1,j)/U(1,1)
    end
    
    for k = 2:n         // Recorremos las filas de U
        // Calculamos el elemento de la diagonal de U
        t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
        if t <= eps then
            printf('Matriz no definida positiva.\n')
            ind = 0
            return
        end
        U(k,k) = sqrt(t)    // Asignamos el elemento diagonal de U
        for j = k+1:n       // Recorremos las columnas de U desde la diagonal en adelante
            U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k) // Calculamos los valores por sobre la diagonal
        end
    end
    ind = 1
    
endfunction



// Resuelve el sistema Ax=b utilizando la factorización de cholesky
// y luego haciendo 2 sustituciones (progresiva y regresiva)
// donde A=U'*U
function x = resolverCholesky(U, b)
    g = resolverTriangularInferior (U', b)
    x = resolverTriangularSuperior (U, g)
endfunction

// Vemos que con un pequeño cambio de filas (1 por 3)
// obtenemos una matriz simétrica que podemos factorizar con cholesky
// debido a que es simétrica.
// De esta manera nos ahorramos hacer cálculos con la factorización LU
// ya que solamente trabajamos con la parte superior de la matriz A.
// Además, luego podemos aplicar eficientemente una solución 
// con sustituciones progresivas y regresivas
A = [4 2 1 0 0; 2 5 1 1 0; 1 1 6 2 -1; 0 1 2 4 1; 0 0 -1 1 3]
[U, ind] = cholesky(A)

// Almacenaremos en x todas las soluciones de Ax=b
// donde cada fila nos representará la solución
// partiendo de c=-10 hasta c=10
x = zeros(21, 5)
for c = -10:10
    // Para mantener la coherencia con el sistema original,
    // también haremos el cambio de filas en b
    b = [c 1 -c 1 1]'
    x(c + 11, :) = resolverCholesky(U, b)
end

// Vemos las soluciones de x para los c
disp(x)

// Gráficamos x1 en función de c
plot(-10:10, x(:, 1)')

// Resultado de Scilab:

// column 1 to 4
// -4.1428571   1.6785714   3.2142857  -2.3214286
// -3.7391304   1.5326087   2.8913043  -2.076087 
// -3.3354037   1.386646    2.568323   -1.8307453
// -2.931677    1.2406832   2.2453416  -1.5854037
// -2.5279503   1.0947205   1.9223602  -1.3400621
// -2.1242236   0.9487578   1.5993789  -1.0947205
// -1.7204969   0.802795    1.2763975  -0.8493789
// -1.3167702   0.6568323   0.9534161  -0.6040373
// -0.9130435   0.5108696   0.6304348  -0.3586957
// -0.5093168   0.3649068   0.3074534  -0.113354 
// -0.1055901   0.2189441  -0.015528    0.1319876
//  0.2981366   0.0729814  -0.3385093   0.3773292
//  0.7018634  -0.0729814  -0.6614907   0.6226708
//  1.1055901  -0.2189441  -0.984472    0.8680124
//  1.5093168  -0.3649068  -1.3074534   1.113354 
//  1.9130435  -0.5108696  -1.6304348   1.3586957
//  2.3167702  -0.6568323  -1.9534161   1.6040373
//  2.7204969  -0.802795   -2.2763975   1.8493789
//  3.1242236  -0.9487578  -2.5993789   2.0947205
//  3.5279503  -1.0947205  -2.9223602   2.3400621
//  3.931677   -1.2406832  -3.2453416   2.5854037
//        column 5
//  2.1785714
//  1.9891304
//  1.7996894
//  1.6102484
//  1.4208075
//  1.2313665
//  1.0419255
//  0.8524845
//  0.6630435
//  0.4736025
//  0.2841615
//  0.0947205
// -0.0947205
// -0.2841615
// -0.4736025
// -0.6630435
// -0.8524845
// -1.0419255
// -1.2313665
// -1.4208075
// -1.6102484
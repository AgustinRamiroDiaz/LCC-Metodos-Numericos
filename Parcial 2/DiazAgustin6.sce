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

// La función mixto resuelve el sistema Ax=b
// partiendo de un x0 inicial aplicando el método SOR con parámetro w
// hasta alcanzar una tolerancia de 10^(-beta/2).
// Luego aplica el método de Gauss-Seidel hasta alcanzar
// una tolerancia de 10^(-beta)
function x = mixto(A, b, x0, beta, w)
    x = SOR(A, b, w, x0, 10^(-beta/2), 1000)
    x = gauss_seidel(A, b, x, 10^(-beta), 1000)
endfunction

n = 9
A = diag(ones(1, n) * 2) + diag(ones(1, n-1), -1) + diag(ones(1, n-1), 1)
b = [1 2 3 4 5 4 3 2 1]'

m = mixto(A, b, zeros(9, 1), 8, 1.5)
disp(m)

g = gauss_seidel(A, b, zeros(9, 1), 1e-8, 1000)
disp(g)


// Respuesta de Scilab:

// 24.

// 73.

// 0.5      
// 2.867D-08
// 1.5      
// 4.195D-08
// 2.5      
// 3.795D-08
// 1.5      
// 2.121D-08
// 0.5      

// 170.

// 0.5      
// 2.858D-08
// 1.5      
// 4.183D-08
// 2.5      
// 3.784D-08
// 1.5      
// 2.115D-08
// 0.5      


// Vemos que las soluciones son idénticas
// pero que el método mixto realizó 97 iteraciones:
//     24 iteraciones con SOR
//     73 iteraciones con Gauss Seidel
// mientras que el método de Gauss Seidel
// realizó 170 iteraciones
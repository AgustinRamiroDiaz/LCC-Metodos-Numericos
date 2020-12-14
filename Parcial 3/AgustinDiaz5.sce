// Pregunta 5


// Aproxima la funcion que pasa por los puntos
// x con imagen y
// por el metodo de los minimos cuadrados
// a partir de los 2 polinomios:
//      x 
//      x^2
// retorna los coeficientes correspondientes
function a = minimosCuadrados(x, y)
    [_, n] = size(x)
    A = zeros(n, 2)
    for i = 1:n
        A(i, 1) = x(i)
        A(i, 2) = x(i)^2
    end
    // Calculamos los coeficientes utilizando el Teorema 2
    a = (A' * A) \ (A' * y')
endfunction

x = 0:0.4:8
y = [10.61 14.89 11.50 -4.58 13.65 23.47 24.29 23.33 22.80 45.03 25.85 46.46 40.84 65.50 75.96 79.91 93.19 126.21 147.40 158.78 193.23]


// f(0) = 8 => a = 8
// obtenemos los coeficientes b y c
coeficientes = minimosCuadrados(x, y-8)
a = 8
b = coeficientes(1)
c = coeficientes(2)

p = poly([a b c], 'x', 'coeff')

plot2d(x, horner(p, x), style = color('magenta'))   // En magenta la aproximacion por minimos cuadrados
plot2d(x, y, style = color('blue'))                 // En azul los datos

// Aunque parezca tener error en varios casos
// hay que recordar que esta tabla de valores tiene errores de medicion
// Y que una curva suave muchas veces tiende a ser una funcion aproximada
// que describe mejor la realidad
// Ademas el error parece ser aleatorio, lo cual es una propiedad del error de medicion

// Por estas razones creo que es una buena aproximacion
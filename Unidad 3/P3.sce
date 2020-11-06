// Ejercicio 1
deff('y=f1(x)', 'y=1 + (cos(x) * ((exp(x) + exp(-x))/ 2))')
function y = dibujar(f, x)
  plot(x,f)
  a=gca();
  a.x_location = "origin";
  a.y_location = "origin";
  xgrid(2)
  y = 0
endfunction
// --> dibujar(f1, [0:0.1:8])
// Podemos ver que las primeras 3 raices positivas son aproximadamente:
// x1 = 2
// x2 = 4.5
// x3 = 8

// metodo_biseccion :: (Float -> Float) Float Float Float -> Float
// Dado una funcion y dos puntos (a y b / f(a) * f(b) < 0)
// calcula la raiz de f en el interbalo [a, b] con un error máximo de eps
function v = metodo_biseccion(f, a, b, eps)
  if (f(a) * f(b)) >= 0 then
    v = %nan
  else
    c = (a + b) / 2
    while b - c > eps then
      if (f(b) * f(c)) <= 0 then
        a = c
      else
        b = c
      end
      c = (a + b) / 2
    end

    v = c
  end
end

// Ejercicio 2

// a)

// --> deff('y=g(x)', 'y = sin(x) - ((x^2)/2)')
// --> dibujar(g, [-1:0.1:2])
// --> metodo_biseccion(g, -1, 1, 1e-2)
//  ans  =
//    0.0078125
// --> metodo_biseccion(g, 1, 2, 1e-2)
//  ans  =
//    1.3984375

// b)
// --> deff('y=g(x)', 'y = %e^(-x) - (x^4)')
// --> dibujar(g, [-10:0.1:2])
// --> metodo_biseccion(g, -10, -5, 10^-2)
//  ans  =
//    -8.603515625
// --> metodo_biseccion(g, -5, 0, 1e-2)
//  ans  =
//    -1.435546875
// --> metodo_biseccion(g, 0, 2, 10^-2)
//  ans  =
//    0.8203125

// c)

// Utilizando logaritmo natural
// --> deff('y=g(x)', 'y = x - 1 - log(x)')
// --> dibujar(g, [0:0.1:5])
// --> metodo_biseccion(g, 0.5, 1.5, 1e-2)
// ans  =
//    Nan
// Vemos que no se puede calcular con el método de la bisección
// ya que la función es no negativa

// Utilizando logaritmo base 10
// --> deff('y=g(x)', 'y = x - 1 - log10(x)')
// --> dibujar(g, [0:0.1:5])
// --> metodo_biseccion(g, 0.5, 1.5, 1e-2)
//  ans  =
//    1.0078125

// --> metodo_biseccion(g, 0, .5, 1e-2)
//  ans  =
//    0.1328125


// Ejercicio 3
// Calcula por el método de la secante la raíz de f
// En el intervalo [anterior, actual] con una tolerancia epsilon
// (supone f continua)
// Convergencia mas rapida que lineal
// Convergencia no asegurada
// Tratar de que F'(a) != 0
function output = secante(f, anterior, actual, epsilon)
    while abs (f (actual) - f (anterior)) >= epsilon
        siguiente = actual - f (actual) * (actual - anterior) / (f (actual) - f (anterior))
        anterior = actual        
        actual = siguiente 
    end
    output = actual
endfunction

// --> deff('y = g(x)', 'y = x * x / 4 - sin(x)')
// --> dibujar(g, [0:.1:5]);
// --> secante(g, 1.5, 2, 1e-12)
//  ans  =
//    1.933753762827021




// Ejercicio 4
// Usamos el Teorema 2 con g = cos, y obtenemos que existe una única solución de g(x) = x.
// Por lo tanto, al aplicar reiteradas veces la función coseno obtenemos su punto fijo

// Aplica el método de punto fijo con la función f
// partiendo desde x con una tolerancia de epsilon
function output = metodo_punto_fijo(f, x, epsilon)
    while abs (f (x) - x) > epsilon
        x = f (x)
    end
    output = x
endfunction


//--> metodo_punto_fijo(cos, 100, 1e-12)
//   ans  =
//     0.739085133214757


// Ejercicio 5
// Usamos el Teorema 2
// g (x) = 2 ^ (x - 1)
// Buscamos que g' (x) < 1 => x < 1.53
// Definimos a = -inf, b = 1.53
// Corroboramos que si -inf <= x <= 1.53 => -inf <= g (x) <= 1.53
// Por lo tanto, converge si partimos de un x < 1.53
// Como la solución es única y sabemos que 1 es solución, entonces 1 es la única solución
// Por lo tanto converge a 1

// --> deff ('y = g(x)', 'y = 2^(x-1)')
// --> metodo_punto_fijo(g, 1.5, 1e-12)
//  ans  =
//    1.000000000002622

// --> metodo_punto_fijo(g, -100, 1e-12)
//  ans  =
//    0.999999999997722



// Ejercicio 6
// REVISAR HOJAS
// Utilizaremos el Corolario 1 con alpha = -srqt(5)
// Tenemos que g y g' son continuas en R, en particular alrededor de alpha
// además 
//      g'(x) = 1 + 2*c*x => g'(-sqrt(5)) = 1 - 2 * c * sqrt(5)
// por lo tanto
//      |g'(-sqrt(5))| < 1   <=>    0 < c < 1 / sqrt(5)
// Además notemos que cuando c = 1/(2 * sqrt(5)), g'(-sqrt(5)) = 0 
// dandonos el valor óptimo para la convergencia

// c = 1/(2 * sqrt(5))
// deff('y = g(x)', 'y = x + c * (x^2 - 5)')
// --> metodo_punto_fijo(g, 0, 1e-12)
//  ans  =
//   -2.23606797749979
// --> -sqrt(5)
//   ans  =
//    -2.23606797749979


// Aplica el método de newton con la función f
// partiendo desde x con una tolerancia de epsilon 
// y un cantidad máxima de iteraciones iter
function v = metodo_newton(f, x0, e, iter)
  i = 0
  xn = x0 - (f(x0) / numderivative(f, x0))
  while (abs(xn - x0) > e) && (i < iter) then
    x0 = xn
    xn = x0 - (f(x0) / numderivative(f, x0))
    i = i + 1
  end

  v = xn
endfunction



// Ejercicio 7

function output = longitudDeOnda(d)
    output = 4 * %pi ^ 2 / (5 ^ 2 * 9.8 * tanh (4 * d)) 
endfunction

// a)
// --> metodo_punto_fijo(longitudDeOnda, 1, 1e-1)
//  ans  =
//    0.283551260367738


// b)
// TODO




// Ejercicio 8

// --> deff("y=g1(x)", "y=%e^x / 3")
// --> dibujar(g1, [-1:0.1:2]);
// Podemos ver que 0 <= x <= 1 => 0 <= g1(x) <= 1
// Y como |g1'(x)| < 1 en [0, 1]
// Podemos obtener la solucion con:
// --> metodo_punto_fijo(g1, 0, 1e-4)
//  ans  =
//    0.6188108

// --> deff("y=g2(x)", "y=(exp(x) - x) / 2")
// --> dibujar(g2, [-1:0.1:2])
// Podemos ver que 0 <= x <= 1 => 0 <= g1(x) <= 1
// Y como |g2'(x)| < 1 en [0, 1]
// Podemos obtener la solucion con:
// --> metodo_punto_fijo(g2, .5, 1e-4)
//  ans  =
//    0.618952

// --> deff("y=g3(x)", "y=log(3*x)")
// --> dibujar(g3, [0:0.1:5]);
// Podemos ver que 1.5 <= x <= 2 => 1.5 <= g3(x) <= 2
// Entonces podemos obtener la solucion con:
// Y como |g3'(x)| < 1 en [1.5, 2]
// --> metodo_punto_fijo(g3, 1.5, 10^-4)
//  ans  =
//    1.5119381

// --> deff("y=g4(x)", "y= %e^x - 2 * x")
// --> dibujar(g4, [0:0.1:2]);
// Vemos que 0 <= x <= 1 => 0 <= g4(x) <= 1
// Y como |g4'(x)| < 1 en [0, 1]
// Obtenemos el resultado con:
// --> metodo_punto_fijo(g4, .5, 10^-4)
//  ans  =
//    0.6189904

// Todas las funciones nos fueron útiles
// En particular, g3 nos permitió encontrar el punto fijo cercano a 1.5
// mientras que el resto nos permitieron entcontrar el cercano a 0.6



// Dado un sitema de ecuaciones F y un punto inicial x,
// obtiene una solución aproximada al sistema F(v) = 0 con el método de newton multivariable
// con una tolerancia de epsilon y un máximo iter de iteraciones
function v = metodo_newton_multivariable(F, x, eps, iter)
    for i = 1: iter
        Jinv = inv(numderivative(F, x))
        x = x - (Jinv * F(x))
        if norm(F(x)) < eps then
            v = x; return
        end
    end
    v = %nan
endfunction


// Ejercicio 9
function y = F9(x)
  y(1) = 1 + x(1)^2 - x(2)^ 2 + exp(x(1)) * cos(x(2))
  y(2) = (2 * x(1) * x(2)) + (exp(x(1)) * sin(x(2)))
endfunction

// --> metodo_newton_multivariable(F9, [-1,4]', 1e-8, 5)
//  ans  =
//   -0.2931627
//    1.1726598


// Ejercicio 10
function y = F10(x)
  y(1) = x(1)^2 + (x(1) * x(2)^3) - 9
  y(2) = (3 * x(1)^2 * x(2)) - 4 - x(2)^3
endfunction

// a)
// --> metodo_newton_multivariable(F10, [1.2, 2.5]', 1e-8, 100)
//  ans  =
//    1.3363554
//    1.7542352
// --> F10(ans)
//  ans  =
//   -1.776D-15
//   -1.776D-15


// b)
// --> metodo_newton_multivariable(F10, [-2, 2.5]', 1e-8, 100)
//  ans  =
//   -0.9012662
//   -2.0865876
// --> F10(ans)
//   ans  =
//   7.726D-10
//   1.821D-10


// c)
// --> metodo_newton_multivariable(F10, [-1.2, -2.5]', 1e-8, 100)
//  ans  =
//   -0.9012662
//   -2.0865876
// --> F10(ans)
//   ans  =
//     1.776D-15
//     3.553D-15


// d)
// --> metodo_newton_multivariable(F10, [2, -2.5]', 1e-8, 100)
//  ans  =
//   -3.0016249
//    0.148108 
// --> F10(ans)
//    ans  =
//     -4.086D-14
//      5.381D-13
  


// Ejercicio 11
function y = f(x)
  y = 2 * x(1) + 3 * x(2)^2 + exp((2 * x(1)^2) + x(2)^2)
endfunction

function y = dF(x)
  y(1) = 2 + exp(2*x(1)^2 + x(2)^2) * 4*x(1)
  y(2) = 6 * x(2) + exp(2*x(1)^2 + x(2)^2) * 2*x(2)
endfunction

// --> [v, tipo] = minimo_maximos_locales(dF, [1,1]', 10^-12)
//  tipo  = 
//    1.
//  v  = 
//   -0.3765446
//    2.640D-16

// --> f(v)
//  ans  =
//    0.5747748



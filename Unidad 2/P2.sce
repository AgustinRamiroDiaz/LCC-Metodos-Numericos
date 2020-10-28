// Ejercicio 1

//Calcula las raíces de forma robusta para un polinomio de grado 2
function x = raices(p)
    a = coeff(p, 2)
    b = coeff(p, 1)
    c = coeff(p, 0)
    // Veo que el discriminante no sea negativo
    if b*b - 4 * a * c < 0 then x(1) = %nan; x(2) = %nan; return x end

    if b < 0 
        then
        x1 = 2*c / (-b + sqrt(b ^ 2 - 4 * a * c));
        x2 = -b + sqrt(b ^ 2 - 4 * a * c);
        else
            if b == 0
                then 
                    x1 = sqrt(-c / a)
                    x2 = - x1
                else
                    x1 = -b - sqrt(b ^ 2 - 4 * a * c);
                    x2 = 2*c / (-b - sqrt(b ^ 2 - 4 * a * c));
            end
    end
    x(1) = x1
    x(2) = x2
endfunction

// Ejercicio 3
// a
function v = hornerNew(p, x)
    v = hornerAux (p, x, 0)
    
endfunction

function v = hornerAux(p, x, n)
    if n == degree(p) 
        then v = coeff(p, n)
        else v = coeff(p, n) + x * hornerNewAux(p, x, n+1)
    end
endfunction

function b = horner1(p, x)
    b = coeff(p, degree(p)) 
    for (i = degree(p) - 1: -1: 0)
        b = coeff(p, i) + x * b
    end
endfunction

//d


function x = horner2(p, x)
    b = coeff(p, degree(p)) 
    c = b * x ^ (degree(p) - 1)
    for (i = degree(p) - 1: -1: 1)
        b = coeff(p, i) + x * b
        c = c + b * x ^ (i - 1)
    end
    b = coeff(p, 0) + x * b
    x(1) = b
    x(2) = c
endfunction


//Ejercicio 4
function r = derivar(f, v, n, epsilon)
    if n == 0 
        then r = f (v)
        else r = (derivar (f, v + epsilon, n-1) - derivar (f, v, n-1)) / epsilon
    end
endfunction


function x = derivarIterativo(f,v,n,h)
    deff("y=D0F(x)", "y="+f)
    for i = 1:1:n-1
        deff("y=D"+string(i)+"F(x)", "y=(D"+string(i-1)+"F(x+h)-D"+string(i-1)+"F(x))/h")
    end
    deff ("y=DnF(x)", "y=(D"+string(n-1)+"F(x+h)-D"+string(n-1)+"F(x))/h")
    x = DnF(v)
endfunction

// a
// Estamos aproximando la función derviada a partir de un epsilon
// También tenemos errores de redondeo

//b
// dividimos n veces por epsilon, lo cual hace que dividimos por un número muy cercano a 0

//5
// crea un polinomio de Taylor de la función f alrededor del punto a
// de grado n, aproximando sus derivadas con epsilon
function r = Taylor(f, v, n, epsilon, a)
    r = 0
    for (i = 0:n) 
        r = r + derivar (f, a, i, epsilon) * (v - a) ^ i / factorial (i)
    end
endfunction

function f = factorial(n)
    f = 1
    for (i = 1:n)
        f = f * i
    end
endfunction

//6

//a
// poly([1 1 1/2 1/6 1/factorial(4) 1/factorial(5) 1/factorial(6) 1/factorial(7) 1/factorial(8) 1/factorial(9) 1/factorial(10)], "x", "c")
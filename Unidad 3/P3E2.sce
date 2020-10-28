//2

//Recibe una funcion continua en un intervalo y dos puntos del intervalo
//Convergencia lineal asegurada
function output = biseccion(f, a, b, epsilon)
    // Retorna una raíz de f entre a y b con un error de a lo sumo epsilon
    if f(a)*f(b) > 0 then output = %nan; return end
    c = (a + b) / 2
    if b - c <= epsilon then 
        output = c
    else 
        if f (b) * f (c) <= 0 then 
            output = biseccion (f, c, b, epsilon)
        else 
            output = biseccion (f, a, c, epsilon)
        end
    end
endfunction


// a
function output = unoa(x)
    output = sin (x) - x * x / 2
endfunction

biseccion(unoa, 0, 1, 10^-2) 
biseccion(unoa, 1, 2, 10^-2) 
// 0.0078125    1.3984375


// b
function output = unob(x)
    output = %e ^ (-x) - x ^ 4
endfunction

biseccion(unob, -2, 0, 10^-2) 
biseccion(unob, 0, 1, 10^-2) 
// -1.4296875   0.8203125


// c
function output = unoc(x)
    output = log (x) - x + 1
endfunction

// 
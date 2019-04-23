# solving nonlinear equations

# solving : f(x) = 0
# need:
#   > bisection or any other algo to find approximation x_0
#   > iterative algorithms for actual solving

def bisect(f, a, b, eps = 0.01):
    if abs(b - a) <= eps:
        return [(a, b)]
    
    x_0 = (a + b) / 2
    
    f_a, f_x, f_b = f(a), f(x_0), f(b)
    
    roots = []
    
    if f_x == 0:
        roots += [(x_0, x_0)]
    else:
        if f_a * f_x < 0:
            roots += bisect(f, a, x_0, eps)
        if f_x * f_b < 0:
            roots += bisect(f, x_0, b, eps)
    
    return roots


def simple_iteration_base(f, x, kernel, eps = 0.00001, MAX_ITER = 100000):
    x_next = x
    iters = 0
    delta = 1e5
    while abs(delta) >= eps and iters <= MAX_ITER:
        x = x_next
        x_next = kernel(f, x)
        print(f'#it : {iters}-> ', x, x_next)
        
        delta = x_next - x
        iters += 1
        
        
    return x_next


def aitken_iteration_base(f, x, kernel, eps = 0.00001, MAX_ITER = 100000):
    x_next = x
    x_next_next = x
    iters = 0
    x_star = x
    delta = 1e5
    
    while abs(delta) >= eps and iters <= MAX_ITER:
        x = x_next
        x_next = kernel(f, x)
        x_next_next = kernel(f, x_next)
        x_star = x_next
        
        if (x_next_next - 2*x_next + x) == 0:
            print('oops')
            return x_star
        
        if iters >= 3:
            x_star = x_next_next - (x_next_next - x_next)**2 / (x_next_next - 2*x_next + x)
            
        print(f'#it : {iters}-> ', x, x_next, x_next_next, x_star)
            
        delta = x - x_next_next    
        iters += 1
        
    return x_star


def s_tau(tau):
    def _tau(f, x):
        nonlocal tau
        return x - tau * f(x)
    return _tau


def simple_iter_tau(f, x, tau):
    s = s_tau(tau)
    return simple_iteration_base(f, x, s)


def approx_deriv(f, x, eps = 0.0001):
    return (f(x+eps) - f(x)) / eps

def secants_method(f, x_0, eps = 0.00001, MAX_ITER = 100000):
    iters = 0
    delta = 1e4
    x = {'n-1' : x_0, 'n' : approx_deriv(f, x_0), 'n+1' : None}
    while abs(delta) >= eps and iters <= MAX_ITER:
        x['n+1'] = x['n'] - (x['n'] - x['n-1']) * f(x['n']) / (f(x['n']) - f(x['n-1']))
        delta = x['n+1'] - x['n']
        iters += 1
        
    return x['n+1']
    
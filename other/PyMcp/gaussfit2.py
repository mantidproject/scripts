## Parametric function: 'v' is the parameter vector, 'x' the independent varible

fp = lambda v, x: (v[0]**2)*exp(-(x-v[1])**2/(2*v[2]**2)) 


## Error function
e = lambda v, x, y: (fp(v,x)-y)


## Initial parameter value guess from moments

cen = sum(x*y)/sum(y)
width = sqrt(abs(sum((x-cen)**2*y)/sum(y)))
max = y.max()

v0 = [max, cen, width]

## Fitting
v, success = leastsq(e, v0, args=(x,y), maxfev=10000)

plot(x, fp(v,x))


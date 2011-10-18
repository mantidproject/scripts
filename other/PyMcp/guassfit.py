h=.0010
cen=782
fw=10
gaussian = lambda x: h*exp(-(cen-x)**2/fw)

data = gaussian(arange(100))

plot(data)

X = arange(data.size)
x = sum(X*data)/sum(data)
width = sqrt(abs(sum((X-x)**2*data)/sum(data)))

max = data.max()

fit = lambda t : max*exp(-(t-x)**2/(2*width**2))

plot(fit(X))

show()
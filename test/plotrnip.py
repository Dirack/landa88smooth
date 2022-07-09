import matplotlib.pyplot as plt

f =  open("vel.txt","r")
xstr = filter(None,f.read().split('\n'))
x = map(lambda i: float(i),xstr)

f =  open("rnip.txt","r")
zstr = filter(None,f.read().split('\n'))
z = map(lambda i: float(i),zstr)

plt.plot(list(x),list(z))
plt.show()

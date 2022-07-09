import matplotlib.pyplot as plt

f =  open("vel.txt","r")
xstr = filter(None,f.read().split('\n'))
x = map(lambda i: float(i),xstr)

f =  open("semb.txt","r")
ystr = filter(None,f.read().split('\n'))
y = map(lambda i: float(i),ystr)

plt.plot(list(x),list(y))
plt.show()

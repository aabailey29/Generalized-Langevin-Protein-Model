import numpy as np
import matplotlib.pyplot as plt

# a = 0.0002787
# gamma = 0.02

a = 0.05178716240853
gamma = 0.2


timev = np.zeros(51)
data = np.zeros(51)
f = open(r"C:\Users\alana\OneDrive\Documents\Thesis\VACF_10ps.txt", "r")
line = f.readline()
timev[0] = float(line.split()[0])
data[0] = float(line.split()[1])
# print(line)
i = 1
while line and i < 51:
    line = f.readline()
    str = line
    #print(str)
    spl = str.split()
    #print(spl)
    timev[i] = float(spl[0])
    data[i] = float(spl[1])
    i += 1
          
f.close()
timev = timev


t_not = 0
t_end = 500
steps = 1000
R = 100
kBT = 0.16707718278387142
h = (t_end - t_not)/steps
timevec = np.linspace(t_not, t_end, steps+1)
theta =np.exp(-h*gamma)
alpha =np.sqrt(np.square(1-theta)/h)
alpha1 = alpha*np.sqrt(2*kBT*a*gamma)
vx = np.zeros(steps+1)
vy = np.zeros(steps+1)
x = np.zeros(steps+1)
y = np.zeros(steps+1)
svecx = np.zeros(steps+1)
svecy = np.zeros(steps+1)
vt100 = np.zeros(R, dtype=float)
svecx[0] = 0
svecy[0] = 0
VACF = np.zeros(steps+1)
v0v0 = 0


for j in range (R):
    vx[0] = 0.1*np.random.normal(loc = 0.0, scale = 1.0)
    vy[0] = 0.1*np.random.normal(loc = 0.0, scale = 1.0)
    x[0] = 0.1*np.random.normal(loc = 0.0, scale = 1.0)
    y[0] = 0.1*np.random.normal(loc = 0.0, scale = 1.0)
    for i in range (1, steps+1):        
        vhalfx = vx[i-1] + (h/2)*svecx[i-1]
        vhalfy = vy[i-1] + (h/2)*svecy[i-1]
        x[i] = x[i-1] + h*vhalfx 
        y[i] = y[i-1] + h*vhalfy 
        svecx[i] = theta*svecx[i-1] - (1-theta)*(a/gamma)*vhalfx + alpha1*np.random.normal()
        svecy[i] = theta*svecy[i-1] - (1-theta)*(a/gamma)*vhalfy + alpha1*np.random.normal()
        vx[i] = vhalfx + (h/2)*svecx[i]
        vy[i] = vhalfy + (h/2)*svecy[i]
    vt100[j] = vx[steps - 1]
    vt100[j] = vy[steps - 1]
    VACF += vx[0]*vx + vy[0]*vy
    v0v0 += vx[0]*vx[0] + vy[0]*vy[0]
VACF = VACF/v0v0
#    plt.plot(x,y)
#print(vt100)    
#plt.hist(vt100, bins=40, weights=np.ones(len(vt100)) / len(vt100))
#plt.xlabel("velocity at t = 100")
#plt.ylabel("percent")

# plt.xlabel("time")
# plt.ylabel("VACF")
# plt.title("VACF Comparison")
# plt.grid()
# plt.plot(data, label = "MD VACF",)
# plt.plot(VACF, label = "GLE VACF")
# plt.legend()
# plt.show()

# plt.xlabel("time")
# plt.ylabel("position x")
# plt.plot(timevec, x)
# plt.show()

# plt.xlabel("time")
# plt.ylabel("position y")
# plt.plot(timevec, y)
# plt.show()

plt.title("Position of GLE Simulated Protein")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x,y)
plt.show()

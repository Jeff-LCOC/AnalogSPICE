import matplotlib.pyplot as plt

filename = input("Enter the output file name (without .lis extension): ")
filename = './tests/' + filename + '.lis'
nodename = input("Enter the node to be measured: ")
currentT = []
node = []

with open(filename, 'r') as file:
    for line in file:
        if line.startswith('currentT:'):
            currentT.append(float(line.split(':')[1]))
        elif line.startswith('node ' + nodename):
            node.append(float(line.split(':')[1]))

plt.plot(currentT, node)
plt.xlabel('currentT')
plt.ylabel('node')
plt.title('Graph')
plt.show()

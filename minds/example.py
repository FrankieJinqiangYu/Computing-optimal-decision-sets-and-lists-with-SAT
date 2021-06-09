from gurobipy import *

import sys
import os
# Create a model.
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
model = Model()
sys.stdout = old_stdout

# model.Params.OutputFlag = 0


# Add by row (easier)
p2 = model.addVar(obj=2)
p3 = model.addVar(obj=2)
p7 = model.addVar(obj=3)

e1 = model.addConstr(p7 >= 1)
e2 = model.addConstr(p3 >= 1)
e3 = model.addConstr(p2 >= 1)
e4 = model.addConstr(p3 >= 1)
e5 = model.addConstr(p7 >= 1)



# Add by column
# e1 = model.addConstr(0 >= 1)
# e2 = model.addConstr(0 >= 1)
# e3 = model.addConstr(0 >= 1)
# e4 = model.addConstr(0 >= 1)
# e5 = model.addConstr(0 >= 1)

# col = Column([0,1,1,1,0],[e1,e2,e3,e4,e5])
# p4 = model.addVar(obj=2, column=col)








model.optimize()
print('--------------------------------------')

print(f'{p2.X}')
print(f'{p3.X}')
print(f'{p7.X}')

# Retrieve dual value of the rows (constraints)
u1 = e1.Pi
u2 = e2.Pi
u3 = e3.Pi
u4 = e4.Pi
u5 = e5.Pi

print(f'Dual value of row e1: u1 = {u1}')
print(f'Dual value of row e2: u2 = {u2}')
print(f'Dual value of row e3: u3 = {u3}')
print(f'Dual value of row e4: u4 = {u4}')
print(f'Dual value of row e5: u5 = {u5}')
print('--------------------------------------')



# Find any negative c or (better) minimize c and add all negative c
# Replace this with maxsat oracle

print(f'P4 oracle: c = {2 - u2 - u3 - u4}')
# if c < 0: do this
col = Column([0,1,1,1,0],[e1,e2,e3,e4,e5])
p4 = model.addVar(obj=2, column=col)





print(f'P5 oracle: c = {2 - u1 - u5}')
print(f'P1 oracle: c = {2 - u4}')
print(f'P6 oracle: c = {2 - u1}')













model.optimize()
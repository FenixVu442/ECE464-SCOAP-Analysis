# SCOAP method for finding Controlability
# This file contains SCOAP for Boolean Operation:
# NOT , BUF
# AND , OR
# NAND , NOR
# XOR, XNOR
# note: listInputs is a list of { [c0, c1] }
def scoapNOT(listInputs):
  c0 = listInputs[0][1] + 1
  c1 = listInputs[0][0] + 1
  return [c0, c1]

def scoapAND(listInputs):
  c0 = min([pair[0] for pair in listInputs]) + 1
  c1 = sum([pair[1] for pair in listInputs]) + 1
  return [c0, c1]

def scoapOR(listInputs):
  c0 = sum([pair[0] for pair in listInputs]) + 1
  c1 = min([pair[1] for pair in listInputs]) + 1
  return [c0, c1]

def _scoapXOR2(pair1, pair2):
  p00 = pair1[0] + pair2[0]
  p11 = pair1[1] + pair2[1]
  c0 = min(p00, p11)
  p01 = pair1[0] + pair2[1]
  p10 = pair1[1] + pair2[0]
  c1 = min(p01, p10)
  return [c0, c1]

def scoapXOR(listInputs):
  output = listInputs[0]
  for i in range(1, len(listInputs)):
    output = _scoapXOR2(output, listInputs[i])
  c0 = output[0] + 1
  c1 = output[1] + 1
  return [c0, c1]

def scoapNAND(listInputs):
  c0 = sum([pair[1] for pair in listInputs]) + 1
  c1 = min([pair[0] for pair in listInputs]) + 1
  return [c0, c1]

def scoapNOR(listInputs):
  c0 = min([pair[1] for pair in listInputs]) + 1
  c1 = sum([pair[0] for pair in listInputs]) + 1
  return [c0, c1]

def _scoapXNOR2(pair1, pair2):
  p00 = pair1[0] + pair2[0]
  p11 = pair1[1] + pair2[1]
  c1 = min(p00, p11) + 1
  p01 = pair1[0] + pair2[1]
  p10 = pair1[1] + pair2[0]
  c0 = min(p01, p10) + 1
  return [c0, c1]

def scoapXNOR(listInputs):
  output = listInputs[0]
  for i in range(1, len(listInputs)):
    output = _scoapXNOR2(output, listInputs[i])
  c0 = output[0] + 1
  c1 = output[1] + 1
  return [c0, c1]

def scoapBUFF(listInputs):
  c0 = listInputs[0][0] # + 1
  c1 = listInputs[0][1] # + 1
  return [c0, c1]

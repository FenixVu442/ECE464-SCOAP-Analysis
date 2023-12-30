#######################
# Author: Manh Vu Bui
# Class: ECE 464
# Project 2: SCOAP vs. MC simulation
# Prof. Rao Wenjing
# 11/12/2023
#######################
import time
from random import seed, randint
from Dalgebra import *
from scoap import *


class Parser:
  # Constructor
  def __init__(self):
    self.varIndex = {}
    self.varMap = {}
    self.nodeLevel = {}
    self.sortedNode = []
    self.faultList = []
    self.inputList = []
    self.outputList = []
    self.fault_count = 0  # number of faults for a given bench
    self.size = 0  # The total number of nodes
    self.faultNode = ''
    # varIndex (key, value) = (nodeName, nodeIndex), possibily not used at all
    # varMap (key, value) = (nodeName, list)
    #   ***********************************************
    #   * list[0]| list[1]| list[2]| list[3]| list[4]| list[5]
    #   * index  | name   | type   | value  | gate   | listInputs[]
    #   *   added: list[6] = list of associated faults
    #   *   added: list[7] = controlability, [-1, -1] for unknown
    #   ***********************************************
    #   Ex:
    #   INPUT: 'a': [1, 'a', 'INPUT', 'u', '___', []]
    #   OUTPUT: 'y': [8, 'y', 'OUTPUT', 'u', 'OR', [14, 15]]
    #   in-network: 'd': [10, 'd', '___', 'u', 'AND', [1, 12]]
    # nodeLevel (key, value) = (nodeName, nodeLevel)
    # sortedNode (list) = list of all nodes
    #                  in the order of ascending level values
    # faultList (list) = list of all Faults
    # inputList (list) = list of nodes with type='INPUT'
    # outputList (list) = list of nodes with type='OUTPUT'

  # Read file
  def readFile(self, fileName):
    f = open(fileName, "r")
    lines = f.readlines()
    lines = [i for i in lines if i != '\n' and i[0] != '#']
    # Get rid of line break
    for i in range(len(lines)):
      lines[i] = lines[i].replace("\n", "").replace(" ", "")
    # Load inputs
    self.__loadCode(lines)
    return

  # Analyze Bench Code
  def __loadCode(self, lines):
    counter = 1
    for code in lines:
      if code.strip() == '':
        continue
      # INPUT / OUTPUT
      elif code.casefold().find("=") == -1:
        var = code[code.index("(") + 1:code.index(")")]
        type = "INPUT" if (code.casefold().find("input") != -1) else "OUTPUT"
        if var in self.varMap:
          # Special case of a node being both INPUT and OUTPUT
          if type == "OUTPUT" and self.varMap[var][2] == "INPUT":
            newVar = var + 'o'
            self.varMap[newVar] = [
              counter, newVar, type, "_", "BUFF", [var], [], [1, 1]
            ]
            self.outputList.append(var)
            counter = counter + 1
            self.size = self.size + 1
            continue
            # End special
          self.varMap[var][2] = type
          if type == 'OUTPUT' and type not in self.outputList:
            self.outputList.append(var)
        else:
          self.varMap[var] = [counter, var, type, "_", "___", [], [], [-1, -1]]
          self.varIndex[var] = counter
          if type == 'INPUT':
            self.inputList.append(var)
            self.varMap[var][7] = [1, 1]
          else:
            self.outputList.append(var)
          counter = counter + 1
          self.size = self.size + 1
      else:
        # GATE
        var = code[:code.index("=")]
        gate = code[code.index("=") + 1:code.index("(")]
        inputs = code[code.index("(") + 1:code.index(")")].split(",")
        if var in self.varMap:
          self.varMap[var][4] = gate
          self.varMap[var][5] = inputs
        else:
          self.varMap[var] = [
            counter, var, "___", "_", gate, inputs, [], [-1, -1]
          ]
          self.varIndex[var] = counter
          counter = counter + 1
          self.size = self.size + 1
    # self.__loadCode2()
    self.__cirLevelization()
    self.__sortByLevel()
    self.__genSCOAP()
    return

  ######### Part 2 - Circuit Levelization ##########
  # For a node, determine whether
  #  all of its associate inputs have finite level
  def __nodeIsLevelable(self, node):
    return all(inp in self.nodeLevel for inp in self.varMap[node][5])

  # Find the maximum among inputs' levels
  def __maxOfInp(self, inpLst):
    return max(self.nodeLevel[inp] for inp in inpLst)

  # Create a list of sorted node by level value
  def __nodeLv(self, node):
    return self.nodeLevel[node]

  def __sortByLevel(self):
    self.sortedNode = sorted(self.nodeLevel.keys(), key=self.__nodeLv)
    self.inputList = sorted(self.inputList, key=self.__nodeLv)
    return

  # Circuit Levelization
  def __cirLevelization(self):
    isUpdate = True
    allAssigned = False
    while (isUpdate and not allAssigned):
      isUpdate = False  # Assume nothing was updated
      for line in self.varMap.values():
        node = line[1]
        if node in self.nodeLevel:
          continue  # node level is already assigned
        # node is an input
        elif line[2] == "INPUT":
          self.nodeLevel[node] = 0
          isUpdate = True
        # node whose inputs' level is finite
        elif self.__nodeIsLevelable(node):
          max = self.__maxOfInp(line[5])
          self.nodeLevel[node] = 1 + max
          isUpdate = True
          if len(self.nodeLevel) == self.size:
            allAssigned = True
            # all nodes' level are assigned, stop loop
        # If no update whatsoever, isUpdate remains false and stop loop
    return

  # SCOAP generator for each node
  def __genSCOAP(self):
    for var in self.sortedNode:
      # Input node already have (c0, c1) identified as (1, 1)
      if self.varMap[var][2] == 'INPUT':
        continue
      gate = self.varMap[var][4]
      inputs = self.varMap[var][5]
      self.varMap[var][7] = self.__operateSCOAP(var, gate, inputs)
    return



  # Function to print the expected number of bits for input and output
  # Part B.1
  def printInOutBits(self):
    print("Number of Input bits :", len(self.inputList))
    print("Number of Output bits : ", len(self.outputList))
    return

  # This Function reset all node values to '_'
  def __clearValue(self): 
    for var in self.varMap:
      self.varMap[var][3] = '_'
    return

  # Function for good circuit simulation
  # Given an input vector
  # Return/print the expected values at outputs
  # Part B.2
  def circuitSimulation(self, inputVector):
    # Clear all values at the beginning of the simulation
    self.__clearValue()
    # Initialize inputs
    self.__initInput(inputVector)
    for var in self.sortedNode:
        if var in self.inputList:
            continue
        gate = self.varMap[var][4]
        inputs = self.varMap[var][5]
        result = self.__operate(var, gate, inputs)
        # Assign the result to the node
        self.varMap[var][3] = result
    ## Display the system detail table after simulation
    #print("System Detail Table (After Simulation):")
    #self.printCode()

    ## Display the expected values at outputs
    #print("Expected output values:")
    #for var in self.outputList:
    #    print(f"{var}: {self.varMap[var][3]}")
    return

  # This function is to initialie the inputs
  # in the system from given inputVector
  def __initInput(self, inputVector):
    for i in range(len(self.inputList)):
      node = self.inputList[i]
      value = inputVector[i]
      self.varMap[node][3] = value
    return


  # This function is for the convenience of
  # find gate output
  # nodeOut: name of node at the gate output
  # gate: name of gate
  # inputNodes: list of {nodes}
  # Important: this function assumes all inputs have valid value
  def __operate(self, nodeOut, gate, inputNodes):
    inputs = []  # list of {0, 1, D, D'}
    for node in inputNodes:
      canonical = nodeOut + '-' + node
      inputName = canonical if canonical in self.varMap else node
      inputs.append(self.varMap[inputName][3])
    if gate == 'NOT':
      return NOT(inputs)
    elif gate == 'AND':
      return AND(inputs)
    elif gate == 'OR':
      return OR(inputs)
    elif gate == 'NAND':
      return NAND(inputs)
    elif gate == 'NOR':
      return NOR(inputs)
    elif gate == 'XOR':
      return XOR(inputs)
    elif gate == 'XNOR':
      return XNOR(inputs)
    elif gate == 'BUFF':
      return BUFF(inputs)

  # This function is for the convenience of
  # find a node controlability
  # nodeOut: name of node at the gate output
  # gate: name of gate
  # inputNodes: list of {nodes}
  # Important: this function assumes all inputs have known (c0, c1)
  def __operateSCOAP(self, nodeOut, gate, inputNodes):
    inputs = []  # list of { [c0, c1] } for each input node
    for node in inputNodes:
      canonical = nodeOut + '-' + node
      inputName = canonical if canonical in self.varMap else node
      inputs.append(self.varMap[inputName][7])
    if gate == 'NOT':
      return scoapNOT(inputs)
    elif gate == 'AND':
      return scoapAND(inputs)
    elif gate == 'OR':
      return scoapOR(inputs)
    elif gate == 'NAND':
      return scoapNAND(inputs)
    elif gate == 'NOR':
      return scoapNOR(inputs)
    elif gate == 'XOR':
      return scoapXOR(inputs)
    elif gate == 'XNOR':
      return scoapXNOR(inputs)
    elif gate == 'BUFF':
      return scoapBUFF(inputs)

  # Print System Detail in nice format
  def printCode(self):
    print('\n')
    print('System Detail: ')
    print("|     |      |        |    |     |ctrl-ablt|              |")
    label = "|{:>5}|{:>6}|{:>8}|{:>4}|{:>5}|{:^9}|{:>14}|".format(
        "node",
        "level",
        "type",
        "val",
        "gate",
        "(c0,c1)",
        "input",
    )
    print(label)
    print("|-----|------|--------|----|-----|---------|--------------|")
    for key in self.sortedNode:
      line = self.varMap[key]
      # num = line[0] Not needed
      node = line[1]
      level = "inf" if node not in self.nodeLevel else self.nodeLevel[node]
      type = line[2]
      val = line[3]
      gate = line[4]
      ctrl = ("(" + ','.join([str(key) for key in line[7]]) + ")") if line[7] != [-1, -1] else "(u,u)"
      input = ','.join([str(key) for key in line[5]]) if line[5] else "___"
      # fault = ','.join([str(key) for key in line[6]]) if line[6] else "___"
      printLine = "|{:{f}>5}|{:{f}>6}|{:{f}>8}|{:{f}>4}|{:{f}>5}|{:{f}^9}|{:{f}>14}|".format(
          node, level, type, val, gate, ctrl, input, f='_')
      print(printLine)
    return

  # Print Controlability Graph in nice format
  def printMC(self, mcResult):
    print('\n')
    print('System Detail: ')
    print("|     |      |        |    |ctrl-ablt|              |")
    label = "|{:>5}|{:>6}|{:>8}|{:>4}|{:^9}|{:^14}|".format(
        "node",
        "level",
        "type",
        "gate",
        "(c0,c1)",
        "(m0, m1)",
    )
    print(label)
    print("|-----|------|--------|----|---------|--------------|")
    for var in self.outputList:
      line = self.varMap[var]
      node = line[1]
      level = "inf" if node not in self.nodeLevel else self.nodeLevel[node]
      type = line[2]
      gate = line[4]
      ctrl = ("(" + ','.join([str(key) for key in line[7]]) + ")") if line[7] != [-1, -1] else "(u,u)"
      simOutput = ("(" + ','.join([str(key) for key in mcResult[var]]) + ")")
      printLine = "|{:{f}>5}|{:{f}>6}|{:{f}>8}|{:{f}>4}|{:{f}^9}|{:{f}^14}|".format(
        node, level, type, gate, ctrl, simOutput, f='_')
      print(printLine)
    return

  def getValue(self, node):
    return self.varMap[node][3]
##### End of Class Parser ###


######### Outside Class #########
# Menu Option 4
def goodCircuitSim(prog):
  print("Good Circuit Simulation:")
  print("1. All-0 input")
  print("2. All-1 input")
  print("3. Custom input")

  sub_choice = input("Choose an option (1-3): ")
  if sub_choice == "1":
    custom_input = ['0' for _ in prog.inputList]
    prog.circuitSimulation(custom_input)
  elif sub_choice == "2":
    custom_input = ['1' for _ in prog.inputList]
    prog.circuitSimulation(custom_input)
  elif sub_choice == "3":
    custom_input = input("Enter custom input vector (e.g., 0101): ")
    prog.circuitSimulation(list(custom_input))
  else:
    print("Invalid sub-choice. Returning to Main Menu.")
  prog.printCode()
  return

# Menu Option 5
def monteCarloSim(prog):
  print("Monte-Carlo Simulation: for each of 1000 random input patterns with")
  print("equal probability, we will collect the number of times each output")
  print("node has value 0, or 1, appeared in the follwing table as (m0, m1).")

  mcResult = {}

  # seed random number generator
  # By default the random number generator uses the current system time.
  seed(24)
  while True:
    for node in prog.sortedNode:
      mcResult[node] = [0, 0]
    for _ in range(1000):
      inputVector = [str(randint(0, 1)) for _ in prog.inputList]
      prog.circuitSimulation(inputVector)
      for node in prog.outputList:
        mcResult[node][int(prog.getValue(node))] += 1  
    prog.printMC(mcResult)
    yn = input("Run MC-Sim one more time (Y/N)? ")
    if yn.casefold() == "n":
      break
  return

######## Main ########
def main():
  prog = Parser()
  fileName = input("Enter bench file name: ")
  prog.readFile(fileName)
  prog.printCode()

  while True:
    print("\nMenu:")
    print("1. Print System Detail Table")
    print("2. Circuit Analysis")
    print("3. Good Circuit Simulation")
    print("4. Monte-Carlo Simulation")
    print("5. Exit")

    choice = input("Enter your choice (1-7): ")
    print("\n")

    if choice == "1":
      prog.printCode()
    elif choice == "2":
      prog.printInOutBits()
    elif choice == "3":
      goodCircuitSim(prog)
    elif choice == "4":
      monteCarloSim(prog)
    elif choice == "5":
      print("Exiting the program.")
      break
    else:
      print("Invalid choice. Please enter a valid option.")


main()

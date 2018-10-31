import sys
import copy

########################### NODE CLASS ################################
class Node:
    
    def __init__(self, parent, state, hValue = 0, gValue = 0):
        self.parent = parent
        self.state = state
        self.hValue = hValue
        self.gValue = gValue     
        
    def isEqual(self, otherNode):
        if self.state == otherNode.state:
            return True
        else:
            return False
            
    def changeParent(self,otherNode):
        self.parentName = otherNode.parentName[:]
        self.hValue = otherNode.hValue
        self.gValue = otherNode.gValue
            
    def printNode(self):
        i = 0 
        while i<len(self.state):
            print ' '.join(self.state[i])
            i = i+1
        print ""
            
############################ TREE CLASS #####################################

# Every puzzle is represented by a Tree
class Tree:

    def __init__(self, root, goalStates, typeOfH, numberOfBlocks, rows, cols):
        # startState is a Root object
        self.root = root
        # OPEN is a set of Node
        self.OPEN = []
        self.OPEN.append(root)
        # CLOSED is a set of Node
        self.CLOSED = []
        # goalStates is a list of goal states
        self.goalStates = goalStates[:]
        # type of heuristic function that will be used during processing this tree
        self.typeOfH = int(typeOfH)
        # number of blocks
        self.numberOfBlocks = int(numberOfBlocks)
        # number of rows
        self.rows = int(rows)
        # number of columns
        self.cols = int(cols)

    def isGoal(self, node):
        i = 0
        while i < len(self.goalStates):
            if node.state == self.goalStates[i]:
                return True
            i = i+1
        return False
        
    def isInOPEN(self,node) :
        i = 0
        while i < len(self.OPEN):
            if node.state == self.OPEN[i].state:
                return True
            i = i+1
        return False

    def isInCLOSED(self,node) :
        i = 0
        while i < len(self.CLOSED):
            if node.state == self.CLOSED[i].state:
                return True
            i = i+1
        return False
        
    def getFromOPEN(self,node):
        i = 0
        while i < len(self.OPEN):
            if node.state == self.OPEN[i].state:
                return self.OPEN[i]
            i = i+1
            
    def getFromCLOSED(self,node):
        i = 0
        while i < len(self.CLOSED):
            if node.state == self.CLOSED[i].state:
                return self.CLOSED[i]
            i = i+1
        
    def printTree(self):
        root.printRoot()
        
    def printOPEN(self):
        if len(self.OPEN) == 0:
            print "No node in OPEN"
        i = 0
        while i < len(self.OPEN) : 
            self.OPEN[i].printNode()
            print " "
            i = i+1
            
    def printPath(self,node):
        if node.parent == "root":
            node.printNode()
            
        else :
            self.printPath(node.parent)
            node.printNode()
    
    
############### MOVE FUNCTIONS ############### 

    def moveRight(self, blockID, state) :
        
        indexesToChange = findIndexes(blockID, state, self.rows, self.cols)
        result = copy.deepcopy(state)
        i = 0
        while i  < len(indexesToChange) :
            x = indexesToChange[i][0]
            y = indexesToChange[i][1]
            result[x][y+1] = str(blockID)
            result[x][y] = '0' 
            i = i+1
        return copy.deepcopy(result)
        
    def moveLeft(self, blockID, state):
        
        indexesToChange = findIndexes(blockID, state, self.rows, self.cols)
        result = copy.deepcopy(state)
        i = 0
        while i  < len(indexesToChange) :
            x = indexesToChange[i][0]
            y = indexesToChange[i][1]
            result[x][y-1] = str(blockID)
            result[x][y] = '0' 
            i = i+1
        return copy.deepcopy(result)
                
    def moveUp(self,blockID, state):
        indexesToChange = findIndexes(blockID, state, self.rows, self.cols)
        result = copy.deepcopy(state)
        i = 0
        while i < len(indexesToChange) : 
            x = indexesToChange[i][0]
            y = indexesToChange[i][1]
            result[x-1][y] = str(blockID)
            result[x][y] = '0'
            i = i+1
        return copy.deepcopy(result)
    
    def moveDown(self,blockID,state):
        indexesToChange = findIndexes(blockID, state, self.rows, self.cols)
        result = copy.deepcopy(state)
        i = 0
        while i < len(indexesToChange) : 
            x = indexesToChange[i][0]
            y = indexesToChange[i][1]
            result[x+1][y] = str(blockID)
            result[x][y] = '0'
            i = i+1
        return copy.deepcopy(result)
    
############## canGo FUNCTIONS ###############
    def canGoUp(self, indexesOfComponents,state):
        # take the indexes at the uppermost part
        indexesAtTop = getUppermost(indexesOfComponents,self.rows)

        if indexesAtTop[0][0] <= 0 :
            return False
        else :
            i = 0
            while i < len(indexesAtTop):
                x = indexesAtTop[i][0]
                y = indexesAtTop[i][1]
                
                # there is an object at top, return False
                if int(state[x-1][y]) != 0 :
                    return False
                i = i+1
            # there is no obstacle to go up , return True
        return True
        
    def canGoDown(self, indexesOfComponents,state):
        # take the indexes at the undermost part
        indexesAtBottom = getUndermost(indexesOfComponents)

        # if the block at the bottom already, return False
        if indexesAtBottom[0][0] >= (self.rows-1):
            return False
        else :
            i=0
            while i < len(indexesAtBottom) :
                x = indexesAtBottom[i][0]
                y = indexesAtBottom[i][1]
                # there is another object under , return False
                if int(state[x+1][y]) != 0:
                    return False
                i = i+1
            # there is no obstacle to go down, return True
        return True

    
    def canGoLeft(self, indexesOfComponents, state):
    	# take the indexes at the leftmost 
        indexesAtLeft = getLeftmost(indexesOfComponents, self.cols)
    	# if the block at the leftmost place already, return False
    	if indexesAtLeft[0][1] <= 0 : 
            return False
    	else :
            i = 0
            while i < len(indexesAtLeft) :
                x = indexesAtLeft[i][0]
                y = indexesAtLeft[i][1]
    			# if there is another block at left , return False
                if state[x][y-1] != '0' :
                    return False			
                i = i+1
    		# there is no obstacle to go left, return True
        return True
    
    
    def canGoRight(self, indexesOfComponents, state):
        
        # take the indexes at the rightmost 
        indexesAtRight = getRightmost(indexesOfComponents)
        # if the block is at the rightmost place already, return False
        if indexesAtRight[0][1]  >=  int(self.cols) -1 :
            return False
        else : 
            i = 0
            while i  < len(indexesAtRight):
                x = indexesAtRight[i][0]
                y = indexesAtRight[i][1]
                # if there is another block at right, return False
                if int(state[x][y+1])  !=  0 :
                    return False
                i = i+1
        # there is no obstacle for going right, return True
        return True
       
############################ CORE FUNCTIONS ##################################
        
    # extends the given node and creates a list of extended nodes
    def extend(self,node):
        blockID = 1
        result = []
        while blockID <= self.numberOfBlocks:
            
            indexesOfBlockA = findIndexes(blockID,node.state,self.rows,self.cols)
            
            goRight = self.canGoRight(indexesOfBlockA , node.state) 
            goLeft = self.canGoLeft(indexesOfBlockA, node.state)
            goUp = self.canGoUp(indexesOfBlockA,node.state)
            goDown = self.canGoDown(indexesOfBlockA,node.state)
           
            # use Manhattan Distance as heuristic
            if self.typeOfH == 0 :
                
                if goRight == True : 
                    result.append(Node(node,copy.deepcopy(self.moveRight(blockID,node.state)) , manhattanDistance(self.goalStates, copy.deepcopy(self.moveRight(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
                if goLeft == True :
                    result.append(Node(node,copy.deepcopy(self.moveLeft(blockID,node.state)) , manhattanDistance(self.goalStates, copy.deepcopy(self.moveLeft(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
                if goUp == True :
                    result.append(Node(node,copy.deepcopy(self.moveUp(blockID,node.state)) , manhattanDistance(self.goalStates, copy.deepcopy(self.moveUp(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
                if goDown == True :
                    result.append(Node(node,copy.deepcopy(self.moveDown(blockID,node.state)) , manhattanDistance(self.goalStates, copy.deepcopy(self.moveDown(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
            
            # use my Heuristic as heuristic
            else :
                
                if goRight == True : 
                    result.append(Node(node,copy.deepcopy(self.moveRight(blockID,node.state)) , myHeuristic(self.goalStates, copy.deepcopy(self.moveRight(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
                if goLeft == True :
                    result.append(Node(node,copy.deepcopy(self.moveLeft(blockID,node.state)) , myHeuristic(self.goalStates, copy.deepcopy(self.moveLeft(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
                if goUp == True :
                    result.append(Node(node,copy.deepcopy(self.moveUp(blockID,node.state)) , myHeuristic(self.goalStates, copy.deepcopy(self.moveUp(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
                if goDown == True :
                   result.append(Node(node,copy.deepcopy(self.moveDown(blockID,node.state)) , myHeuristic(self.goalStates, copy.deepcopy(self.moveDown(blockID,node.state)), self.numberOfBlocks, self.rows, self.cols)  , (node.gValue +1 )))
            
            # loop variable
            blockID = blockID+1
        
        return copy.deepcopy(result)
    
    def findLeastValueInOPEN(self):
        i = 0
        minimumValue = 99999
        while i < len(self.OPEN):
            if minimumValue > (self.OPEN[i].hValue + self.OPEN[i].gValue):
                minimumValue = self.OPEN[i].hValue + self.OPEN[i].gValue
            i=i+1
        j=0
        while j < len(self.OPEN) :
            if minimumValue == (self.OPEN[ j ].hValue + self.OPEN[ j ].gValue) :
                return self.OPEN[ j ]
            j = j+1  
    
    def process(self):
        ## START A* HERE
        
        ## loop over OPEN list until it is empty
        while len(self.OPEN) > 0 :
            
            ## search a node N with least h value in OPEN  list, returns reference
            N = self.findLeastValueInOPEN()
            
            ##  if N is a goal state return the path 
            if self.isGoal(N) == True:
                self.printPath(N)
                return
            else:
                ## delete N from OPEN , since it is a reference we can directlr remove it 
                self.OPEN.remove(N)

                self.CLOSED.append(N)
                
                # extend N 
                nodeList = self.extend(N)                      
                i = 0
                # for each node in nodeList
                while i < len(nodeList) :
                   #  if node.state  is not in OPEN and is not in CLOSED, insert it into OPEN
                    if self.isInOPEN(nodeList[i]) == False and self.isInCLOSED(nodeList[i]) == False :
                        self.OPEN.append(nodeList[i])
           
            ## else if , node.state is in OPEN but its gValue is less than the one's in the OPEN, change parent of the one in the OPEN list
                    elif self.isInOPEN(nodeList[i]) == True :
                        # returns reference
                        o = self.getFromOPEN(nodeList[i])
                        if o.gValue > nodeList[i].gValue:
                            o.gValue = nodeList[i].gValue
                            o.parent = nodeList[i].parent
            
            ## else if , node.state is in CLOSED, and gValue is less than the one's in the CLOSED, change parent of it , delete it from CLOSED, append it to OPEN   
                    elif self.isInCLOSED(nodeList[i]) == True:
                        # returns reference
                        c = self.getFromCLOSED(nodeList[i])
                        if c.gValue > nodeList[i].gValue:
                            c.gValue = nodeList[i].gValue
                            c.parent = nodeList[i].parent
                            self.CLOSED.remove(c)
                            self.OPEN.append(c) 
                    i=i+1
        
        print "There is no solution"
        
        
################# INDEX FUNCTIONS ###########################

# returns rightmost indexes in indexes list
def getRightmost(indexes):
    i = 0
    maximum = -1
    result = []
    while i < len(indexes):
        if indexes[i][1] > maximum:
            maximum = indexes[i][1]
        i = i+1
    j = 0
    while j < len(indexes):
        if indexes[j][1] == maximum:
            result.append(indexes[j])
        j = j+1
    return result[:]

# returns the leftmost indexes in indexes list
def getLeftmost(indexes,cols):
    i = 0
    minimum = cols+1
    result = []
    while i < len(indexes):
        if indexes[i][1] < minimum:
            minimum = indexes[i][1]
        i = i+1
    j = 0
    while j < len(indexes):
        if indexes[j][1] == minimum:
            result.append(indexes[j])
        j = j+1
    return result[:]

def getUppermost(indexes,rows) : 
    i = 0
    result = []
    minimum = rows+1
    while i < len(indexes) :
        if indexes[i][0] < minimum:
            minimum = indexes[i][0]
        i = i+1
    j=0
    while j < len(indexes):
        if indexes[j][0] == minimum:
            result.append(indexes[j])
        j = j+1
    return result[:]

def getUndermost(indexes) : 
    maximum = 0
    i = 0
    result = []
    while i < len(indexes) :
        if indexes[i][0] > maximum:
            maximum =  indexes[i][0]
        i = i+1
    j = 0
    while j < len(indexes):
        if indexes[j][0] == maximum:
            result.append(indexes[j])
        j = j + 1
             
    return result[:]


# returns indexes of given block in the state as a tuple form
def findIndexes(blockID, state, rowNumber, colNumber):
    result = []
    i=0
    while i < int(rowNumber) :
        j = 0
        while j < int(colNumber) :
            if int(state[i][j]) == int(blockID) :
                result.append((i,j))
            j = j+1
        i = i+1    
    return result

########################## HEURISTICS ############################

# returns the minimum manhattan distance for given block 
def manhattanDistanceOfOneBlock( goalStates, currentState, blockID,rows,cols):
    distances = []
    (a,b) = findIndexes(blockID,currentState,rows,cols)[0]
    j = 0
    while j < len(goalStates):
        (x,y) = findIndexes(blockID,goalStates[j],rows,cols)[0]
        distances.append(abs(a-x)+abs(b-y))
        j = j + 1
    return min(distances)

# returns the heuristic value of current state
def manhattanDistance(goalStates, currentState, pieces, rows,cols):
    i = 1
    result = 0
    while i <= int(pieces):
        result = result +  manhattanDistanceOfOneBlock(goalStates,currentState,i,rows,cols)
        i = i+1
    return result

# helper for myHeuristic
def isInRightPlace(blockID,currentState,goalState, rows,cols) :
    indexesInCurrent = findIndexes(blockID,currentState,rows,cols)
    indexesInGoal = findIndexes(blockID,goalState,rows,cols)    
    if indexesInGoal == indexesInCurrent :
        return True
    else : 
        return False 

#  implementation of "number of misplaced tiles" 
def myHeuristic(goalStates, currentState,pieces,rows,cols):
    result = []
    i = 0
    while i < len(goalStates) :
        temp = 0
        j = 1
        while j <= pieces:
             if isInRightPlace(j,currentState,goalStates[i],rows,cols) == False :
                temp = temp+1            
             j = j+1
        result.append(temp)     
        i = i+1
    return min(result)
        
    
################################# MAIN CALLS ###################################

file_name = sys.argv[1]

read_file = open(file_name,"r")

numOfTrees = (read_file.readline())[:-1:]

t = 0
while t < int(numOfTrees):

    # wait for enter 
    while raw_input() == "\n" :
        print "waiting for enter"
    
    typeOfHeuristic = (read_file.readline())[:-1:]
    
    inputInfo = (read_file.readline())[:-1:]
    rows = inputInfo[0]
    columns = inputInfo[2]
    pieces = inputInfo[4]
    numOfGoalStates = inputInfo[6]

    # to read S symbol and startingState
    symbol = read_file.readline()

    if symbol[:-1:] == 'S':
        startingState = []
        n = 0
        while n < int(rows):
            temp = read_file.readline()
            startingState.append(list(temp[:-1:2]))
            n = n + 1
            
    # to read F and goalStates, there can be more than one goalState
    g = 0
    goalStates = []
    while g < int(numOfGoalStates):
        symbol = read_file.readline()
        if symbol[:-1:] == 'F':
            goalState = []
            m = 0
            while m < int(rows):
                temp = read_file.readline()
                goalState.append(list(temp[:-1:2]))
                m = m + 1
            goalStates.append(goalState)
        g = g + 1
        
    root = Node('root', copy.deepcopy(startingState))
    tree = Tree(copy.deepcopy(root), copy.deepcopy(goalStates), copy.deepcopy(typeOfHeuristic), copy.deepcopy(pieces), rows, columns)
    tree.process()
    t = t + 1



        
    

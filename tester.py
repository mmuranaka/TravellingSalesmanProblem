# python tester.py CAMPUS.CSV 1 1000 0.01
# python tester.py CAMPUS.CSV 2 1000 0.01
import sys
import csv
import math
import random
import copy
import timeit
import matplotlib.pyplot as plt
# a Point object is a state. it contains info on the label, x, and y coordinate
class Point:
    def __init__(self, label, x, y):
        self.label = label
        self.x = float(x)
        self.y = float(y)

    def __str__(self):
        out = "Label: " + str(self.label) + ", x: " + str(self.x) + ", y: " + str(self.y)
        return out
    
    def getLabel(self):
        return self.label
    
    def distance(self, other):
        differenceX = self.x-other.x
        differenceY = self.y-other.y
        distance = math.sqrt((differenceX**2) + (differenceY**2))
        return distance
    
    def changeX(self, newX):
        self.x = newX    
        
# a Space object is a list of Point objects that include all of the states in the csv
class Space:
    def __init__(self, list):
        self.statespace = list
    
    def __str__(self):
        out = ""
        for i in self.statespace:
            out += str(i) + "\n"
        return out
    
    def __len__(self):
        return len(self.statespace)+1
    
    def route(self, initial):
        out = [[initial.label]]
        for i in self.statespace:
            out.append([i.label])
            # print(i.label)
        return out
    
    def printPath(self, initial):
        out = initial.getLabel() + ", "
        for i in range(0, len(self.statespace)):
            if i==len(self.statespace)-1:
                out += self.statespace[i].getLabel()
            else:
                out += self.statespace[i].getLabel() + ", "
        return out
    
    def shuffle(self):
        new = copy.deepcopy(self)                                                   # we use a deep copy becasue we dont want to lose our current config
        random.shuffle(new.statespace)
        return new

#performs a 2-edge swap
    def twoEdgeSwap(self):
        point1 = random.randint(0, len(self.statespace)-1)
        point2 = random.randint(0, len(self.statespace)-1)
        while point1==point2:
            point2 = random.randint(0, len(self.statespace)-1)
        temp = self.statespace[point1]
        self.statespace[point1] = self.statespace[point2]
        self.statespace[point2] = temp
        return self



    def totalDistance(self, initial):
        total = 0
        total += initial.distance(self.statespace[0])                               #distance from initial to next     
        for index in range(0, len(self.statespace)-1):
            total += self.statespace[index].distance(self.statespace[index+1])
        total += initial.distance(self.statespace[len(self.statespace)-1])
        return total

# computes two new states (paths) where p1 and p2 are the crossover points
    def crossover(self, other):
        point1 = random.randint(0, len(self.statespace)-2)
        point2 = random.randint(point1+1, len(self.statespace)-1)
        return Space(self.statespace[0:point1] + other.statespace[point1:point2] + self.statespace[point2:len(self.statespace)])

    
# we implement the inversion mutation mechanism
    def mutate(self, probability):
        result = copy.copy(self)
        if random.random()<probability:
            point1 = random.randint(0, len(self.statespace)-2)
            point2 = random.randint(point1+1, len(self.statespace)-1)
            result.statespace[point1:point2] = result.statespace[point1:point2][::-1]
        return result


    def validRoute(self):
        myset = set(self.statespace)
        if (len(myset)==len(self.statespace)):
            return True
        else:
            return False
        

# states is a Space object, initial is a Point object
# initial is always first and last
def SimulatedAnnealing(states, initial, temperature, cooling):
    currentTemp = temperature
    currentSolution = copy.copy(states)
    newSolution = copy.copy(currentSolution)
    iteration = 1
    totalDistances = []
    while currentTemp>0:
        newSolution = newSolution.twoEdgeSwap()                                                          #perform two-edge swap
        if currentSolution.totalDistance(initial) > newSolution.totalDistance(initial):
# 100% acceptance rate
            currentSolution = copy.copy(newSolution)
        else:
# worse solution
            negchange = currentSolution.totalDistance(initial) - newSolution.totalDistance(initial)
            probability = math.exp(negchange/currentTemp)
            decision = random.random()
            if decision<probability:
                currentSolution=copy.copy(newSolution)
        total_distance = currentSolution.totalDistance(initial)
        totalDistances.append(total_distance)
        currentTemp = currentTemp * (math.exp(-1 * iteration * cooling))
        iteration +=1
    return (currentSolution, iteration)

# representation:   each state along the path is a chromosome 
# fitness function: cost of all paths
# population pool: 10?
# selection mechanism: roulette wheel
#       since higher fitness, worse solution, it is fitness ratios below probability

# pool is the size of the initial population
# two point crossover, 2/10 of population are elitists
#       half of elite population will be untouched (directly moved to next generation)
#       other half will undergo mutation
# other 8/10 of population will be considered for selection
# Mutation: two positions are randomly selected states between those positions are inverted
def GeneticAlgorithm(states, initial, iterations, mutation, pool=10):
    eliteAmount = pool//10                      # one tenth of population
    initialSolution = copy.copy(states)
    population = [[initialSolution, initialSolution.totalDistance(initial)]]
    # population is tuple of routes and their cost (and eventually fitness)
    # (route, cost, fitness)
    sum = 0                                                               # for calculating fitness
    for i in range(0, pool-1):
        new = states.shuffle()
        population.append([new, new.totalDistance(initial)])                                            # initialize population with 10 paths
        sum += new.totalDistance(initial)                                 # used for fitness ratio
# at this point we have pool amount of inidividuals in population list
    for i in population:
        fitnessRatio = i[1]/sum
        # print(fitnessRatio)
        i.append(fitnessRatio)
# at this point we have a pool amount of different individuals and their percentage of being selected (smaller percentage better)
    population.sort(key=lambda list: list[1])                 #POSSIBLY CHANGE TO COPY 
    previousGeneration = copy.deepcopy(population)   
    # print(type(previousGeneration[0][0]))
    # print(str(previousGeneration[0]))

    for _ in range(0, iterations):
        children = []
        sum = 0
        for i in range(0, eliteAmount):                                     # these are kept the same
            # print(len(previousGeneration))
            child = previousGeneration.pop(0)
            children.append([child[0], child[0].totalDistance(initial)])
            sum += child[0].totalDistance(initial)
            # print(type(child))
        for i in range(eliteAmount, eliteAmount):                         # these are only mutated
            # print(len(previousGeneration))
            child = previousGeneration.pop(0)
            child = child[0].mutate(mutation)
            children.append(child)
            # print(type(child))
            sum += child.totalDistance(initial)
        selected = []                                                   # contains Space objects
        #Selection process (Roulette Wheel)
        decision = random.random()
        for i in previousGeneration:
            if decision>i[2]:                   #MIGHT BE WRONG
                selected.append(i[0])
            else:
                break
        # perform crossover
        for i in range(0, len(selected)):
            if (i == len(selected)-1):
                child = selected[i].crossover(selected[0])
            else:
                child = selected[i].crossover(selected[i+1])
            child = child.mutate(mutation)
            if (child.validRoute()):
                children.append([child, child.totalDistance(initial)])
                sum += child.totalDistance(initial)
            else:
                children.append([child, child.totalDistance(initial)**2])      #large emphasis against valid routes
                sum += child[0].totalDistance(initial)
        # at this point we have all of our children
        # I gotta figure out fitness and sort in preparation for next iteration
        for i in population:
            i.append(i[1]/sum)
        population.sort(key=lambda list: list[1])
        previousGeneration = copy.deepcopy(population)
        # print(previousGeneration)
        # print(type(previousGeneration))
    return (previousGeneration[0][0])

def runTests(filename):
    try:
        with open(filename, mode ='r')as file:
            csvFile = csv.reader(file)
            states = []
            lineNumber = 0
            for index, i in enumerate(csvFile):
                if index==0:
                    initial = Point(i[0], i[1], i[2])
                else:
                    states.append(Point(i[0], i[1], i[2]))
    except FileNotFoundError:
        print("ERROR: Not enough/too many/illegal input arguments.")
        exit()
    
    
    header = ['Min Path cost','Max Path cost', 'Average Path cost', 'Min search time in seconds', 'Max search time in seconds', 'Average search time in seconds']

    print("ENTERING TESTING MODE")
    # print("Enter Algo to Test: ")
    # algo=str(input())
    stateSpace = Space(states)
    attempts = 5
    p1 = 10
    p2 = 0.1
    csvrows = []
    for algo in range(1,2):
        for k in range(1, 2):
            for j in range(2,3):
                times = []
                costs = []
                for _ in range(0, attempts):
                    timeStart = timeit.default_timer()
                    if (algo == 0):
                        (solution, _) = SimulatedAnnealing(stateSpace, initial, p1**j, p2**k)
                    elif (algo==1):
                        solution = GeneticAlgorithm(stateSpace, initial, int(p1**j), p2**k)
                    else:
                        exit()
                    timeEnd = timeit.default_timer()
                    elapsedTimeInSec = timeEnd - timeStart
                    times.append(elapsedTimeInSec)
                    costs.append(solution.totalDistance(initial))
                iterations = range(1, 5)
                plt.plot(iterations, costs)
                plt.xlabel('Iteration Number')
                plt.ylabel('Distance Traveled')
                plt.title('Simulated Annealing Trials')
                plt.show()

                times.sort()
                averageTimes = 0
                averageCosts = 0
                print('p1: ' + str(k) + '    p2: ' + str(j))
                for i in range(0, 5):
                    averageTimes += times[i]
                    averageCosts += costs[i]
                averageTimes = averageTimes/attempts
                averageCosts = averageCosts/attempts
                csvrows.append([costs[0], costs[4], averageCosts, times[0], times[4], averageTimes])
        print('finished SA')
    with open("algotests.csv", 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(header)
        csvwriter.writerows(csvrows)
    exit()
        
        
            
    
    



if (len(sys.argv) != 5):
    if (sys.argv[2] == '3'):
        runTests(sys.argv[1])
    else: 
        print("ERROR: Not enough or too many input arguments.")
        exit()
filename = sys.argv[1]
    # csv file
algo = sys.argv[2]
    # 1 - Simulated Annealing
    # 2 - Genetic Algorithm
p1 = float(sys.argv[3])
    # temperature/iterations
p2 = float(sys.argv[4])
    # cooling/mutation probability
try:
    with open(filename, mode ='r')as file:
        csvFile = csv.reader(file)
        states = []
        lineNumber = 0
        for index, i in enumerate(csvFile):
            if index==0:
                initial = Point(i[0], i[1], i[2])
            else:
                states.append(Point(i[0], i[1], i[2]))
except FileNotFoundError:
    print("ERROR: Not enough/too many/illegal input arguments.")
    exit()
# we made this into a Space object for some ease of calculations
# print(str(stateSpace))
if algo=='1':
    print("Initial state: " + initial.getLabel())
    print("\nSimulated Annealing:")
    print("Command Line Parameters:    Temperature: " + str(p1) + ", Cooling: " + str(p2))
    # initial is the assumed starting state (made into Point object
    # states is an array containing Point objects of all the states (except initial) in the given csv
    stateSpace = Space(states)
    # we made this into a Space object for some ease of calculations
    print("Initial Solution: " + stateSpace.printPath(initial))
    timeStart = timeit.default_timer()
    (solution, iterations) = SimulatedAnnealing(stateSpace, initial, p1, p2)
    # print(type(solution))
    timeEnd = timeit.default_timer()
    elapsedTimeInSec = timeEnd - timeStart
    print("Final Solution: " + solution.printPath(initial))
    print("Number of Iterations: " + str(iterations))
    print("Execution Time: " + str(elapsedTimeInSec))
    print("Complete Path Cost: " + str(solution.totalDistance(initial)))

    filewrite = filename[:-4]
    filewrite += "_SOLUTION_SA.csv"
    with open(filewrite, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([str(solution.totalDistance(initial))])
        path = solution.route(initial)
        csvwriter.writerows(path)

elif algo=='2':
    print("Initial state: " + initial.getLabel())
    print("\nGenetic Algorithm:")
    print("Command Line Parameters:    Iterations: " + str(p1) + ", Mutation Probability: " + str(p2))
    # initial is the assumed starting state (made into Point object
    # states is an array containing Point objects of all the states (except initial) in the given csv
    stateSpace = Space(states)
    # we made this into a Space object for some ease of calculations
    print("Initial Solution: " + stateSpace.printPath(initial))
    timeStart = timeit.default_timer()
    solution = GeneticAlgorithm(stateSpace, initial, int(p1), p2)
    timeEnd = timeit.default_timer()
    elapsedTimeInSec = timeEnd - timeStart
    print("Final Solution: " + solution.printPath(initial))
    print("Number of Iterations: " + str(p1))
    print("Execution Time: " + str(elapsedTimeInSec))
    print("Complete Path Cost: " + str(solution.totalDistance(initial)))

    filewrite = filename[:-4]
    filewrite += "_SOLUTION_GA.csv"
    with open(filewrite, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([str(solution.totalDistance(initial))])
        path = solution.route(initial)
        csvwriter.writerows(path)
# elif algo=='3':
    # header = ['p1','p2','Min Path cost','Max Path cost', 'Average Path cost', 'Min search time in seconds', 'Max search time in seconds', 'Average search time in seconds']

    # print("ENTERING TESTING MODE")
    # print("Enter Algo to Test: ")
    # algo=str(input())
    # stateSpace = Space(states)
    # attempts = 5
    # p1 = 10
    # p2 = 0.1
    # csvrows = []
    # for i in range(1, 4):
    #     for j in range(2,5):
    #         times = []
    #         costs = []
    #         for _ in range(0, attempts):
    #             timeStart = timeit.default_timer()
    #             if (algo == '1'):
    #                 (solution, _) = SimulatedAnnealing(stateSpace, initial, p1**j, p2**i)
    #             elif (algo=='2'):
    #                 solution = GeneticAlgorithm(stateSpace, initial, int(p1**j), p2**i)
    #             else:
    #                 exit()
    #             timeEnd = timeit.default_timer()
    #             elapsedTimeInSec = timeEnd - timeStart
    #             times.append(elapsedTimeInSec)
    #             costs.append(solution.totalDistance(initial))
    #         times.sort()
    #         averageTimes = 0
    #         averageCosts = 0
    #         for i in range(0, 5):
    #             averageTimes += times[i]
    #             averageCosts += costs[i]
    #         averageTimes = averageTimes/attempts
    #         averageCosts = averageCosts/attempts
    #         csvrows.append([p1**j, p2**i, costs[0], costs[4], averageCosts, times[0], times[4], averageTimes])
    # with open("algotests.csv", 'w') as csvfile:
    #     csvwriter = csv.writer(csvfile)
    #     csvwriter.writerow(header)
    #     csvwriter.writerows(csvrows)

    
    
    
    # print("Times: Min: " + str(times[0]) + "    Max: " + str(times[4]) + "    Average:  " + str(averageTimes))
    # costs.sort()
    # averageCosts = 0
    # for i in costs:
    #     averageCosts +=i
    # averageCosts = averageCosts/5
    # print("Costs: Min: " + str(costs[0]) + "    Max: " + str(costs[4]) + "    Average:  " + str(averageCosts))
        


else:
    print("ERROR: Invalid algorithm argument.")
    exit()



# testInitial = Point("it", 7,6)
# test = [Point("sb", 0, 3), Point("mt", 7.2, 3.4), Point("qt", 4, 4.5)]
# testSpace = Space(test)
# print(str(testSpace))
# print("test total distance: ")
# print(testSpace.totalDistance(testInitial))


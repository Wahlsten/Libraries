import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from NeuralNetwork_Object import NeuralNetworkClass

def CreateDomain(type):
    """
    Creates a domain
    """

    if type == 'Circle':
        def upper_left(x):  return np.sqrt(np.abs(11*11 - x**2))
        def upper_right(x): return np.sqrt(np.abs(8*8 - x**2))
        def lower_left(x):  return - np.sqrt(np.abs(11*11 - x**2))
        def lower_right(x): return - np.sqrt(np.abs(8*8 - x**2))
        domain = [[-11, 11, upper_left], [-8, 8, upper_right], [-11, 11, lower_left], [-8, 8, lower_right]]

    return domain

def PlotDomain(domain, points, angle):
    """
    Plots the domain
    """

    global domain_g
    global x_g
    global y_g
    global angle_g

    domain_g = domain
    x_g = points[0, :]
    y_g = points[1, :]
    angle_g = angle[0]

    def getDomain():
        global domain_g

        numPoints = 1000
        x = []
        y = []

        for k in range(0, np.shape(domain_g)[0]):
            x.append(np.linspace(domain_g[k][0], domain_g[k][1], numPoints))
            y.append(domain[k][2](x[k]))

        #plt.plot(x[0], y[0], 'b', x[1], y[1], 'b', x[2], y[2], 'b', x[3], y[3], 'b', points[0, :], points[1, :], 'r')
        '''
        visionpoint_right = [[points[0][0], vision[0][0]], [points[1][0], vision[0][1]]]
        visionpoint_left  = [[points[0][0], vision[1][0]], [points[1][0], vision[1][1]]]
        visionpoint_front = [[points[0][0], vision[2][0]], [points[1][0], vision[2][1]]]
        plt.plot(x[0], y[0], 'b', x[1], y[1], 'b', x[2], y[2], 'b', x[3], y[3], 'b', points[0, :], points[1, :], 'r', visionpoint_right[0], visionpoint_right[1], 'g', visionpoint_left[0], visionpoint_left[1], 'g', visionpoint_front[0], visionpoint_front[1], 'g')
        '''
        #plt.show()

        return x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3]

    def drawCar(x, y, angle):
        '''
        Draw a car with center at (x, y)
        '''

        #print('----- Draw Car -----')
        car_height = 2.
        car_width = 1.

        front_left  = np.array([ car_height / 2.,  car_width / 2.])
        front_right = np.array([ car_height / 2., -car_width / 2.])
        rear_left   = np.array([-car_height / 2.,  car_width / 2.])
        rear_right  = np.array([-car_height / 2., -car_width / 2.])
        rotation_matrix = np.matrix([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

        rot_front_left  = np.dot(rotation_matrix, front_left)  + np.array([x, y])
        rot_front_right = np.dot(rotation_matrix, front_right) + np.array([x, y])
        rot_rear_left   = np.dot(rotation_matrix, rear_left)   + np.array([x, y])
        rot_rear_right = np.dot(rotation_matrix, rear_right)  + np.array([x, y])

        car_x = [rot_front_left[0,0], rot_front_right[0,0], rot_rear_right[0,0], rot_rear_left[0,0], rot_front_left[0,0]]
        car_y = [rot_front_left[0,1], rot_front_right[0,1], rot_rear_right[0,1], rot_rear_left[0,1], rot_front_left[0,1]]

        return car_x, car_y

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(-20, 20), ylim=(-20, 20))
    line2, = ax.plot([], [], 'b-', lw=2)
    line3, = ax.plot([], [], 'b-', lw=2)
    line4, = ax.plot([], [], 'b-', lw=2)
    line5, = ax.plot([], [], 'b-', lw=2)
    line6, = ax.plot([], [], 'k-', lw=2)

    # initialization function: plot the background of each frame
    def init():
        #line6.set_data(x0, y0)
        line2.set_data([], [])
        line3.set_data([], [])
        line4.set_data([], [])
        line5.set_data([], [])
        line6.set_data([], [])
        return line2, line3, line4, line5, line6

    # animation function. This is called sequentially
    def animate(i):
        global x_g
        global y_g
        global angle_g

        x0, y0, x1, y1, x2, y2, x3, y3 = getDomain()
        line2.set_data(x0, y0)
        line3.set_data(x1, y1)
        line4.set_data(x2, y2)
        line5.set_data(x3, y3)
        #theta_1, theta_2 = PendulumPosition(i)
        car_x, car_y = drawCar(x_g[i], y_g[i], angle_g[i])
        #x = np.linspace(0, 2, 1000)
        #y = np.sin(2 * np.pi * (x - 0.01 * i))

        #line2.set_data(0,0)
        #line3.set_data(0,0)
        #line4.set_data(0,0)
        #line5.set_data(0,0)
        line6.set_data(car_x, car_y)
        return line2, line3, line4, line5, line6,

    #print('plot')
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=100, interval=20, blit=True)
    plt.show()
def ComputeDistance(points, rotation, domain):
    """
    Compute the distances to the walls
    """
    angle = np.pi/180. * 30.
    rotation_right = rotation + angle
    rotation_left = rotation - angle
    rotation_front = rotation
    vec_right = np.array([np.cos(rotation_right), np.sin(rotation_right)])
    vec_left = np.array([np.cos(rotation_left), np.sin(rotation_left)])
    vec_front = np.array([np.cos(rotation_front), np.sin(rotation_front)])

    len_right_low = 0
    len_right_high = 10

    len_left_low = 0
    len_left_high = 10

    len_front_low = 0
    len_front_high = 10

    for k in range(0, 10):

        len_right_mid = (len_right_low + len_right_high) / 2.
        len_left_mid = (len_left_low + len_left_high) / 2.
        len_front_mid = (len_front_low + len_front_high) / 2.

        point_right = points + len_right_mid * vec_right
        point_left  = points + len_left_mid * vec_left
        point_front = points + len_front_mid * vec_front

        outOfBounds_right = CheckOutOfBounds(domain, point_right)
        outOfBounds_left  = CheckOutOfBounds(domain, point_left)
        outOfBounds_front = CheckOutOfBounds(domain, point_front)

        if outOfBounds_right == True:
            len_right_high = len_right_mid
        else:
            len_right_low = len_right_mid

        if outOfBounds_left == True:
            len_left_high = len_left_mid
        else:
            len_left_low = len_left_mid

        if outOfBounds_front == True:
            len_front_high = len_front_mid
        else:
            len_front_low = len_front_mid


        #print(len_right_mid, len_left_mid, len_front_mid)

    dist = [len_right_mid, len_left_mid, len_front_mid]
    #point = [point_right, point_left, point_front]

    return dist

def InitializeModels(sizes, numberOfANNs):
    """
    Initialize ANNs
    """

    ANNs = []
    initialData_chromosome = []

    for nANNs in range(numberOfANNs):

        initialData_tmp, initialData_chromosome_tmp = ComputeInitialData(sizes)
        ANN = NeuralNetworkClass(sizes, initialData_tmp)
        ANNs.append(ANN)
        initialData_chromosome.append(initialData_chromosome_tmp)

    return ANNs, initialData_chromosome

def DecodeChromosome(chromosome):
    """
    Decode chromosome
    """

    decimal = 0
    dimension = len(np.shape(chromosome))
    if dimension == 2:
        len_chrome = np.shape(chromosome)[1]

        for k in range(1, len_chrome):

            decimal = decimal + chromosome[:, k] * 2**(len_chrome - k - len_chrome/2. - 1)

        decimal = decimal * (-1)**chromosome[:, 0]

    else:
        len_chrome = np.shape(chromosome)[2]
        for k in range(1, len_chrome):

            decimal = decimal + chromosome[:, :, k] * 2**(len_chrome - k - len_chrome/2. - 1)

        decimal = decimal * (-1)**chromosome[:, :, 0]

    return decimal

def ComputeInitialData(sizes):
    """
    Compute initial data
    """

    numLayers = len(sizes)
    chromosome_size = 8
    chromosome_bias = []
    chromosome_weight = []
    biases = []
    weights = []

    for k in range(0, numLayers - 1):

        bias   = np.random.randint(0, 2, [sizes[k + 1], chromosome_size])
        weight = np.random.randint(0, 2, [sizes[k + 1], sizes[k], chromosome_size])
        chromosome_bias.append(bias)
        chromosome_weight.append(weight)
        bias_decimal = DecodeChromosome(bias)
        weight_decimal = DecodeChromosome(weight)

        biases.append(bias_decimal)
        weights.append(weight_decimal)

    initialData = [biases, weights]
    initialData_chromosome = [chromosome_bias, chromosome_weight]

    return initialData, initialData_chromosome

def EvaluateNN(ANNs, dist, speed, rotation):
    """
    Evaluate artificial neural network
    """

    inputs = [dist[0], dist[1], dist[2], speed]

    outputs = ANNs.Feedforward(inputs)
    d_rotation = (outputs[0] - 0.5)/30.
    d_speed = (outputs[1] - 0.5)/150.
    rotation = rotation + d_rotation
    speed = max(speed + d_speed, 0.)

    return rotation, speed

def UpdatePoints(points, rotation):
    """
    Upodate location of points
    """

    points[0] = points[0] + np.cos(rotation) * speed
    points[1] = points[1] + np.sin(rotation) * speed

    return points

def CheckOutOfBounds(domain, points):
    """
    Check if position is out of bounds
    """

    outOfBounds = False
    if points[1] > domain[0][2](points[0]):
        outOfBounds = True
    elif points[1] < domain[2][2](points[0]):
        outOfBounds = True
    elif points[1] < domain[1][2](points[0]) and points[1] > domain[3][2](points[0]) and points[0] < 8. and points[0] > -8.:
        outOfBounds = True
    elif points[0] > 11.:
        outOfBounds = True
    elif points[0] < -11.:
        outOfBounds = True


    return outOfBounds

def UpdateANNs(ANNs, data_chromosome, fitness, sizes):
    """
    Update ANNs through crossovers and mutation
    """

    #print('############ Update ANNs ############')
    numLayers = len(sizes)

    # Selection
    selected_chromosomes = Selection('Tournament', data_chromosome, fitness)

    # Crossover
    crossover_chromosomes = Crossover(data_chromosome, selected_chromosomes)

    # Mutation
    chromosomes_new = Mutation(crossover_chromosomes)

    ANNs_new = []
    #print('np.shape(chromosomes_new[k][0][l]) = ', np.shape(chromosomes_new[0][0]))
    for k in range(len(ANNs)):

        bias_vec = []
        weight_vec = []
        print(np.shape(chromosomes_new[k][0]))
        for l in range(0, numLayers - 1):
            print('l = ', l)
            bias           = chromosomes_new[k][0][l]
            weight         = chromosomes_new[k][1][l]
            bias_decimal   = DecodeChromosome(bias)
            weight_decimal = DecodeChromosome(weight)
            bias_vec.append(bias_decimal)
            weight_vec.append(weight_decimal)

        ANNs_tmp = NeuralNetworkClass(sizes, [bias_vec, weight_vec])
        ANNs_new.append(ANNs_tmp)

    return ANNs_new, chromosomes_new

def EncodeChromosomes(ANNs):
    """
    Encode Chromosome
    """
    chromosomes = []
    lenANN = 1
    for k in range(0, lenANN):

        chromosomes[k] = CodeChromosome(ANNs[k])

    return chromosomes

def CodeChromosome(ANN):
    """

    """

    return chromosome
def Selection(method, chromosomes, fitness):
    """

    """

    fitnessLen = len(fitness)
    #print('fitness = ', fitness)
    if method == 'Tournament':
        # Normalize fitness vector
        fitness = fitness/sum(fitness)

        # Sort in descending order
        fitness_sorted_ind = [x for _,x in sorted(zip(fitness, range(fitnessLen)))]
        fitness_sorted = sorted(fitness)

        # Cumulative sum
        fitness_sorted = np.cumsum(fitness_sorted)

        # Draw random numbers in [0, 1]
        randNum = np.random.uniform(0, 1, fitnessLen * 2)

        # Find parents
        parents = FindParents(fitness_sorted, fitness_sorted_ind, randNum)

    return parents

def FindParents(fitnessSorted, fitnessSortedInd, randomNumbers):
    """
    Find parents through tournament selection
    """
    parentsVec = np.zeros(len(randomNumbers))
    for k in range(len(randomNumbers)):

        for l in range(len(fitnessSorted)):

            if randomNumbers[k] < fitnessSorted[l]:
                parentsVec[k] = fitnessSortedInd[l]
                break

    return parentsVec

def Crossover(chromosomes, parents):
    """

    """

    chromosomes_new = []
    for k in range(int(len(parents)/2)):

        chromosomes_tmp = CrossParents(chromosomes[int(parents[k*2])], chromosomes[int(parents[k*2 + 1])])
        chromosomes_new.append(chromosomes_tmp)

    return chromosomes_new

def CrossParents(parent1, parent2):
    """
    Cross parents
    """

    bias_parent1   = parent1[0]
    bias_parent2   = parent2[0]
    weight_parent1 = parent1[1]
    weight_parent2 = parent2[1]
    num_layers = len(parent1)

    lower_cutoff = 6
    upper_cutoff = 2
    bias_vec = []
    weight_vec = []

    for k in range(0, num_layers):

        bias_new = bias_parent1[k]
        bias_new[upper_cutoff:lower_cutoff] =  bias_parent2[k][upper_cutoff:lower_cutoff]
        bias_vec.append(bias_new)

        weight_new = weight_parent1[k]
        weight_new[:, upper_cutoff:lower_cutoff] =  weight_parent2[k][:, upper_cutoff:lower_cutoff]
        weight_vec.append(weight_new)

    parents = [bias_vec, weight_vec]

    return parents

def Mutation(chromosomes):
    """
    Mutate with probability mutate_prob
    """
    mutate_prob = 0.01
    
    num_chromosomes = np.shape(chromosomes)[0]
    for m in range(0, num_chromosomes):

        num_layers = np.shape(chromosomes[m][0])[0]

        for n in range(0, num_layers):
            
            num_neurons_layer = np.shape(chromosomes[m][0][n])[0]

            for k in range(0, num_neurons_layer):
               
                num_chrom_size = np.size(chromosomes[m][0][n][k])

                for l in range(0, num_chrom_size):

                    randNum = np.random.uniform(0, 1)
                    
                    if randNum < mutate_prob:
                        if chromosomes[m][0][n][k][l] == 1:
                            chromosomes[m][0][n][k][l] = 0
                        else:
                            chromosomes[m][0][n][k][l] = 1
                    

    return chromosomes

def DecodeGenes(chromosomes):
    """

    """
    return chromosomes

def ComputeTrackDistance(points):
    """
    Compute track distance
    """

    startingPoint = [-1.0, 9.5]
    u = startingPoint
    u_n = u/(np.sqrt(u[0]**2 + u[1]**2))
    v = points/(np.sqrt(u[0]**2 + u[1]**2))
    v_n = v/(np.sqrt(v[0]**2 + v[1]**2))
    argument = (u_n[0]*v_n[0] + u_n[1]*v_n[1])
    angle = np.arccos(argument)

    return angle

# Initialization
numberOfANNs = 200
numberOfGenerations = 20
time_samples = 500
time_range = np.linspace(0, 1, time_samples)
sizes = [4, 10, 2, 2]
ANNs, data_chromosome = InitializeModels(sizes, numberOfANNs)
domain = CreateDomain('Circle')


for nGens in range(numberOfGenerations):
    trackDistance = np.zeros(numberOfANNs)
    bestDistance = 0
    bestANN = 0
    bestPath = []
    totalDistanceGen = 0
    for nANNs in range(numberOfANNs):
        #print('######### Next ANN #########')
        #print('nANNs = ', nANNs)

        points = [0, 9.5]
        rotation = 0 #np.pi - 0.01
        speed = 1/100.
        pointsVec = np.zeros((2, time_samples + 1))
        rotationVec = np.zeros((1, time_samples + 1))
        pointsVec[:, 0] = points
        rotationVec[:, 0] = rotation

        k = 0
        outOfBounds = False
        totalDistance = 0

        for t in time_range:
            k = k + 1

            # Compute the distances to the walls
            dist = ComputeDistance(points, rotation, domain)

            # Compute speed and rotation by Evaluating the model
            rotation, speed = EvaluateNN(ANNs[nANNs], dist, speed, rotation)
            #print('speed', speed)
            # Update points
            points = UpdatePoints(points, rotation)
            pointsVec[:, k] = points
            rotationVec[:, k] = rotation
            outOfBounds = CheckOutOfBounds(domain, points)

            
            if outOfBounds == True:
                totalDistance = k
                break

        if outOfBounds == False:
            totalDistance = k

        distance = ComputeTrackDistance(pointsVec[:, totalDistance])
        trackDistance[nANNs] = distance
    
        totalDistanceGen = totalDistanceGen + distance
        
        if distance > bestDistance:
            bestDistance = distance
            bestANN = nANNs
            #bestPath = pointsVec[:,:totalDistance+1]
            #bestAngle = rotationVec[:totalDistance+1]
            bestPath = pointsVec[:,:]
            bestAngle = rotationVec[:]

        #print(pointsVec[:, k])
        #print(trackDistance)
    
    print('total distance = ', totalDistanceGen)
    print('best distance =', bestDistance)
    ANNs, data_chromosome = UpdateANNs(ANNs, data_chromosome, trackDistance, sizes)
    #PlotDomain(domain, bestPath, bestAngle)

#PlotDomain(domain, pointsVec[:,:totalDistance+1])
PlotDomain(domain, bestPath, bestAngle)

import numpy as np
import matplotlib.pyplot as plt

def GenerateData(N, dist_min, dist_max):

    x1 = np.random.uniform(dist_min, dist_max, N)
    x2 = np.random.uniform(dist_min, dist_max, N)
    Y  = np.sqrt(x1**2 + x2**2) < 0.5

    X = [x1, x2]
    X = np.reshape(np.ravel(X), (N, -1), order='F')

    return X, Y

def CreateSineData(n_data, mean_data, std_data):

    x_size = 10
    x_vals = np.array(np.linspace(0, 1, x_size))

    rand_on  = np.random.normal(1, 1, (int(n_data/2), x_size))
    rand_off = np.random.normal(1, 1, (int(n_data/2), x_size))

    x_train_on  = mean_data[0] * np.sin(x_vals * np.pi) + rand_on
    x_train_off = mean_data[1] * np.sin(x_vals * np.pi) + rand_off

    y_train_on  = np.ones(int(n_data/2))
    y_train_off = np.zeros(int(n_data/2))
    
    X = np.concatenate((x_train_on, x_train_off))
    Y = np.concatenate((y_train_on, y_train_off))

    return X, Y

def GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max):

    X_train, Y_train = CreateSineData(n_train_points, dist_min, dist_max)
    X_eval, Y_eval   = CreateSineData(n_eval_points, dist_min, dist_max)

    return X_train, Y_train, X_eval, Y_eval

def VisualizeData(X, Y):
    plt.plot(X[Y < 0.5, 0], X[Y < 0.5, 1], 'b*')
    plt.plot(X[Y > 0.5, 0], X[Y > 0.5, 1], 'r*')
    plt.show()

def VisualizeMLData(X_train, Y_train, X_eval, Y_predicted):

    VisualizeData(X_train, Y_train)
    VisualizeData(X_eval,  Y_predicted)

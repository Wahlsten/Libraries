import tensorflow as tf
from keras.layers import Dense, Input
from keras.models import Sequential
from sklearn import metrics
from DataGeneration import GenerateMLData
import numpy as np
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def ComputeANNClassifier(X_train, Y_train, X_eval):

    n_size_in  = np.shape(X_train)
    n_inputs   = n_size_in[1]
    n_classes  = len(np.unique(Y_train))

    model = Sequential()
    model.add(Input(shape=(n_inputs,)))
    model.add(Dense(20, activation='relu'))
    model.add(Dense(10, activation = 'sigmoid'))
    model.add(Dense(n_classes, activation = 'softmax'))

    model.compile(optimizer='adam',
                  loss = tf.keras.losses.SparseCategoricalCrossentropy(),
                  metrics = ['accuracy'])

    model.fit(X_train, Y_train, epochs=10, verbose=0)
    
    Y_predicted = model.predict(X_eval, verbose=0)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 10000
    n_eval_points  = 10000
    dist_min = -1
    dist_max =  1

    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    Y_predicted = ComputeANNClassifier(X_train, Y_train, X_eval)

    #VisualizeMLData(X_train, Y_train, X_eval, Y_predicted[:, 1])
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5))
import tensorflow as tf
import numpy as np
import os
from keras.layers import Dense, LSTM
from keras.models import Sequential
from sklearn import metrics
from DataGeneration import GenerateMLData

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def ComputeLSTMClassifier(X_train, Y_train, X_eval):

    n_size_in = np.shape(X_train)
    n_inputs  = n_size_in[1]
    n_classes = len(np.unique(Y_train))

    model = Sequential()
    model.add(LSTM( 10, input_shape = (n_inputs, 1)))
    model.add(Dense(20, activation  = 'relu'))
    model.add(Dense(n_classes, activation = 'softmax'))

    model.compile(optimizer = 'adam',
                  loss      = tf.keras.losses.SparseCategoricalCrossentropy(),
                  metrics   = ['accuracy'])

    model.fit(X_train, Y_train, epochs=10, verbose=0)

    Y_predicted = model.predict(X_eval, verbose=0)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 10000
    n_eval_points  = 10000
    dist_min       = -1
    dist_max       =  1

    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    Y_predicted                      = ComputeLSTMClassifier(X_train, Y_train,       X_eval)

    #VisualizeMLData(X_train, Y_train, X_eval, Y_predicted[:, 1])
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5))

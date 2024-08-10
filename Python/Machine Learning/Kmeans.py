from sklearn.cluster import KMeans
from sklearn import metrics
from DataGeneration import GenerateMLData, VisualizeMLData
import numpy as np

def ComputeKmeansClassifier(X_train, Y_train, X_eval):

    # Create a svm Classifier
    k_means = KMeans(n_clusters=2, random_state=0, n_init="auto")

    # Train the model using the training sets
    k_means.fit(X_train, Y_train)

    # Predict the response for test dataset
    Y_predicted = k_means.predict(X_eval)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 10000
    n_eval_points  = 10000
    dist_min = -1
    dist_max =  1

    mean_data = [8, 6]
    std_data  = [2, 2]

    #X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, mean_data, std_data)
    Y_predicted = ComputeKmeansClassifier(X_train, Y_train, X_eval)

    #VisualizeMLData(X_train, Y_train.astype(float), X_eval, Y_predicted)
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted))
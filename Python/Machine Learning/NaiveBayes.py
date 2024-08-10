from sklearn.naive_bayes import GaussianNB
from sklearn import metrics
from DataGeneration import GenerateMLData, VisualizeMLData

def ComputeNaiveBayesClassifier(X_train, Y_train, X_eval):

    # Create a svm Classifier
    naive_bayes_model = GaussianNB()

    # Train the model using the training sets
    naive_bayes_model.fit(X_train, Y_train)

    # Predict the response for test dataset
    Y_predicted = naive_bayes_model.predict(X_eval)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 100000
    n_eval_points  = 100000
    #dist_min = -1
    #dist_max =  1

    dist_min = [8, 6]
    dist_max  = [2, 2]

    #X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    #X_train, Y_train = CreateData(n_train_points, mean_data, std_data)
    #X_eval, Y_eval   = CreateData(n_train_points, mean_data, std_data)
    Y_predicted = ComputeNaiveBayesClassifier(X_train, Y_train, X_eval)

    #VisualizeMLData(X_train, Y_train.astype(float), X_eval, Y_predicted)
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted))
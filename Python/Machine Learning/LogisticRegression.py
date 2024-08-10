from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from DataGeneration import GenerateMLData, VisualizeMLData

def ComputeLogisticRegressionClassifier(X_train, Y_train, X_eval):

    # Create a svm Classifier
    log_reg_model = LogisticRegression(solver='lbfgs')

    # Train the model using the training sets
    log_reg_model.fit(X_train, Y_train)

    # Predict the response for test dataset
    Y_predicted = log_reg_model.predict(X_eval)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 100
    n_eval_points  = 100
    dist_min = -1
    dist_max =  1

    mean_data = [9, 6]
    std_data  = [2, 2]

    #X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, mean_data, std_data)
    #X_train, Y_train = CreateData(n_train_points, mean_data, std_data)
    #X_eval, Y_eval   = CreateData(n_train_points, mean_data, std_data)
    Y_predicted = ComputeLogisticRegressionClassifier(X_train, Y_train, X_eval)

    #VisualizeMLData(X_train, Y_train.astype(float), X_eval, Y_predicted)
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted))
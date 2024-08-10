from sklearn import svm
from sklearn import metrics
from DataGeneration import GenerateMLData

def ComputeSVMClassifier(X_train, Y_train, X_eval):

    #Create a svm Classifier
    svm_classifier = svm.SVC(kernel='rbf') # Linear Kernel

    #Train the model using the training sets
    svm_classifier.fit(X_train, Y_train)

    #Predict the response for test dataset
    Y_predicted = svm_classifier.predict(X_eval)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 10000
    n_eval_points  = 10000
    dist_min = -1
    dist_max =  1

    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    Y_predicted = ComputeSVMClassifier(X_train, Y_train, X_eval)

    #VisualizeMLData(X_train, Y_train, X_eval, Y_predicted[:, 1])
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted > 0.5))
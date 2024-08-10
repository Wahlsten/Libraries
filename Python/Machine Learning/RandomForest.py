from sklearn import ensemble
from sklearn import metrics
from DataGeneration import GenerateMLData

def ComputeDecisionRandomForest(X_train, Y_train, X_eval):
    clf = ensemble.RandomForestClassifier(criterion = "entropy",
                                          max_depth = 10,
                                          min_samples_leaf = 5)
    clf = clf.fit(X_train, Y_train)

    Y_predicted = clf.predict_proba(X_eval)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 100000
    n_eval_points  = 100000
    dist_min = -1
    dist_max =  1

    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)
    Y_predicted = ComputeDecisionRandomForest(X_train, Y_train, X_eval)

    #VisualizeMLData(X_train, Y_train, X_eval, Y_predicted[:, 1])
    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5))

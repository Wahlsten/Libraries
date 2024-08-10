from sklearn import tree
from sklearn import metrics
from sklearn.tree import export_graphviz
from six import StringIO
from IPython.display import Image  
import pydotplus
from DataGeneration import GenerateMLData

def SaveTree(clf):
    dot_data = StringIO()
    export_graphviz(clf, out_file = dot_data,  
                    filled = True, rounded = True,
                    special_characters = True, feature_names = ['x1', 'x2'], class_names=['0','1'])
    graph = pydotplus.graph_from_dot_data(dot_data.getvalue())

    graph.write_png('tree.png')
    Image(graph.create_png())

def ComputeDecisionTree(X_train, Y_train, X_eval):

    model_tree  = tree.DecisionTreeClassifier(criterion="entropy", 
                                              max_depth = 10,
                                              min_samples_leaf = 5)
    model_tree  = model_tree.fit(X_train, Y_train)
    Y_predicted = model_tree.predict_proba(X_eval)

    return Y_predicted

if __name__ == '__main__':

    n_train_points = 100000
    n_eval_points  = 100000
    dist_min = -1
    dist_max =  1

    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)

    Y_predicted = ComputeDecisionTree(X_train, Y_train, X_eval)

    print("Accuracy:", metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5))

    #VisualizeMLData(X_train, Y_train, X_eval, Y_predicted[:, 1])
    #SaveTree(model_tree)
from DataGeneration import GenerateMLData
from DecisionTree import ComputeDecisionTree
from RandomForest import ComputeDecisionRandomForest
from SupportVectorMachine import ComputeSVMClassifier
from ArtificialNeuralNetwork import ComputeANNClassifier
from ConvolutionalNeuralNetwork import ComputeCNNClassifier
from LongShortTermMemory import ComputeLSTMClassifier
from LogisticRegression import ComputeLogisticRegressionClassifier
from NaiveBayes import ComputeNaiveBayesClassifier
from Kmeans import ComputeKmeansClassifier
from sklearn import metrics
import numpy as np

def test_Accuracy(method):
    
    n_train_points = 10000
    n_eval_points  = 10000
    #dist_min = -1
    #dist_max =  1

    dist_min = [8, 6]
    dist_max = [2, 2]

    X_train, Y_train, X_eval, Y_eval = GenerateMLData(n_train_points, n_eval_points, dist_min, dist_max)

    if method == 'DecisionTree':
        Y_predicted = ComputeDecisionTree(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5)
    elif method == 'RandomForest':
        Y_predicted = ComputeDecisionRandomForest(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5)
    elif method == 'SupportVectorMachine':
        Y_predicted = ComputeSVMClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted > 0.5)
    elif method == 'ArtificialNeuralNetwork':
        Y_predicted = ComputeANNClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5)
    elif method == 'ConvolutionalNeuralNetwork':
        Y_predicted = ComputeCNNClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5)
    elif method == 'LongShortTermMemory':
        Y_predicted = ComputeLSTMClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted[:, 1] > 0.5)
    elif method == 'LogisticRegression':
        Y_predicted = ComputeLogisticRegressionClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted > 0.5)
    elif method == 'NaiveBayes':
        Y_predicted = ComputeNaiveBayesClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted > 0.5)
    elif method == 'Kmeans':
        Y_predicted = ComputeKmeansClassifier(X_train, Y_train, X_eval)
        accuracy = metrics.accuracy_score(Y_eval, Y_predicted > 0.5)

    #print('Accuracy:', accuracy)
    if accuracy < 0.95:
        result = False
    else:
        result = True

    return result

def print_Tests(tests, res):

    print('---------------- Test Results ----------------')

    total_length = 40

    for k in range(len(res)):
        tmpLen = len(tests[k])
        tmpString = tests[k] + (total_length - tmpLen)*'.'
        if res[k] == 1.0:
            result = 'PASS'
        elif res[k] == 0.0:
            result = 'FAIL'
        else:
            result = 'SKIPPED'
        print(tmpString, result)

def almost_Equal(given_value, correct_value, tolerance):
    if given_value > correct_value + tolerance or given_value < correct_value - tolerance:
        print('Given value = ', given_value, ' while the tolerated value = ', correct_value, ' +- ', tolerance)
        return False
    else:
        return True

def test_act_01():
    ''' Testing ...'''
    print('Running test 01... ', end='', flush=True)

    method = 'DecisionTree'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_02():
    ''' Testing ...'''
    print('Running test 02... ', end='', flush=True)

    method = 'RandomForest'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_03():
    ''' Testing ...'''
    print('Running test 03... ', end='', flush=True)

    method = 'SupportVectorMachine'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_04():
    ''' Testing ...'''
    print('Running test 04... ', end='', flush=True)

    method = 'ArtificialNeuralNetwork'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_05():
    ''' Testing ...'''
    print('Running test 05... ', end='', flush=True)

    method = 'ConvolutionalNeuralNetwork'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_06():
    ''' Testing ...'''
    print('Running test 06... ', end='', flush=True)

    method = 'LongShortTermMemory'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_07():
    ''' Testing ...'''
    print('Running test 07... ', end='', flush=True)

    method = 'LogisticRegression'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_08():
    ''' Testing ...'''
    print('Running test 08... ', end='', flush=True)

    method = 'NaiveBayes'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_09():
    ''' Testing ...'''
    print('Running test 09... ', end='', flush=True)

    method = 'Kmeans'

    result = test_Accuracy(method)

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def run_tests(tests):

    test_description = []

    test_description.append('Decision tree')
    test_description.append('Random forest')
    test_description.append('Support vector machine')
    test_description.append('Artificial neural network')
    test_description.append('Convolutional neural network')
    test_description.append('Long short term memory')
    test_description.append('Logistic Regression')
    test_description.append('Naive Bayes')
    test_description.append('K means')
    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run Machine learning tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

if __name__ == '__main__':

    tests_to_run = [0]

    run_tests(tests_to_run)
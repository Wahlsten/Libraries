from RootFinding_class import Solve
import numpy as np

def test_Convergence(method):
    def f (x): return x**3 - 27
    def fp(x): return 3*x**2

    x_n = Solve(1, 10, 10, method, f, fp, 1e-4)

    if almost_Equal(x_n[-1], 3, 0.5) is False:
        result = False
    else:
        result = True

    return result

def test_accuracy(method, accuracy):
    ''' Testing the accuracy '''


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

def test_Accuracy(method):
    def f (x): return (x-2)**4 - 1
    def fp(x): return 4*(x-2)**3

def almost_Equal(given_value, correct_value, tolerance):
    if given_value > correct_value + tolerance or given_value < correct_value - tolerance:
        print('Given value = ', given_value, ' while the tolerated value = ', correct_value, ' +- ', tolerance)
        return False
    else:
        return True

def test_act_01():
    ''' Testing ...'''
    print('Running test 01... ', end='', flush=True)
    #print('Installing XXX...      ', end='', flush=True)

    method = 'Bisection'

    result = True

    if test_Convergence(method) is False:
        result = False

    if result:
        print('\b OK')
    else:
        print('\b FAIL')
    return result

def test_act_02():
    ''' Testing ...'''
    print('Running test 02... ', end='', flush=True)

    method = 'Brent'

    result = True

    if test_Convergence(method) is False:
        result = False

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_03():
    ''' Testing ...'''
    print('Running test 03... ', end='', flush=True)

    method = 'NewtonRaphson'

    result = True

    if test_Convergence(method) is False:
        result = False

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_04():
    ''' Testing ...'''
    print('Running test 04... ', end='', flush=True)

    method = 'RegulaFalsi'

    result = True

    if test_Convergence(method) is False:
        result = False

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_05():
    ''' Testing ...'''
    print('Running test 05... ', end='', flush=True)

    method = 'Secant'

    result = True

    if test_Convergence(method) is False:
        result = False

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def run_tests(tests):

    test_description = []

    test_description.append('Bisection')
    test_description.append('Brent')
    test_description.append('Newton Raphson')
    test_description.append('Regula Falsi')
    test_description.append('Secant')
    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run root finding tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

from Optimization import Optimize
import numpy as np

def test_Convergence(method):
    def f (x): return (x-2)**4 - 1
    def fp(x): return 4*(x-2)**3

    [min, xmin] = Optimize([0, 2.1], 1e-2, 20, method, f, fp)

    if almost_Equal(min, -1, 0.01) is False or \
       almost_Equal(xmin, 2, 0.01) is False:
        result = False
    else:
        result = True

    return result

def test_2DConvergence(method):
    def f (x, y):  return (x-2)*(x-2) + (y + 1)*(y + 1)
    def fp (x, y): return (x-2)*(x-2) + (y + 1)*(y + 1)

    [min, xmin] = Optimize([0, 2.1], 1e-4, 50, method, f, fp)

    if almost_Equal(min,      0, 0.01) == False or \
       almost_Equal(xmin[0],  2, 0.01) == False or \
       almost_Equal(xmin[1], -1, 0.01) == False:
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

    method = 'FibonacciOptimization'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result == True:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_02():
    ''' Testing ...'''
    print('Running test 02... ', end='', flush=True)

    method = 'GoldenSearch'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result == True:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_03():
    ''' Testing ...'''
    print('Running test 03... ', end='', flush=True)

    method = 'NelderMead'

    if test_2DConvergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_04():
    ''' Testing ...'''
    print('Running test 04... ', end='', flush=True)

    method = 'QuadraticInterpolation'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_05():
    ''' Testing ...'''
    print('Running test 05... ', end='', flush=True)

    method = 'SteepestDescent'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def run_tests(tests):

    test_description = []

    test_description.append('Fibonacci')
    test_description.append('Golden search')
    test_description.append('Nelder-Mead')
    test_description.append('Quadratic interpolation')
    test_description.append('Steepest descent')
    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run optimization tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

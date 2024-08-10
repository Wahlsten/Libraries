from ODESolver import ODESolve
import numpy as np

def test_CheckConvergence(method):
    def y_a(t):   return np.sin(t)
    def ft(t, y): return y - np.sin(t) + np.cos(t)
    N = 100

    [t, y] = ODESolve(0, 1, N, method, ft, y_a(0))
    
    error = abs(y[-1] - y_a(1))
    
    if almost_Equal(error, 0, 0.01) == False:
        result = False
    else:
        result = True

    return result

def test_CheckConvergenceRate(method, real_convergence_rate):
    def y_a(t):  return np.sin(t)
    def ft(t, y): return y - np.sin(t) + np.cos(t)
    N = 100

    [t, y]     = ODESolve(0, 1, N,   method, ft, y_a(0))
    [t_2, y_2] = ODESolve(0, 1, 2*N, method, ft, y_a(0))
    
    normed_error   = abs(y[-1]   - y_a(1))
    normed_error_2 = abs(y_2[-1] - y_a(1))
    
    convergence_rate = np.sqrt(normed_error / normed_error_2)

    if almost_Equal(convergence_rate, real_convergence_rate, 0.5) == False:
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

    method = 'EulerForward'
    convergence_rate = 1
    if test_CheckConvergence(method)                       is False or \
       test_CheckConvergenceRate(method, convergence_rate) is False:
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

    method = 'Heun'
    convergence_rate = 2
    if test_CheckConvergence(method)                       is False or \
       test_CheckConvergenceRate(method, convergence_rate) is False:
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

    method = 'RungeKutta4'
    convergence_rate = 4
    if test_CheckConvergence(method)                       is False or \
       test_CheckConvergenceRate(method, convergence_rate) is False:
        result = False
    else:
        result = True

    if result == True:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def run_tests(tests):

    test_description = []

    test_description.append('EulerForward')
    test_description.append('Heun')
    test_description.append('RungeKutta4')

    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run ordinary differential equation solver tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

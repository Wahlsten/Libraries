from Integration_class import IntegrationClass
import numpy as np

print()
def test_Convergence(method):
    def f (x): return np.sin(x)

    Int = IntegrationClass(0, np.pi, 5, method, f, 0.001)
    Int.integrate()

    result = True

    if almost_Equal(Int.I, 2, 0.1) is False:
        result = False

    return result

def test_Convergence2D(method):
    def f (x, y): return (x-2)**4 - 1 + (y + 1)**2
    def fp(x, y): return 4*(x-2)**3 + 2*(y + 1)

    Opt = OptimizationClass([0, 2.1], 1e-4, 20, method, f, fp)
    Opt.optimize()

    result = True

    if almost_Equal(Opt.min, -1, 0.1) is False:
        result = False

    if almost_Equal(Opt.xmin[0], 2, 0.1) is False:
        result = False

    if almost_Equal(Opt.xmin[1], -1, 0.1) is False:
        result = False

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

    method = 'BoolesRule'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_02():
    ''' Testing ...'''
    print('Running test 02... ', end='', flush=True)

    method = 'GaussLegendre'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_03():
    ''' Testing ...'''
    print('Running test 03... ', end='', flush=True)

    method = 'GaussLobatto'

    if test_Convergence(method) is False:
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

    method = 'MidpointRule'

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

    method = 'Romberg'

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_06():
    ''' Testing ...'''
    print('Running test 06... ', end='', flush=True)

    method = 'SimpsonsRule'

    result = True

    if test_Convergence(method) is False:
        result = False
    else:
        result = True

    if result:
        print('\b OK')
    else:
        print('\b FAIL')

    return result

def test_act_07():
    ''' Testing ...'''
    print('Running test 07... ', end='', flush=True)

    method = 'TrapeziodalRule'

    result = True

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

    test_description.append('Booles rule')
    test_description.append('Gauss-Legendre')
    test_description.append('Gauss-Lobatto')
    test_description.append('Midpoint rule')
    test_description.append('Romberg')
    test_description.append('Simpsons rule')
    test_description.append('Trapeziodal rule')
    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run numerical integration tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

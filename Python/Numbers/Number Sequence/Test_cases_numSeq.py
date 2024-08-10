from NumberSequence import GetSequence
import numpy as np

def test_CheckSequence(method):
    N = 10
    sequence = GetSequence(method, N)
    
    if method == 'Fibonacci':
        correct_sequence = np.array([1, 1, 2, 3, 5, 8, 13, 21, 34, 55])
    else:
        correct_sequence = np.array([0.0, 0.1, 0.1111111111111111, 0.125, 0.14285714285714285, \
                                    0.16666666666666666, 0.2, 0.2222222222222222, 0.25, \
                                    0.2857142857142857, 0.3, 0.3333333333333333, 0.375, \
                                    0.4, 0.42857142857142855, 0.4444444444444444, 0.5, \
                                    0.5555555555555556, 0.5714285714285714, 0.6, 0.625, \
                                    0.6666666666666666, 0.7, 0.7142857142857143, 0.75, \
                                    0.7777777777777778, 0.8, 0.8333333333333334, \
                                    0.8571428571428571, 0.875, 0.8888888888888888, 0.9, 1.0])
        
    n_sequence         = len(sequence)
    n_correct_sequence = len(correct_sequence)

    if not (n_sequence == n_correct_sequence):
        result = False
    else:
        result = True
        for k in range(n_sequence):
            if not (sequence[k] == correct_sequence[k]):
                result = False
                break

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

    method = 'Fibonacci'

    if test_CheckSequence(method) is False:
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

    method = 'Farey'

    if test_CheckSequence(method) is False:
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

    test_description.append('Fibonacci')
    test_description.append('Farey')
    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run number sequence tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

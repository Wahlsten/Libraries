from NumberAlg import CheckSequence
import numpy as np

def test_CheckSequence(method):

    check_true  = CheckSequence(7905125121, method)
    check_false = CheckSequence(7905125122, method)

    if check_true == False or \
       check_false == True:
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

    method = 'LuhnAlgorithm'

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

    test_description.append('Luhn algorithm')
    numTests = len(test_description)
    res = -np.ones(numTests)

    testString = 'test_act_0'

    if tests == [0]:
        tests = range(1, numTests + 1)

    print('---------------- Run number algorithm tests ('+ str(len(tests)) + '/' + str(numTests) + ')' + '----------------')

    for k in tests:

        res[k-1] = eval(testString + str(k) + '()')

    print_Tests(test_description, res)

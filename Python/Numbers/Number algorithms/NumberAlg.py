import numpy as np
from LuhnAlgorithm import LuhnAlgorithm

def CheckSequence(number, method):

    if method == 'LuhnAlgorithm':
        check = LuhnAlgorithm(number)

    return check

if __name__ == '__main__':
    check = CheckSequence(7905125121, 'LuhnAlgorithm')

    print(check)

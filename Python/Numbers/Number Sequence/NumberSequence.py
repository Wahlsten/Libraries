import numpy as np
from Fibonacci import Fibonacci
from Farey import Farey
    
def GetSequence(method, N):

    if method == 'Fibonacci':
        sequence = Fibonacci(N)
    elif method == 'Farey':
        sequence = Farey(N)
    
    return sequence

if __name__ == '__main__':
    sequence = GetSequence('Fibonacci', 10)

    print('Sequence: ', sequence)
    print('Length: ', len(sequence))

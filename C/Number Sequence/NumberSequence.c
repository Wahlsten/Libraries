#include <stdio.h>

int * Fibonacci(int n);

int main(void)
{
    int n = 10;
    int *numbers;

    numbers = Fibonacci(n);
	
    for (int i = 0; i < n; i++)
    {
	printf("Fibonacci(%i) = %i \n", i, numbers[i]);
    }
}

int * Fibonacci(int n) {
    
    static int sequence[10];

    if (n >= 1) 
    {
	sequence[0] = 1;
    }

    if (n >= 2)
    {
	sequence[1] = 1;
    }

    if (n > 2) 
    {

	for (int i = 2; i<n; i++)
	{
	    sequence[i] = sequence[i-1] + sequence[i-2];
	}
    }
	
    return sequence;

}

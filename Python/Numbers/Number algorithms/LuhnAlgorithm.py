import numpy as np

def LuhnAlgorithm(number):
   '''Luhn method'''

   strNumber = str(number)
   numberSum = 0
   k = 0
   for n in reversed(strNumber[:-1]):
      k = k + 1
      if np.mod(k,2) == 1:
         num = str(int(n)*2)
         for l in num:
            numberSum = numberSum + int(l)
      else:
         numberSum = numberSum + int(n)

   numberSum = numberSum + int(strNumber[-1])

   if np.mod(numberSum, 10) == 0:
      check = True
   else:
      check = False

   return check

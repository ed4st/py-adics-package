class Number:
  INFINITY = 2147483647
  p = 0 #prime p
  __m = 0  #minimum positive index where the sum can start 
  __M = 0  #maximum positive index where the sum finish 
  n = 0 #negative number such that the norm of the number is greater that p^n
  N = 0 #positive number such that the norm of the number is lower that p^N
  digits = [] #contains the digits of a p-adic number
  def __init__(self,p ,n ,N, digits):
    self.p = p 
    self.n = n 
    self.N = N
    self.__m = self.N
    self.__M = -self.n 
    self.digits = digits 

  #---------------Getters and setters-------------------------

  def getm(self):
    return self.__m
  def setn(self, m):
    self.__m = m

  def getM(self):
    return self.__M
  def setM(self, M):
    self.__M = M

  #-------------------------------------------------------
  def show(self):#shows the digits of current number
    print(self.digits)
  
  def order(self): #computes the order of a p-adic number
                   #by taking the index of the first nonzero digit
    
    c = 0
    for i in self.digits:
      if(i == 0):
        c += 1
    if(c == len(self.digits)):
      return self.INFINITY
    digits_reverse = self.digits[::-1]
    
    i = 0
    while(digits_reverse[i] == 0 and(i+1 !=  self.__m+self.__M+1)):     
      i += 1
    digits_reverse = self.digits[::-1]
    return -int(self.__m) + i
  
  def len(self): #quantity of digits
        return self.__M + self.__m +1
      
  def norm(self): #computes the norm of the p-adic number
    c = 0
    for i in self.digits:
      if(i == 0):
        c += 1
    if(c == len(self.digits)):
      return 0

    return self.p**(-self.order())  

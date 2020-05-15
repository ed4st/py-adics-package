#NOTE: We suppose that all Numbers in GmMp have the same Number of digits
#therefore, following algorithms follow this fact 
from padics.Number import Number
import random
from random import randint
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import networkx as nx
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import gif
from matplotlib import cm
from numpy import linspace

class GpnN:
  p = 0 #prime p
  __m = 0  #minimum positive index where the sum can start 
  __M = 0  #maximum positive index where the sum finish 
  n = 0 #negative number such that the norm of the number is greater that p^n
  N = 0 #positive number such that the norm of the number is lower that p^N
  numbers = []


  __all_solutions_set = set()
  __all_solutions_dict = dict() # a dictionary that contains sorted solutions as keys 
                                # and a colors as values

  def __init__(self,p ,n ,N):
    self.p = p 
    self.n = n 
    self.N = N
    self.__m = self.N
    self.__M = -self.n 


#---------------Getters and setters-------------------------

  def getm(self):
    return self.__m
  def setm(self, m):
    self.__m = m

  def getM(self):
    return self.__M
  def setM(self, M):
    self.__M = M 
#-------------------------numbers_initialization-------------------------- 
  def generate_numbers(self):
    to_permute  =  [i for i in range(self.p)]
    perms  =  product(to_permute,repeat = self.__m + self.__M+1)
    i  =  0
    numbers_aux  =  []
    for p in perms:
      arr  =  []
      arr.extend(p) 
          
      numbers_aux.append(Number(self.p,self.n,self.N,arr[::-1])) 
          
      self.numbers  =  numbers_aux
      i  =  i + 1
#-------------------------GmMp_export-------------------------- 
  def GmMp_export(self, name):
    f  =  open(name + ".txt","w+")
    for num in self.numbers:
      f.write(" ".join(str(num.digits)) + "   "+str(num.norm()) +"\n")
    
    f.close()
#-------------------------console_printing--------------------------
  def console_printing(self):
    for i in self.numbers:
      i.show()
#-------------------------sum---------------------------
  #following function computes the sum of two p-adic 
  #numbers and truncates the result to the first M-m digits
  def p_sum(self, n1: Number, n2 : Number):
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]

    psum = []
    carry  =  0
    for i in range (len(n1.digits)):
      sum = n1.digits[i]+n2.digits[i]+ carry
      r = sum %n1.p
      psum.append(r)
      carry  =  int((sum)/n1.p)
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]
    return Number(n1.p, n1.n, n1.N, psum[::-1] )
#-------------------------substraction-----------------------
  #following function computes the substraction of two
  #p-adic numbers and truncates the result to the first M-m digits
  def p_sub(self,n1: Number, n2 : Number):
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  
    result = []#[0 for i in range(n1.digits)][1,0]
    carry  =  0
    for i in range (len(n1.digits)):
      sub  =  n1.digits[i]-n2.digits[i]+carry
      #sub = 0 - 1 + 0 = -1 
       
      r  =  sub%n1.p
      #r = -1%2 = 1 
      
      result.append(r)
       
      carry  =  int((sub-r)/n1.p)
      # int(-1/2)=0

    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  
    return Number(n1.p, n1.n, n1.N, result[::-1] )
#-------------------------product-------------------------
  #following function computes the multiplication of two
  #p-adic numbers and truncates the result to the first M-m digits
  def p_mul(self, n1: Number, n2 : Number):

    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  

    resultados  =  []
    for i in range(len(n2.digits)):
      resultado  =  []
      carry  =  0
      for j in range(len(n2.digits)):
        prod  =  n2.digits[j] * n1.digits[i] + carry
        resultado.append(prod % n2.p)
        carry  =  int(prod / n2.p)
      resultados.append(resultado)
      
    aux  =  Number(n1.p,n1.n,n1.N,resultados[0])
    for i in range(1,len(resultados),1):
      for j in range(i):
        resultados[i].pop(len(resultados[i])-1)
        resultados[i].insert(0,0)
      aux  = self.p_sum(aux,Number(n1.p,n1.n,n1.N,resultados[i]))
    aux.digits = aux.digits[::-1]
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  
    return(aux)
#-------------------------division---------------------------
  #following function returns the solution of congruence equation ax = b(mod p)
  #where a,b are the first nonzero digit of divisor and dividend
  #respectively. This solution can be given because p is prime

  def first_dividend_digit_anihilator(self, dividend : Number,divisor : Number):
    dr  =  divisor.digits[0]
    dn  =  dividend.digits[0]
    if(dr == dn and(dr != 0 and dn != 0)):
      return 1
    elif(dn == 0 ):
      return 0
    else:          
      for i in range(dividend.p):
        if((dr*i-dn)and((dr*i-dn)%dividend.p)  ==  0):
          return i
        if(dr*i-dn == 0):
          return i
  #following funtion returns p-adic division when p is prime. 
  # In other case we cannot guarantee the correct result of 
  # first_dividend_digit_anihilator function
  def p_div(self,dividend: Number, divisor : Number):
    dividend_copy = dividend
    divisor_copy = divisor
    dividend_copy .digits = dividend_copy .digits[::-1]  
    divisor_copy .digits = divisor_copy .digits[::-1]  
    i = 0
    while(divisor_copy .digits[i] == 0):
      divisor_copy .digits.pop(0)
      divisor_copy .digits.append(0)
    resul = []
    for i in range(len(dividend_copy .digits)):
      r_temp  =  self.first_dividend_digit_anihilator(dividend_copy ,divisor_copy )
      c_tem_digits  =  [0 for i in range(len(dividend_copy .digits))]
      c_tem_digits[len(dividend_copy .digits)-1]  =  r_temp
      
      
      c_temp_number  =  Number(dividend_copy .p, dividend_copy .n, dividend_copy .N,c_tem_digits)
      divisor_copy_reverse = Number(divisor_copy .p,divisor_copy .n,divisor_copy .N,divisor_copy .digits[::-1])

      
      dividend_copy .digits = dividend_copy .digits[::-1] 
      product_aux = self.p_mul(c_temp_number,divisor_copy_reverse)
      substraction_aux = self.p_sub(dividend_copy ,product_aux)
      substraction_aux.digits.pop(len(substraction_aux.digits)-1)
      substraction_aux.digits.insert(0,0)

      dividend_copy  = substraction_aux
      dividend_copy .digits = dividend_copy .digits[::-1]
      resul.append(r_temp)
    
    return Number(dividend_copy .p,dividend_copy .n,dividend_copy .N,resul[::-1])     
#-------------------------inverse------------------------
  def p_inverse(self,n : Number):
    identity_digits  =  [0 for i in range(len(n.digits))]
    identity_digits[self.__M] = 1
    identity  =  Number(n.p,n.n,n.N,identity_digits)
    return self.p_div(identity,n)
#-------------------------GmMp_representation_tree---------------------------
  def representation_tree(self):
    try:
        import pygraphviz
        from networkx.drawing.nx_agraph import graphviz_layout
    except ImportError:
        try:
            import pydot
            from networkx.drawing.nx_pydot import graphviz_layout
        except ImportError:
            raise ImportError("This example needs Graphviz and either "
                              "PyGraphviz or pydot")
    
    G  =  nx.balanced_tree(self.p,self.__m + self.__M + 1)
    for u,v,d in G.edges(data = True):
      d['weight']  =  (v%self.p)/self.p #coloring edges               
    leaf_nodes  =  []
    for node in G.nodes():
        if nx.degree(G,node)  ==  1:
            leaf_nodes.append(node)


    
    norm  =  0.0
    nx.set_node_attributes(G,norm,'norm')
    for i in range(len(G.nodes)):
          G.nodes[i]['norm']  =  0.0
    if(len(self.numbers) == 0):
      print("GPnN have no numbers initialized")
    for i in range(len(leaf_nodes)):
      G.nodes[leaf_nodes[i]]['norm']  =  self.numbers[i].norm()
      
    norms  =  []
    edges,weights  =  zip(*nx.get_edge_attributes(G,'weight').items())
    nodes,norms  =  zip(*nx.get_node_attributes(G,'norm').items())
    
    
    #Creating a list of p colors
    basic_colors = ['k','r','b','gray','green','c','y','m']
    if(self.p>len(basic_colors)-1):
      
      cm_linspace = linspace(0.0, 1.0, self.p-len(basic_colors))
      basic_colors = basic_colors + [cm.brg(x) for x in cm_linspace] 
    color_of_edges = []

    for i in range(0,len(edges),3):
      for j in range(self.p):
        color_of_edges.append(basic_colors[j])

    pos  =  graphviz_layout(G, prog = 'twopi', args = '')
    tam = self.__M + self.__m + self.p
    plt.figure(figsize = (int(tam**1.1),int(tam**1.1)))

    name = 'G' + str(self.p) + '_' + str(abs(self.n)) + str(self.N)
    plt.title(name)
    
    norms_set = set(norms)
    unique_norms = []
    for i in norms_set: 
      unique_norms.append(i)
    #creating a list of colors based 
    # on how many diferents norms are
    start = 0.0
    stop = 1.0
    cm_subsection = linspace(start, stop, len(norms_set)) 
    norms_color_set = [ cm.jet(x) for x in cm_subsection ]
    
    unique_color = []  
    for i in norms_color_set:
      unique_color.append(i)
    
    unique_norms.sort()
    color_norms = []
    for i in norms:
      if(i == 0.0):
        color_norms.append('white')
      else:
        index = unique_norms.index(i)
        color_norms.append(unique_color[index])
    
    cm_aux = ListedColormap(unique_color)#creating a color map
                        
    nx.draw(G, pos, node_size = int(1000/(tam**1.5)), alpha = 0.7, node_color  =  color_norms, edge_color  = color_of_edges,with_labels = False)
    plt.axis('equal')

    #vertical colorbar
    sm = plt.cm.ScalarMappable(cmap = cm_aux)
    sm._A = []

    mn = min(unique_norms)      # colorbar min value
    mx = max(unique_norms)       # colorbar max value
    
    m1 = (mx-mn)/4
    m2 = (mx-mn)/2
    m3 = (mx-mn)*3/4  
    tks = linspace(0,1,2*len(unique_norms)+1)

    cbar = plt.colorbar(sm, ticks = tks)
    
    labels = []

    labels.append('')
    labels.append('$0$')

    for i in range(1,len(unique_norms)):
      labels.append('')
      labels.append('$'+str(self.p)+'^{' + str(i-self.__m-1)+'}$')
    labels.append('')
    #cbar.set_ticks([mn,md,mx])
    #cbar.set_ticklabels([mn,md,mx])
    cbar.ax.set_yticklabels(labels)
    #print(labels)
    cbar.set_label('Norm', rotation=270)
    #plt.show()
    plt.savefig(name + '.png')    
#-------------------------Monna Map---------------------------
  '''Monna map: takes a p-adic number that will be mapped 
  into a positive real number. Following function is going to return a vector
  which entries are the evaluation of such map on every number in the numbers
  atributte of current class'''
  def monna_map(self):
    real_numbers = []
    for pnumber in self.numbers:
      j = pnumber.len() - 1
      real_number = 0
      for i in range(-self.__m, self.__M, 1):
        real_number += (pnumber.p**(-i-1))*pnumber.digits[j] 
        j -= 1
      real_numbers.append(real_number)

    return real_numbers
#--------------------------Parisi_Matrix--------------------------
  #following function returns the (i,j)-th value of Parisi Matrix
  def fij(self, alpha, constant, i: Number, j: Number):
    return constant/(self.p_sub(i,j).norm()**alpha + 1)

  def matrix(self):
    alpha = 2
    constant = 3
    W = []
    for i in self.numbers:
      aux = []
      for j in self.numbers:
        aux.append(self.fij(alpha,constant,i,j))
      W.append(aux)

    sum_row1 = sum(W[0]) 
    for i in range(len(W)):
      W[i] = [j*(sum_row1**(-1)) for j in (W[i])]
    #upgrading the diagonal of W
    row_copy = []
    row_copy = W[0].copy()
    row_copy.pop(0)
    w0 =  - sum(row_copy)
    for i in range(len(W)):  
      W[i][i] = w0 
    
    '''
    #matrix visualization
    fig, ax = plt.subplots()
    ax.matshow(W, cmap = plt.cm.get_cmap("jet"))
    for i in range(len(self.numbers)):
      for j in range(len(self.numbers)):
        c = W[i][j]
        #ax.text(i, j, f"{c:.2f}", va='center', ha='center')#text(i, j, f"{c:.2f}", va='center', ha='center')
    
    #colorbar
    sm = plt.cm.ScalarMappable(cmap = cm.jet)
    sm._A = []

    mn = min(W[0])      # colorbar min value
    mx = max(W[0])       # colorbar max value
    
    md = (mx-mn)/2
    tks = linspace(0,1,3)#2*(self.__m+self.__M+1))

    cbar = plt.colorbar(sm, ticks = tks)
    cbar.ax.set_yticklabels([f"{mn:.2f}",f"{md:.2f}" , f"{mx:.2f}"])
    #plt.show()
    name = 'G'+str(self.p)+'_'+str(-self.n)+str(self.N)
    plt.title('Transition Matrix of ' + name + ' with $\\alpha='+str(alpha)+'$ and $C='+str(constant)+'$')
    plt.savefig('matrix'+name+'.png')'''


    return W
#-------------------------Solving the Master Equation System---------------------------
  def __model(self, u, t):
    W = self.matrix()
    dudts = []
    for i in range(len(W)):
      auxeq = 0
      for j in range(len(W)):        
        auxeq += W[i][j]*u[j]
      dudts.append(auxeq)
    return dudts
  
  def ODESols(self, type_ic = None):
    #time points
    t = np.linspace(0,10,4)

    #initial condition
    u0 = []
    
    if(type_ic != None):
      #random initial condition
      if(type_ic == 'random'):
        u0 = [random.random() for i in range(len(self.numbers))]
        plt.plot(u0)
        plt.show()
      #ones initial condition
      elif(type_ic == 'ones'):
        u0 = [1 for i in range(len(self.numbers))]
        plt.plot(u0)
        plt.show()
    #gauss bell like initial condition by default
    else:
      #gauss bell centered at GpnN order divided by 2
      center = len(self.numbers)/2
      #creating a partition on the interval [center-1,center+1] 
      #of size |GpnN|  
      partition = linspace(center-0.5 ,center+0.5,len(self.numbers))
      #here we evaluate the partition created above in
      #gauss bell function centered at GpnN order divided by 2 
      u0= np.exp((partition-center)**2)/self.p
      #plt.plot(u0)
      #plt.show()

    u = odeint(self.__model,u0,t)
    for u_i in u:
      self.__all_solutions_set.update(set(u_i))
    self.__color_solutions()
    return [u,t]

  def __color_solutions(self):
    all_sols_lenght =  len(self.__all_solutions_set)
    cm_subsection = linspace(0, 1, all_sols_lenght) 
    all_solutions_list = sorted(list(self.__all_solutions_set))
    all_colors = [ cm.jet(x) for x in cm_subsection ]
    
    for i in range(all_sols_lenght):
      key = all_solutions_list[i]
      value = all_colors[i]
      self.__all_solutions_dict[key] = value


  #--------------------------Creating image transitions---------------
  
  @gif.frame
  def animate(self, boundary, time, type_ic = None):
    try:
        import pygraphviz
        from networkx.drawing.nx_agraph import graphviz_layout
    except ImportError:
        try:
            import pydot
            from networkx.drawing.nx_pydot import graphviz_layout
        except ImportError:
            raise ImportError("This example needs Graphviz and either "
                              "PyGraphviz or pydot")
    
    G  =  nx.balanced_tree(self.p,self.__m + self.__M + 1)
    for u,v,d in G.edges(data = True):
      d['weight']  =  (v%self.p)/self.p #coloring edges               
    leaf_nodes  =  []
    for node in G.nodes():
        if nx.degree(G,node)  ==  1:
            leaf_nodes.append(node)
    
    norm  =  0.0
    nx.set_node_attributes(G,norm,'norm')
    for i in range(len(G.nodes)):
          G.nodes[i]['norm']  =  0.0
    
    for i in range(len(leaf_nodes)):
      G.nodes[leaf_nodes[i]]['norm']  =  boundary[i]
      
    norms  =  []
    edges,weights  =  zip(*nx.get_edge_attributes(G,'weight').items())
    nodes,norms  =  zip(*nx.get_node_attributes(G,'norm').items())
    
    
    #Creating a list of p colors
    basic_colors = ['k','r','b','gray','green','c','y','m']
    if(self.p>len(basic_colors)-1):
      cm_linspace = linspace(0.0, 1.0, self.p-len(basic_colors))
      basic_colors = basic_colors + [cm.brg(x) for x in cm_linspace] 
    color_of_edges = []

    for i in range(0,len(edges),3):
      for j in range(self.p):
        color_of_edges.append(basic_colors[j])

    pos  =  graphviz_layout(G, prog = 'twopi', args = '')
    tam = self.__M + self.__m + self.p
    plt.figure(figsize = (int(tam**1.1),int(tam**1.1)))


    name = 'G' + str(self.p) + '_' + str(abs(self.n)) + str(self.N)
    
    if(type_ic != None):
      #random initial condition
      if(type_ic == 'random'):
        plt.title(name + ' difussion model with randomly distributed initial condition vector')
      #ones initial condition
      elif(type_ic == 'ones'):
        plt.title(name + ' difussion model with $\\vec{1}$ initial condition vector')
    #normal initial condition by default
    else:
        #median = 0, variance = 1
        plt.title(name + ' difussion model with normally distributed initial condition vector')
    plt.text(2, 6, '$t =' + f"{time:.2f}" + '$', fontsize=15)

    
    

    boundary_colors = []
    for i in norms:
      if(i == 0.0):
        boundary_colors.append('white')
      else:
        boundary_colors.append(self.__all_solutions_dict[i])
    
                        
    nx.draw(G, pos, node_size = int(1000/(tam**1.5)), alpha = 0.7, node_color  =  boundary_colors , edge_color  = color_of_edges,with_labels = False)
    plt.axis('equal')

    #vertical colorbar
    sm = plt.cm.ScalarMappable(cmap = cm.jet)#cm_aux)
    sm._A = []
    mn = min(self.__all_solutions_set)
    mx = max(self.__all_solutions_set)
    md = (mn+mx)/2
    tks = linspace(0,1,3)
    cbar = plt.colorbar(sm, ticks = tks )
    cbar.ax.set_yticklabels([f"{mn:.3f}", f"{md:.3f}",f"{mx:.3f}"])
    cbar.set_label('solution values', rotation=270)

    #plt.show()
  '''
  following function returns a gif showing
  the behaviour of ultrametric difussion equations: 
  type_ic: type of initial condition(normal,random,ones)
           default = normal
           establishes the values of initial condition vector.
           It can be distributed normally with mean 0 and variance 1,
           randomically and constant(ones in every entry of initial
           condition vector).  
  '''
  def export_gif(self, type_ic = None ):
    frames = []
    u = self.ODESols(type_ic)
    
    for u_i, t_i in zip(u[0],u[1]):
      frame = self.animate(u_i,t_i,type_ic)
      frames.append(frame)

    name = 'G' + str(self.p) + '_' + str(abs(self.n)) + str(self.N)
    if(type_ic!=None):
      gif.save(frames,name + "_" +type_ic +"_difussion.gif",duration=1000)
    else:
      gif.save(frames,"gifs/"+name + "_normal_difussion.gif",duration=1000)
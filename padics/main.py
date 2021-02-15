from Number import Number
from GpnN import GpnN

import matplotlib.pyplot as plt
import time

if __name__ == '__main__':
    start_time = time.time()
    G2_22 = GpnN(2,-1,3) #initialization
    x = Number(7,-2,4,[2,4,6,3,2,2,6])
    G2_22.generate_numbers()
    #print(G2_22.monna_map())
    #G2_22.ODESols()
    #G2_22.matrix()
    G2_22.representation_tree()
    #G2_22.console_printing()

    #G2_22.representation_tree()

    #G2_22.matrix()

    #G2_22.ODESols()
    #G2_22.export_gif()

    print("--- %s seconds ---" % (time.time() - start_time))
  

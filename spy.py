import matplotlib.pyplot as plt

def spy(data):
  plt.close('all')
  f,ax = plt.subplots(1,1)
  ax.spy(data, precision=0.1, markersize=5)
  plt.show()
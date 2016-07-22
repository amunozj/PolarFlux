from zaw_coord import CRD
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


x = CRD('MDI\\fd_M_96m_01d.1222.0005.fits')
area = x.eoa(x.im_raw.data).value.flatten()
area = area[~np.isnan(area)]

loc = np.arange(0, 5e17, 5e17/999)
n, bins, patches = plt.hist(area, loc, log=True)
print(n.max())

plt.axis([0,5e17, 0, 60000])
l = plt.plot(bins, 'r--', linewidth=1)
plt.show()

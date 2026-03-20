import numpy as np
import matplotlib.pyplot as plt


wavelength = np.linspace(1395, 1410, 1000)

def gaussian(x, amp, center, sigma):
    return amp * np.exp(-(x - center)**2 / (2 * sigma**2))


center = 1402.77
amplitude = 1.0
sigma = 0.2

flux = gaussian(wavelength, amplitude, center, sigma)

# Plot
plt.figure(figsize=(8, 5))

plt.plot(wavelength, flux, color="red")

plt.xlabel("Wavelength (Å)")
plt.ylabel("Intensity")
plt.title("Simulated Si IV 1403 Emission Line")

plt.grid(alpha=0.3)

plt.show()

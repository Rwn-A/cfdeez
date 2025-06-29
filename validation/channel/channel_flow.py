import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from scipy.interpolate import griddata
import matplotlib.ticker as ticker

directory = "./validation/channel"

subprocess.run(["./.build/cfdeez.exe", f"{directory}/laminar_flow.fml"])

data = pd.read_csv(f"{directory}/laminar_channel_flow_1.csv")
points = data[["x", "y"]].values
vel_x = data["velocity.x"].values
vel_y = data["velocity.y"].values
vel_mag = np.sqrt(vel_x**2 + vel_y**2)

x_target = 4.5

y_min, y_max = data["y"].min(), data["y"].max()
y_interp = np.linspace(y_min, y_max, 200)
query_points = np.column_stack((np.full_like(y_interp, x_target), y_interp))
vel_mag_interp = griddata(points, vel_mag, query_points, method='linear')

plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(0.05)) 

plt.plot(y_interp, vel_mag_interp, label='Velocity Magnitude')
plt.xlabel('position y (m)')
plt.ylabel('velocity (m / s)')
plt.title(f'velocity profile at x = {x_target}')
plt.grid(True)
plt.savefig(f"{directory}/laminar_channel_velocity_profile{x_target}.png")
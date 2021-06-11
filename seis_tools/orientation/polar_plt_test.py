import matplotlib.pyplot as plt

r = [1, 1,]
theta = [2, 4]

fig = plt.figure(figsize=(5, 5))
ax = plt.subplot(111, projection = 'polar')

for i in range(len(r)):
    ax.plot([0, theta[i]], [0, r[i]])


# Go clockwise
ax.set_theta_direction(-1)
# Start from the top
ax.set_theta_offset(1.570796327)


plt.show()
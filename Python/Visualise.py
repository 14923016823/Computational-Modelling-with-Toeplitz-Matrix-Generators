import matplotlib.pyplot as plt
import numpy as np

with open("b.txt") as f:
    Nx, Ny, Nz = list(map(int, f.readline().split()[1:]))
b = np.loadtxt("b.txt")

class IndexTracker:
    # 1. Pass the desired 'axis' (X, Y, or Z) during initialization
    def __init__(self, ax, B, axis):
        self.index = 0
        self.B = B
        self.ax = ax
        self.axis = axis
        
        # Determine the dimension index (0, 1, or 2) to slice along
        if axis == 'X':
            self.dim = 0
            x_label, y_label = 'Y', 'Z'
        elif axis == 'Y':
            self.dim = 1
            x_label, y_label = 'X', 'Z'
        elif axis == 'Z':
            self.dim = 2
            x_label, y_label = 'X', 'Y'
        else:
            raise ValueError("Axis must be 'X', 'Y', or 'Z'")

        # Set up the initial slicing based on the chosen dimension
        slice_tuple = [slice(None)] * 3
        slice_tuple[self.dim] = self.index
        
        # Display the initial slice and set labels (using correct Matplotlib methods)
        self.im = ax.imshow(self.B[tuple(slice_tuple)], origin='lower')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        self.update()

    def on_scroll(self, event):
        increment = 1 if event.button == 'up' else -1
        
        # 2. Max index must be determined by the selected dimension (self.dim)
        max_index = self.B.shape[self.dim] - 1
        
        # Update and clip the index
        self.index = np.clip(self.index + increment, 0, max_index)
        self.update()

    def update(self):
        # 3. Create a dynamic slicing tuple based on self.dim
        slice_tuple = [slice(None)] * 3
        slice_tuple[self.dim] = self.index
        
        # Set the data using the dynamic slice
        self.im.set_data(self.B[tuple(slice_tuple)])
        
        # Update the title
        self.ax.set_title(
            f'Slice along {self.axis} axis\n{self.axis} Index = {self.index}')
        self.im.axes.figure.canvas.draw()


B = np.reshape(b, (Nz, Ny, Nx))

fig, ax = plt.subplots()
# create an IndexTracker and make sure it lives during the whole
# lifetime of the figure by assigning it to a variable
tracker = IndexTracker(ax, B, input('Axis to scroll along: '))

fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
plt.show()
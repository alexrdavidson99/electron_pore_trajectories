import os
import numpy as np

def append_data_with_header(filename, header, data):
    # Convert data to a 2D array if it's 1D
    if data.ndim == 1:
        data = np.expand_dims(data, axis=1)
    
    # Open the file in append mode
    with open(filename, 'a') as f:
        # Write the header as a comment line
        f.write(f"# {header}\n")
        # Append the new data
        np.savetxt(f, data, delimiter=',', fmt='%s')


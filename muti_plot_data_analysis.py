import numpy as np
import matplotlib.pyplot as plt

def read_data_sections(filename):
    sections = []
    current_section = []
    labels = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                if current_section:
                    sections.append(current_section)
                    current_section = []
                labels.append(line[1:].strip())
            else:
                if line:
                    current_section.append(float(line))
                    
    if current_section:
        sections.append(current_section)
    
    return sections, labels

def plot_histograms(sections, labels):
    plt.figure()
    for i,  (section, label) in enumerate(zip(sections, labels)):
        plt.hist(section, bins=100, alpha=0.75, label=f'{label}', density=True, range=(0,50))
    plt.title(f'Histogram for number of hits with different length diameter ratios')
    plt.xlabel('number of hits')
    plt.ylabel('Frequency')
    
    plt.legend()


# hits 
filename = 'Hits_data'
sections, labels = read_data_sections(filename)
plot_histograms(sections, labels)
plt.show()

# energy
filename = 'energy_data.txt'
sections, labels = read_data_sections(filename)
plot_histograms(sections, labels)
plt.show()

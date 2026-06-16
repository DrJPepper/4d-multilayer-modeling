# A script to plot results for inclusion in dissertation

import sys
import pandas as pd
import matplotlib.pyplot as plt

# Match LaTeX styling
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "font.size": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20,
})

def plot_csv(file_path):
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    if len(df.columns) < 2:
        print("Error: The CSV file must contain at least two columns.")
        return

    x_label = df.columns[0]
    y_label = df.columns[1]
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Scatter plot
    ax.plot(x, y, marker='o', linestyle='-', color='#1f77b4', 
            linewidth=1.5, markersize=6, markerfacecolor='white', 
            markeredgewidth=1.5, zorder=3)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if len(sys.argv) > 2:
        plt.title(sys.argv[2])
    else:
        plt.title(f'{y_label} vs {x_label}')
    
    ax.grid(True, linestyle=':', color='#d3d3d3', zorder=0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.spines['left'].set_color('#333333')
    ax.spines['bottom'].set_color('#333333')

    plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    csv_filename = sys.argv[1]
    plot_csv(csv_filename)

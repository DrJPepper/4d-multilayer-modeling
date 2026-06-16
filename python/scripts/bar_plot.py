# A script to plot results for inclusion in dissertation

import sys
import pandas as pd
import matplotlib.pyplot as plt

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

def plot_csv_bar_with_error(file_path):
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    if len(df.columns) < 3:
        print("Error: The CSV file must contain at least three columns (Label, Value, Standard Deviation).")
        return

    x_name = df.columns[0]
    y_name = df.columns[1]

    x_labels = df.iloc[:, 0]
    y_values = df.iloc[:, 1]
    std_devs = df.iloc[:, 2]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Bar chart
    ax.bar(x_labels, y_values, yerr=std_devs, capsize=5, color='#1f77b4', 
           edgecolor='black', linewidth=1, alpha=0.8, zorder=3)
    
    ax.set_xlabel(x_name, labelpad=20)
    ax.set_ylabel(y_name)
    ax.set_title(f'{y_name} by {x_name}')
    
    ax.grid(True, axis='y', linestyle=':', color='#d3d3d3', zorder=0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.spines['left'].set_color('#333333')
    ax.spines['bottom'].set_color('#333333')

    plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    csv_filename = sys.argv[1]
    plot_csv_bar_with_error(csv_filename)

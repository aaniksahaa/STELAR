import pandas as pd
import matplotlib.pyplot as plt
import os

# Create a directory to save plots
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# Load the CSV files
rf_data = pd.read_csv("average_rfs.csv")
diff_data = pd.read_csv("average_diffs.csv")

# Function to create bar chart with subplots for each inner folder
def plot_combined_subplots(data, title, filename, value_col):
    # Get unique inner folders
    inner_folders = data['inner_folder'].unique()
    num_folders = len(inner_folders)
    
    # Create a subplot grid (adjust rows and cols based on number of inner folders)
    cols = 2  # Number of columns in the subplot grid
    rows = (num_folders + 1) // cols  # Calculate rows needed
    fig, axes = plt.subplots(rows, cols, figsize=(12, rows * 4), sharey=True)
    axes = axes.flatten()  # Flatten for easier iteration

    # Plot each inner folder as a subplot
    for idx, inner_folder in enumerate(inner_folders):
        ax = axes[idx]
        group = data[data['inner_folder'] == inner_folder]
        ax.bar(group['summary_new_method'], group[value_col], 
               color=['#FF9999', '#66B2FF', '#99FF99', '#FFCC99'])
        ax.set_title(inner_folder, fontsize=10)
        ax.set_xlabel("Method")
        ax.set_ylabel(value_col.replace("_", " ").title())
        ax.tick_params(axis='x', rotation=45)

    # Hide unused subplots if any
    for idx in range(num_folders, len(axes)):
        fig.delaxes(axes[idx])

    # Adjust layout and save
    plt.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename), bbox_inches='tight')
    plt.close()

# Function to create average bar chart (same as before)
def plot_avg_bar_chart(data, title, filename, value_col):
    avg_data = data.groupby('summary_new_method')[value_col].mean().reset_index()
    plt.figure(figsize=(12, 6))
    plt.bar(avg_data['summary_new_method'], avg_data[value_col], 
            color=['#FF9999', '#66B2FF', '#99FF99', '#FFCC99'])
    plt.xticks(rotation=45, ha='right')
    plt.title(title)
    plt.xlabel("Method")
    plt.ylabel(value_col.replace("_", " ").title())
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

# Create combined subplot charts
plot_combined_subplots(rf_data, "RF Distance by Inner Folder", "combined_rf_distance_subplots.png", "average_rf")
plot_combined_subplots(diff_data, "Average Differences by Inner Folder", "combined_average_diffs_subplots.png", "average_diff")

# Create average bar charts (same as before)
plot_avg_bar_chart(rf_data, "Average RF Distance Across Inner Folders", "avg_rf_distance_bar.png", "average_rf")
plot_avg_bar_chart(diff_data, "Average Differences Across Inner Folders", "avg_diffs_bar.png", "average_diff")

print(f"Plots saved in {output_dir}")